library(tidyverse)
library(lme4)
library(ggplot2)
library(lmerTest)

if (startsWith(getwd(), "/home/lakrids")) {
  path_prefix <- "/home/lakrids/GenomeDK"
} else {path_prefix <- "/faststorage/project"}

output_dir <- paste0(path_prefix, "/megaFauna/sa_megafauna/results/shared/modeling")
dir.create(output_dir, showWarnings = FALSE)

df <- read_csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/modeling_data.csv")) %>%
  rename(temperature = mean_temp, npp = mean_npp, precipitation = mean_prec) %>%
  arrange(species, time_kya) %>%
  group_by(species) %>%
  mutate(
    lag_temperature   = lag(temperature),
    temperature_s     = as.numeric(scale(temperature)),
    lag_temperature_s = as.numeric(scale(lag_temperature))
  ) %>%
  ungroup() %>%
  mutate(window_size = window_end - time_kya)


# =============================================================================
# FIGURE 1: Ne trajectories for all species
# =============================================================================

fig1 <- df %>%
  filter(!is.na(Ne)) %>%
  ggplot(aes(x = log10(time_kya * 1000), y = log10(Ne), colour = species)) +
  geom_step() +
  scale_x_continuous(
    labels = function(x) parse(text = paste0("10^", x)),
    name   = "Years ago"
  ) +
  scale_y_continuous(
    labels = function(x) parse(text = paste0("10^", x)),
    name   = expression(N[e])
  ) +
  scale_colour_brewer(palette = "Paired", name = "Species") +
  theme_classic(base_size = 12) +
  theme(legend.text = element_text(face = "italic"))

ggsave(paste0(output_dir, "/fig1_ne_trajectories.pdf"), fig1, width = 10, height = 6)
ggsave(paste0(output_dir, "/fig1_ne_trajectories.png"), fig1, width = 10, height = 6, dpi = 300)


# =============================================================================
# MODEL FITTING
# =============================================================================

m_climate <- lmer(
  log(Ne) ~ temperature_s +
    (1 | species),
  data    = df,
  weights = window_size
)

m_human <- lmer(
  log(Ne) ~ p_human +
    (1 | species),
  data    = df,
  weights = window_size
)

m_combined <- lmer(
  log(Ne) ~ temperature_s + p_human +
    (1 | species),
  data    = df,
  weights = window_size
)

m_mass <- lmer(
  log(Ne) ~ temperature_s * log(adult_mass) + p_human * log(adult_mass) +
    (1 | species),
  data    = df,
  weights = window_size
)


# =============================================================================
# TABLE 1: Model comparison
# =============================================================================

model_comparison <- AIC(m_climate, m_human, m_combined, m_mass) %>%
  as.data.frame() %>%
  rownames_to_column("Model") %>%
  mutate(
    Model = recode(Model,
                   m_climate  = "Climate only (temperature)",
                   m_human    = "Human only (p_human)",
                   m_combined = "Combined (temperature + p_human)",
                   m_mass     = "Combined + body mass interactions"
    ),
    delta_AIC = AIC - min(AIC)
  ) %>%
  arrange(AIC) %>%
  rename(
    `Num. parameters` = df,
    `Delta AIC`       = delta_AIC
  )

write_csv(model_comparison, paste0(output_dir, "/table1_model_comparison.csv"))
print(model_comparison)


# =============================================================================
# FIGURE 2: Species-level slopes (dot plot)
# =============================================================================

species_slopes <- coef(m_combined)$species %>%
  rownames_to_column("species") %>%
  rename(
    intercept   = `(Intercept)`,
    temp_slope  = temperature_s,
    human_slope = p_human
  ) %>%
  pivot_longer(
    cols      = c(temp_slope, human_slope),
    names_to  = "predictor",
    values_to = "slope"
  ) %>%
  mutate(
    predictor = recode(predictor,
                       temp_slope  = "Temperature",
                       human_slope = "Human pressure"
    ),
    species = gsub("_", " ", species)
  )

fig2 <- species_slopes %>%
  ggplot(aes(x = slope, y = reorder(species, slope), colour = slope > 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  scale_colour_manual(values = c("TRUE" = "#E07B54", "FALSE" = "#5B8DB8"), guide = "none") +
  facet_wrap(~predictor, scales = "free_x") +
  labs(
    x = "Species-level slope",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(face = "italic"))

ggsave(paste0(output_dir, "/fig2_species_slopes.pdf"), fig2, width = 10, height = 6)
ggsave(paste0(output_dir, "/fig2_species_slopes.png"), fig2, width = 10, height = 6, dpi = 300)


# =============================================================================
# FIGURE 3: Observed vs predicted Ne per species
# =============================================================================

df_pred <- df %>%
  filter(!is.na(Ne)) %>%
  mutate(
    Ne_pred_climate  = exp(predict(m_climate,  newdata = .)),
    Ne_pred_human    = exp(predict(m_human,    newdata = .)),
    Ne_pred_combined = exp(predict(m_combined, newdata = .))
  ) %>%
  select(species, time_kya, Ne, Ne_pred_climate, Ne_pred_human, Ne_pred_combined) %>%
  pivot_longer(
    cols      = c(Ne, Ne_pred_climate, Ne_pred_human, Ne_pred_combined),
    names_to  = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = recode(series,
                    Ne               = "Observed",
                    Ne_pred_climate  = "Climate only",
                    Ne_pred_human    = "Human only",
                    Ne_pred_combined = "Combined"
    ),
    series  = factor(series, levels = c("Observed", "Climate only", "Human only", "Combined")),
    species = gsub("_", " ", species)
  )

fig3 <- df_pred %>%
  ggplot(aes(x = log10(time_kya * 1000), y = log10(value),
             colour = series, linetype = series, linewidth = series)) +
  geom_step() +
  scale_colour_manual(values = c(
    "Observed"     = "black",
    "Climate only" = "#E07B54",
    "Human only"   = "#5B8DB8",
    "Combined"     = "#6BAE75"
  )) +
  scale_linetype_manual(values = c(
    "Observed"     = "solid",
    "Climate only" = "dashed",
    "Human only"   = "dashed",
    "Combined"     = "dashed"
  )) +
  scale_linewidth_manual(values = c(
    "Observed"     = 0.9,
    "Climate only" = 0.6,
    "Human only"   = 0.6,
    "Combined"     = 0.6
  )) +
  scale_x_continuous(labels = function(x) parse(text = paste0("10^", x)), name = "Years ago") +
  scale_y_continuous(labels = function(x) parse(text = paste0("10^", x)), name = expression(N[e])) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  labs(colour = NULL, linetype = NULL, linewidth = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.text      = element_text(face = "italic"),
    legend.position = "bottom"
  )

ggsave(paste0(output_dir, "/fig3_observed_vs_predicted.pdf"), fig3, width = 12, height = 14)
ggsave(paste0(output_dir, "/fig3_observed_vs_predicted.png"), fig3, width = 12, height = 14, dpi = 300)


# =============================================================================
# TABLE 2: Fixed effects summary for combined model
# =============================================================================

fixed_effects <- summary(m_combined)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  mutate(
    Term = recode(Term,
                  `(Intercept)`   = "Intercept",
                  `temperature_s` = "Temperature (scaled)",
                  `p_human`       = "Human pressure"
    )
  ) %>%
  rename(
    `Std. Error` = `Std. Error`,
    `t value`    = `t value`
  )

write_csv(fixed_effects, paste0(output_dir, "/table2_fixed_effects.csv"))
print(fixed_effects)

message("All outputs saved to: ", output_dir)

