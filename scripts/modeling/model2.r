library(tidyverse)
library(lme4)

if (startsWith(getwd(), "/home/lakrids")) {
  path_prefix <- "/home/lakrids/GenomeDK"
} else {path_prefix <- "/faststorage/project"}

df <- read_csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/modeling_data.csv"))

df <- df %>% rename(temperature = mean_temp, npp = mean_npp, precipitation = mean_prec) %>%
  arrange(species, time_kya) %>%
  group_by(species) %>%
  mutate(lag_temperature = lag(temperature), temperature_s = scale(temperature)) %>%
  ungroup() 

df <- df %>%
  mutate(window_size = window_end - time_kya)

# distributions
df %>% 
  pivot_longer(c(temperature, precipitation, npp, Ne)) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")

# correlations
df %>%
  select(temperature, precipitation, npp, p_human, adult_mass, generation_time) %>%
  cor(use = "complete.obs")

m0 <- lmer(log(Ne) ~ 1 + (1 | species), data = df)

m1 <- lmer(
  log(Ne) ~ temperature + precipitation + npp + p_human +
    (1 + temperature + p_human | species),
  data = df,
  weights = window_size
)

m2 <- lmer(
  log(Ne) ~ temperature + lag_temperature + p_human +
    (1 + temperature + lag_temperature + p_human | species),
  data = df,
  weights = window_size
)

m3 <- lmer(
  log(Ne) ~ temperature * log(adult_mass) +
    p_human * log(adult_mass) +
    (1 + temperature + p_human | species),
  data = df,
  weights = window_size
)

m4 <- lmer(
  log(Ne) ~ temperature * log(generation_time) +
    p_human * log(generation_time) +
    (1 + temperature + p_human | species),
  data = df,
  weights = window_size
)

m5 <- lmer(
  log(Ne) ~ temperature * log(adult_mass) +
    p_human * log(adult_mass) +
    temperature * log(generation_time) +
    p_human * log(generation_time) +
    (1 + temperature + p_human | species),
  data = df,
  weights = window_size
)

model_list <- list(m0, m1, m2, m3, m4, m5)
AICtab <- AIC(m0, m1, m2, m3, m4, m5)
AICtab

ggplot(df, aes(x = log10(time_kya), y = log(Ne), colour = species)) +
  geom_step() +
  labs(x = "Log10 Time (kya)", y = "Log Ne")

library(broom.mixed)

fixef_df <- broom.mixed::tidy(m5, effects = "fixed")

ggplot(fixef_df, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = estimate - std.error*1.96,
                     xmax = estimate + std.error*1.96)) +
  geom_vline(xintercept = 0, linetype = "dashed")

ranef_df <- coef(m1)$species %>%
  as.data.frame() %>%
  rownames_to_column("species")

ggplot(ranef_df, aes(x = temperature, y = `(Intercept)`)) +
  geom_point()

library(ggeffects)

pred <- ggpredict(m3, terms = c("temperature", "adult_mass [low,mean,high]"))

plot(pred)

glimpse(df)
