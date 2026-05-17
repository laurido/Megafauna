library(tidyverse)
library(reticulate)

species = "Panthera_leo"
if (startsWith(getwd(), "/home/lakrids")) {
  path_prefix <- "/home/lakrids/GenomeDK"
} else {path_prefix <- "/faststorage/project"}

pickle <- import("pickle")
builtins <- import_builtins()
f <- builtins$open(paste0(path_prefix, "/megaFauna/sa_megafauna/results/", species, "/parameters_", species, ".pkl"), "rb")
params <- pickle$load(f)
f$close()
mu <- as.numeric(params$mu)
gen_time <- as.numeric(params$generation)

df_gone <- read.csv(Sys.glob(paste0(path_prefix, "/megaFauna/sa_megafauna/results/", species, "/GONE/GONE_pop0_18_Panthera_leo_GONE2_Ne")), sep = "\t") %>%
  rename(Ne = Ne_diploids) %>%
  mutate(time_bp_ka = (-Generation * gen_time) / 1000)

df_smc <- read.csv(Sys.glob(paste0(path_prefix, "/megaFauna/sa_megafauna/results/", species, "/smcpp/smcpp_plot_*_*_", species, ".csv"))) %>%
  rename(time_gen_a = x, Ne = y, population = label) %>%
  select(-c(plot_type, plot_num)) %>%
  mutate(time_bp_ka = (-time_gen_a * gen_time) / 1000)
end_time <- min(df_smc$time_bp_ka)

df_ne <- bind_rows(
  gone = df_gone %>% select(time_bp_ka, Ne),
  smc  = df_smc  %>% select(time_bp_ka, Ne),
  .id = "source")

df_ne <- df_ne %>% mutate(log_time = ifelse(time_bp_ka != 0, -log10(-time_bp_ka), NA),
                          log_ne = ifelse(Ne != 0, log10(Ne), NA))
df_ne %>% ggplot(aes(x = log_time, y = Ne, colour = source)) + geom_line()

df_npp <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/npp.csv")) %>%
  filter(Binomial.1.2 == species, time_bp_ka >= ceiling(end_time))

df_npp %>% ggplot(aes(x = time_bp_ka, y = mean)) + geom_point()


df_bioclim01 <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/bioclim01.csv")) %>%
  filter(Binomial.1.2 == species, time_bp_ka >= ceiling(end_time))

df_bioclim12 <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/bioclim12.csv")) %>%
  filter(Binomial.1.2 == species, time_bp_ka >= ceiling(end_time))

df_ne_agg <- df_ne %>%
  mutate(time_bp_ka = ceiling(time_bp_ka)) %>%  # assigns each Ne point to its 1ka window
  group_by(time_bp_ka) %>%
  summarise(Ne_mean = mean(Ne),
            Ne_median = median(Ne))

df_joined <- df_ne_agg %>%
  left_join(df_npp %>% select(-c(pn.range, range.size, Binomial.1.2)) %>%
            rename_with(~paste0(.x, "_npp"), -time_bp_ka), 
            by = c("time_bp_ka" = "time_bp_ka")) %>%
  left_join(df_bioclim01 %>% select(-c(pn.range, range.size, Binomial.1.2)) %>%
            rename_with(~paste0(.x, "_temp"), -time_bp_ka), 
            by = c("time_bp_ka" = "time_bp_ka")) %>%
  left_join(df_bioclim12 %>% select(-c(pn.range, range.size, Binomial.1.2)) %>%
            rename_with(~paste0(.x, "_pre"), -time_bp_ka), 
            by = c("time_bp_ka" = "time_bp_ka"))

df_joined <- df_ne_agg %>%
  left_join(df_npp, by = c("time_bp_ka" = "time_bp_ka"))

df_ne_agg %>% ggplot(aes(x = log10(-time_bp_ka + 1), y = Ne_mean)) + geom_point()

model <- lm(log(Ne_mean) ~ df_npp$mean, data = df_ne_agg)
summary(model)

df_joined <- df_joined %>%
  mutate(Ne_predicted = exp(predict(model, df_joined)))

df_joined %>%
  ggplot(aes(x = time_bp_ka)) +
  geom_line(aes(y = Ne_mean, colour = "observed")) +
  geom_line(aes(y = Ne_predicted, colour = "predicted")) +
  scale_x_reverse() +
  labs(x = "Time (ka BP)", y = "Ne", colour = "")

m_npp   <- lm((Ne_mean) ~ mean_npp,   data = df_joined)
m_temp <- lm((Ne_mean) ~ mean_temp, data = df_joined)
m_pre <- lm((Ne_mean) ~ mean_pre, data = df_joined)

df_joined <- df_joined %>%
  mutate(
    Ne_pred_npp   = (predict(m_npp,   df_joined)),
    Ne_pred_temp = (predict(m_temp, df_joined)),
    Ne_pred_pre = (predict(m_pre, df_joined))
  )

df_joined %>%
  pivot_longer(cols = c(Ne_mean, Ne_pred_npp, Ne_pred_temp, Ne_pred_pre),
               names_to = "series", values_to = "Ne") %>%
  ggplot(aes(x = log10(-time_bp_ka + 1), y = log10(Ne), colour = series)) +
  geom_line() +
  scale_x_reverse() +
  labs(x = "Time (ka BP)", y = "Ne", colour = "")

df_combine <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/trait_data_reported.csv")) %>%
  filter(phylacine_binomial == "Panthera leo")

biogeo <- strsplit(df_combine$biogeographical_realm[1], ", ")[[1]][1]
