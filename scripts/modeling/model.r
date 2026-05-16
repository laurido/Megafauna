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
generation <- as.numeric(params$generation)

df_smc <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/results/shared/smcpp/", "smcpp_plot_pop0_18_5.42_Panthera_leo.csv")) %>%
  rename(time_bp_ka = x, Ne = y, population = label) %>%
  select(-c(plot_type, plot_num)) %>%
  mutate(time_bp_ka = (-time_bp_ka) / 1000)
end_time <- min(df_smc$time_bp_ka)

df_npp <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/npp.csv")) %>%
  filter(Binomial.1.2 == species, time_bp_ka >= floor(end_time))

df_bioclim01 <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/bioclim01.csv")) %>%
  filter(Binomial.1.2 == species, time_bp_ka >= floor(end_time))

df_bioclim12 <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/bioclim12.csv")) %>%
  filter(Binomial.1.2 == species, time_bp_ka >= floor(end_time))

