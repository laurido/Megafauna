library(tidyverse)
library(reticulate)
library(dplyr)
library(purrr)

if (startsWith(getwd(), "/home/lakrids")) {
  path_prefix <- "/home/lakrids/GenomeDK"
} else {path_prefix <- "/faststorage/project"}

df_all <- data.frame()

genus_list <- c("Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros")

# Read and combine sample files
data <- map_dfr(genus_list, ~read.delim(paste0(path_prefix, "/megaFauna/sa_megafauna/metadata/samples_", .x, ".txt"),
                                        colClasses = "character"))

ref_folders <- sort(unique(data$REFERENCE_FOLDER))

# Read and filter references
references <- read.delim(paste0(path_prefix, "/megaFauna/sa_megafauna/metadata/references.txt")) %>%
  filter(REFERENCE_FOLDER %in% ref_folders)

# Make species_and_refs dataframe
species_and_refs <- data %>%
  mutate(GVCF_FOLDER = paste(GENUS, SPECIES, sep = "_")) %>%
  select(FOLDER, REFERENCE_FOLDER, GVCF_FOLDER) %>%
  distinct() %>%
  left_join(references, by = "REFERENCE_FOLDER")

# Loop
for (i in seq_len(nrow(species_and_refs))) {
  group      <- species_and_refs$FOLDER[i]
  ref_folder <- species_and_refs$REFERENCE_FOLDER[i]
  species_name <- gsub("_", " ", group)
  ref_name <-gsub("_", " ", ref_folder)
  print(group)
    # Import species parameters
  pickle <- import("pickle")
  builtins <- import_builtins()
  f <- builtins$open(paste0(path_prefix, "/megaFauna/sa_megafauna/results/", group, "/parameters_", group, ".pkl"), "rb")
  params <- pickle$load(f)
  f$close()
  mu <- as.numeric(params$mu)
  biogeo <- params$biogeography
  gen_time <- as.numeric(params$generation)
  
  
  # COMBINE life history traits
  df_traits <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/trait_data_reported.csv"))
  df_traits <- if (any(df_traits$phylacine_binomial == species_name)) {
    df_traits %>% filter(phylacine_binomial == species_name)} else {
      df_traits %>% filter(phylacine_binomial == ref_name)}
  
  if (any(df_traits$phylacine_binomial == species_name)) {
    gen_time <- as.numeric(df_traits$generation_length_d)/365}
  print(gen_time)
  
  # GONE data
  df_gone <- read.csv(Sys.glob(paste0(path_prefix, "/megaFauna/sa_megafauna/results/", group, "/GONE/*GONE2_Ne")), sep = "\t") %>%
    rename(Ne = Ne_diploids) %>%
    mutate(time_ya = Generation * gen_time)
  
  # SMC data GONE takes over when possible
  df_smc <- read.csv(Sys.glob(paste0(path_prefix, "/megaFauna/sa_megafauna/results/", group, "/smcpp/smcpp_plot_*_*_", group, ".csv"))) %>%
    rename(generations_ago = x, Ne = y, population = label) %>%
    select(-c(plot_type, plot_num)) %>%
    mutate(time_ya = generations_ago * gen_time, Ne = if_else(row_number() == n(), NA, Ne)) %>%
    filter(generations_ago > max(df_gone$Generation))
  
  end_time <- min(max(df_smc$time_ya)/1000, 800)
  print(end_time)
  
  # Join GONE and SMC
  df_ne <- bind_rows(
    gone = df_gone %>% select(time_ya, Ne),
    smc  = df_smc  %>% select(time_ya, Ne),
    .id = "source")
  df_ne <- df_ne %>% mutate(log_time = ifelse(time_ya != 0, log10(time_ya), NA),
                            log_ne = ifelse(Ne != 0, log10(Ne), NA),
                            time_kya = time_ya / 1000)
  
  # Plot GONE and SMC
  df_ne %>% ggplot(aes(x = log_time, y = Ne)) +
    geom_line() + scale_x_continuous(labels = function(x) parse(text = paste0("10^", x))) + xlab("Years ago")
  
  # Make minimum 1000 year windows for merging with other data
  df_ne_windows <- df_ne %>%
    arrange(time_kya) %>%
    mutate(window_end  = lead(time_kya), 
           window_end  = pmin(window_end, end_time),
           window_size = window_end - time_kya) %>%
    filter(!is.na(window_end), time_kya <= end_time) %>%
    # add a final row to close the last window
    bind_rows(tibble(time_kya   = max(.$window_end),
                     window_end  = max(.$window_end),
                     window_size = 0,
                     Ne          = last(.$Ne)))
  
  df_ne_agg <- df_ne_windows %>% rowwise() %>%
    do({row <- .
    bins <- seq(floor(row$time_kya), floor(row$window_end))
    tibble(
      time_kya    = bins,
      Ne     = row$Ne,
      weight = pmin(bins + 1, row$window_end) - pmax(bins, row$time_kya))
    }) %>% ungroup() %>% group_by(time_kya) %>% summarise(Ne_mean = weighted.mean(Ne, weight))
  
  # Collapse windows so only the first one is kept and the very last one in the dataset
  df_ne_collapsed <- df_ne_agg %>%
    filter(row_number() == 1 | Ne_mean != lag(Ne_mean) | row_number() == n()) %>%
    mutate(Ne_mean = if_else(row_number() == n(), NA, Ne_mean))
  
  df_ne_collapsed %>% ggplot(aes(x = log10(time_kya), y = Ne_mean)) +
    geom_line() + scale_x_continuous(labels = function(x) parse(text = paste0("10^", x)))
  
  # Helper function to filter with fallback
  filter_bioclim <- function(df, group, ref, time_filter) {
    if (any(df$Binomial.1.2 == group)) {
      df %>% filter(Binomial.1.2 == group)
    } else {
      df %>% filter(Binomial.1.2 == ref_folder)
    }}
  
  # Net primary production
  df_npp <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/npp.csv")) %>%
    filter_bioclim(group, ref_folder)
  df_npp %>% ggplot(aes(x = time_bp_ka, y = mean)) + geom_line()
  
  # Temperature
  df_temp <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/bioclim01.csv")) %>%
    filter_bioclim(group, ref_folder, ceiling(end_time))
  
  # Precipitation
  df_precipitation <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/bioclim12.csv")) %>%
    filter_bioclim(group, ref_folder, ceiling(end_time))
  
  join_climate <- function(df_ne, df_climate, time_col, suffix) {
    df_ne %>% rowwise() %>%
      mutate(
        !!paste0("mean_", suffix) := mean(df_climate$mean[df_climate[[time_col]] >= -window_end & df_climate[[time_col]] < -time_kya]),
        !!paste0("var_",  suffix) := mean(df_climate$var[ df_climate[[time_col]] >= -window_end & df_climate[[time_col]] < -time_kya])    ) %>%
      ungroup()}
  
  df_joined <- df_ne_collapsed %>%
    arrange(time_kya) %>%
    mutate(window_end = lead(time_kya, default = max(time_kya))) %>%
    join_climate(df_npp,           "time_bp_ka", "npp") %>%
    join_climate(df_temp,          "time_bp_ka", "temp") %>%
    join_climate(df_precipitation, "time_bp_ka", "prec") %>%
    rowwise() %>%
    mutate(
      pn.range   = mean(df_npp$pn.range[  df_npp$time_bp_ka >= -window_end & df_npp$time_bp_ka < -time_kya]),
      range.size = mean(df_npp$range.size[df_npp$time_bp_ka >= -window_end & df_npp$time_bp_ka < -time_kya]),
      Binomial.1.2       = group) %>%
    ungroup() %>% rename(species = Binomial.1.2)#%>% select(-window_end)
  
  df_human <- read.table(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/biogeo.txt"), header = TRUE)
  
  arrival <- df_human %>%
    filter(Biogeographic_realm == biogeo)
  
  arrival_mean <- mean(c(arrival$al, arrival$au))
  
  df_joined <- df_joined %>%
    mutate(human_pressure = case_when(
      time_kya >= arrival_mean  ~ 0,                                        # window entirely before humans
      window_end <= arrival_mean ~ 1,                                        # window entirely after humans
      TRUE ~ (arrival_mean - time_kya) / (window_end - time_kya))) # partial overlap
  
    # Add biogeography column
    df_joined <- df_joined %>% mutate(biogeography = biogeo, generation_time = gen_time, adult_mass = df_traits$adult_mass_g)
    
    df_joined <- df_joined %>% rename(Ne = Ne_mean, p_human = human_pressure)
  
  df_all <- bind_rows(df_all, df_joined)
}
