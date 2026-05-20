library(tidyverse)
library(reticulate)
library(dplyr)
library(purrr)
if (startsWith(getwd(), "/home/lakrids")) {
  path_prefix <- "/home/lakrids/GenomeDK"
} else {path_prefix <- "/faststorage/project"}

df_all <- data.frame()
genus_list <- c("Loxodonta", "Elephas", "Boselaphus", "Panthera", "Rhinoceros", "Ceratotherium", "Diceros")

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
print(species_and_refs)
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

  # Join GONE and SMC
  df_ne <- bind_rows(
    gone = df_gone %>% select(time_ya, Ne),
    smc  = df_smc  %>% select(time_ya, Ne),
    .id = "source") %>% mutate(log_time = ifelse(time_ya != 0, log10(time_ya), NA),
                            log_ne = ifelse(Ne != 0, log10(Ne), NA),
                            time_kya = time_ya / 1000)
  df_ne <- df_ne %>%
    arrange(time_kya) %>%
    filter(row_number() == 1 | Ne != lag(Ne) | row_number() == n())
  
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
    bins <- seq(floor(row$time_kya), ceiling(row$window_end))
    tibble(
      time_kya = bins,
      Ne       = row$Ne,
      weight   = pmax(0, pmin(bins + 1, row$window_end) - pmax(bins, row$time_kya))
    )
    }) %>% ungroup() %>% group_by(time_kya) %>% summarise(Ne_mean = weighted.mean(Ne, weight))
  
  # Collapse windows so only the first one is kept and the very last one in the dataset
    df_ne_collapsed <- df_ne_agg %>%
    filter(row_number() == 1 | !is.na(lag(Ne_mean)) & Ne_mean != lag(Ne_mean) | row_number() == n()) %>%
    mutate(Ne_mean = if_else(row_number() == n(), NA, Ne_mean))
  
  df_ne_collapsed %>% ggplot(aes(x = log10(time_kya), y = Ne_mean)) +
    geom_line() + scale_x_continuous(labels = function(x) parse(text = paste0("10^", x)))
  
  # Helper function to filter with fallback
  filter_bioclim <- function(df, group, ref, time_filter) {
    if (any(df$Binomial.1.2 == group)) {
      df %>% filter(Binomial.1.2 == group)
    } else {
      df %>% filter(Binomial.1.2 == ref)
    }}
  
  # Net primary production
  df_npp <- read.csv(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/npp.csv")) %>%
    filter_bioclim(group, ref_folder, ceiling(end_time))
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
    df_joined <- df_joined %>% mutate(biogeography = biogeo, generation_time = gen_time, adult_mass = df_traits$adult_mass_g, mutation_rate = mu)
    
    df_joined <- df_joined %>% rename(Ne = Ne_mean, p_human = human_pressure)
  
  df_all <- bind_rows(df_all, df_joined)
}

write.csv(df_all, paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/modeling_data.csv"), row.names = FALSE)

summary_table <- data %>%
  mutate(GVCF_FOLDER = paste(GENUS, SPECIES, sep = "_")) %>%
  group_by(GVCF_FOLDER) %>%
  summarise(total_count = n_distinct(IND_ID)) %>%  # replace SAMPLE_ID with your actual column
  left_join(species_and_refs %>% select(GVCF_FOLDER, REFERENCE_FOLDER, ref_genome_name), 
            by = "GVCF_FOLDER") %>%
  rename(species        = GVCF_FOLDER,
         reference_used = REFERENCE_FOLDER,
         reference_name = ref_genome_name)

library(gt)
Sys.setenv(CHROMOTE_CHROME = "/usr/bin/brave-browser")
summary_table %>%
  mutate(
    species        = gsub("_", " ", species),
    reference_used = gsub("_", " ", reference_used)
  ) %>%
  gt() %>%
  tab_header(title = "Species Summary") %>%
  cols_label(
    species        = "Species",
    total_count    = "Total Count",
    reference_used = "Reference Used",
    reference_name = "Reference Name"
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = c(species, reference_used))
  ) %>%
  gtsave(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/summary_table.png"))

library(gt)

summary_table2 <- data %>%
  mutate(GVCF_FOLDER = paste(GENUS, SPECIES, sep = "_")) %>%
  group_by(GVCF_FOLDER, DATASET) %>%
  summarise(project_count = n_distinct(IND_ID), .groups = "drop") %>%
  group_by(GVCF_FOLDER) %>%
  mutate(total_count = sum(project_count)) %>%
  ungroup() %>%
  rename(species = GVCF_FOLDER)

library(gt)

summary_table2 <- data %>%
  mutate(GVCF_FOLDER = paste(GENUS, SPECIES, sep = "_")) %>%
  group_by(GVCF_FOLDER, DATASET) %>%
  summarise(project_count = n_distinct(SAMPLE_ID), .groups = "drop") %>%
  group_by(GVCF_FOLDER) %>%
  mutate(total_count = sum(project_count)) %>%
  ungroup() %>%
  rename(species = GVCF_FOLDER)

summary_table2 %>%
  mutate(species = gsub("_", " ", species)) %>%
  gt(groupname_col = "species") %>%  # groups rows by species
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_row_groups()
  ) %>%
  cols_label(
    total_count   = "Total Count",
    DATASET       = "Project",
    project_count = "Project Count"
  ) %>%
  gtsave(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/summary_table2.png"))

library(flextable)
library(officer)
summary_table2 %>%
  arrange(species) %>%
  mutate(species = gsub("_", " ", species)) %>%
  flextable() %>%
  merge_v(j = c("species", "total_count")) %>%
  valign(j = c("species", "total_count"), valign = "center") %>%
  italic(j = "species") %>%
  hline(i = cumsum(table(summary_table2$species))[-length(table(summary_table2$species))],
        border = fp_border(color = "gray40", width = 2)) %>%
  save_as_html(path = paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/summary_table2.html"))

# Split species into two halves
species_list <- unique(summary_table2$species)
half <- ceiling(length(species_list) / 2)

make_ft <- function(species_subset) {
  df <- summary_table2 %>%
    arrange(species) %>%
    mutate(species = gsub("_", " ", species)) %>%
    filter(species %in% species_subset)
  
  flextable(df) %>%
    merge_v(j = c("species", "total_count")) %>%
    valign(j = c("species", "total_count"), valign = "center") %>%
    italic(j = "species") %>%
    hline(i = cumsum(table(df$species))[-length(table(df$species))],
          border = fp_border(color = "gray40", width = 2))
}

ft1 <- make_ft(gsub("_", " ", species_list[1:half+1]))
ft2 <- make_ft(gsub("_", " ", species_list[(half+2):length(species_list)]))

# Save as HTML with two columns
html <- htmltools::tagList(
  htmltools::tags$div(
    style = "display: flex; gap: 40px;",
    htmltools::tags$div(htmltools::HTML(knit_print(ft1))),
    htmltools::tags$div(htmltools::HTML(knit_print(ft2)))
  )
)

htmltools::save_html(html, paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/summary_table2.html"))

species_info <- df_all %>%
  group_by(species) %>%
  slice(1) %>%
  ungroup() %>%
  select(species, generation_time, mutation_rate)  # replace with actual column names

summary_table <- data %>%
  mutate(GVCF_FOLDER = paste(GENUS, SPECIES, sep = "_")) %>%
  group_by(GVCF_FOLDER) %>%
  summarise(total_count = n_distinct(IND_ID)) %>%
  left_join(species_and_refs %>% select(GVCF_FOLDER, REFERENCE_FOLDER, ref_genome_name),
            by = "GVCF_FOLDER") %>%
  rename(species        = GVCF_FOLDER,
         reference_used = REFERENCE_FOLDER,
         reference_name = ref_genome_name) %>%
  left_join(species_info, by = "species") %>%
  select(-total_count) %>%
  mutate(
    species        = gsub("_", " ", species),
    reference_used = gsub("_", " ", reference_used)
  ) %>%
  gt() %>%
  fmt_number(columns = generation_time, decimals = 1) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = c(species, reference_used))
  ) %>%
  cols_label(
    species         = "Species",
    reference_used  = "Reference Used",
    reference_name  = "Reference Name",
    generation_time = "Generation Time",
    mutation_rate   = "Mutation Rate"
  ) %>%
  gtsave(paste0(path_prefix, "/megaFauna/sa_megafauna/data/bioclim/summary_table.png"))
