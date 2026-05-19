library(dplyr)
library(purrr)

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

# Loop
for (i in seq_len(nrow(species_and_refs))) {
  group      <- species_and_refs$FOLDER[i]
  ref_folder <- species_and_refs$REFERENCE_FOLDER[i]
  species <- gsub("_", " ", group)
  print(species)
  # rest of your code here
}



m_npp   <- lm(Ne_mean ~ mean_npp,   data = df_joined)
m_temp <- lm(Ne_mean ~ mean_temp, data = df_joined)
m_prec <- lm(Ne_mean ~ mean_prec, data = df_joined)
m_range <- lm(Ne_mean ~ (range.size/pn.range), data = df_joined)
summary(m_range)

df_joined <- df_joined %>%
  mutate(
    Ne_pred_npp   = predict(m_npp,   df_joined),
    Ne_pred_temp = predict(m_temp, df_joined),
    Ne_pred_prec = predict(m_prec, df_joined),
    Ne_pred_range = predict(m_range, df_joined)
  )

df_joined %>%
  pivot_longer(cols = c(Ne_mean, Ne_pred_npp, Ne_pred_temp, Ne_pred_prec, Ne_pred_range),
               names_to = "series", values_to = "Ne") %>%
  ggplot(aes(x = log10(time_kya), y = Ne, colour = series)) +
  geom_line() +
  labs(x = "Time (ka BP)", y = "Ne", colour = "")

