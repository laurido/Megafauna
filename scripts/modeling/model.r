library(dplyr)
library(purrr)

genus_list <- c("Panthera")

# Read and combine sample files
data <- map_dfr(genus_list, ~read.delim(paste0(path_prefix, "/megaFauna/sa_megafauna/metadata/samples_", .x, ".txt")))

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
