########################################################################################################################################################################
###################################################------------------------------------------------------------------###################################################
########################################################################################################################################################################
library(dplyr)
genus = "Elephas"
species = "maximus"
if (startsWith(getwd(), "/home/lakrids")) {
  path_prefix <- "/home/lakrids/GenomeDK"
  } else {path_prefix <- "/faststorage/project"}

#### read in sra file
sra_result <- read.csv(paste0(path_prefix, "/megaFauna/people/laurids/new_samples/", genus, "_all.csv"))
sra_result <- sra_result %>% filter(
  any(startsWith(scientific_name, paste0(genus, " ", species))), library_strategy == "WGS" & library_source == "GENOMIC" & library_selection == "RANDOM" & library_layout == "PAIRED")
study_titles_and_bioproject_ids <- unique(data.frame(Study.Title=sra_result$study_title, bioproject_id=sra_result$bioproject_id))
dim(study_titles_and_bioproject_ids)[1]
for(i in 1:dim(study_titles_and_bioproject_ids)[1]){
  print(paste(study_titles_and_bioproject_ids$Study.Title[i], "---", study_titles_and_bioproject_ids$bioproject_id[i], "--- Exclusion_criteria: NA --- Inclusion_criteria: NA"))
}

selected_BioProject_accession <- c(rep("PRJEB98188", each=23),
                                   rep("PRJDB9017", each=4),
                                   rep("PRJNA1013751", each=28),
                                   rep("PRJNA247496", each=4),
                                   rep("PRJNA840935", each=7),
                                   rep("PRJNA622303", each=5),
                                   rep("PRJNA512907", each=1),
                                   rep("PRJEB4905", each=2))

organism_name <-                 c(rep("Elephas_maximus", each=23),
                                   rep("Elephas_maximus_indicus", each=4),
                                   rep("Elephas_maximus", each=28), 
                                   rep("Elephas_maximus", each=4),
                                   rep("Elephas_maximus", each=7),
                                   rep("Elephas_maximus", each=5),
                                   rep("Elephas_maximus", each=1),
                                   rep("Elephas_maximus", each=2))

#### manual inspection
prjna = "PRJDB9017"
dim(subset(sra_result, bioproject_id == prjna))
table(subset(sra_result, bioproject_id == prjna)$scientific_name)

#### make dataframe of selected projects
selected_data <- data.frame(Organism.Name = organism_name, BioProject.Accession = selected_BioProject_accession)
print(sort(unique(selected_data$BioProject.Accession)))


#### subset the sra_result file
sra_result_subset <- subset(sra_result, bioproject_id %in% selected_data$BioProject.Accession)
table(sra_result_subset$scientific_name)

#### size in tb
sum(sra_result_subset$gb_size, na.rm = T)/1000
length(unique(sra_result_subset$sample_accession))

#### make dataframe needed to extract ftp links and extract
sample_accessions_biop_ids <- unique(data.frame(sample_accession = sra_result_subset$sample_accession, bioproject_id = sra_result_subset$bioproject_id, scientific_name = sra_result_subset$scientific_name,
                                                sex = sra_result_subset$sex, biosample_id = sra_result_subset$biosample_id))
ENA.SRA_SAMPLE_ID <- c(); run_accession <- c(); R1_or_R2 <- c(); fastq_ftp <- c(); fastq_bytes <- c(); fastq_md5 <- c()
IND_IDs <- c(); DATASETs <- c(); SEXes<-c(); FOLDERs <- c()
GENUSes <- c(); SPECIESes <- c(); SUBSPECIESes <- c()
SAMPLE_ACCESSIONs <- c()
error_accessions <- c()
for (i in 1:dim(sample_accessions_biop_ids)[1]){
  err <- read.csv(paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",sample_accessions_biop_ids$sample_accession[i],"&result=read_run"), sep = "\t")
  for(j in 1:dim(err)[1]){
    if (is.na(err$fastq_ftp[j]) || err$fastq_ftp[j] == "") next
    fastq_ftps <- unlist(strsplit((err$fastq_ftp[j]), split = ";"))
    for(k in fastq_ftps){
      if(length(unlist(strsplit(k, "_"))) == 2 & unlist(strsplit(k, "_"))[2] == "1.fastq.gz"){
        # find if there is read 2
        if(paste0(unlist(strsplit(k, "_"))[1], "_2.fastq.gz") %in% fastq_ftps){
          index_read_1 = match(k, fastq_ftps)
          index_read_2 = match(paste0(unlist(strsplit(k, "_"))[1], "_2.fastq.gz"), fastq_ftps)
          
          fastq_ftp <- c(fastq_ftp, c(k, paste0(unlist(strsplit(k, "_"))[1], "_2.fastq.gz")))
          ENA.SRA_SAMPLE_ID <- c(ENA.SRA_SAMPLE_ID, rep(sample_accessions_biop_ids$sample_accession[i], 2))
          name_in_list <- unlist(strsplit(sample_accessions_biop_ids$scientific_name[i], " "))
          IND_IDs <- c(IND_IDs, rep(sample_accessions_biop_ids$biosample_id[i], 2))
          DATASETs <- c(DATASETs, rep(sample_accessions_biop_ids$bioproject_id[i], 2))
          SEXes <- c(SEXes, rep(sample_accessions_biop_ids$sex[i], 2))
          if(length(name_in_list)==2){
            FOLDERs <- c(FOLDERs, rep(paste(c(name_in_list, "ssp"), collapse = "_"), 2))
            GENUSes <- c(GENUSes, rep(name_in_list[1],2)); SPECIESes <- c(SPECIESes, rep(name_in_list[2],2)); SUBSPECIESes <- c(SUBSPECIESes, rep("ssp",2))
          }else if(length(name_in_list)==3){
            FOLDERs <- c(FOLDERs, rep(paste(name_in_list, collapse = "_"), 2))
            GENUSes <- c(GENUSes, rep(name_in_list[1],2)); SPECIESes <- c(SPECIESes, rep(name_in_list[2],2)); SUBSPECIESes <- c(SUBSPECIESes, rep(name_in_list[3],2))
          }else if(length(name_in_list)==1){
            FOLDERs <- c(FOLDERs, rep(paste(c(name_in_list, "sp", "ssp"), collapse = "_"), 2))
            GENUSes <- c(GENUSes, rep(name_in_list[1],2)); SPECIESes <- c(SPECIESes, rep("sp",2)); SUBSPECIESes <- c(SUBSPECIESes, rep("ssp",2))
          }else if(sample_accessions_biop_ids$scientific_name[i] == "Papio cynocephalus x Papio anubis"){
            FOLDERs <- c(FOLDERs, rep("Papio_cynocephalus_anubishybrid", 2))
            GENUSes <- c(GENUSes, rep(name_in_list[1],2)); SPECIESes <- c(SPECIESes, rep(name_in_list[2],2)); SUBSPECIESes <- c(SUBSPECIESes, rep("anubishybrid",2))
          }else{
            FOLDERs <- c(FOLDERs, rep(NA, 2))
            GENUSes <- c(GENUSes, rep(NA,2)); SPECIESes <- c(SPECIESes, rep(NA,2)); SUBSPECIESes <- c(SUBSPECIESes, rep(NA,2))
            error_accessions <- c(error_accessions, sample_accessions_biop_ids$sample_accession[i])
            print(paste0("Error: ", sample_accessions_biop_ids$sample_accession[i]))
          }
          run_accession <- c(run_accession, rep(err$run_accession[j], each = 2))
          SAMPLE_ACCESSIONs <- c(SAMPLE_ACCESSIONs, rep(sample_accessions_biop_ids$sample_accession[i], 2))
          R1_or_R2 <- c(R1_or_R2, c("R1", "R2"))
          fastq_bytes <- c(fastq_bytes, unlist(strsplit((err$fastq_bytes[j]), split = ";"))[c(index_read_1, index_read_2)])
          fastq_md5 <- c(fastq_md5, unlist(strsplit((err$fastq_md5[j]), split = ";"))[c(index_read_1, index_read_2)])
        }
      }
    }
  }
}
print(length(error_accessions))

## take only data that is needed to run the workflow
data_stripped_pre = data.frame(IND_ID = IND_IDs, DATASET = DATASETs,
                               FOLDER = FOLDERs,
                               #GENUS = data_merged$GENUS, SPECIES = data_merged$SPECIES, SUBSPECIES = data_merged$SUBSPECIES,
                               GENUS = GENUSes, SPECIES = SPECIESes, SUBSPECIES = SUBSPECIESes,
                               SEX = SEXes,
                               REFERENCE_FOLDER = NA,
                               # run_accession = data_merged$run_accession, R1_or_R2 = data_merged$R1_or_R2,
                               # fastq_ftp = data_merged$fastq_ftp, fastq_bytes = data_merged$fastq_bytes, fastq_md5 = data_merged$fastq_md5)
                               run_accession = run_accession, sample_accession = SAMPLE_ACCESSIONs, R1_or_R2 = R1_or_R2,
                               fastq_ftp = fastq_ftp, fastq_bytes = fastq_bytes, fastq_md5 = fastq_md5)

# Control that there is a maximum of around 40x coverage per sample
accession_list <- character()
for (sample in unique(data_stripped_pre$IND_ID)){
  accessions <- unique(filter(arrange(data_stripped_pre, desc(fastq_bytes)), IND_ID == sample)$run_accession)
  gbases <- 0
  for (accession in accessions){
    gbases <- gbases + sum(filter(sra_result, run_accession == accession)$total_bases)
    if (gbases <= 0) next
    accession_list <- append(accession_list, accession)
    if (gbases >= 120 * 1000000000) break }}

data_stripped_pre <- filter(data_stripped_pre, run_accession %in% accession_list)
sum(as.numeric(data_stripped_pre$fastq_bytes))/1000000000

#### see if there are run_accessions that map to multiple biosamples
multi_run_accessions <- c()
for(ra in unique(run_accession)){
  ss <- subset(data_stripped_pre, run_accession == ra)
  ss <- unique(data.frame(IND_ID = ss$IND_ID, run_accession = ss$run_accession))
  if(dim(ss)[1] > 1){
    multi_run_accessions <- c(multi_run_accessions, ra)
  }
}
print(multi_run_accessions)

#### include only Panthera tigris tigris and change folder name
data_stripped_pre <- subset(data_stripped_pre, grepl("^Elephas_maximus", FOLDER))
data_stripped_pre$FOLDER <- "Elephas_maximus"

#### exclude multi_run_accessions
length(unique(data_stripped_pre$IND_ID))
hist(as.numeric(data_stripped_pre$fastq_bytes)/1000000000, xlab = "Size in GB")
summary(as.numeric(data_stripped_pre$fastq_bytes)/1000000000)
sum(as.numeric(data_stripped_pre$fastq_bytes)/1000000000)

## fix sexes
data_stripped <- data_stripped_pre
unique(data_stripped$SEX)
data_stripped$SEX[data_stripped$SEX == "male"] <- "M"
data_stripped$SEX[data_stripped$SEX == "female"] <- "F"
data_stripped$SEX[data_stripped$SEX == "not collected"] <- "Unknown"
data_stripped$SEX[data_stripped$SEX == "missing"] <- "Unknown"
data_stripped$SEX[data_stripped$SEX == "not recorded"] <- "Unknown"
data_stripped$SEX[data_stripped$SEX == "not applicable"] <- "Unknown"
data_stripped$SEX[data_stripped$SEX == "Not determined"] <- "Unknown"
data_stripped$SEX[data_stripped$SEX == "unspecified"] <- "Unknown"
unique(data_stripped$SEX)

## set reference folders
unique(data_stripped$SPECIES)
data_stripped$REFERENCE_FOLDER[data_stripped$GENUS == "Elephas"] <- "Elephas_maximus"
unique(data_stripped$FOLDER)
unique(data_stripped$REFERENCE_FOLDER)

### write subset file
data_stripped <- data_stripped[order(data_stripped$FOLDER, data_stripped$IND_ID, data_stripped$run_accession, data_stripped$R1_or_R2),] %>%
  filter(run_accession != "SRR12877368")
write.table(data_stripped, paste0(path_prefix, "/megaFauna/people/laurids/new_samples/samples_new_",genus,".txt"), quote = FALSE, sep = "\t", row.names = F)

df1 <- read.table(paste0(path_prefix, "/megaFauna/people/laurids/new_samples/samples_new_",genus,".txt"), sep="\t", header=TRUE)
df2 <- read.table(paste0(path_prefix, "/megaFauna/sa_megafauna/metadata/samples_",genus,"_all.txt"), sep="\t", header=TRUE)

setdiff(df1$run_accession, df2$run_accession)
overlap <- intersect(df1$sample_accession, df2$sample_accession)

library(dplyr)
combined <- dplyr::bind_rows(df1, df2)
combined_unique <- combined %>%
  distinct(sample_accession, run_accession, R1_or_R2, .keep_all = TRUE)

bad_runs <- combined_unique %>%
  count(run_accession) %>%
  filter(n != 2)

if (nrow(bad_runs) > 0) {
  message("Warning: some runs do not have both R1 and R2:")
  print(bad_runs)
}


write.table(combined_unique, paste0(path_prefix, "/megaFauna/people/laurids/new_samples/samples_",genus,".txt"), quote = FALSE, sep = "\t", row.names = F)

