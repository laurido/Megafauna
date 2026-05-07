library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(GenomicRanges)

### 1. Define genome size first
genome_size <- 3401247148  # mEleMax1 primary haplotype

### 2. Load RepeatMasker BED
df <- read_delim(
  "GenomeDK/megaFauna/sa_megafauna/data/Elephas_maximus/masking/Repeat_mask_asianelephant",
  delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE
)
colnames(df) <- c("chrom", "start", "end", "name", "score", "strand",
                  "m_start", "m_end", "rgb")

### 3. Classify repeats
df <- df %>%
  separate(name, into = c("family", "class"), sep = "#", remove = FALSE) %>%
  mutate(
    major = case_when(
      grepl("LINE",           class, ignore.case = TRUE) ~ "LINE",
      grepl("SINE",           class, ignore.case = TRUE) ~ "SINE",
      grepl("LTR",            class, ignore.case = TRUE) ~ "LTR",
      grepl("DNA",            class, ignore.case = TRUE) ~ "DNA",
      grepl("Simple",         class, ignore.case = TRUE) ~ "Simple_repeat",
      grepl("Low_complexity", class, ignore.case = TRUE) ~ "Low_complexity",
      grepl("Satellite",      class, ignore.case = TRUE) ~ "Satellite",
      TRUE ~ "Other"
    )
  )

### 4. Build GRanges
gr <- GRanges(
  seqnames = df$chrom,
  ranges   = IRanges(start = df$start + 1, end = df$end),
  type     = df$major
)

### 5. Sanity check: total unique repeat space
gr_all <- reduce(gr)
total_repeat_bp <- sum(width(gr_all))
cat(sprintf("Total unique repeat bp: %.0f (%.2f%% of genome)\n",
            total_repeat_bp, 100 * total_repeat_bp / genome_size))
### 6. Split by type
gr_list <- split(gr, gr$type)

### 7. Priority-based non-overlapping bp per type
# Each base is assigned to the first type that claims it (highest priority wins).
# Priority: LINE > SINE > LTR > DNA > Satellite > Simple_repeat > Low_complexity > Other
priority <- c("LINE", "SINE", "LTR", "DNA", "Satellite",
              "Simple_repeat", "Low_complexity", "Other")

claimed  <- GRanges()   # tracks already-assigned bases
type_bp  <- setNames(numeric(length(priority)), priority)

for (type in priority) {
  if (!type %in% names(gr_list)) next
  
  # Reduce within this type, then remove already-claimed bases
  this_reduced  <- reduce(gr_list[[type]])
  this_unique   <- setdiff(this_reduced, claimed)
  type_bp[type] <- sum(width(this_unique))
  
  # Mark these bases as claimed
  claimed <- reduce(c(claimed, this_unique))
}

### 8. Compute percentages and sort
percent_summary <- round(100 * type_bp / genome_size, 2)
percent_summary <- sort(percent_summary, decreasing = TRUE)

cat(sprintf("\nTotal accounted for: %.2f%% of genome\n\n",
            sum(percent_summary)))

### 9. Pretty print
for (i in seq_along(percent_summary)) {
  cat(sprintf("%-15s ~%.2f%%\n",
              names(percent_summary)[i],
              percent_summary[i]))
}

