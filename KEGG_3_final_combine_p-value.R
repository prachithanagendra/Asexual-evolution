#!/usr/bin/env Rscript

# ---- Load libraries ----
library(dplyr)
library(purrr)
library(metap)   # for Fisher’s method
library(stringr)

# ---- User Input ----
input_dir <- "/CNV_GO/KEGG/Glu_Del-Gal_Dup/g_del_G_dup_combined"
output_file <- file.path(input_dir, "g_del_G_dup_combined_final_combined_replicates.csv")

# ---- 1. Read all replicate CSVs ----
csv_files <- list.files(input_dir, pattern = "^rep[0-9]+_combined\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No replicate combined CSV files found in: ", input_dir)
}

all_data <- lapply(csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  df$source_file <- basename(f)   # track which replicate
  return(df)
}) %>% bind_rows()

cat("Loaded", nrow(all_data), "rows from", length(csv_files), "replicate CSVs.\n")

# ---- 2. Fisher’s method helper ----
combine_fisher <- function(pvals) {
  pvals <- as.numeric(pvals)
  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0) {
    return(NA)
  } else if (length(pvals) == 1) {
    return(pvals)   # keep the single value
  } else {
    return(sumlog(pvals)$p)
  }
}

# ---- 3. Aggregate across replicates ----
final_combined <- all_data %>%
  group_by(ID, Description) %>%
  summarise(
    bins_present     = sum(bins_present, na.rm = TRUE),
    combined_p       = combine_fisher(combined_p),
    combined_padjust = combine_fisher(combined_padjust),
    combined_q       = combine_fisher(combined_q),
    geneID           = paste(unique(unlist(strsplit(paste(geneID, collapse = "; "), "; "))), collapse = "; "),
    total_count      = sum(as.numeric(total_count), na.rm = TRUE),
    GeneRatio        = paste(unique(na.omit(GeneRatio)), collapse = "; "),
    BgRatio          = paste(unique(na.omit(BgRatio)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(combined_p)

# ---- 4. Save results ----
write.csv(final_combined, output_file, row.names = FALSE)
cat("Final combined results saved to:", output_file, "\n")
