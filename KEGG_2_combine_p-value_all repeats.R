#!/usr/bin/env Rscript


# ---- Load libraries ----
library(dplyr)
library(purrr)
library(metap)   # for Fisher’s method
library(stringr)

# ---- User Input ----
# Parent folder containing rep1, rep2, ..., repN
input_parent <- "/CNV_GO/KEGG/Glu_Del-Gal_Dup"
output_parent <- file.path(input_parent, "g_del_G_dup_combined")

# Create output folder if missing
if (!dir.exists(output_parent)) dir.create(output_parent, recursive = TRUE)

# ---- 1. Fisher’s method ----
combine_fisher <- function(pvals) {
    pvals <- as.numeric(pvals)
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) == 0) {
        return(NA)
    } else if (length(pvals) == 1) {
        return(pvals)   # keep the single value if only one
    } else {
        return(sumlog(pvals)$p)
    }
}

# ---- 2. Column name mapping ----
colmap <- c(
  "ID"         = "ID",
  "Description"= "Description",
  "GeneRatio"  = "GeneRatio",
  "BgRatio"    = "BgRatio",
  "pvalue"     = "pvalue",
  "p.adjust"   = "p.adjust",
  "qvalue"     = "qvalue",
  "geneID"     = "geneID",
  "Count"      = "Count"
)

# ---- 3. Loop over all reps ----
rep_dirs <- list.dirs(input_parent, recursive = FALSE, full.names = TRUE)
rep_dirs <- rep_dirs[grepl("rep[0-9]+$", rep_dirs)]  # keep only rep1, rep2, etc.

cat("Found", length(rep_dirs), "replicate folders.\n")

for (rep_dir in rep_dirs) {
    rep_name <- basename(rep_dir)
    output_file <- file.path(output_parent, paste0(rep_name, "_combined.csv"))

    # --- Read all CSVs inside this rep ---
    csv_files <- list.files(rep_dir, pattern = "KEGG_bin_.*\\.csv$", full.names = TRUE)
    if (length(csv_files) == 0) {
        cat("Skipping", rep_name, "-> No CSV files found.\n")
        next
    }

    all_results <- lapply(csv_files, function(f) {
        df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
        colnames(df) <- trimws(colnames(df))

        # Create aligned output with only mapped columns
        out <- data.frame(matrix(ncol = length(colmap), nrow = nrow(df)))
        colnames(out) <- names(colmap)

        for (nm in names(colmap)) {
            src <- colmap[[nm]]
            if (src %in% colnames(df)) {
                out[[nm]] <- df[[src]]
            } else {
                out[[nm]] <- NA
            }
        }
        return(out)
    }) %>% bind_rows(.id = "bin_id")

    cat(rep_name, "-> Normalized input:", nrow(all_results),
        "rows across", length(csv_files), "files.\n")

    # --- Aggregate across bins ---
    combined <- all_results %>%
      group_by(ID, Description) %>%
      summarise(
        bins_present     = n(),
        combined_p       = combine_fisher(pvalue),
        combined_padjust = combine_fisher(p.adjust),
        combined_q       = combine_fisher(qvalue),
        geneID           = paste(unique(na.omit(geneID)), collapse = "; "),
        total_count      = sum(as.numeric(Count), na.rm = TRUE),
        GeneRatio        = paste(unique(na.omit(GeneRatio)), collapse = "; "),
        BgRatio          = paste(unique(na.omit(BgRatio)), collapse = "; "),
        .groups = "drop"
      ) %>%
      arrange(combined_p)

    # --- Save ---
    write.csv(combined, output_file, row.names = FALSE)
    cat("Saved:", output_file, "\n")
}
