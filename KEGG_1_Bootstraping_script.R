#!/usr/bin/env Rscript

# ============================================
# KEGG Enrichment with Bootstrapped Binning
# ============================================

library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)

# ---- USER INPUT ----
input_file  <- "/CNV_GO/KEGG/Glu_Del-Gal_Dup/heatmap_genes.txt"
output_dir  <- "/CNV_GO/KEGG/Glu_Del-Gal_Dup"
bin_size    <- 100   # genes per bin
n_repeats   <- 3     # number of bootstrap replicates

# ---- SETUP ----
dir.create(output_dir, showWarnings = FALSE)

# ---- 1. Load gene list ----
gene_list <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
gene_list <- unique(gene_list[gene_list != ""])  # remove blanks & duplicates
total_genes <- length(gene_list)
cat("Total genes:", total_genes, "\n")

# ---- 2. Bootstrap repeats ----
for (rep in 1:n_repeats) {
    cat("\n--- Bootstrap replicate", rep, "---\n")
    
    rep_dir <- file.path(output_dir, paste0("rep", rep))
    dir.create(rep_dir, showWarnings = FALSE)

    # Shuffle genes
    set.seed(100 + rep)  # different seed each time
    shuffled <- sample(gene_list)

    # Split into bins of bin_size
    num_bins <- ceiling(total_genes / bin_size)
    gene_bins <- split(shuffled, ceiling(seq_along(shuffled) / bin_size))
    cat("Bin sizes for replicate", rep, ":", sapply(gene_bins, length), "\n")

    # ---- 3. Process each bin ----
    for (i in seq_along(gene_bins)) {
        bin_genes <- gene_bins[[i]]
        cat("Processing replicate", rep, "bin", i, "with", length(bin_genes), "genes...\n")

        # ---- Map GENENAME -> ENTREZ ID ----
        mapped <- bitr(bin_genes,
                       fromType = "GENENAME",
                       toType   = "ENTREZID",
                       OrgDb    = org.Sc.sgd.db)
        entrez_genes <- unique(mapped$ENTREZID)

        if (length(entrez_genes) == 0) {
            cat("No valid ENTREZ IDs for bin", i, "in replicate", rep, "\n")
            next
        }

        # ---- Run KEGG enrichment ----
        kegg_res <- enrichKEGG(
            gene         = entrez_genes,
            organism     = "sce",
            keyType      = "ncbi-geneid",
            pvalueCutoff = 1
        )

        if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
            cat("No enrichment found for bin", i, "in replicate", rep, "\n")
            next
        }

        # ---- Shorten descriptions ----
        kegg_res@result$Description <- str_replace(
            kegg_res@result$Description,
            " - Saccharomyces cerevisiae \\(budding yeast\\)", ""
        )

        # ---- Map Entrez IDs in results to GENENAME ----
        all_entrez <- unique(unlist(strsplit(kegg_res@result$geneID, "/")))
        entrez2gene <- AnnotationDbi::select(
            org.Sc.sgd.db,
            keys = all_entrez,
            columns = c("GENENAME"),
            keytype = "ENTREZID"
        )
        lookup <- setNames(entrez2gene$GENENAME, entrez2gene$ENTREZID)
        kegg_res@result$geneID <- sapply(kegg_res@result$geneID, function(ids) {
            genes <- unlist(strsplit(ids, "/"))
            paste(lookup[genes], collapse = "/")
        })

        # ---- Save CSV ----
        csv_file <- file.path(rep_dir, paste0("KEGG_bin_", i, ".csv"))
        write.csv(as.data.frame(kegg_res), csv_file, row.names = FALSE)

        cat("Results saved:", csv_file, "\n")
    }
}

cat("\nAll bootstrap replicates complete.\n")
