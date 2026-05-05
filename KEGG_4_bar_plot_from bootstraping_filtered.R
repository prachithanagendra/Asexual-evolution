#!/usr/bin/env Rscript

# ---- Load libraries ----
library(ggplot2)
library(dplyr)
library(stringr)
library(igraph)
library(ggraph)
library(tidyr)

# ---- User Input ----
input_file  <- "/CNV_GO/KEGG/Glu_Del-Gal_Dup/g_del_G_dup_combined/g_del_G_dup_combined_final_combined_replicates.csv"
output_dir  <- "/CNV_GO/KEGG/Glu_Del-Gal_Dup/g_del_G_dup_combined/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 1. Load data ----
df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

# ---- 1a. Filter rows with bins_present >= 2 ----
df <- df %>% filter(bins_present >= 2)

# ---- 2. Bar plot (Top 10 by combined_padjust, cutoff < 0.05) ----
df_bar <- df %>%
  filter(combined_padjust < 0.05) %>%
  arrange(combined_padjust) %>%
  head(10) %>%
  # Calculate unique gene count from geneID column
  rowwise() %>%
  mutate(gene_count = length(unique(unlist(str_split(geneID, "/|;\\s*"))))) %>%
  ungroup()

# Reorder factor for decreasing order (most significant on top)
df_bar$Description <- factor(df_bar$Description, levels = rev(df_bar$Description[order(df_bar$combined_padjust)]))

# ---- 3. Plot ----
# ---- 3. Plot ----
p <- ggplot(df_bar, aes(x = Description, y = gene_count, fill = combined_padjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(
    colours = c("#FF6666", "#FFCC66", "#66CC99", "#6699CC"),
    values = scales::rescale(c(0, 0.01, 0.03, 0.05)),
    limits = c(0, 0.05),
    breaks = c(0, 0.01, 0.02, 0.05),
    labels = c("0", "0.01", "0.02", "0.05"),
    name = "Adjusted p-value"
  ) +
  coord_flip() +
  labs(x = "Pathway", y = "Unique Gene Count") +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black")
  )

# Save PNG
png(file.path(output_dir, "g_del_G_dup_kegg_barplot.png"), width = 8, height = 6, units = "in", res = 300)
print(p)
dev.off()
