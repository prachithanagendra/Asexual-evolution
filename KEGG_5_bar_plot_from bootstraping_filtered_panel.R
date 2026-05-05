#!/usr/bin/env Rscript


# ---- Load libraries ----
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(grid)

# ---- Function to create individual barplots ----
make_barplot <- function(df, title) {
  df_bar <- df %>%
    filter(combined_padjust < 0.05) %>%
    arrange(combined_padjust) %>%
    head(10) %>%
    rowwise() %>%
    mutate(gene_count = length(unique(unlist(str_split(geneID, "/|;\\s*"))))) %>%
    ungroup()
  
  # Wrap pathway names
  df_bar$Description <- factor(
    str_wrap(df_bar$Description, width = 30),
    levels = rev(str_wrap(df_bar$Description, width = 30)[order(df_bar$combined_padjust)])
  )
  
  ggplot(df_bar, aes(x = Description, y = gene_count, fill = combined_padjust)) +
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
    scale_y_continuous(expand = c(0, 2), limits = c(0, 70)) +
    labs(x = NULL, y = NULL, title = title) +   # remove local axis labels
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y  = element_text(size = 11, color = "black"),
      axis.text.x  = element_text(size = 11, color = "black"),
      legend.position = "bottom",
      legend.key.width = unit(1.0, "cm"),
      legend.key.height = unit(1.5, "cm"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 11),
      axis.title   = element_blank()   # suppress titles inside each plot
    )
}

# ---- Load your four datasets ----
df_300  <- read.csv("/CNV_GO/KEGG/Gal_dup/300/3G_dup_combined/3G-dup_final_combined_replicates.csv",  stringsAsFactors = FALSE, check.names = FALSE)
df_600  <- read.csv("/CNV_GO/KEGG/Gal_dup/600/6G_dup_combined/6G-dup_final_combined_replicates.csv",  stringsAsFactors = FALSE, check.names = FALSE)
df_900  <- read.csv("/CNV_GO/KEGG/Gal_dup/900/9G_dup_combined/9G-dup_final_combined_replicates.csv",  stringsAsFactors = FALSE, check.names = FALSE)
df_1200 <- read.csv("/CNV_GO/KEGG/Gal_dup/1200/12G_dup_combined/12G-dup_final_combined_replicates.csv", stringsAsFactors = FALSE, check.names = FALSE)

# ---- Make individual plots ----
p300  <- make_barplot(df_300,  "300 generations")
p600  <- make_barplot(df_600,  "600 generations")
p900  <- make_barplot(df_900,  "900 generations")
p1200 <- make_barplot(df_1200, "1200 generations")

# ---- Combine into 2x2 panel with SINGLE legend ----
panel <- (p300 | p600) / (p900 | p1200) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ---- Ensure room for outer labels ----
panel <- panel + plot_annotation(
  theme = theme(plot.margin = margin(t = 20, r = 20, b = 60, l = 80))
)

# draw the panel
print(panel)

# add common horizontal label (centered at bottom)
grid.text("Unique Gene Count", x = 0.5, y = unit(0.05, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))

# add common vertical label on the left
grid.text("Pathway", x = unit(0.05, "npc"), y = 0.5,
          rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))

dev.off()

# ---- Save PNG with manual bottom + left labels ----
out_png <- "C:/PhD/APS/5th APS/NGF/Fig3/CNV_GO/KEGG/Gal_dup/gal_dup_KEGG_panel.png"
png(out_png, width = 14, height = 10, units = "in", res = 600)

print(panel)

grid.text("Unique Gene Count", x = 0.5, y = unit(0.05, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("Pathway", x = unit(0.05, "npc"), y = 0.5,
          rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))

dev.off()
