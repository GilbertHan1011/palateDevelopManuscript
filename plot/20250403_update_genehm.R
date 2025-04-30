#!/usr/bin/env Rscript

# Script to create gene trajectory heatmap from pre-computed data
# This script loads pre-computed GAM model, gene clusters, and scaled matrix
# then generates and saves the final heatmap visualization

# Load required libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(1234)

# Load pre-computed data
gam <- readRDS("../2024.4_scATAC/processed_data/trajectory/20250107_gam_threelineage.Rds")
gene_clusters <- read.csv("descriptive_results/rna_trajectory_gene/20250117_gene_clusters.csv")
scaled_mat <- read.csv("descriptive_results/rna_trajectory_gene/20250117_gene_scaled_mat.csv",row.names = NULL)

pseudotime1 <- read.csv("../2024.4_scATAC/processed_data/trajectory/5.1_rna_celloracle.pseudotime_run2.csv",row.names = 1)
seurat$pseudo <- pseudotime1$Pseudotime
# Create pseudotime metadata
# Load relevant metadata
# (You'll need to adjust this based on where your metadata is stored)
seurat <- dietPalate # Adjust path as needed
meta <- seurat@meta.data[c("level1_anno", "group")]

# Get row names for scaled_mat (assuming they're in the first column)
if(!is.matrix(scaled_mat)) {
  row_names <- scaled_mat[,1]
  scaled_mat <- scaled_mat[,-1]
  scaled_mat <- as.matrix(scaled_mat)
  rownames(scaled_mat) <- row_names
}

# Make sure gene_clusters has the gene names as the first column
if("Gene" %in% colnames(gene_clusters)) {
  rownames(gene_clusters) <- gene_clusters$Gene
}

# Create column split factor for the heatmap
column_split <- rep(c("stem", "K14"), each = 100)
column_split <- factor(column_split, levels = c("stem", "K14"))

# Function to create a color mapping that includes gray for NA
create_color_map <- function(values, palette_name = "Set1") {
  unique_values <- unique(values[!is.na(values)])
  n_colors <- length(unique_values)

  if (n_colors == 0) {
    return(c("NA" = "gray"))
  } else if (n_colors == 1) {
    colors <- brewer.pal(3, palette_name)[1]
  } else if (n_colors == 2) {
    colors <- brewer.pal(3, palette_name)[1:2]
  } else if (n_colors <= 8) {
    colors <- brewer.pal(n_colors, palette_name)
  } else {
    colors <- colorRampPalette(brewer.pal(8, palette_name))(n_colors)
  }

  color_map <- setNames(colors, unique_values)
  if (any(is.na(values))) {
    color_map <- c(color_map, "NA" = "gray")
  }
  return(color_map)
}

# Prepare gene clusters
# Define the desired order and new names
cluster_order <- c(4, 6, 1, 2, 5, 3, 8, 7)
new_cluster_names <- paste0("GC", 1:8)

# Create a named vector for renaming
rename_vector <- setNames(new_cluster_names, cluster_order)

# Extract the cluster assignments
sigGene <- gene_clusters$Cluster

# Reorder and rename the cluster vector
sigGene_ordered <- factor(sigGene, levels = cluster_order)
levels(sigGene_ordered) <- rename_vector[levels(sigGene_ordered)]

# Process metadata for annotations
# You need to provide these files or extract this information from your Seurat object
# This is just a placeholder - you'll need to adjust this based on your actual data
if(exists("seurat")) {
  # Extract metadata from cells if using actual Seurat object
  # This is a simplified version - adjust according to your actual data structure
  lineage1 <- c("stem cells", "Transit", "K6+ cells")
  lineage2 <- c("K14(-)", "K6+ cells")
  lineage1Cell <- colnames(seurat)[seurat$level1_anno%in%lineage1]
  lineage2Cell <- colnames(seurat)[seurat$level1_anno%in%lineage2]


  meta1 <- meta[lineage1Cell,]
  meta1$pseudo <- seurat@meta.data[lineage1Cell, "pseudo"]

  meta2 <- meta[lineage2Cell,]
  meta2$pseudo <- seurat@meta.data[lineage2Cell, "pseudo"]

  # Create time grid
  time_grid <- seq(0, 1, length.out = 100)

  # Process metadata
  # Define a simplified process_metadata function
  process_metadata <- function(meta_df, time_grid) {
    # Convert inputs to appropriate types and add error checking
    meta_df <- as.data.frame(meta_df)
    time_grid <- as.numeric(time_grid)

    # Verify required columns exist
    required_cols <- c("pseudo", "level1_anno", "group")
    if (!all(required_cols %in% colnames(meta_df))) {
      stop("Missing required columns in meta_df: ",
           paste(setdiff(required_cols, colnames(meta_df)), collapse = ", "))
    }

    # Step 1: Create time intervals
    meta_df$time_interval <- findInterval(meta_df$pseudo, time_grid)

    # Step 2: Filter valid intervals
    meta_df <- meta_df[meta_df$time_interval > 0, ]

    # Step 3: Group and summarize - split into separate operations
    level1_anno_summary <- tapply(meta_df$level1_anno, meta_df$time_interval,
                                  function(x) names(which.max(table(x))))
    group_summary <- tapply(meta_df$group, meta_df$time_interval,
                            function(x) names(which.max(table(x))))
    pseudo_time_summary <- tapply(meta_df$pseudo, meta_df$time_interval, mean)

    # Combine the summaries
    new_meta <- data.frame(
      time_interval = as.numeric(names(level1_anno_summary)),
      max_level1_anno = unname(level1_anno_summary),
      max_group = unname(group_summary),
      pseudo_time = unname(pseudo_time_summary)
    )

    # Step 4: Create complete intervals
    complete_intervals <- data.frame(time_interval = 1:100)

    # Step 5: Merge with complete intervals
    final_meta <- merge(complete_intervals, new_meta,
                        by = "time_interval", all.x = TRUE)

    # Step 6: Sort by time interval
    final_meta <- final_meta[order(final_meta$time_interval), ]

    return(final_meta)
  }

  new_meta1 <- process_metadata(meta1, time_grid)
  new_meta2 <- process_metadata(meta2, time_grid)
  new_meta <- rbind(new_meta1, new_meta2)

  # Factor levels for annotations
  new_meta$max_level1_anno <- factor(new_meta$max_level1_anno)

  # Create color mappings for annotations
  level1_anno_colors <- create_color_map(new_meta$max_level1_anno, "Set1")
  group_colors <- create_color_map(new_meta$max_group, "Set2")

  # Create annotation
  top_anno <- HeatmapAnnotation(
    level1_anno = new_meta$max_level1_anno,
    group = new_meta$max_group,
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
      level1_anno = list(title = "Level 1 Anno"),
      group = list(title = "Group")
    ),
    col = list(
      level1_anno = level1_anno_colors,
      group = group_colors
    )
  )
} else {
  # Create an empty annotation if metadata is not available
  top_anno <- NULL
}

# Create and save the heatmap
#pdf("descriptive_results/rna_trajectory_gene/gene_trajectory_heatmap.pdf", width = 6, height = 8)


gene_list <- c("Krt5", "Krt14", "Krt6a", "Wnt3",
               "Grhl3", "Cldn3", "Elf5", "Klf4", "Sprr1a")


# Find genes that exist in our dataset
labelsT <- gene_list
rownames(scaled_mat) <- gene_clusters$Gene
position <- sapply(labelsT, function(gene) which(rownames(scaled_mat) == gene))

# Create row annotation with gene labels
har <- rowAnnotation(link = anno_mark(at = position, labels = labelsT,
                                      labels_gp = gpar(fontsize = 9, fontface = "bold"),
                                      link_width = unit(1.5, "cm")))

hmAll <- Heatmap(scaled_mat,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 cluster_row_slices = FALSE,
                 show_row_dend = FALSE,
                 row_split = gene_clusters$Cluster,
                 column_split = column_split,
                 col = colorRamp2(c(-2, 0, 4), c("Deepskyblue3", "white", "red")),
                 border = TRUE,
                 na_col = "grey",
                 row_title_gp = gpar(fontsize = 10),
                 row_title_rot = 0,  # Horizontal text
                 row_gap = unit(1, "mm"),  # Add some gap between row clusters
                 top_annotation = top_anno,right_annotation = har)

# Draw the heatmap
hmAll_drawn <- draw(hmAll)
pdf("descriptive_results/rna_trajecotry/20250403_rna_gene_hm.pdf",width = 8,height = 8)
hmAll_drawn
dev.off()

# Save the heatmap object
saveRDS(hmAll_drawn, "process/trajectory/gene_trajectory_heatmap.Rds")

cat("Gene trajectory heatmap created and saved.\n")
