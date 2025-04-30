# plot heatmap
hmdata <- read.csv("process/SNP/SNP2ell/20250215_heatmap_res_df_5round.csv",row.names = 1)
colnames(hmdata) <- c("stem cells","K5(-)", "Transit","Mature K5(-) ", 
                      "K6+ cells",  "K14(-)","Shh(+)")

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(grid)

# Filter data for boxplot - only include rows where max value > 2
boxplot_data <- hmdata[apply(hmdata, 1, function(x) max(x) > 2), ]
cat("Number of genes/peaks with max z-score > 2:", nrow(boxplot_data), "\n")

# Identify top 10 genes with highest scores
row_max_scores <- apply(hmdata, 1, max)
top_genes_idx <- order(row_max_scores, decreasing = TRUE)[1:5]
top_genes <- rownames(hmdata)[top_genes_idx]

# Create right annotation to mark top genes
right_annotation = rowAnnotation(
  link = anno_mark(
    at = top_genes_idx,
    labels = top_genes,
    labels_gp = gpar(fontsize = 10, fontface = "bold", family = "sans"),
    link_width = unit(1, "cm"),
    link_gp = gpar(lwd = 1, col = "black"),
    padding = unit(0.5, "mm")
  ),
  width = unit(2.5, "cm")
)

# Create boxplot annotation with filtered data
top_annotation = HeatmapAnnotation(
  boxplot = anno_boxplot(
    boxplot_data,  # Using filtered data for boxplot
    axis = TRUE,
    border = TRUE,
    box_width = 0.6,
    pch = 16,  # Solid circles for points
    size = unit(1.5, "mm"),  # Smaller points
    gp = gpar(fill = "#92C5DE"), # Color of the boxes
    outlier_gp = gpar(col = "grey30", fill = "grey30", alpha = 0.5), # Outlier styling
    height = unit(3, "cm"),
    ylim = c(-3, 18),  # Set fixed y-axis range from -3 to 5
    outline = FALSE,  # Remove outliers for cleaner look
    whisker = 0.1  # Set whiskers at 0.1 the interquartile range
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_rot = 0,
  gap = unit(1, "mm"),
  show_legend = FALSE
)

# Create and draw the heatmap with enhanced aesthetics
pdf("descriptive_results/snp/20250226_snp_hm.pdf", 
    width = 6, height = 8)  # Increased width to accommodate gene labels

# Draw the heatmap with the boxplot annotation and gene labels
draw(Heatmap(hmdata,  # Still use full data for heatmap
        cluster_columns = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        show_row_names = FALSE,
        
        # Add boxplot annotation at top
        top_annotation = top_annotation,
        
        # Add gene labels on right
        right_annotation = right_annotation,
        
        # Enhanced color scheme
        col = colorRamp2(c(-2, 1, 4, 7, 10), 
                         c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")),
        
        # Improved column styling
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
        
        # Enhanced titles
        column_title = "SNP2Cell score",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        name = "SNP2Cell score",
        
        # Better spacing and borders
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 9),
          title_position = "topcenter",
          legend_height = unit(3, "cm")
        ),
        
        # Add border
        border = TRUE
))

dev.off()
