# This script is used for plot atac trajectory
library(Seurat)
library(RColorBrewer)
atac_traj_red <- read.csv("process/framework/reduction/20250118_snap_harmony_drawgraph.csv",row.names = 1)
colnames(atac_traj_red) <- c(1,2)
prefix_map <- c(
  "C2" = "E12_2",
  "C1" = "E12_1",
  "B2" = "E14_2",
  "B1" = "B1"
  # Add more mappings if needed
)


# Function to replace prefix
# Vectorized function to replace prefixes
# Modified function to remove the original prefix
replace_prefix <- function(x, prefix_map) {
  # Extract the prefix (C1, C2, etc.)
  old_prefixes <- sub("^([^_]+)_.*", "\\1", x)

  # Get the barcode part (everything after the underscore, before -1)
  barcodes <- sub("^[^_]+_(.*)-1", "\\1", x)

  # Get new prefixes
  new_prefixes <- prefix_map[old_prefixes]

  # Combine with barcode
  new_strings <- paste0(new_prefixes, "#", barcodes, "-1")
  return(new_strings)
}

newid = replace_prefix(rownames(atac_traj_red), prefix_map)

rownames(atac_traj_red) <- newid
drawReduction <- CreateDimReducObject(embeddings = as.matrix(atac_traj_red),
                                      key = "FA_", assay = DefaultAssay(motifSeurat))
motifSeurat@reductions$drawGraph <- drawReduction
color_palate = c(`K14(-)` = "#282828", `K5(-)` = "#F2C9D5", `K6+ cells` = "#B43E44",
                 `Mature K5(-) ` = "#8B6D9C", `Shh(+)` = "#496496", `stem cells` = "#904869",
                 Transit = "#FADF92")
#Idents(seurat) <- seurat$level1_anno
#pseudotime1 <- read.csv("process/trajectory//20250121_atac_celloracle.pseudotime_.csv",row.names = 1)


DimPlot(motifSeurat,group.by = "level1_anno",reduction = "drawGraph",cols = color_palate)
ggsave("descriptive_results/atac_trajectory/20250124_drawgraph_level1.pdf",width = 7,height = 5)

#== plot pseudotime--------------------

motifSeurat$pseudo_K6 <- pseudotime1$Pseudotime_Lineage_K6
motifSeurat$pseudo_K14 <- pseudotime1$Pseudotime_Lineage_K14

FeaturePlot(motifSeurat,"pseudo_K6",reduction = "drawGraph")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/rna_trajecotry/20250124_atac_drawgraph_drawgraph_K6.pdf",width = 7,height = 5)
FeaturePlot(motifSeurat,"pseudo_K14",reduction = "drawGraph")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/rna_trajecotry/20250124_atac_drawgraph_drawgraph_K14.pdf",width = 7,height = 5)
