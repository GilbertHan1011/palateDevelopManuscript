library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(clusterProfiler)
draw_graph <- read.csv("../2024.4_scATAC/processed_data/embedding/RNA_scvi_20neibor_drawgraph.csv",row.names = 1)
colnames(draw_graph) <- c(1,2)
drawReduction <- CreateDimReducObject(embeddings = as.matrix(draw_graph), key = "FA_", assay = DefaultAssay(seurat))
seurat@reductions$drawGraph <- drawReduction
Idents(seurat) <- seurat$level1_anno


colorPalate3 <- c("#282828","#F2C9D5", "#B43E44", "#8B6D9C", "#496496", "#904869","#FADF92")
names(colorPalate3) <- c("K14(-)",  "K5(-)", "K6+ cells", "Mature K5(-) ","Shh(+)", "stem cells", "Transit")
DimPlot(seurat,group.by = "level1_anno",reduction = "drawGraph",cols = colorPalate3)
ggsave("descriptive_results/rna_trajecotry/20250124_rna_drawgraph.pdf",width = 6,height = 5)

pseudotime_rna <- read.csv("../2024.4_scATAC/processed_data/trajectory/5.1_rna_celloracle.pseudotime_run2.csv",row.names = 1)
seurat$K6_pseudo <- pseudotime_rna$Pseudotime_Lineage_K6
seurat$K14 <- pseudotime_rna$Pseudotime_Lineage_K14
seurat$K6_pseudo[(seurat$K6_pseudo>0.75) & (seurat$level1_anno=="Transit")] <- NA
FeaturePlot(seurat,"K6_pseudo",reduction = "scviumap")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
FeaturePlot(seurat,"K6_pseudo",reduction = "drawGraph")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/atac/20250124_rna_drawgraph_drawgraph.pdf",width = 6,height = 5)

FeaturePlot(seurat,"K14",reduction = "drawGraph")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/rna_trajecotry/20250124_rna_drawgraph_drawgraph_K14.pdf",width = 7,height = 5)


#== clusterprofiler of trajectory gene clusters-----------------

rna_genes <- read.csv("descriptive_results/rna_trajectory_gene/20250117_gene_clusters.csv")
rna_genes_list <- split(rna_genes$Gene, rna_genes$Cluster)
ck <- clusterProfiler::compareCluster(rna_genes_list,fun = enrichGO,OrgDb='org.Mm.eg.db',ont="BP",keyType="SYMBOL")

dotplot(ck)
ggsave("descriptive_results/20250124_gene_enrich.pdf",width = 8,height = 10)
write.csv(ck@compareClusterResult,"process/trajectory/20250124_rna_tradeseq_enrich.csv")
saveRDS(ck,"process/trajectory/20250124_rna_tradeseq_enrich.Rds")
