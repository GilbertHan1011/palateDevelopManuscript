rna <- readRDS("process/framework/obj/RNA_E12E14_diet.Rds")
rna_full <- readRDS("../2024.4_scATAC/data/palate_epi_MAGIC_periderm_MAGIC_harmony_renamed.rds")

FeaturePlot(rna,"Arid3a",reduction = "scviumap")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
FeaturePlot(rna_full,"Arid3a",reduction = "scviumap")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/regulon_arid3a/20250225_arid3a_umap.pdf",width = 6,height = 5)
#rna_full@reductions$scviumap

drawgraph <- read.csv("../2024.4_scATAC/processed_data/embedding/RNA_scvi_20neibor_drawgraph.csv",row.names = 1)
colnames(drawgraph) <- c(1,2)
rna_full@reductions$draw_graph <-  CreateDimReducObject(embeddings = as.matrix(drawgraph),key = "DA_", assay = "RNA")
FeaturePlot(rna_full,"Arid3a",reduction = "draw_graph")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/regulon_arid3a/20250225_arid3a_drawgraph.pdf",width = 6,height = 5)
FeaturePlot(rna_full,"Grhl1",reduction = "draw_graph")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/regulon_arid3a/20250225_grhl1_drawgraph.pdf",width = 6,height = 5)
