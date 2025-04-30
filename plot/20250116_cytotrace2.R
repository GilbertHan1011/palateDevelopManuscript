library(CytoTRACE2)
library(RColorBrewer)
seurat <- readRDS("../2024_7.ATAC_E12/process/framework/obj/RNA_E12E14_diet.Rds")

cytotrace2_result <- cytotrace2(seurat@assays$RNA@data, is_seurat = FALSE,  species = 'mouse')
seurat$ct_score <- cytotrace2_result$CytoTRACE2_Score
annotation <- data.frame(phenotype = seurat@meta.data$level1_anno) %>% set_rownames(., colnames(seurat))
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = FALSE,expression_data = seurat@assays$RNA@data)
plots$CytoTRACE2_Potency_UMAP
seurat$ct_score_relative <- cytotrace2_result$CytoTRACE2_Relative
seurat$preKNN_score <- cytotrace2_result$preKNN_CytoTRACE2_Score
saveRDS(cytotrace2_result, "process//trajectory//20250117_cytotrace.Rds" )

FeaturePlot(seurat,"preKNN_score",reduction = "scviumap")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))


FeaturePlot(seurat,"ct_score",reduction = "scviumap")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/rna_trajecotry/20250116_cytotrace_score.pdf",width = 6,height = 6)
# FeaturePlot(seurat,"ct_score_relative",reduction = "scviumap")+
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Boxplot_byPheno

cytotrace2_result2 <- cytotrace2(seurat@assays$RNA@counts, is_seurat = FALSE,  species = 'mouse')
