library(Seurat)

atac = readRDS("process/framework/obj/combined_E12_E14_signac.Rds")
#atacGlue <- zellkonverter::readH5AD("process/framework/obj/atac_scglue_prepare.h5ad")
peakUmap <- read.csv("process/framework/reduction/7.17_combined_umap_full.csv",row.names = 1)
colnames(peakUmap) <- c("umapPeak_1","umapPeak_2")
atac@reductions$peakvi_umap <-  CreateDimReducObject(embeddings = as.matrix(peakUmap),key = "umapPeak_", assay = "ATAC")

glueUmap <- read.csv("process/framework/reduction/8.17_scglue_combined_umap.csv",row.names = 1)
colnames(glueUmap) <- c("umap_1","umap_2")
atac@reductions$glue_umap <-  CreateDimReducObject(embeddings = as.matrix(glueUmap),key = "umap_", assay = "ATAC")

cluster <- read.csv("process/framework/cluster/20241211_knn_label.csv",row.names = 1)
atac$knn_label <- cluster$knn_label


colorPalate1 <- c("#E6E2F1", "#F8E8E3", "#F6D481", "#F5A556",
                             "#FB9FB9", "#C94D6C",
                             "#4F4E46")
colorPalate2 <-  c("#462D2E", "#BA8F7D", "#AB6B7A",
                   "#CB6B6E", "#E09DBD", "#47324E", "#FFD1BF")
colorPalate3 <- c("#282828","#F2C9D5", "#B43E44", "#8B6D9C", "#496496", "#904869","#FADF92" )
colorPalate4 <- colorRampPalette(colorPalate3)(7)
DimPlot(atac,reduction = "glue_umap",group.by = "knn_label",cols = colorPalate3)


rna <- readRDS("process/framework/obj/RNA_E12E14_diet.Rds")

rna@reductions$glue_umap <-  CreateDimReducObject(embeddings = as.matrix(glueUmap),key = "umap_", assay = "RNA")
DimPlot(rna,reduction = "glue_umap",group.by = "level1_anno",cols = colorPalate3)

combined_all = zellkonverter::readH5AD("process/framework/obj/20241211_combined_scglue.h5ad")
combined_seurat <- as.Seurat(combined_all,assay = "X",counts = NULL,data = "data")
combined_seurat <- as.Seurat(combined_all,counts = NULL,data = "X")

combined_seurat$knn_label <- NULL
combined_seurat@meta.data[rownames(cluster),"knn_label"] <- cluster$knn_label
combined_seurat@meta.data[colnames(rna),"knn_label"] <- rna$level1_anno
DimPlot(combined_seurat,reduction = "X_umap",group.by = "knn_label",cols = colorPalate3)
ggsave("descriptive_results/rna_atac_integration/20241211_scglue_umap.pdf",width = 6,height = 6)

DimPlot(combined_seurat,reduction = "X_umap",group.by = "modal",cols = c("#FADF92","#496496"))

ggsave("descriptive_results/rna_atac_integration/20241211_scglue_umap_modal.pdf",width = 6,height = 6)


p1 <- DimPlot(atac,reduction = "peakvi_umap",group.by = "knn_label",cols = colorPalate3)+ggtitle("ATAC")
p2 <- DimPlot(rna,reduction = "scviumap",group.by = "level1_anno",cols = colorPalate3)+ggtitle("RNA")
p1|p2
ggsave("descriptive_results/rna_atac_integration/20241211_split_cluster_umap.pdf",width = 10,height = 4)
