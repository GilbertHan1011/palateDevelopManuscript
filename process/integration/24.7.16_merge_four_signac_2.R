# merge four obj------------------------------
setwd("../2024_7.ATAC_E12/")
peak10xC1 <- read.table("data/E12.5_1/peaks.bed",col.names = c("chr", "start", "end"))
peak10xC2 <- read.table("data/E12.5_2/peaks.bed",col.names = c("chr", "start", "end"))
grC1 <- makeGRangesFromDataFrame(peak10xC1)
grC2 <- makeGRangesFromDataFrame(peak10xC2)
combined.peaks <- reduce(x = c(grB1,grB2,grC1, grC2))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
cellsB1 <- B1_filtered %>% colnames()
cellsB2 <- B2_filtered %>% colnames()
dir.create("../2024.4_scATAC/processed_data/framework/filtered_cells")
write.csv(cellsB1,"../2024.4_scATAC/processed_data/framework/filtered_cells/B1_cells_Epi")
write.csv(cellsB2,"../2024.4_scATAC/processed_data/framework/filtered_cells/B2_cells_Epi")
cellsC1 <- read.csv("process/framework/filtered_cells/C1_cells_Epi",row.names = 1) %>% unlist
cellsC2 <- read.csv("process/framework/filtered_cells/C2_cells_Epi",row.names = 1) %>% unlist
names(cellsC1) <- NULL
names(cellsC2) <- NULL

frags.B1 <- CreateFragmentObject(
  path = "../2024.4_scATAC/data/raw/B1/atac_fragments.tsv.gz",
  cells = cellsB1
)
frags.B2 <- CreateFragmentObject(
  path = "../2024.4_scATAC/data/raw/B2/fragments.tsv.gz",
  cells = cellsB2
)


frags.C1 <- CreateFragmentObject(
  path = "data/E12.5_1/fragments.tsv.gz",
  cells = cellsC1
)
frags.C2 <- CreateFragmentObject(
  path = "data/E12.5_2//fragments.tsv.gz",
  cells = cellsC2
)

B1.counts <- FeatureMatrix(
  fragments = frags.B1,
  features = combined.peaks,
  cells = cellsB1
)

B2.counts <- FeatureMatrix(
  fragments = frags.B2,
  features = combined.peaks,
  cells = cellsB2
)

C1.counts <- FeatureMatrix(
  fragments = frags.C1,
  features = combined.peaks,
  cells = cellsC1
)

C2.counts <- FeatureMatrix(
  fragments = frags.C2,
  features = combined.peaks,
  cells = cellsC2
)


B1_assay <- CreateChromatinAssay(B1.counts, fragments = frags.B1)
B1_Obj <- CreateSeuratObject(B1_assay, assay = "ATAC", meta.data=B1_filtered@meta.data)

B2_assay <- CreateChromatinAssay(B2.counts, fragments = frags.B2)
B2_Obj <- CreateSeuratObject(B2_assay, assay = "ATAC", meta.data=B2_filtered@meta.data)


C1_assay <- CreateChromatinAssay(C1.counts, fragments = frags.C1)
C1_Obj <- CreateSeuratObject(C1_assay, assay = "ATAC")

C2_assay <- CreateChromatinAssay(C2.counts, fragments = frags.C2)
C2_Obj <- CreateSeuratObject(C2_assay, assay = "ATAC")



B1_Obj$batch <- "B1"
B2_Obj$batch <- "B2"

C1_Obj$batch <- "C1"
C2_Obj$batch <- "C2"

combined4 <- merge(
  x = B1_Obj,
  y = list(B2_Obj,C1_Obj,C2_Obj),
  add.cell.ids = c("B1", "B2","C1","C2")
)
combined4 <- RunTFIDF(combined4)
combined4 <- FindTopFeatures(combined4, min.cutoff = 20)
combined4 <- RunSVD(combined4)
combined4 <- RunUMAP(combined4, dims = 2:50, reduction = 'lsi')
DepthCor(combined4)
DimPlot(combined4,group.by = "batch")
ggsave("result/intergrate/combine4_before_harmony.pdf",width = 5,height = 5)
saveRDS(combined4,"process/framework/obj/combined_E12_E14_signac.Rds")

combined4@reductions$umap_before <- combined4@reductions$umap

combined4 <- RunHarmony(
  object = combined4,
  group.by.vars = 'batch',
  reduction.use = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE,
  dims.use = 2:50
)
combined4 <- RunUMAP(combined4, dims = 1:20, reduction = 'harmony')
p2 <- DimPlot(combined4,group.by = "batch")+ggtitle("after harmony")
p1 <- DimPlot(combined4,group.by = "batch",reduction = "umap_before")+ggtitle("before harmony")
p1|p2
ggsave("result/intergrate/combine4_before_and_after_harmony.pdf",width = 10,height = 5)

combined4@reductions$umap_after <- combined4@reductions$umap
combined4Sce <- as.SingleCellExperiment(combined4,assay = "ATAC")
writeH5AD(combined4Sce,"process/framework/obj/combined_E12_E14_signc_zell.h5ad")
SaveH5Seurat(combined4, filename = "process/framework/obj/combined_E12_E14_signac.h5Seurat")
Convert("process/framework/obj/combined_E12_E14_signac.h5Seurat", dest = "h5ad")
combined4Sce@assays@data$logcounts <- combined4Sce@assays@data$counts
library(reticulate)
use_condaenv('py310')
sceasy::convertFormat(combined4Sce, from="sce", to="anndata",
                      outFile='process/framework/obj/combined_E12_E14_signc_sceasy.h5ad')
combined4 <- FindNeighbors(object = combined4, reduction = 'harmony', dims = 1:20)
combined4 <- FindClusters(object = combined4, verbose = FALSE, algorithm = 3,resolution = 0.05)
DimPlot(combined4)
write.csv(combined4$seurat_clusters,"process/framework/cluster/combine_harmony.csv")
