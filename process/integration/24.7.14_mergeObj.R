#== merge two obj with signac-------------------
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(harmony)
plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
DefaultAssay(C1_filter) <- "peaks"
peak <- C1_filter %>% rownames()
peak10xC1 <- read.table("data/E12.5_1/peaks.bed",col.names = c("chr", "start", "end"))
peak10xC2 <- read.table("data/E12.5_2/peaks.bed",col.names = c("chr", "start", "end"))
grC1 <- makeGRangesFromDataFrame(peak10xC1)
grC2 <- makeGRangesFromDataFrame(peak10xC2)
combined.peaks <- reduce(x = c(grC1, grC2))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

cellsC1 <- C1_filter %>% colnames()
cellsC2 <- C2_filtered %>% colnames()
write.csv(cellsC1,"process/framework/filtered_cells/C1_cells_Epi")
write.csv(cellsC2,"process/framework/filtered_cells/C2_cells_Epi")
frags.C1 <- CreateFragmentObject(
  path = "data/E12.5_1/fragments.tsv.gz",
  cells = cellsC1
)
frags.C2 <- CreateFragmentObject(
  path = "data/E12.5_2//fragments.tsv.gz",
  cells = cellsC2
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

C1_assay <- CreateChromatinAssay(C1.counts, fragments = frags.C1)
C1_Obj <- CreateSeuratObject(C1_assay, assay = "ATAC", meta.data=C1_filter@meta.data)

C2_assay <- CreateChromatinAssay(C2.counts, fragments = frags.C2)
C2_Obj <- CreateSeuratObject(C2_assay, assay = "ATAC", meta.data=C2_filtered@meta.data)

C1_Obj$batch <- "C1"
C2_Obj$batch <- "C2"

combined <- merge(
  x = C1_Obj,
  y = list(C2_Obj),
  add.cell.ids = c("C1", "C2")
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DepthCor(combined)
DimPlot(combined,group.by = "batch")
ggsave("result/intergrate/7.15_before_hamony.pdf",width = 6,height = 6)
combined@reductions$umap_sep <- combined@reductions$umap

p1 <- DimPlot(combined,group.by = "batch") + ggtitle("after harmony")
p2 <- DimPlot(combined,group.by = "batch",reduction = "umap_sep") + ggtitle("before harmony")
p2|p1

ggsave("result/intergrate/7.16_after_harmony.pdf")

combined <- RunHarmony(
  object = combined,
  group.by.vars = 'batch',
  reduction.use = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE,
  dims.use = 2:50
)
combined <- RunUMAP(combined, dims = 1:20, reduction = 'harmony')
saveRDS(combined,"process/framework/obj/C1_C2_combined.Rds")
DefaultAssay(C1_filter) <- "RNA"
DefaultAssay(C2_filtered) <- "RNA"
mergeObj <- merge(C1_filter,C2_filtered,add.cell.ids = c("C1", "C2"))

RNA1 <- CreateSeuratObject(C1_filter@assays$RNA)
RNA2 <- CreateSeuratObject(C2_filtered@assays$RNA)
RNAcombined <- merge(RNA1,list(RNA2),
                     add.cell.ids = c("C1", "C2"))
#mergeActivity <- merge(C1_filter@assays$RNA,C2_filtered@assays$RNA,
#                       add.cell.ids = c("C1", "C2"))



combined@assays$activity <- RNAcombined@assays$RNA
#colnames(mergeActivity) <- colnames(combined)
#mergeActivity
combined@assays$activity <- mergeActivity
FeaturePlot(combined,c("Krt6a","Krt40","Krtdap","Dsg1b"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))


FeaturePlot(combined,c("Krt6a","Krt14","Krt5"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))

saveRDS(RNAcombined,"process/framework/RNA_activity/C1_C2_signac.Rds")



