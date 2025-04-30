# quality control---------------------------------
library(Signac)
library(Seurat)
library(plotly)
library(GGally)


counts_1 <- Read10X_h5(filename = "data/E12.5_1/raw_peak_bc_matrix.h5")

#count_mat <- Read10X("../rawdata/24.4_ATAC/whfs-xs-231165_202404190900/whfs-xs-231165_delivery_20240419/B1/filtered_peak_bc_matrix/")

chrom_assay_1 <- CreateChromatinAssay(
  counts = counts_1,
  sep = c(":", "-"),
  min.cells = 10,
  min.features = 200,
  fragments = "data/E12.5_1/fragments.tsv.gz"
)

C1 <- CreateSeuratObject(
  counts = chrom_assay_1,
  assay = "peaks",
)

C1[["peaks"]]
#granges(B1)
library(EnsDb.Mmusculus.v79)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(C1) <- annotations

C1 <- NucleosomeSignal(object = C1)
C1$nucleosome_group <- ifelse(C1$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = C1, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
ggsave("result/QC//fragment_signal_C1_hist_signac.pdf",width = 6,height = 4)
C1 <- TSSEnrichment(C1, fast = FALSE)
C1$high.tss <- ifelse(C1$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(C1, group.by = 'high.tss') + NoLegend()
ggsave("result/QC//C1__tssplot_signac.pdf",width = 5,height = 5)

# C1_meta <- read.csv("data/E12.5_1/",row.names = 1)
# B1@meta.data[colnames(B1_meta)] = B1_meta[colnames(B1),]
# 
# C1_meta$pct_reads_in_peaks <- C1_meta$? / C1_meta$passed_filters * 100
# C1_meta$blacklist_ratio <- C1_meta$blacklist_region_fragments / C1_meta$peak_region_fragments
countDF <- CountFragments("data/E12.5_1/fragments.tsv.gz",cells = colnames(C1))
C1@meta.data[colnames(countDF)] <- countDF
C1 <- FRiP(C1, "peaks", "frequency_count", col.name = "FRiP", verbose = TRUE)

DensityScatter(C1, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
VlnPlot(
  object = C1,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.01,
  ncol = 5,raster = T
)

qcMetrics <- C1@meta.data[c("nCount_peaks", "nFeature_peaks", "nucleosome_signal", 
                             "nucleosome_percentile", "TSS.enrichment", 
                              "FRiP")]
write.csv(qcMetrics,"process/QC/C1_qc_metrics.csv")

qcATAC <- function(obj,fragmentFile,outDir){
  obj <- NucleosomeSignal(object = obj)
  obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  #FragmentHistogram(object = obj, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
  #ggsave("result/QC//fragment_signal_C1_hist_signac.pdf",width = 6,height = 4)
  obj <- TSSEnrichment(obj, fast = FALSE)
  obj$high.tss <- ifelse(obj$TSS.enrichment > 2, 'High', 'Low')
  TSSPlot(obj, group.by = 'high.tss') + NoLegend()
  #ggsave("result/QC//C1__tssplot_signac.pdf",width = 5,height = 5)
  
  # C1_meta <- read.csv("data/E12.5_1/",row.names = 1)
  # B1@meta.data[colnames(B1_meta)] = B1_meta[colnames(B1),]
  # 
  # C1_meta$pct_reads_in_peaks <- C1_meta$? / C1_meta$passed_filters * 100
  # C1_meta$blacklist_ratio <- C1_meta$blacklist_region_fragments / C1_meta$peak_region_fragments
  countDF <- CountFragments(fragmentFile,cells = colnames(obj))
  obj@meta.data[colnames(countDF)] <- countDF
  obj <- FRiP(obj, "peaks", "frequency_count", col.name = "FRiP", verbose = TRUE)
  qcMetrics <- obj@meta.data[c("nCount_peaks", "nFeature_peaks", "nucleosome_signal", 
                              "nucleosome_percentile", "TSS.enrichment", 
                              "FRiP")]
  write.csv(qcMetrics,paste0(outDir,"/",obj@project.name,"_qc_metrics.csv"))
  return(obj)
} 




#load sample 2--------------------
counts_2 <- Read10X_h5(filename = "data/E12.5_2/raw_peak_bc_matrix.h5")

#count_mat <- Read10X("../rawdata/24.4_ATAC/whfs-xs-231165_202404190900/whfs-xs-231165_delivery_20240419/B1/filtered_peak_bc_matrix/")

chrom_assay_2 <- CreateChromatinAssay(
  counts = counts_2,
  sep = c(":", "-"),
  min.cells = 10,
  min.features = 200,
  fragments = "data/E12.5_2/fragments.tsv.gz"
)

C2 <- CreateSeuratObject(
  counts = chrom_assay_2,
  assay = "peaks"
)
C2@project.name <- "C2"
C1@project.name <- "C1"
Annotation(C2) <- annotations

C2 <- qcATAC(C2,fragmentFile = "data/E12.5_2/fragments.tsv.gz",outDir = "process/QC/")

Annotation(C2) <- annotations
DensityScatter(C2, x = 'nCount_peaks', y = 'FRiP', log_x = TRUE, quantiles = TRUE)




qcC1 <- C1@meta.data
qcC2 <- C2@meta.data
# 
# ggplot(longAnimateMeta,aes(x = pct_reads_in_peaks, y = n_fragment,color = value))+
#   geom_point(size= 0.01)+
#   geom_density_2d()+theme_bw()
# 


qcMetricsC1 <- C1@meta.data[c("nCount_peaks",  "nucleosome_signal","TSS.enrichment", 
                               "FRiP")]
qcMetricsC2 <- C2@meta.data[c("nCount_peaks",  "nucleosome_signal","TSS.enrichment", 
                              "FRiP")]
ggpairs(qcMetricsC2)

#== try to use Kmeans----------------------
km4 <- kmeans(qcMetricsC1, centers = 4, nstart = 25)
km4_2 <- kmeans(qcMetricsC2, centers = 4, nstart = 25)

qcMetricsC1$km_4 <- km4$cluster %>% as.character()
ggpairs(qcMetricsC1, mapping = aes(color = km_4), columns =c("nCount_peaks",  "nucleosome_signal","TSS.enrichment", 
                                                             "FRiP"))
qcMetricsC2$km_4 <- km4_2$cluster %>% as.character()
ggpairs(qcMetricsC2, mapping = aes(color = km_4), columns =c("nCount_peaks",  "nucleosome_signal","TSS.enrichment", 
                                                             "FRiP"))


#== use arbitary threshold--------------
VlnPlot(
  object = C1,
  features = c("nCount_peaks",  "nucleosome_signal","TSS.enrichment", 
               "FRiP"),
  pt.size = 0.1,
  ncol = 5
)
C1_cells <- C1@meta.data %>% 
  filter(nCount_peaks > 3000 &
           FRiP > 0.4 &
           nucleosome_signal < 4 &
           TSS.enrichment > 3)%>%
  rownames()

qcMetricsC1$group <- "low"
qcMetricsC1$group[colnames(C1)%in%C1_cells] <- "high"
qcMetricsC1$logPeak <- log(qcMetricsC1$nCount_peaks)
ggpairs(qcMetricsC1, mapping = aes(color = group), columns =c("logPeak",  "nucleosome_signal","TSS.enrichment", 
                                                             "FRiP"))
ggsave("result/QC/C1_metric_ggpair.pdf",width = 8,height = 8)

C2_cells <- C2@meta.data %>% 
  filter(nCount_peaks > 3000 &
           FRiP > 0.4 &
           nucleosome_signal < 4 &
           TSS.enrichment > 3)%>%
  rownames()

qcMetricsC2$group <- "low"
qcMetricsC2$group[colnames(C2)%in%C2_cells] <- "high"
qcMetricsC2$logPeak <- log(qcMetricsC2$nCount_peaks)
ggpairs(qcMetricsC2, mapping = aes(color = group), columns =c("logPeak",  "nucleosome_signal","TSS.enrichment", 
                                                              "FRiP"))
ggsave("result/QC/C2_metric_ggpair.pdf",width = 8,height = 8)

write.table(C2_cells,"process/framework/filtered_cells/C2_cells",quote = F)
write.table(C1_cells,"process/framework/filtered_cells/C1_cells",quote = F)

C1 <- C1[,C1_cells]
C2_filtered <- C2[,C2_cells]
Annotation(C2_filtered) <- annotations
saveRDS(C1,"process/framework/obj/C1_filtered_7.13.Rds")
saveRDS(C2_filtered,"process/framework/obj/C2_filtered_7.13.Rds")

#== doublets---------
library(scDblFinder)
repeats =  import('../2024.4_scATAC/data/peak_info/mm10_repeated.bed')
otherChroms <- GRanges(c("chrM","chrX","chrY","MT"),IRanges(1L,width=10^8)) # check which chromosome notation you are using c("M", "X", "Y", "MT")
toExclude <- suppressWarnings(c(repeats, otherChroms))
frag_path="data/E12.5_1/fragments.tsv.gz"
#future::plan("multisession", workers = 10) # do parallel
res1 <- amulet(frag_path, regionsToExclude=toExclude)
str(res)
write.csv(res1,"process/QC/7.14_amulet_B1.csv")

frag_path2="data/E12.5_2/fragments.tsv.gz"
#future::plan("multisession", workers = 10) # do parallel
res2 <- amulet(frag_path2, regionsToExclude=toExclude)
write.csv(res2,"process/QC/7.14_amulet_B2.csv")


C1$amulet <-  -log(res1[colnames(C1),]$p.value)
DefaultAssay(C1) <- "RNA"
FeaturePlot(C1,"Krt5")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
FeaturePlot(C1,c("Krt5","Krt14","Krt6a"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
FeaturePlot(C2_filtered,c("Krt5","Krt14","Krt6a"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
FeaturePlot(C2_filtered,c("Tie1","S100a14","Ly6d","Krt40","Krtdap","Zfp872","Sprr4","Klk7"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
FeaturePlot(C1,c("Tie1","S100a14","Ly6d","Krt40","Krtdap","Zfp872","Sprr4","Klk7"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))

#=========== scdbl finder------------
sce <- as.SingleCellExperiment(C1)
sce <- scDblFinder(sce,  aggregateFeatures=TRUE, nfeatures=25, 
                   processing="normFeatures")

dbl_score <- sce$scDblFinder.score
C1$scDblFinderScore <- dbl_score
FeaturePlot(C1,"scDblFinderScore")
B1_dbl <- as.data.frame(sce@colData[c("scDblFinder.class", 
                                      "scDblFinder.score", "scDblFinder.weighted", "scDblFinder.cxds_score")])
write.csv(B1_dbl,"process//QC/7.15_C1_filter_scDblFinder.csv")
B1_dbl$amulet_score <- C1$amulet 

ggplot(B1_dbl,aes(x = -log(1-scDblFinder.score), y = amulet_score))+
  geom_point()+
  theme_bw()



B1_dbl$scDblFinder.p <- 1-B1_dbl$scDblFinder.score
B1_dbl$amulet_p <- res1[rownames(B1_dbl),]$p.value
B1_dbl$combined <- apply(B1_dbl[,c("scDblFinder.p", "amulet_p")], 1, FUN=function(x){
  x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
  suppressWarnings(aggregation::fisher(x))
})

ggplot(B1_dbl,aes(x = amulet_score, y = log(combined)))+
  geom_point()+
  theme_bw()

ggplot(B1_dbl,aes(x = -log(1-scDblFinder.score), y = amulet_score))+
  geom_point()+
  theme_bw()

write.csv(B1_dbl,"process//QC/7.15_C1_filter_scDblFinder.csv")

C1$doublet <- "single"
C1$doublet[B1_dbl$combined<0.005] <- "doublets"


DimPlot(C1,group.by = "doublet")
ggsave("result/QC/C1_doublet_umap.pdf")

#DimPlot(C1)

#== signac reduction--------------
C1 <- RunTFIDF(C1)
C1 <- FindTopFeatures(C1, min.cutoff = 'q0')
C1 <- RunSVD(C1)
DepthCor(C1)
C1 <- RunUMAP(object = C1, reduction = 'lsi', dims = 2:30)
C1 <- FindNeighbors(object = C1, reduction = 'lsi', dims = 2:30)
C1 <- FindClusters(object = C1, verbose = FALSE, algorithm = 3)
DimPlot(object = C1, label = TRUE) + NoLegend()

VlnPlot(C1,"nCount_peaks")
FeaturePlot(C1,"nCount_peaks")

runReduction <- function(obj){
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = 'q0')
  obj <- RunSVD(obj)
  DepthCor(obj)
  obj <- RunUMAP(object = obj, reduction = 'lsi', dims = 2:30)
  obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30)
  obj <- FindClusters(object = obj, verbose = FALSE, algorithm = 3)
  DimPlot(object = obj, label = TRUE) + NoLegend()
  return(obj)
}

C2_filtered <- runReduction(C2_filtered)
DimPlot(object = C2_filtered, label = TRUE) + NoLegend()
gene.activities_C1 <- GeneActivity(C1)
C1[['RNA']] <- CreateAssayObject(counts = gene.activities_C1)
C1 <- NormalizeData(
  object = C1,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(C1$nCount_RNA)
)

marker_7 <- FindMarkers(C1,ident.1 = "7",assay = "RNA")
marker_19 <- FindMarkers(C1,ident.1 = "19",assay = "RNA")
marker_18 <- FindMarkers(C1,ident.1 = "18",assay = "RNA")

write.csv(marker_7,"process/annotation/C1_marker7.csv")  # mesenchyme
write.csv(marker_19,"process/annotation/C1_marker19.csv")
write.csv(marker_18,"process/annotation/C1_marker18.csv") # neuron

gene.activities_C2 <- GeneActivity(C2_filtered)

C2_filtered[['RNA']] <- CreateAssayObject(counts = gene.activities_C2)
C2_filtered <- NormalizeData(
  object = C2_filtered,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(C2_filtered$nCount_RNA)
)

DimPlot(object = C2_filtered, label = TRUE) + NoLegend()
C2_marker_17 <- FindMarkers(C2_filtered,ident.1 = "17",assay = "RNA") # blood
C2_marker_20 <- FindMarkers(C2_filtered,ident.1 = "20",assay = "RNA")
C2_marker_19 <- FindMarkers(C2_filtered,ident.1 = "19",assay = "RNA")
C2_marker_14 <- FindMarkers(C2_filtered,ident.1 = "14",assay = "RNA") # mesenchyme

saveRDS(C1,"process/framework/obj/C1_signac_processed.Rds")
saveRDS(C2_filtered,"process/framework/obj/C2_signac_processed.Rds")

write.csv(C2_marker_17,"process/annotation/C2_marker17.csv")  # blood
write.csv(C2_marker_20,"process/annotation/C2_marker20.csv")# mesenchyme
write.csv(C2_marker_19,"process/annotation/C2_marker19.csv") 
write.csv(C2_marker_14,"process/annotation/C2_marker14.csv") 

#== annotation--------------------------

FeaturePlot(C1,"Krt14",label = T)
newLabel <- c("Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "mesenchyme", "Epi", "Epi", "Epi", "Epi", 
              "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "neuron", "Msn+")
names(newLabel) <- levels(C1)
C1 <- RenameIdents(C1,newLabel)
C1$coarse_label <- Idents(C1)

write.csv(C1$coarse_label,"process/framework/cluster/C1_coarse.csv")

newlabel2 <- c("Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", 
               "Epi", "Epi", "mesenchyme", "Epi", "Epi", "blood", "Epi", "mesenchyme", "mesenchyme", "Epi")

names(newlabel2) <- levels(C2_filtered)
C2_filtered <- RenameIdents(C2_filtered,newlabel2)
C2_filtered$coarse_label <- Idents(C2_filtered)

write.csv(C2_filtered$coarse_label,"process/framework/cluster/C2_coarse.csv")


#FeaturePlot(C2_filtered,"Krt14",label = T)

C1_filter <- C1[,C1$coarse_label=="Epi"]

C2_filtered <- C2_filtered[,C2_filtered$coarse_label=="Epi"]


