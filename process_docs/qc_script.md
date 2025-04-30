
## Meta
Date : 2024.7.16

Sample : B2

Filename : 24.7.16_B2_signac.R

## Input and output

### Input
- Raw data files:
  - `data/raw/B2/raw_peak_bc_matrix.h5`: Peak count matrix
  - `data/raw/B2/fragments.tsv.gz`: Fragment file
  - `data/raw/B2/singlecell.csv`: Cell metadata

### Output
- QC metrics and filtered cells:
  - `processed_data/QC/7.16_B2_qc_metric.csv`: QC metrics for all cells
  - `processed_data/cluster/7.16_B2_coarse_label.csv`: Cell type labels
  - `processed_data/markerGene/7.16_B2_cluster9_mesencyme.csv`: Marker genes for mesenchymal cluster
  - `processed_data/framework/B2_signac.Rds`: Final filtered Seurat object with epithelial cells only

```R
library(Signac)
library(GGally)


setwd("../2024.4_scATAC/")

counts_2 <- Read10X_h5(filename = "data/raw/B2/raw_peak_bc_matrix.h5")

#count_mat <- Read10X("../rawdata/24.4_ATAC/whfs-xs-231165_202404190900/whfs-xs-231165_delivery_20240419/B1/filtered_peak_bc_matrix/")

chrom_assay_2 <- CreateChromatinAssay(
  counts = counts_2,
  sep = c(":", "-"),
  min.cells = 10,
  min.features = 200,
  fragments = "data/raw/B2/fragments.tsv.gz"
)

C2 <- CreateSeuratObject(
  counts = chrom_assay_2,
  assay = "peaks",
)

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
  # countDF <- CountFragments(fragmentFile,cells = colnames(obj))
  # obj@meta.data[colnames(countDF)] <- countDF
  # obj <- FRiP(obj, "peaks", "frequency_count", col.name = "FRiP", verbose = TRUE)
  # qcMetrics <- obj@meta.data[c("nCount_peaks", "nFeature_peaks", "nucleosome_signal", 
  #                              "nucleosome_percentile", "TSS.enrichment", 
  #                              "FRiP")]
  # write.csv(qcMetrics,paste0(outDir,"/",obj@project.name,"_qc_metrics.csv"))
  return(obj)
} 

library(EnsDb.Mmusculus.v79)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
Annotation(B2) <- annotations

B2 <- qcATAC(B2,fragmentFile = "data/raw/B2/fragments.tsv.gz",outDir = "processed_data//QC/")
metadata <- read.csv(
  file = "data/raw/B2/singlecell.csv",
  header = TRUE,
  row.names = 1
)
B2 <- AddMetaData(B2,metadata = metadata)
B2$FRIP <- B2$peak_region_fragments / B2$passed_filters * 100
B2@meta.data$logCount <- log(B2@meta.data$nCount_peaks)
ggpairs(B2@meta.data, columns =c("logCount",  "nucleosome_signal","TSS.enrichment", 
                                                             "FRIP"))
B2_cells <- B2@meta.data %>% 
  dplyr::filter(nCount_peaks > 3000 &
           FRIP > 15 &
           nucleosome_signal < 4 &
           TSS.enrichment > 3)%>%
  rownames()
B2@meta.data$group <- "low"
B2@meta.data$group[colnames(B2)%in%B2_cells] <- "high"
#qcMetricsC2$logPeak <- log(qcMetricsC2$nCount_peaks)
ggpairs(B2@meta.data, mapping = aes(color = group), columns =c("logCount",  "nucleosome_signal","TSS.enrichment", 
                                                              "FRIP"))
B2_filtered <- B2[,B2_cells]

B2_filtered <- RunTFIDF(B2_filtered)
B2_filtered <- FindTopFeatures(B2_filtered, min.cutoff = 'q0')
B2_filtered <- RunSVD(B2_filtered)
DepthCor(B2_filtered)
B2_filtered <- RunUMAP(object = B2_filtered, reduction = 'lsi', dims = 2:30)
B2_filtered <- FindNeighbors(object = B2_filtered, reduction = 'lsi', dims = 2:30)
B2_filtered <- FindClusters(object = B2_filtered, verbose = FALSE, algorithm = 3)
DimPlot(object = B2_filtered, label = TRUE) + NoLegend()
gene.activities_B2 <- GeneActivity(B2_filtered)


B2_filtered['RNA']('RNA') <- CreateAssayObject(counts = gene.activities_B2)
B2_filtered <- NormalizeData(
  object = B2_filtered,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(B2_filtered$nCount_RNA)
)
newLabel <- c("Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Epi", "Mesenchyme", "Epi", "Epi"
)
names(newLabel) <- levels(B2_filtered)
B2_filtered <- RenameIdents(B2_filtered,newLabel)
B2_filtered$coarse_label <- Idents(B2_filtered)
write.csv(B2_filtered$coarse_label,"processed_data/cluster/7.16_B2_coarse_label.csv")
write.csv(B2_filtered@meta.data,"processed_data/QC/7.16_B2_qc_metric.csv")
DimPlot(B2_filtered)


FeaturePlot(B2_filtered,c("Krt14","Krt5","Krt10","Krt6a"))
FeaturePlot(B2_filtered,c("Krt6a","Krt40","Krtdap","Dsg1b"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
B2_marker_9 <- FindMarkers(B2_filtered,ident.1 = "9",assay = "RNA") # mesenchyme with Prrx1
write.csv(B2_marker_9,"processed_data/markerGene/7.16_B2_cluster9_mesencyme.csv")
FeaturePlot(B2_filtered,c("Krt14","Krt5","Krt10"))

B2_filtered <- B2_filtered[,B2_filtered$coarse_label=="Epi"]
saveRDS(B2_filtered,"processed_data/framework/B2_signac.Rds")
```
