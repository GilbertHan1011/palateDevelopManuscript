```R
setwd("~/lh/2024.4_scATAC/")

library(Seurat)

palate <- readRDS("data/palate_epi_MAGIC_periderm_MAGIC_harmony_renamed.rds")

#palatesc <- zellkonverter::readH5AD("processed_data/framework/4.29_scvi.h5ad")

label <- read.csv("processed_data/cluster/rna_annotation_scvi_mod.csv",row.names = 1)

  
  
  
  

palate@meta.data[colnames(label)] <- label

reduction <- read.csv("processed_data/embedding/RNA_scvi_umap.csv",row.names = 1)

colnames(reduction) <- c("umap_1","umap_2")

palate@reductions$scvi_umap = CreateDimReducObject(embeddings = as.matrix(reduction),key = "umap_", assay = DefaultAssay(palate))

FeaturePlot(palate,"Krt6a",reduction = "scvi_umap")

FeaturePlot(palate,"Krt6a",reduction = "umap")

DimPlot(palate,group.by = "level1_anno")

DimPlot(palate,group.by = "level2_anno",reduction = "scvi_umap")

DimPlot(palate,group.by = "level2_anno",reduction = "umap")

FeaturePlot(palate,"Krt6a")

DefaultAssay(palate) <- "RNA"

Idents(palate) <- palate$level1_anno

markersAll <- FindAllMarkers(palate,only.pos = T)

  
  

#== marker-----

  

markersTb <- markersAll%>%dplyr::filter(p_val_adj<0.01)

dir.create("processed_data/differential_gene")

dir.create("result/differential_gene")

diffGenes <- markersTb$gene%>%unique()

  

palate$level1_anno%>%unique

palate$level1_anno[palate$level1_anno=="Mature K5(-) "] <- "K5(-)"

palate$anno_sample <- paste0(palate$orig.ident,"_",palate$level1_anno)

  

ident_keep <- table(palate$anno_sample)%>%names()%>%.[table(palate$anno_sample)>30]

  

palate <- palate[,palate$anno_sample%in%ident_keep]

palatesce <- as.SingleCellExperiment(palate)

library(dreamlet)

pb <- aggregateToPseudoBulk(palatesce,

assay = "counts",

cluster_id = "cell",

sample_id = "anno_sample",

verbose = FALSE

)

  

pb <- aggregateData(palatesce,

assay = "counts", fun = "sum",

by = c("anno_sample"))

pbDf <- assay(pb)

library(edgeR)

pbDf <- cpm(pbDf)

  

pbDf <- pbDf[diffGenes,]

  

split_data <- strsplit(colnames(pbDf), "_")

  

# Extract the first two parts of each element

time <- vapply(split_data, function(x) paste(x[1], sep = "_"), character(1))

  

celltype <- vapply(split_data, function(x) paste(x[3], sep = "_"), character(1))

  

pbDf <- t(scale(t(pbDf)))

  

haAll = HeatmapAnnotation(

ident = factor(timeHa),

col = list(

ident=my_color

)

)

  
  
  

library(RColorBrewer)

my_color_time<-brewer.pal(2,'Set1')[1:2]

names(my_color_time) <- c("E125","E145")

my_color_celltype<-brewer.pal(length(unique(celltype)),'Set3')

names(my_color_celltype) <- unique(celltype)

  
  

haAll = HeatmapAnnotation(

time = factor(time,levels = c("E125","E145")),

celltype=factor(celltype),

col = list(

time=my_color_time,

celltype=my_color_celltype

)

)

  
  

hm <- Heatmap(pbDf,km = 8,show_row_names = F,show_row_dend = F,top_annotation=haAll,use_raster = F)

  

hmOrderAll <- row_order(hm)

hmGeneAll <- rownames(hm@matrix)

hmGeneListAll <- lapply(hmOrderAll, function(x) hmGeneAll[x])

library(stringi)

resAll <- as.data.frame((stri_list2matrix(hmGeneListAll)))

colnames(resAll) <- names(hmGeneListAll)

  

set.seed(124)

hm2 <- Heatmap(pbDf,km = 9,show_row_names = F,show_row_dend = F,top_annotation=haAll,use_raster = F)

hm2 <- draw(hm2)

hmOrderAll <- row_order(hm2)

hmGeneAll <- rownames(hm2@ht_list$matrix_15@matrix)

hmGeneListAll <- lapply(hmOrderAll, function(x) hmGeneAll[x])

library(stringi)

resAll <- as.data.frame((stri_list2matrix(hmGeneListAll)))

colnames(resAll) <- names(hmGeneListAll)

  

#7,2,6,1,4,5,8,9,3

  

geneOrdered <- c(hmOrderAll["7"][[1]],hmOrderAll["2"][[1]],hmOrderAll["6"][[1]],

hmOrderAll["1"][[1]],hmOrderAll["4"][[1]],hmOrderAll["5"][[1]],

hmOrderAll["8"][[1]],hmOrderAll["9"][[1]],hmOrderAll["3"][[1]])

  
  

rowSplit <- c(rep("E12.5", length(hmOrderAll["7"][[1]])),rep("E14.5", length(hmOrderAll["2"][[1]])),

rep("stem cell", length(hmOrderAll["6"][[1]])),

rep("K5(-)", length( hmOrderAll["1"][[1]])),rep("E12.5 Shh(+)", length(hmOrderAll["4"][[1]])),

rep("E14.5 Shh(+)", length(hmOrderAll["5"][[1]])),

rep("K6(+)", length(hmOrderAll["8"][[1]])),rep("K14(-)", length(hmOrderAll["9"][[1]])),

rep("matureK6(+)", length(hmOrderAll["3"][[1]])))

  

hm_3 <- Heatmap(pbDf[geneOrdered,],show_row_names = F,col=colorRamp2(c(-2, 0, 4), c("DeepSkyBlue3", "white", "red")),

show_row_dend = F,top_annotation=haAll,use_raster = F,row_split = rowSplit,

cluster_rows = F,cluster_row_slices = F)

  

pdf("result/differential_gene/hm3.pdf",width = 6,height = 10)

draw(hm_3)

dev.off()

  

#length(unique(markersTb$gene))

rowSplit <- factor(rowSplit,levels = unique(rowSplit))

Heatmap(pbDf[geneOrdered,],show_column_names = F,cluster_rows = F,

show_row_names = F,col=colorRamp2(c(-2, 0, 4), c("DeepSkyBlue3", "white", "red")),top_annotation = haAll,

row_split = rowSplit,,border=T,)

  

Heatmap(pbDf[geneOrdered,],show_column_names = F,cluster_rows = F,cluster_row_slices = F,

show_row_names = F,col=colorRamp2(c(-2, 0, 4), c("DeepSkyBlue3", "white", "red")),top_annotation = haAll,

row_split = rowSplit,,border=T)

  
  

genelist <- cbind(hmGeneAll[geneOrdered],as.character(rowSplit))

colnames(genelist) <- c("gene","cluster")

write.csv(genelist,"processed_data/5.24_de_gene_hm_km9.csv")

  
  

pheatmap::pheatmap(pbDf[geneOrdered,],cluster_rows = F,show_rownames = F)

  

write.csv(pbDf[geneOrdered,],"processed_data/differential_gene/5.30_hm_matrix.csv")
```