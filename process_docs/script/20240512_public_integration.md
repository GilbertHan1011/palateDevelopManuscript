```R
rnaRef <- zellkonverter::readH5AD("processed_data/framework/4.29_scvi.h5ad")

rnaSeurat <- rnaRef%>%as.Seurat()

epiSubset <- DietSeurat(epiCells)

epiSubset <- epiSubset[intersect(rownames(epiSubset),rownames(rnaSeurat)),]

rnaSeurat@assays$RNA <- rnaSeurat@assays$originalexp

rnaSeurat@assays$originalexp <- NULL

#rnaSeurat <- rnaSeurat[intersect(rownames(epiSubset),rownames(rnaSeurat)),]

mergeObj <- merge(rnaSeurat,epiSubset)

mergeObj <- NormalizeData(mergeObj, normalization.method = "LogNormalize", scale.factor = 10000)

mergeObj <- ScaleData(mergeObj)

mergeObj <- RunPCA(mergeObj, features = rownames(mergeObj))

mergeObj <- RunHarmony(mergeObj, "orig.ident",reduction.use = "pca")

mergeObj <- FindNeighbors(mergeObj, dims = 1:30,reduction = "harmony")

mergeObj <- FindClusters(mergeObj, resolution = 0.5)

mergeObj <- FindClusters(mergeObj, resolution = 1)

mergeObj <- FindClusters(mergeObj, resolution = 1.5)

mergeObj <- RunUMAP(mergeObj, dims = 1:30,reduction = "harmony")

FeaturePlot(mergeObj,c("Krt6a","Krt14","Krt5"),ncol = 3)

DimPlot(mergeObj,group.by = "orig.ident")

DimPlot(mergeObj,group.by = "celltype")

DimPlot(mergeObj)

  
  

dayDict <- c("E121_h5" = "E12.5", "E122_h5" = "E12.5",

"E132_h5" = "E13.5", "E133_h5" = "E13.5",

"E1402_h5" = "E14.0", "E140_h5" = "E14.0",

"E1451_h5" = "E14.5", "E1452_h5" = "E14.5",

"E125_1" = "E12.5", "E125_2" = "E12.5",

"E145_1" = "E14.5", "E145_2" = "E14.5")

  

# Create the new column 'age'

mergeObj$age <- dayDict[match(mergeObj$orig.ident, names(dayDict))]

projDict <- c("E121_h5" = "public", "E122_h5" = "public",

"E132_h5" = "public", "E133_h5" = "public",

"E1402_h5" = "public", "E140_h5" = "public",

"E1451_h5" = "public", "E1452_h5" = "public",

"E125_1" = "lh", "E125_2" = "lh",

"E145_1" = "lh", "E145_2" = "lh")

mergeObj$proj <- projDict[match(mergeObj$orig.ident, names(projDict))]

FeaturePlot(mergeObj,c("Tgfb2","Shh","Maf","Sgk1","Krt6a","Gal"),ncol = 3)

DimPlot(mergeObj,group.by = "proj")

  

label = read.csv("processed_data/cluster/rna_annotation_scvi_mod.csv",row.names = 1)

mergeObj$label <- NaN

mergeObj$label[mergeObj$proj=="lh"] <- label$level1_anno

DimPlot(mergeObj,group.by = "label")

reduction_harmony1_run1 <- mergeObj@reductions

harmonyObj <- mergeObj[,mergeObj$proj=="public"]

FeaturePlot(harmonyObj,c("Tgfb2","Shh","Maf","Sgk1","Krt6a","Gal"),ncol = 3)

  
  

FeaturePlot(harmonyObj,c("Krt6a","Krt5","Krt14"),ncol = 3)

harmonyObj$label <- label$cell_type

DimPlot(harmonyObj,group.by = "label")

DimPlot(harmonyObj,group.by = "RNA_snn_res.1")

harmonyObj$RNA_snn_res.1 <- mergeObj$RNA_snn_res.1[mergeObj$proj=="public"]

  

newname <- c("K5(-)", "Transit", "K6+ cells", "Transit", "Mature K5(-)",

"K6+ cells", "Transit")

  

names(newname) <- c("K5(-)", "Transit", "K6+ cells", "stem cells", "Mature K5(-) ",

"K14(-)", "Shh(+)")

  

Idents(harmonyObj) <- harmonyObj$label

harmonyObj <- RenameIdents(harmonyObj,newname)

DimPlot(harmonyObj,group.by = "confi_label")

FeaturePlot(harmonyObj,"Krt6a")

harmonyObj$confi_label <- Idents(harmonyObj)

harmonyObj$confi_label[harmonyObj$RNA_snn_res.1%in%c("15","16")] = "K6+ cells"

  

epiCells$confi_label <- harmonyObj$confi_label

DimPlot(epiCells,group.by = "confi_label")

meta <- harmonyObj@meta.data[,c("age", "confi_label")]

grouped <- meta %>%

group_by(age, confi_label) %>%

count() %>%

ungroup()

  

# Create the bar plot

ggplot(grouped, aes(x = confi_label, y = freq, fill = age)) +

geom_bar(stat = "identity", position = "dodge") +

xlab("Cell Type") +

ylab("Count") +

ggtitle("Count of Cell Types by Age") +

theme_minimal()

epiCells@reductions$umap_harmony <- harmonyObj@reductions$umap

confi_label <- epiCells$confi_label%>%as.data.frame()

write.csv(confi_label,"processed_data/cluster/public_confident_label.csv")

reductionEpi <- epiCells@reductions

epiCells <- DietSeurat(epiCells)

epiCells@reductions <- reductionEpi

saveRDS(epiCells,"processed_data/framework/5.12_public_integration_scvi_harmony.Rds")

  

#== use seurat reference mapping -------

  

#

# rnaSeurat <- SCTransform(rnaSeurat, verbose = FALSE)

# epiSubset <- SCTransform(epiSubset, verbose = FALSE)

# rnaSeurat <- ScaleData(rnaSeurat)

# rnaSeurat <- RunPCA(rnaSeurat, assay = 'RNA',features = rownames(rnaSeurat))

# anchors <- FindTransferAnchors(

# reference = rnaSeurat,

# query = epiSubset,

# normalization.method = "SCT",

# reference.reduction = "PCA",

# dims = 1:50

# )

#

# anchors <- FindTransferAnchors(

# reference = reference,

# query = pbmc3k,

# normalization.method = "SCT",

# reference.reduction = "spca",

# dims = 1:50

# )
```