```R
library(slingshot)

library(destiny)

library(Seurat)

library(lattice)

  

coembed = zellkonverter::readH5AD("processed_data/framework/5.10_rna_atac_14.5_combine_drawgraph.h5ad")

find_sigmas(coembed@int_colData$reducedDims$X_glue)

dm <- run_diffMap(t(coembed@int_colData$reducedDims$X_glue), colData(coembed)$level2_anno, sigma = 2)

plot_eigenVal(dm = dm)

  

coembedSeurat <- as.Seurat(coembed,counts = "X",data = "X")

  

coembedSeurat@reductions$dm = CreateDimReducObject(embeddings = as.matrix(dm@eigenvectors),key = "DM_", assay = DefaultAssay(coembedSeurat))

DimPlot(coembedSeurat,reduction = "dm",group.by = "level2_anno")

  

#plot_eigenVal(dm = fibro_dm_int)

mycolor <- colorRampPalette(brewer.pal(7,'Set1'))(length(unique(colData(coembed)$level2_anno)))

names(mycolor) <- unique(colData(coembed)$level2_anno)

splom(~dm@eigenvectors[, c(1:6)], groups = colData(coembed)$level2_anno, col = mycolor, main = "Lineage",

key = list(space="right", points = list(pch = 19, col = mycolor), text = list(names(mycolor))))

  

plot(dm@eigenvectors[, c(1,5)], col = mycolor[as.character(colData(coembed)$level2_anno)],

pch=16,

asp = 1,

cex = .5)

legend("left",legend = names(mycolor),

fill = mycolor)

  

two_lineage <- slingshot(dm@eigenvectors[, c(1,2,5)],

clusterLabels = colData(coembed)$level2_anno,

start.clus = 'Transit 1',end.clus=c("Mature K5(-) ","K6(+)1"),allow.breaks=FALSE,

maxit = 1000, shrink.method = "density", thresh = 0.001, extend = "n")

rownames(dm@eigenvectors) <- colnames(coembed)

colData(coembed)$level2_anno <- factor(colData(coembed)$level2_anno,levels = unique(colData(coembed)$level2_anno))

dmVal <- dm@eigenvectors[, c(1,5)]

dmVal <- dmVal[colData(coembed)$level2_anno%in%c("Shh(+)", "K5(-)", "Transit 2", "Transit 3",

"Transit 1", "Transit 5", "K6(+)2", "K6(+)1", "Mature K5(-) "),]

colMeta <- colData(coembed)

colMeta <- colMeta[colData(coembed)$level2_anno%in%c("Shh(+)", "K5(-)", "Transit 2", "Transit 3",

"Transit 1", "Transit 5", "K6(+)2", "K6(+)1", "Mature K5(-) "),]

  

# 1,2,5 dim

two_lineage <- slingshot(dmVal,

clusterLabels = colMeta$level2_anno,

start.clus = 'Transit 1',end.clus=c("Mature K5(-) ","K6(+)1"),allow.breaks=FALSE,

maxit = 1000, shrink.method = "density", thresh = 0.001, extend = "n")

# 1,2 dim

two_lineage2 <- slingshot(dmVal,

clusterLabels = colMeta$level2_anno,

start.clus = 'Transit 1',end.clus=c("Mature K5(-) ","K6(+)1"),allow.breaks=FALSE,

maxit = 15, shrink.method = "density", thresh = 0.001, extend = "n")

#1,5 dim

two_lineage2 <- slingshot(dmVal,

clusterLabels = colMeta$level2_anno,

start.clus = 'Transit 1',end.clus=c("Mature K5(-) ","K6(+)1"),allow.breaks=FALSE,

maxit = 15, shrink.method = "density", thresh = 0.001, extend = "n")

  

pseudotime <- two_lineage2@assays@data$pseudotime

  

saveRDS(dm,"processed_data/trajectory/5.10_diffusionmap.Rds")

saveRDS(two_lineage,"processed_data/trajectory/5.10_slingshot_lineage.Rds")

saveRDS(two_lineage2,"processed_data/trajectory/5.10_slingshot_lineage2.Rds")

  

two_lineage@metadata$lineages

  

pseudotime <- two_lineage@assays@data$pseudotime

slingshotKrt6<- pseudotime[,4]

slingshotKrt6 <- slingshotKrt6[colnames(scviObjSeurat)]

scviObjSeurat$k6_pseudotime_slingthot<- slingshotKrt6

FeaturePlot(scviObjSeurat,features = "k6_pseudotime_slingthot",reduction = "scglue_umap")+

scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
```