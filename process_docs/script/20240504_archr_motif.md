```R
library(Seurat)

library(ggmotif)

  

devtools::install_local("~/software/chromVARmotifs-master/")

#== change umap-----

reductionUmapDf <- read.csv("processed_data/embedding/scglue_umap_embedding.csv",row.names = 1)

#signacObj <- readRDS("processed_data/framework/B1_signac_filterd_integration.Rds")

cellname <- colnames(signacObj )

reductionUmap <- reductionUmapDf[cellname,]

rownames(reductionUmap) <- paste0("B1#",rownames(reductionUmap))

reductionUmap <- reductionUmap[proj$cellNames,]

colnames(reductionUmap) <- c("IterativeLSI#scglue_UMAP_Dimension_1", "IterativeLSI#scglue_UMAP_Dimension_2")

proj@embeddings$scglueUMAP <- proj@embeddings$UMAP

proj@embeddings$scglueUMAP$df <- reductionUmap

  

p1 <- plotEmbedding(

ArchRProj = proj,

colorBy = "MotifMatrix",

name = sort(markerMotifs),

embedding = "scglueUMAP",

imputeWeights = getImputeWeights(proj)

)

p2 <- plotEmbedding(

ArchRProj = proj,

colorBy = "cellColData",

name = "level1_anno",

embedding = "scglueUMAP",

imputeWeights = getImputeWeights(proj)

)

ggAlignPlots(p1, p2,p3, type = "h")

#ggsave("result/2_archR/Plots/5.5_chomvar_scglue_umap2.pdf")

plotPDF(p1,p2,p3, name = "5.5_chomvar_scglue_umap.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

#p1|p2

  

p3 <- plotEmbedding(

ArchRProj = proj,

colorBy = "cellColData",

name = "level2_anno",

embedding = "scglueUMAP",

imputeWeights = getImputeWeights(proj)

)

  
  
  
  
  

#== marker peak-----------------

  

markersPeaks <- getMarkerFeatures(

ArchRProj = proj,

useMatrix = "PeakMatrix",

groupBy = "level1_anno",

bias = c("TSSEnrichment", "log10(nFrags)"),

testMethod = "wilcoxon"

)

  

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

  

DARlist= list()

for (i in names(markerList)[1:4]){

markerdf <- markerList[[i]]%>%as.data.frame()

markerdf$cluster <- i

DARlist[[i]] <- markerdf

}

DARdf <- do.call(rbind,DARlist)

write.csv(DARdf,"processed_data/DAR/5.4_archR_DAR_0.05_0.5.csv")

  
  
  

#=== peak annotation---------

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

motifRes <- proj@peakAnnotation$Motif

motifRes$motifSummary

arid3aMotif <- motifRes$motifs$Arid3a_7

arid3aMotif

  

motifSummary <- proj@peakAnnotation$Motif$motifSummary%>%as.data.frame()

motif_1 <- proj@peakAnnotation$Motif$motifs

arid3aMotif <- motif_1$Arid3a_7

proj@peakAnnotation$Motif$Positions

positionFile <- readRDS("/home/zhanglab/lh/2024.4_scATAC/result/2_archR//Annotations/Motif-Positions-In-Peaks.rds")

matchFile <- readRDS("/home/zhanglab/lh/2024.4_scATAC/result/2_archR//Annotations/Motif-Matches-In-Peaks.rds")

  

mat <- matchFile@assays@data$matches

peakLogic <- mat[colnames(mat) =="Arid3a_7"]

  

arid3aGene <- matchFile@rowRanges$nearestGene[peakLogic]

  

gr <- matchFile@rowRanges[peakLogic]

library(rtracklayer)

export.bed(gr,con='processed_data/Arid3a_centric/archr_cisbp_predict_peak.bed')

table(gr@ranges@NAMES)

  

matchFile@metadata

  
  

arid3aPos <- positionFile$Arid3a_7

  

motifsUp <- peakAnnoEnrichment(

seMarker = markersPeaks,

ArchRProj = proj,

peakAnnotation = "Motif",

cutOff = "FDR <= 0.1 & Log2FC >= 0.5"

)

  
  

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])

df <- df[order(df$mlog10Padj, decreasing = TRUE),]

df$rank <- seq_len(nrow(df))

  
  

motifPositions <- getPositions(proj)

  

library(ggmotif)

mat <- TFBSTools::as.matrix(arid3aMotif)

probmat <- exp(mat) * matrix(TFBSTools::bg(arid3aMotif), nrow = nrow(mat),

ncol = ncol(mat), byrow = FALSE)

ggmotif::ggmotif_plot(probmat)

  

#=== Chromvar------------

proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(

ArchRProj = proj,

peakAnnotation = "Motif",

force = TRUE

)

  

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

  

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

  

motifs <- c("Arid3a")

markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")

markerMotifs

  

markerMotifs <- grep("z:", markerMotifs, value = TRUE)

markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]

markerMotifs

p <- plotGroups(ArchRProj = proj,

groupBy = "level1_anno",

colorBy = "MotifMatrix",

name = markerMotifs,

imputeWeights = getImputeWeights(proj)

)

plotPDF(p, name = "5.5_level1anno_motif_chromvar.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p <- plotGroups(ArchRProj = proj,

groupBy = "level2_anno",

colorBy = "MotifMatrix",

name = markerMotifs,

imputeWeights = getImputeWeights(proj)

)

plotPDF(p, name = "5.5_level1anno_motif_chromvar_group2.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

  
  
  

#

# p2 <- lapply(seq_along(p), function(x){

# if(x != 1){

# p[[x]] + guides(color = FALSE, fill = FALSE) +

# theme_ArchR(baseSize = 6) +

# theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +

# theme(

# axis.text.y=element_blank(),

# axis.ticks.y=element_blank(),

# axis.title.y=element_blank()

# ) + ylab("")

# }else{

# p[[x]] + guides(color = FALSE, fill = FALSE) +

# theme_ArchR(baseSize = 6) +

# theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +

# theme(

# axis.ticks.y=element_blank(),

# axis.title.y=element_blank()

# ) + ylab("")

# }

# })

# do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

  

p <- plotEmbedding(

ArchRProj = proj,

colorBy = "MotifMatrix",

name = sort(markerMotifs),

embedding = "scglue",

imputeWeights = getImputeWeights(proj)

)

  

p2 <- plotEmbedding(

ArchRProj = proj,

colorBy = "level1_anno",

name = sort(markerMotifs),

embedding = "UMAP",

imputeWeights = getImputeWeights(proj)

)

  

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "level1_anno", embedding = "UMAP")

  

#== footprint----------

motifPositions <- getPositions(proj)

  

motifs <- c("Arid3a")

markerMotifs_fp <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "level1_anno",force = T)

seFoot <- getFootprints(

ArchRProj = proj,

positions = motifPositions[markerMotifs_fp],

groupBy = "level1_anno"

)

  

p_fp1 <- plotFootprints(

seFoot = seFoot,

ArchRProj = proj,

normMethod = "Subtract",

plotName = "Footprints-Subtract-Bias",

addDOC = FALSE,

smoothWindow = 5,

plot = FALSE

)

  
  

p_fp2 <- plotFootprints(

seFoot = seFoot,

ArchRProj = proj,

normMethod = "Divide",

plotName = "Footprints-Divide-Bias",

addDOC = FALSE,

smoothWindow = 5,

plot = FALSE

)

  
  

seFoot2 <- getFootprints(

ArchRProj = proj,

positions = motifPositions["Tcfap2a_1"],

groupBy = "level1_anno"

)

plotFootprints(

seFoot = seFoot2,

ArchRProj = proj,

normMethod = "Subtract",

plotName = "Footprints-Subtract-Bias",

addDOC = FALSE,

smoothWindow = 5

)

test <- plotFootprints(

seFoot = seFoot2,

ArchRProj = proj,

normMethod = "Subtract",

plotName = "Footprints-Subtract-Bias",

addDOC = FALSE,

smoothWindow = 5,

plot = FALSE

)

test <- plotFootprints(

seFoot = seFoot2,

ArchRProj = proj,

normMethod = "Subtract",

plotName = "Footprints-Subtract-Bias",

addDOC = FALSE,

smoothWindow = 5,

plot = TRUE

)

  

grid.newpage() # Open a new plotting page

grid.draw(test$Tcfap2a_1)

  
  

saveArchRProject(ArchRProj = proj, outputDirectory = "processed_data/framework", load = FALSE)
```