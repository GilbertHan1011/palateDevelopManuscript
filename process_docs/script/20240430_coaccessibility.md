```R
proj <- readRDS("processed_data/framework/Save-ArchR-Project.rds")

proj@sampleColData$ArrowFiles <- "/home/zhanglab/lh/2024.4_scATAC/processed_data/framework/ArrowFiles/B1.arrow"

proj@projectMetadata$outputDirectory <- "/home/zhanglab/lh/2024.4_scATAC/result/2_archR/"

#== find peaks-----------------------

cluster <- read.csv("processed_data/cluster/B1_label_transfer_scvianno.csv",row.names = 1)

  

# # Define the patterns and replacements

# patterns <- c("Mature K5- cells ", "K5- cells", "K6+ cells", "Shh+ cells")

# replacements <- c("Mature K5(-)", "K5(-)", "K6(+)", "Shh(+)")

#

# # Substitute values

# for (i in seq_along(patterns)) {

# cluster$level1_anno <- gsub(patterns[i], replacements[i], cluster$level1_anno, fixed = TRUE)

# cluster$level2_anno <- gsub(patterns[i], replacements[i], cluster$level2_anno, fixed = TRUE)

# }

  
  

rownames(cluster) <- paste0("B1#",rownames(cluster))

cluster <- cluster[proj$cellNames,]

proj$level1_anno <- cluster$level1_anno

proj$level2_anno <- cluster$level2_anno

pathToMacs2 <- findMacs2()

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "level1_anno")

proj <- addReproduciblePeakSet(

ArchRProj = proj,

groupBy = "level1_anno",

pathToMacs2 = "/home/zhanglab/mambaforge/envs/py311/bin/macs2"

)

  

saveArchRProject(ArchRProj = proj, outputDirectory = "processed_data/framework", load = FALSE)

  

proj <- addPeakMatrix(proj)

  

scglueDimObj <- zellkonverter::readH5AD("processed_data/scglue/E14_step2_combine.h5ad")

scglueDim <- scglueDimObj@int_colData$reducedDims@listData$X_glue

scglueDim <- as.data.frame(scglueDim)

rownames(scglueDim) <- colnames(scglueDimObj)

scglueDimAtac <- scglueDim[length(proj$cellNames),]

write.csv(scglueDim,"processed_data/embedding/B1_scglue_combineRna.csv")

proj@reducedDims$scglue <- scglueDimAtac

  

names(proj@sampleColData$ArrowFiles) <- "B1"

  

getPeakSet(proj)

proj <- addCoAccessibility(

ArchRProj = proj,

reducedDims = "IterativeLSI"

)

  
  

cA <- getCoAccessibility(

ArchRProj = proj,

corCutOff = 0.5,

resolution = 1,

returnLoops = FALSE

)

metadata(cA)[[1]]

cA <- getCoAccessibility(

ArchRProj = projHeme5,

corCutOff = 0.5,

resolution = 1,

returnLoops = TRUE

)

  

markerGenes <- c(

"Ssh",

"Tgfb2", #Erythroid

"Krt6",

"Krt5", #Monocytes

"Krt10", "Krt14"

)

  

p <- plotBrowserTrack(

ArchRProj = proj,

groupBy = "level1_anno",

geneSymbol = markerGenes,

upstream = 50000,

downstream = 50000,

loops = getCoAccessibility(proj)

)

  

grid::grid.newpage()

grid::grid.draw(p$Tgfb2)

saveArchRProject(ArchRProj = proj, outputDirectory = "processed_data/framework", load = FALSE)
```