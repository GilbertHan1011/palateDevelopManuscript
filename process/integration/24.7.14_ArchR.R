BiocManager::install(c("Seurat","ArchR"))
setwd("lh/2024_7.ATAC_E12/")
proj_C1 <- readRDS("process/framework/archr/Save-ArchR-Project.rds")
setwd
proj_C1@sampleColData$ArrowFiles <- "~/lh/2024_7.ATAC_E12/E12_1.arrow"
proj_C1@projectMetadata$outputDirectory <- "~/lh/2024_7.ATAC_E12/process/framework/archr/"
@C1_cells <- read.table("process/framework/archr/")
C1_cells <- read.table("process/framework/filtered_cells/C1_cells") %>% unlist
C1_cells <- paste0("E12_1#",C1_cells)

cellCol <- proj_C1@cellColData %>% as
proj_C1$quality <- "low"
proj_C1$quality[proj_C1$cellNames %in% C1_cells] <- "high"

#proj_C1[proj_C1]
proj_C1 <- addIterativeLSI(ArchRProj = proj_C1, useMatrix = "TileMatrix", name = "IterativeLSI")
proj_C1 <- addClusters(input = proj_C1, reducedDims = "IterativeLSI")

proj_C1 <- addUMAP(ArchRProj = proj_C1, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj_C1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_C1, colorBy = "cellColData", name = "quality", embedding = "UMAP")
ggsave(plot = p2,filename = "result/QC/C1_qc_umap.pdf")

plotPDF(plotList = p2, 
        name = "Plot-UMAP-low_quality.pdf", 
        ArchRProj = proj_C1, 
        addDOC = FALSE, width = 5, height = 5)

proj_C1 <- proj_C1[proj_C1$cellNames%in%C1_cells,]

proj_C1 <- addIterativeLSI(ArchRProj = proj_C1, useMatrix = "TileMatrix", name = "IterativeLSI_filtered",force = T)
proj_C1 <- addClusters(input = proj_C1, reducedDims = "IterativeLSI_filtered",name = "cluster2_filtered",force = T)

proj_C1 <- addUMAP(ArchRProj = proj_C1, reducedDims = "IterativeLSI_filtered",name = "UMAP_filtered",force = T)

plotEmbedding(ArchRProj = proj_C1, colorBy = "cellColData", name = "cluster2_filtered", embedding = "UMAP_filtered",force = T)
plotEmbedding(ArchRProj = proj_C1, colorBy = "cellColData", name = "cluster2_filtered", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_C1, colorBy = "cellColData", name = "quality", embedding = "UMAP")
ggsave(plot = p2,filename = "result/QC/C1_qc_umap.pdf")

proj_C2 <- readRDS("process/framework/archr_sample2/Save-ArchR-Project.rds")
proj_C2@sampleColData$ArrowFiles <- "~/lh/2024_7.ATAC_E12/E12_2.arrow"
#proj_C2@sampleColData$ArrowFiles <- "~/lh/2024_7.ATAC_E12/E12_1.arrow"
proj_C2@projectMetadata$outputDirectory <- "~/lh/2024_7.ATAC_E12/process/framework/archr_sample2/"
C2_cells <- read.table("process/framework/filtered_cells/C2_cells") %>% unlist

C2_cells <- paste0("E12_2#",C2_cells)
proj_C2 <- proj_C2[proj_C2$cellNames%in%C2_cells,]

proj_C2 <- addIterativeLSI(ArchRProj = proj_C2, useMatrix = "TileMatrix", name = "IterativeLSI")
proj_C2 <- addClusters(input = proj_C2, reducedDims = "IterativeLSI")

proj_C2 <- addUMAP(ArchRProj = proj_C2, reducedDims = "IterativeLSI")
plotEmbedding(ArchRProj = proj_C2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")



plotEmbedding(
  ArchRProj = proj_C2, 
  colorBy = "GeneScoreMatrix", 
  name = "Krt14", 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
plotEmbedding(
  ArchRProj = proj_C2, 
  colorBy = "GeneScoreMatrix", 
  name = "Krt14", 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotEmbedding(
  ArchRProj = proj_C1, 
  colorBy = "GeneScoreMatrix", 
  name = "Krt14", 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
plotEmbedding(
  ArchRProj = proj_C1, 
  colorBy = "GeneScoreMatrix", 
  name = "Krtdap", 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
plotEmbedding(
  ArchRProj = proj_C2, 
  colorBy = "GeneScoreMatrix", 
  name = "Krtdap", 
  embedding = "UMAP",size = 1,plotAs = "point",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
pList2 <- plotEmbedding(
  ArchRProj = proj_C2, 
  colorBy = "GeneScoreMatrix", 
  name = c("Tie1","S100a14","Ly6d","Krt40","Krtdap","Zfp872","Sprr4","Klk7","Dsg1b"), 
  embedding = "UMAP",size = 1,plotAs = "point",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
pList1 <- plotEmbedding(
  ArchRProj = proj_C1, 
  colorBy = "GeneScoreMatrix", 
  name = c("Tie1","S100a14","Ly6d","Krt40","Krtdap","Zfp872","Sprr4","Klk7","Dsg1b"), 
  embedding = "UMAP",size = 2,plotAs = "point",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
