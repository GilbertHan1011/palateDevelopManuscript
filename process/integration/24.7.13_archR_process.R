#rm(list=ls())
setwd("../2024_7.ATAC_E12/")
library(ArchR)
archrFile <- "data/E12.5_1/fragments.tsv.gz"
addArchRGenome("mm10")
addArchRThreads(threads = 16) 
ArrowFiles <- createArrowFiles(
  inputFiles = archrFile,
  sampleNames = "E12_1",
  filterTSS = 2, #Dont set this too high because you can always increase later
  filterFrags = 500, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
proj_sample1 <- ArchRProject(
  ArrowFiles = "process/framework/archr/ArrowFiles/E12_1.arrow", 
  outputDirectory = "process/framework/archr/",
  copyArrows = F #This is recommened so that you maintain an unaltered copy for later usage.
)


archrFile2 <- "data/E12.5_2/fragments.tsv.gz"
ArrowFiles <- createArrowFiles(
  inputFiles = archrFile2,
  sampleNames = "E12_2",
  filterTSS = 2, #Dont set this too high because you can always increase later
  filterFrags = 500, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
proj2 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "process/framework/archr_sample2/",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)



saveArchRProject(ArchRProj = proj_sample1, outputDirectory = "process/framework/archr/", load = FALSE)
saveArchRProject(ArchRProj = proj, outputDirectory = "process/framework/archr_sample2/", load = FALSE)
