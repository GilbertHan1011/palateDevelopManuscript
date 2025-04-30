setwd("../2024_7.ATAC_E12/")
peak10xC1 <- read.table("data/E12.5_1/peaks.bed",col.names = c("chr", "start", "end"))
peak10xC2 <- read.table("data/E12.5_2/peaks.bed",col.names = c("chr", "start", "end"))
peak10xB1 <- fread("../2024.4_scATAC/data/raw/B1/atac_peak_annotation.tsv") %>% as.data.frame()
peak10xB2 <- fread("../2024.4_scATAC/data/raw/B2/peak_annotation.tsv") %>% as.data.frame()
grB1 <- makeGRangesFromDataFrame(peak10xB1)
grB2 <- makeGRangesFromDataFrame(peak10xB2)
grC1 <- makeGRangesFromDataFrame(peak10xC1)
grC2 <- makeGRangesFromDataFrame(peak10xC2)
combined.peaks <- reduce(x = c(grB1,grB2,grC1, grC2))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

cellsB1 <- read.csv("../2024.4_scATAC/processed_data/framework/filtered_cells/B1_cells_Epi",row.names = 1) %>% unlist()
cellsB2 <- read.csv("../2024.4_scATAC/processed_data/framework/filtered_cells/B2_cells_Epi",row.names = 1) %>% unlist()

frags.B1 <- CreateFragmentObject(
  path = "../2024.4_scATAC/data/raw/B1/atac_fragments.tsv.gz",
  cells = cellsB1
)
frags.B2 <- CreateFragmentObject(
  path = "../2024.4_scATAC/data/raw/B2/fragments.tsv.gz",
  cells = cellsB2
)

frags.C1 <- CreateFragmentObject(
  path = "data/E12.5_1/fragments.tsv.gz",
  cells = cellsC1
)
frags.C2 <- CreateFragmentObject(
  path = "data/E12.5_2//fragments.tsv.gz",
  cells = cellsC2
)


B1.counts <- FeatureMatrix(
  fragments = frags.B1,
  features = combined.peaks,
  cells = cellsB1
)

B2.counts <- FeatureMatrix(
  fragments = frags.B2,
  features = combined.peaks,
  cells = cellsB2
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
