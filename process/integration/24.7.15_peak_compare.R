peak10xC1 <- read.table("data/E12.5_1/peaks.bed",col.names = c("chr", "start", "end"))
peak10xC2 <- read.table("data/E12.5_2/peaks.bed",col.names = c("chr", "start", "end"))

grC1 <- makeGRangesFromDataFrame(peak10xC1)
grC2 <- makeGRangesFromDataFrame(peak10xC2)

peak10xB1 <- fread("./../2024.4_scATAC/data/raw/B1/atac_peak_annotation.tsv")
peak10xB2 <- fread("./../2024.4_scATAC/data/raw/B1/atac_peak_annotation.tsv")
grB1 <- makeGRangesFromDataFrame(peak10xB1)
grB2 <- makeGRangesFromDataFrame(peak10xB2)
dir.create("process/E12_E14_compare")
dir.create("process/E12_E14_compare/peaks")
rtracklayer::export(grC1,"process/E12_E14_compare/peaks/peak_E12_1.bed")
rtracklayer::export(grC2,"process/E12_E14_compare/peaks/peak_E12_2.bed")

rtracklayer::export(grB1,"process/E12_E14_compare/peaks/peak_E14_1.bed")
rtracklayer::export(grB2,"process/E12_E14_compare/peaks/peak_E14_2.bed")

opC1 <- findOverlaps(grC1,grC2)
