library(scDblFinder)
library(rtracklayer)
library(GenomicRanges)
library(Signac)
setwd("/home/zhanglab/lh/2024_7.ATAC_E12/")
repeats =  import('../2024.4_scATAC/data/peak_info/mm10_repeated.bed')
otherChroms <- GRanges(c("chrM","chrX","chrY","MT"),IRanges(1L,width=10^8)) # check which chromosome notation you are using c("M", "X", "Y", "MT")
toExclude <- suppressWarnings(c(repeats, otherChroms))
frag_path2="data/E12.5_2/fragments.tsv.gz"
#future::plan("multisession", workers = 10) # do parallel
res2 <- amulet(frag_path2, regionsToExclude=toExclude)
write.csv(res2,"process/QC/7.14_amulet_B2.csv")
