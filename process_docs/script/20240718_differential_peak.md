```R
#== peak_region--------------

library(ChIPseeker)

library(clusterProfiler)

peakDiff <- read.csv("process/differential_peak/7.18_peak_vi_cluster1_vs_3.csv",row.names = 1)

peakDiffFilter <- peakDiff %>% dplyr::filter(is_da_fdr=="True")

cluster1DA <- peakDiffFilter %>% dplyr::filter(effect_size<0)

cluster3DA <- peakDiffFilter %>% dplyr::filter(effect_size>0)

  

Annotation(combined4["ATAC"]("ATAC"))

  

annoATAC <- distanceToNearest(combined4, subject = Annotation(combined4["ATAC"]("ATAC")))

  

peakDiff1 <- cluster1DA %>% rownames()

peakDiff3 <- cluster3DA %>% rownames()

createRange <- function(ranges){

gr <- GRanges(

seqnames = sapply(strsplit(ranges, "-"), `[`, 1),

ranges = IRanges(

start = as.integer(sapply(strsplit(ranges, "-"), function(x) as.numeric(x[2]))),

end = as.integer(sapply(strsplit(ranges, "-"), function(x) as.numeric(x[3])))

)

)

return(gr)

}

peakDiff1Gr <- createRange(peakDiff1)

peakDiff3Gr <- createRange(peakDiff3)

  

#== chipseeker pipeline-----------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakDiff1GrAnno <- annotatePeak(peakDiff1Gr, tssRegion=c(-3000, 3000),

TxDb=txdb, annoDb="org.Mm.eg.db")

peakDiff3GrAnno <- annotatePeak(peakDiff3Gr, tssRegion=c(-3000, 3000),

TxDb=txdb, annoDb="org.Mm.eg.db")

  

library(EnsDb.Mmusculus.v79)

edb <- EnsDb.Mmusculus.v79

seqlevelsStyle(edb) <- "UCSC"

  

peakAnno1.edb <- annotatePeak(peak = peakDiff1Gr, tssRegion=c(-3000, 3000),

TxDb=edb, annoDb="org.Mm.eg.db")

peakAnno3.edb <- annotatePeak(peakDiff3Gr, tssRegion=c(-3000, 3000),

TxDb=edb, annoDb="org.Mm.eg.db")

  

plotAnnoPie(peakDiff1GrAnno)

plotAnnoPie(peakDiff3GrAnno)

  

peakAnno <- list(peakDiff1GrAnno,peakDiff3GrAnno)

names(peakAnno) <- c("Differnetiation","Stem")

plotAnnoBar(peakAnno)

dir.create("result/differential_peak")

ggsave("result/differential_peak/diff_peak_region.pdf")

  
  

allPeak <- combined4@assays$ATAC@ranges

peakAllGrAnno <- annotatePeak(allPeak, tssRegion=c(-3000, 3000),

TxDb=txdb, annoDb="org.Mm.eg.db")

plotAnnoBar(peakAllGrAnno)

peakAnno <- list(peakDiff1GrAnno,peakDiff3GrAnno,peakAllGrAnno)

names(peakAnno) <- c("Differnetiation","Stem","All peak")

plotAnnoBar(peakAnno)

ggsave("result/differential_peak/diff_peak_region_all.pdf")

  

annoDf <- as.data.frame(peakDiff1GrAnno)

annoDf3 <- as.data.frame(peakDiff3GrAnno)

write.csv(annoDf,"process/differential_peak/cluster1_anno.csv")

write.csv(annoDf3,"process/differential_peak/cluster3_anno.csv")

library(ReactomePA)

  

annoGene1 <- annoDf$SYMBOL %>% unique

annoGene3 <- annoDf3$SYMBOL %>% unique

  

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)

head(pathway1, 2)

  

library(clusterProfiler)

geneList <- list(annoGene1,annoGene3)

names(geneList) <- c("Differentiation","Stem")

ck <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb='org.Mm.eg.db',ont="ALL",keyType="SYMBOL")

dotplot(ck,showCategory=30)

  

View(ck@compareClusterResult,)

goCompare <- ck

make_cluster_comparison <- function(goCompare,goPerGroup = 5) {

ckRes <- goCompare@compareClusterResult

ckResSub <- ckRes[c("Description", "Cluster", "p.adjust")]

ckWide <- pivot_wider(ckResSub, names_from = Cluster, values_from = p.adjust)

ckWide <- ckWide %>% column_to_rownames("Description")

ckWideMod <- ckWide %>%

log10() * -1

ckWideMod[is.na(ckWideMod)] <- 0

ckRowList <- list()

for (i in 1:dim(ckWideMod)[2]) {

ckRow <- ckWideMod[i] - rowMeans(ckWideMod[-i])

ckRowList[i](i) <- ckRow

}

ckRowDf <- data.frame(do.call(cbind, ckRowList))

ckSelect <- c()

for (i in 1:dim(ckWideMod)[2]) {

ckSelect <- c(ckSelect, ckRowDf %>% dplyr::arrange(desc(across(i))) %>% rownames %>% .[1:goPerGroup])

}

ckSelect <- ckSelect %>% unique

ckResSelect <- ckRes %>% dplyr::filter(Description %in% ckSelect)

ckResSelect$p.adjust <- log10(ckResSelect$p.adjust)

ckResSelect$Description <- factor(ckResSelect$Description, levels = ckSelect)

ckResSelect <- dplyr::arrange(ckResSelect, Description)

goCompare@compareClusterResult <- ckResSelect

return(goCompare)

}

  

ck_2 <- make_cluster_comparison(ck,goPerGroup = 10)

dotplot(ck_2,showCategory=40)

ggsave("result/differential_peak/7.18_differential_GO.pdf",width = 6,height = 6)
```