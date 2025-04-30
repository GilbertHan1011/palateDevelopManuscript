library(clusterProfiler)
library(tidyverse)

marker <- read.csv("descriptive_results/rna_marker/20240930_level1_anno.csv",row.names = 1)
marker <- marker %>% dplyr::filter(p_val_adj<0.001 & avg_log2FC > 2)
markerList <- split(marker$gene, marker$cluster)

res <- compareCluster(geneCluster = markerList, fun = enrichGO,OrgDb='org.Mm.eg.db',ont="BP",keyType="SYMBOL")
dotplot(res)

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
    ckRowList[[i]] <- ckRow
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
res_select <- make_cluster_comparison(res,goPerGroup = 10)
dotplot(res_select,showCategory = 15)


res2 <- compareCluster(geneCluster = markerList, fun = enrichGO,OrgDb='org.Mm.eg.db',ont="BP",keyType="SYMBOL")
res_select2 <- make_cluster_comparison(res2,goPerGroup = 10)
dotplot(res_select2,showCategory = 10)
ggsave("descriptive_results/rna_marker/20250116_marker_go.pdf",width = 8,height = 8)

View(res2@compareClusterResult)
GOStem <- enrichGO(gene         = markerList$`stem cells`,
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
dotplot(GOStem,showCategory = 15)
dotplot(GOStem,showCategory = 20)
ggsave("descriptive_results/rna_marker/20250116_marker_stem_go.pdf",height = 6,width = 6)
