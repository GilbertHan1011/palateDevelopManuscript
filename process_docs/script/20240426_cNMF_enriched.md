# Gene module enrichment

```R
cnmfGene <- read.csv("processed_data/geneModule/cNMF/424_cNMF/geneModule.csv",row.names = 1)

cnmfGeneList <- list()

  

for (col in seq_along(colnames(cnmfGene))){

cnmfGeneList[col](col) <- cnmfGene[col]%>%sort(decreasing = T)%>%

head(150)%>%rownames()

}

  

cnmfGeneListSorted <- list()

  

for (col in seq_along(colnames(cnmfGene))){

cnmfGeneListSorted[col](col) <- cnmfGene[col]%>%sort(decreasing = T)%>%

rownames()

}

  

names(cnmfGeneList) <- paste0("cNMF",1:9)

  

library(aPEAR)

  

library(clusterProfiler)

library(org.Mm.eg.db)

library(DOSE)

data(geneList)

  

GOCluster1 <- enrichGO(gene = cnmfGeneList$cNMF1,

OrgDb = org.Mm.eg.db,

keyType = 'SYMBOL',

ont = "ALL",

pAdjustMethod = "BH",

pvalueCutoff = 0.01,

qvalueCutoff = 0.05)

  
  

# enrich <- gseGO(geneList, OrgDb = org.Hs.eg.db, ont = 'CC')

p <- enrichmentNetwork(GOCluster1@result, drawEllipses = TRUE, fontSize = 2.5)

ggsave("result/geneModule/enrich/4.26_cluster1_GOnet.pdf")

  
  

library(simplifyEnrichment)

goid <- GOCluster1@result$ID

mat = GO_similarity(goid,ont = "BP")

hmCluster1 <- simplifyGO(mat)

  

goEnrichList <- lapply(cnmfGeneList, enrichGO,

OrgDb = org.Mm.eg.db,

keyType = 'SYMBOL',

ont = "ALL",

pAdjustMethod = "BH",

pvalueCutoff = 0.01,

qvalueCutoff = 0.05)

names(goEnrichList) <- names(cnmfGeneList)

goTermEnrichList <- lapply(goEnrichList, function(x) x@result)

  

dotplot(goEnrichList[7](7))

  

pdf("result/geneModule/enrich/4.26_mutliplelist_GO.pdf",width = 10,height = 6)

multiGO <- simplifyGOFromMultipleLists(goTermEnrichList, padj_cutoff = 0.001,ont="BP")

dev.off()

  
  

ck <- compareCluster(geneCluster = cnmfGeneList, fun = enrichGO,OrgDb='org.Mm.eg.db',ont="BP",keyType = 'SYMBOL')

dotplot(ck)

ggsave("result/geneModule/enrich/compare_cluster_raw.pdf",width = 8,height = 8)
```