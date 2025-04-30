library(clusterProfiler)

originalArid3a <- read.table("process/20241219_scenicplus_outs/20241222_regulon/Arid3a_regulon_pos-neg.csv")
newArid3a <- read.table("process/20241226_scenicplus_outs/Arid3a_regulon/arid3a_regulon.tsv")
newGene <- newArid3a$V2[newArid3a$V7 == "Arid3a"] %>% unique
oldGene <- originalArid3a$V2[originalArid3a$V7 == "Arid3a"] %>% unique
intersect(newGene,oldGene) %>% length()

geneList <- list("unpaired" = oldGene, "paired" = newGene)

ck <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb='org.Mm.eg.db',ont="BP",keyType="SYMBOL")
dotplot(ck)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
ggsave("result/regulation/20241226_scenicplus_paired/20241228_paired_GO.pdf")

ck2 <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb='org.Mm.eg.db',ont="ALL",keyType="SYMBOL")
dotplot(ck2,showCategory = 20)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
ggsave("result/regulation/20241226_scenicplus_paired/20241228_paired_GO2.pdf",width = 6,height = 10)

View(ck2@compareClusterResult)
