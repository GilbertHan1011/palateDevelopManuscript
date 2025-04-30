# ArchR run 2-------------

clusterInfo <- read.csv("process/framework/cluster/20241211_knn_label.csv",row.names = 1) 

cellName <- rownames(clusterInfo)

cellName_new <- gsub("_", "#", cellName) %>% gsub("B2","E14_2",.) %>% 
  gsub("C1","E12_1",.) %>% 
  gsub("C2","E12_2",.)

#proj = proj[intersect(cellName_new,proj$cellNames)]

clusterInfo <- clusterInfo[cellName_new%in%proj$cellNames,]
#cluster <- clusterInfo$rna_cluster
names(clusterInfo) <- cellName_new[cellName_new%in%proj$cellNames]


# add cluster metadata
proj$cluster_knn <- clusterInfo

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "cluster_knn",force = T)


## Gene Score

markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "cluster_knn",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
dflist = list()
for (i in names(markerList)){
  df = markerList[[i]] %>% as.data.frame()
  df[["cluster"]] = i
  dflist[[i]] = df
}
Res = do.call(rbind,dflist)
write.csv(Res,"descriptive_results/atac_genescore/20241211_archR_genescore.csv")
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cluster_knn", 
  geneSymbol = "Mir203", 
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Mir203)

saveArchRProject(ArchRProj = proj)

proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = proj,
  corCutOff = 0.2,
  resolution = 1,
  returnLoops = FALSE
)

peakName <- paste0(proj@peakSet)

caTable <- data.frame( Peak1 = peakName[cA$queryHits], Peak2 = peakName[cA$subjectHits], coaccess = cA$correlation  )
write.csv(caTable,"process/regulation/coaccessibility/20241213_archR_ca_celloracle_prepare.csv")
write.csv(x = peakName, file = "process/regulation/coaccessibility//celloracle_prepare_all_peaks.csv")

cAloop <- getCoAccessibility(
  ArchRProj = proj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = T
)
cA
saveRDS(cA,"process/regulation/coaccessibility/20241213_archR_coaccessibiliy.Rds")
