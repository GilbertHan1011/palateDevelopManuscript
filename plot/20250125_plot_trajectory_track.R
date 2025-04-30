proj <- readRDS("process/ArchR_combined/Save-ArchR-Project.rds")
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cluster_knn_order", 
  geneSymbol = c("Krtdap"),
  #features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Krtdap)
plotPDF(p$Krtdap, name = "20250125_Krtdap_track.pdf", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
p2 <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "cluster_knn_order", 
  geneSymbol = c("Sprr1a","Krt16","Sprr1b","Sprr4","Klf2","Snai2","Vgll2"),
  #features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p2$Sprr1a)
plotPDF(p2, name = "20250125_trajectory_gene_track.pdf", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

#== add reduction to snapatac2--------------------------------------------

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

reductionUmapDf <- read.csv("process/framework/reduction/20250118_snap_harmony_drawgraph.csv",row.names = 1)
rowid <- rownames(reductionUmapDf)
rowid <- sapply(rowid, function(x) {substr(x,3,3) = "#";x})
rownames(reductionUmap) <- paste0("B1#",rownames(reductionUmap))
reductionUmap <- reductionUmap[proj$cellNames,]

prefix_map <- c(
  "C2_" = "E12_2#",
  "C1_"= "E12_1#",
  "B2_"=  "E14_2#",
  "B1_" = "B1#"
  # Add more mappings if needed
)


# Function to replace prefixes
replace_prefix <- function(x, prefix_map) {
  # Find which prefix matches
  for (old_prefix in names(prefix_map)) {
    if (startsWith(x, old_prefix)) {
      # Replace the old prefix with the new one
      return(sub(old_prefix, prefix_map[old_prefix], x))
    }
  }
  return(x)  # Return unchanged if no prefix matches
}

# Apply the replacement to all barcodes
new_barcodes <- sapply(rowid, replace_prefix, prefix_map = prefix_map)

rownames(reductionUmapDf) <- new_barcodes
colnames(reductionUmapDf) <- c("IterativeLSI#snap_UMAP_Dimension_1", "IterativeLSI#snap_UMAP_Dimension_2")
# First I copy its own UMAP, then I replace with the other reduction dataframe
# proj@embeddings$scglueUMAP <- proj@embeddings$
# proj@embeddings$scglueUMAP$df <- reductionUmap

# First I copy its own UMAP, then I replace with the other reduction dataframe
proj@embeddings$scglueUMAP <- proj@embeddings$UMAP
proj@embeddings$scglueUMAP$df <- reductionUmapDf

proj <- addImputeWeights(proj)

p_motif <- plotGroups(ArchRProj = proj, 
                groupBy = "cluster_knn_order", 
                colorBy = "MotifMatrix", 
                name = c("z:Rela_698"),
                imputeWeights = getImputeWeights(proj)
)



p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name =  c("z:Rela_698"), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)
