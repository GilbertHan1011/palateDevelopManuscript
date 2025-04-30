library(Seurat)
geneScoreMat <- readRDS("process/framework/genescore/20250117_genescoreMat.Rds")
geneScoreSeurat <- CreateSeuratObject(geneScoreMat@assays@data@listData$GeneScoreMatrix)
rownames(geneScoreSeurat) <- geneScoreMat@elementMetadata$name
geneScoreSeurat <- FindVariableFeatures(geneScoreSeurat,nfeatures = 5000)
#readRDS()


count_score <- as.matrix(geneScoreSeurat@assays$RNA@layers$counts)
count_score <- count_score[,filterCell]
rownames(count_score) <- rownames(geneScoreSeurat)


gam_score <- fitGAM(counts = count_score,
                    pseudotime =pseudotimeCol,
                    cellWeights = cellsCol,
                    nknots = 5,
                    genes = VariableFeatures(geneScoreSeurat),
                    parallel=TRUE,
                    BPPARAM=BPPRARM,
                    verbose = T)
saveRDS(gam_score,"process/trajectory/20250123_genescore_gam.Rds")
rowData(gam_score)$assocRes <- associationTest(gam_score, lineages = TRUE, l2fc = 1)
assocRes <- rowData(gam_score)$assocRes
#View(assocRes)
line1Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05)
]
line2Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.05)
]
# line3Genes <-  rownames(assocRes)[
#   which(p.adjust(assocRes$pvalue_3, "fdr") <= 0.01)
# ]

UpSetR::upset(UpSetR::fromList(list(stem = line1Genes, K14 = line2Genes)))



yhatSmooth <- lapply(list(line1Genes,line2Genes),function(x)
  predictSmooth(gam_score, gene = x, nPoints = 100, tidy = FALSE))

geneUnion <- unique(c(line1Genes,line2Genes))
yhatSmoothAll <- predictSmooth(gam_score, gene = geneUnion, nPoints = 100, tidy = FALSE)
yhatSmoothAllTidy <- predictSmooth(gam_score, gene = geneUnion, nPoints = 500, tidy = TRUE)


#time_grid <- seq(max(yhatSmoothAllTidy$time), min(yhatSmoothAllTidy$time), length.out = 100)
library(pracma)
# Function to interpolate and average yhat values
time_grid <- seq(min(yhatSmoothAllTidy$time), max(yhatSmoothAllTidy$time)-0.04, length.out = 100)


# Process the data
result <- yhatSmoothAllTidy %>%
  group_by(lineage, gene) %>%
  group_modify(~ data.frame(time = time_grid,
                            yhat = interpolate_with_nan(., time_grid))) %>%
  ungroup()

# Create the 3D array
genes <- unique(yhatSmoothAllTidy$gene)
lineages <- unique(yhatSmoothAllTidy$lineage)

array_3d <- array(NaN, dim = c(100, length(genes), length(lineages)))

for (i in 1:length(lineages)) {
  for (j in 1:length(genes)) {
    subset_data <- result %>%
      filter(lineage == lineages[i], gene == genes[j])
    array_3d[, j, i] <- subset_data$yhat
  }
}

# Name the dimensions
dimnames(array_3d) <- list(
  time = time_grid,
  gene = genes,
  lineage = lineages
)

mat1 <- array_3d[,,1] %>% t
mat2 <- array_3d[,,2] %>% t
#mat3 <- array_3d[,,3] %>% t
mat <- cbind(mat1,mat2)
mat %>% head
str(mat)



# Apply the custom scaling function to mat
scaled_mat <- scale_rows_preserve_nan(mat)


# Check for infinite values
sum(is.infinite(mat))

# Check for NaN values
sum(is.nan(mat))


set.seed(1234)
hmAll<- Heatmap(scaled_mat,cluster_columns = F,show_column_names = F,
                show_row_names = F,
                col=colorRamp2(c(-2, 0, 4), c("blue", "white", "red")),na_col = "grey")

hmAll<- draw(hmAll)

meta <- motifSeurat@meta.data[c("level1_anno","group")]

lineage1Cell <- rownames(cellsCol)[cellsCol$stem==1]
lineage2Cell <- rownames(cellsCol)[cellsCol$K14==1]

meta1 <- meta[lineage1Cell,]
meta1$pseudo <- pseudotimeCol[lineage1Cell,]$stem
meta2 <- meta[lineage2Cell,]
meta2$pseudo <- pseudotimeCol[lineage2Cell,]$K14
# meta3 <- meta[lineage3Cell,]
# meta3$pseudo <- pseudotimeCol[lineage3Cell,]$K14

# First, let's clean the data
meta1_clean <- meta1[!is.infinite(meta1$pseudo) & !is.na(meta1$pseudo) & !is.na(meta1$group), ]

pseudotime1
# Create a boxplot (recommended for group comparisons)
boxplot(pseudotime1$Pseudotime_Lineage_K6 ~ group, data = motifSeurat@meta.data,
        main = "Pseudo Values by Group",
        xlab = "Group",
        ylab = "Pseudo Values")
boxplot(pseudo ~ group, data = motifSeurat@meta.data,
        main = "Pseudo Values by Group",
        xlab = "Group",
        ylab = "Pseudo Values")


boxplot(pseudo ~ group, data = meta1_clean,
        main = "Pseudo Values by Group",
        xlab = "Group",
        ylab = "Pseudo Values")


# Process the metadata
#new_meta3 <- process_metadata(meta3, time_grid)
new_meta2 <- process_metadata(meta2, time_grid)
new_meta1 <- process_metadata(meta1, time_grid)
new_meta <- rbind(new_meta1,new_meta2)


set.seed(1234)

# Remove rows with any NaN values for k-means clustering
scaled_mat_no_nan <- scaled_mat[, colSums(is.na(scaled_mat)) == 0]


# Perform k-means clustering
k <- 5  # Number of clusters, you can adjust this
km <- kmeans(scaled_mat_no_nan, centers = k, nstart = 25)

# Create a vector of cluster assignments for all rows
# (including those with NaN that were excluded from clustering)
cluster_assignments <- km$cluster


# Create column split factor
column_split <- rep(c("stem","K14"), each = 100)
column_split <- factor(column_split,levels = c("stem","K14"))
# Create color mappings for level1_anno and group
level1_anno_colors <- create_color_map(new_meta$max_level1_anno, "Set1")
group_colors <- create_color_map(new_meta$max_group, "Set2")
hmAll <- Heatmap(scaled_mat,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 row_split = cluster_assignments,  # Assuming you have this from k-means
                 column_split = column_split,  # Splitting into 3 groups of 100
                 col = colorRamp2(c(-2, 0, 4), c("Deepskyblue3", "white", "red")),
                 border = T,
                 na_col = "grey",
                 top_annotation = HeatmapAnnotation(
                   level1_anno = new_meta$max_level1_anno,
                   group = new_meta$max_group,
                   annotation_name_gp = gpar(fontsize = 8),
                   annotation_legend_param = list(
                     level1_anno = list(title = "Level 1 Anno"),
                     group = list(title = "Group")
                   ),
                   col = list(
                     level1_anno = level1_anno_colors,
                     group = group_colors
                   )
                 ))

hmAll <- draw(hmAll)


sigGene <- cluster_assignments
# Define the desired order and new names
cluster_order <- c(5,3,1,4,2)
new_cluster_names <- paste0("GC", 1:5)

# Create a named vector for renaming
rename_vector <- setNames(new_cluster_names, cluster_order)

# Reorder and rename the sigGene2 vector
sigGene2_ordered <- factor(sigGene, levels = cluster_order)
levels(sigGene2_ordered) <- rename_vector[levels(sigGene2_ordered)]

new_meta$max_level1_anno <- factor(new_meta$max_level1_anno, c("stem cells", "Transit","K14(-)", "K6+ cells"))
pdf("descriptive_results/rna_trajectory_gene/20250117_genescore_traj.pdf",width = 6,height = 8)
hmAll_score <- Heatmap(scaled_mat,
                  cluster_columns = FALSE,
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  cluster_row_slices = F,
                  show_row_dend = F,
                  row_split = sigGene2_ordered,  # Use the reordered and renamed factor
                  column_split = column_split,  # Splitting into 3 groups of 100
                  col = colorRamp2(c(-2, 0, 4), c("Deepskyblue3", "white", "red")),
                  border = TRUE,
                  na_col = "grey",
                  row_title_gp = gpar(fontsize = 10),
                  row_title_rot = 0,  # Horizontal text
                  row_gap = unit(1, "mm"),  # Add some gap between row clusters
                  top_annotation = HeatmapAnnotation(
                    level1_anno = new_meta$max_level1_anno,
                    group = new_meta$max_group,
                    annotation_name_gp = gpar(fontsize = 8),
                    annotation_legend_param = list(
                      level1_anno = list(title = "Level 1 Anno"),
                      group = list(title = "Group")
                    ),
                    col = list(
                      level1_anno = level1_anno_colors,
                      group = group_colors
                    )
                  ))

# Draw the heatmap
hmAll_score <- draw(hmAll_score)
dev.off()
# Export the table
# First, let's create a data frame with the gene clusters
gene_clusters_score <- data.frame(
  Gene = rownames(scaled_mat),
  Cluster = sigGene2_ordered
)
# Write the table to a CSV file
write.csv(gene_clusters_score, file = "descriptive_results/rna_trajectory_gene/20250123_genescore_clusters.csv", row.names = FALSE)

write.csv(scaled_mat, file = "descriptive_results/rna_trajectory_gene/20250123_motif_scaled_mat.csv", row.names = FALSE)
#scaled_mat
saveRDS(hmAll2,"process/trajectory/20250123_motif_traj_heatmap.Rds")

