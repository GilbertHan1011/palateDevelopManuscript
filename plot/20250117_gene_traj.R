# pseudotime1 <- read.csv("../2024.4_scATAC/processed_data/trajectory/5.1_rna_celloracle.pseudotime_run2.csv",row.names = 1)
# seurat$pseudo <- pseudotime1$Pseudotime
# FeaturePlot(seurat,"pseudo",label = T,reduction = "scviumap")+
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))

library(slingshot)
library(tidyverse)
library(destiny)
library(tradeSeq)
library(circlize)
library(ComplexHeatmap)
source("../function/seurat_utils.R")
scvi <- read.csv("../2024.4_scATAC/processed_data/embedding/RNA_scvi_latent.csv",row.names = 1)
colnames(scvi) <- paste("scvi",1:30,sep = "_")
seurat@reductions$scvi = CreateDimReducObject(embeddings = as.matrix(scvi),
                                              key = "scVI_",global = T)

find_sigmas(seurat@reductions$scvi@cell.embeddings)
dm <- run_diffMap(t(seurat@reductions$scvi@cell.embeddings), seurat$level1_anno, sigma = 2.2)
mycolor3<-colorRampPalette(brewer.pal(8,'Spectral'))(7)
# pdf("descriptive_results/rna_trajectory_gene/20241020_splom_destiny.pdf",width = 8,height = 8)
# splom(~dm@eigenvectors[, 1:5], groups = seurat$level1_anno, col = mycolor3,main = "", 
#       key = list(space="right", points = list(pch = 20,col = mycolor3), 
#                  text = list(c(unique(seurat$level1_anno)))))
# dev.off()
# png("descriptive_results/rna_trajectory_gene/20241020_splom_destiny.png",width = 800,height = 800)
# splom(~dm@eigenvectors[, 1:5], groups = seurat$level1_anno, col = mycolor3,main = "", 
#       key = list(space="right", points = list(pch = 20,col = mycolor3), 
#                  text = list(c(unique(seurat$level1_anno)))))
# dev.off()

lineage1 <- c("stem cells", "Transit", "K6+ cells")
lineage2 <- c("K14(-)", "K6+ cells")
cells1 <- colnames(seurat)[seurat$level1_anno%in%lineage1]
cells2 <- colnames(seurat)[seurat$level1_anno%in%lineage2]

cellsCol <- data.frame(colnames(seurat))
colnames(cellsCol) <- "cells"
cellsCol$stem <- 0
cellsCol$stem[cellsCol$cells %in% cells1] <- 1

cellsCol$K14 <- 0
cellsCol$K14[cellsCol$cells %in% cells2] <- 1
cellsCol <- cellsCol %>% column_to_rownames("cells")


pseudotimeCol <- data.frame(colnames(seurat))
colnames(pseudotimeCol) <- "cells"
pseudotimeCol$stem <- pseudotime1$Pseudotime_Lineage_K6
pseudotimeCol$K14 <- pseudotime1$Pseudotime_Lineage_K14
pseudotimeCol <- pseudotimeCol %>% column_to_rownames("cells")

cellsCol$stem[is.na(pseudotimeCol$stem)] <- 0
cellsCol$K14[is.na(pseudotimeCol$K14)] <- 0

sum(is.na(pseudotimeCol$stem))


# slingshot_res <- slingshot(dm@eigenvectors[, c(1:5)],
#                            clusterLabels = factor(seurat$level1_anno),
#                            start.clus = 'stem cells',end.clus=c("K6+ cells"),allow.breaks=FALSE,
#                            maxit = 1000, shrink.method = "density", thresh = 0.001, extend = "n")
# 
# pseudo <- slingshot_res@assays@data$pseudotime
# lineage <- slingshot_res@assays@data$weights

BPPRARM <-BiocParallel::bpparam()
BPPRARM$workers<- 10

filterCell <- rowSums(cellsCol)!=0
count <- as.matrix(seurat@assays$RNA@counts)
count <- count[,filterCell]
pseudotimeCol <- pseudotimeCol[filterCell,]
cellsCol <- cellsCol[filterCell,]
#undebug(evaluateK)
batch <- factor(seurat$orig.ident)
U <- model.matrix(~batch)
U <- U[filterCell,]
icMat <- evaluateK(counts = count,
                   pseudotime =pseudotimeCol,
                   cellWeights = cellsCol,
                   k = 3:12,
                   parallel=TRUE,
                   nGenes=100,
                   BPPARAM=BPPRARM,
                   verbose = T)

pseudotimeCol[is.na(pseudotimeCol)] <- 0
varGene <- read.csv("../2024.4_scATAC/processed_data/variable_gene/5.1_rna_var_scanpy_2000.csv",row.names = 1) %>% unlist

gam <- fitGAM(counts = count,
              pseudotime =pseudotimeCol,
              cellWeights = cellsCol,
              nknots = 8,
              U=U,
              parallel=TRUE,
              genes=varGene,
              BPPARAM=BPPRARM,
              verbose = T)

saveRDS(gam,"../2024.4_scATAC/processed_data/trajectory/20250107_gam_threelineage.Rds")

rowData(gam)$assocRes <- associationTest(gam, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(gam)$assocRes

line1Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.01)
]
line2Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.01)
]
# line3Genes <-  rownames(assocRes)[
#   which(p.adjust(assocRes$pvalue_3, "fdr") <= 0.01)
# ]

UpSetR::upset(UpSetR::fromList(list(line1 = line1Genes, line2 = line2Genes)))

yhatSmooth <- lapply(list(line1Genes,line2Genes),function(x) 
  predictSmooth(gam, gene = x, nPoints = 100, tidy = FALSE))

geneUnion <- unique(c(line1Genes,line2Genes))
yhatSmoothAll <- predictSmooth(gam, gene = geneUnion, nPoints = 100, tidy = FALSE)
yhatSmoothAllTidy <- predictSmooth(gam, gene = geneUnion, nPoints = 500, tidy = TRUE)


#time_grid <- seq(max(yhatSmoothAllTidy$time), min(yhatSmoothAllTidy$time), length.out = 100)
library(pracma)
# Function to interpolate and average yhat values
time_grid <- seq(min(yhatSmoothAllTidy$time), max(yhatSmoothAllTidy$time), length.out = 100)

# Function to interpolate and handle missing values
interpolate_with_nan <- function(data, new_time) {
  # Find the min and max time for this data
  min_time <- min(data$time)
  max_time <- max(data$time)
  
  # Perform interpolation
  interp_values <- approx(data$time, data$yhat, xout = new_time, method = "linear", rule = 1)$y
  
  # Set values outside the original time range to NaN
  interp_values[new_time < min_time | new_time > max_time] <- NaN
  
  return(interp_values)
}

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

# Custom row scaling function that preserves NaN
scale_rows_preserve_nan <- function(x) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- apply(x, 1, sd, na.rm = TRUE)
  
  # Create a matrix of means and sds
  means_matrix <- matrix(row_means, nrow = nrow(x), ncol = ncol(x))
  sds_matrix <- matrix(row_sds, nrow = nrow(x), ncol = ncol(x))
  
  # Scale the matrix
  scaled_mat <- (x - means_matrix) / sds_matrix
  
  # Preserve original NaN values
  scaled_mat[is.nan(x)] <- NaN
  
  return(scaled_mat)
}

# Apply the custom scaling function to mat
scaled_mat <- scale_rows_preserve_nan(mat)

# Get a summary of the values in mat
#summary(as.vector(mat))

# Check for infinite values
sum(is.infinite(mat))

# Check for NaN values
sum(is.nan(mat))


set.seed(1234)
hmAll<- Heatmap(scaled_mat,cluster_columns = F,show_column_names = F,
                show_row_names = F,
                col=colorRamp2(c(-2, 0, 4), c("blue", "white", "red")),na_col = "grey")
hmAll<- draw(hmAll)


set.seed(1234)

# Remove rows with any NaN values for k-means clustering
scaled_mat_no_nan <- scaled_mat[, colSums(is.na(scaled_mat)) == 0]


# Perform k-means clustering
k <- 8  # Number of clusters, you can adjust this
km <- kmeans(scaled_mat_no_nan, centers = k, nstart = 25)

# Create a vector of cluster assignments for all rows
# (including those with NaN that were excluded from clustering)
cluster_assignments <- km$cluster

# Create column split factor
column_split <- rep(c("stem","K14"), each = 100)
column_split <- factor(column_split,levels = c("stem","K14"))

# Create the heatmap with clustering results and column split
hmAll <- Heatmap(scaled_mat,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 row_split = cluster_assignments,  # Use k-means clusters to split rows
                 column_split = column_split,  # Split columns into 3 groups
                 col = colorRamp2(c(-2, 0, 4), c("Deepskyblue3", "white", "red")),
                 na_col = "grey",border = T)

hmAll <- draw(hmAll)


meta <- seurat@meta.data[c("level1_anno","group")]

lineage1Cell <- rownames(cellsCol)[cellsCol$stem==1]
lineage2Cell <- rownames(cellsCol)[cellsCol$K14==1]

meta1 <- meta[lineage1Cell,]
meta1$pseudo <- pseudotimeCol[lineage1Cell,]$stem
meta2 <- meta[lineage2Cell,]
meta2$pseudo <- pseudotimeCol[lineage2Cell,]$K14
# meta3 <- meta[lineage3Cell,]
# meta3$pseudo <- pseudotimeCol[lineage3Cell,]$K14

find_interval <- function(x, grid) {
  findInterval(x, grid)
}

# Process the metadata

process_metadata <- function(meta_df, time_grid) {
  # Convert inputs to appropriate types and add error checking
  meta_df <- as.data.frame(meta_df)
  time_grid <- as.numeric(time_grid)
  
  # Verify required columns exist
  required_cols <- c("pseudo", "level1_anno", "group")
  if (!all(required_cols %in% colnames(meta_df))) {
    stop("Missing required columns in meta_df: ", 
         paste(setdiff(required_cols, colnames(meta_df)), collapse = ", "))
  }
  
  # Process step by step with explicit data.frame operations
  # Step 1: Create time intervals
  meta_df$time_interval <- findInterval(meta_df$pseudo, time_grid)
  
  # Step 2: Filter valid intervals
  meta_df <- meta_df[meta_df$time_interval > 0, ]
  
  # Step 3: Group and summarize
  new_meta <- aggregate(
    cbind(pseudo = meta_df$pseudo) ~ time_interval, 
    data = meta_df,
    FUN = function(x) {
      data.frame(
        max_level1_anno = names(which.max(table(meta_df$level1_anno[meta_df$pseudo %in% x]))),
        max_group = names(which.max(table(meta_df$group[meta_df$pseudo %in% x]))),
        pseudo_time = mean(x)
      )
    }
  )
  
  # Step 4: Create complete intervals
  complete_intervals <- data.frame(time_interval = 1:100)
  
  # Step 5: Merge with complete intervals
  final_meta <- merge(complete_intervals, new_meta, 
                      by = "time_interval", all.x = TRUE)
  
  # Step 6: Sort by time interval
  final_meta <- final_meta[order(final_meta$time_interval), ]
  
  return(final_meta)
}
process_metadata <- function(meta_df, time_grid) {
  # Convert inputs to appropriate types and add error checking
  meta_df <- as.data.frame(meta_df)
  time_grid <- as.numeric(time_grid)
  
  # Verify required columns exist
  required_cols <- c("pseudo", "level1_anno", "group")
  if (!all(required_cols %in% colnames(meta_df))) {
    stop("Missing required columns in meta_df: ", 
         paste(setdiff(required_cols, colnames(meta_df)), collapse = ", "))
  }
  
  # Step 1: Create time intervals
  meta_df$time_interval <- findInterval(meta_df$pseudo, time_grid)
  
  # Step 2: Filter valid intervals
  meta_df <- meta_df[meta_df$time_interval > 0, ]
  
  # Step 3: Group and summarize - split into separate operations
  level1_anno_summary <- tapply(meta_df$level1_anno, meta_df$time_interval, 
                                function(x) names(which.max(table(x))))
  group_summary <- tapply(meta_df$group, meta_df$time_interval, 
                          function(x) names(which.max(table(x))))
  pseudo_time_summary <- tapply(meta_df$pseudo, meta_df$time_interval, mean)
  
  # Combine the summaries
  new_meta <- data.frame(
    time_interval = as.numeric(names(level1_anno_summary)),
    max_level1_anno = unname(level1_anno_summary),
    max_group = unname(group_summary),
    pseudo_time = unname(pseudo_time_summary)
  )
  
  # Step 4: Create complete intervals
  complete_intervals <- data.frame(time_interval = 1:100)
  
  # Step 5: Merge with complete intervals
  final_meta <- merge(complete_intervals, new_meta, 
                      by = "time_interval", all.x = TRUE)
  
  # Step 6: Sort by time interval
  final_meta <- final_meta[order(final_meta$time_interval), ]
  
  return(final_meta)
}
# 
# process_metadata <- function(meta_df, time_grid) {
#   # First ensure meta_df is a data frame/tibble
#   meta_df <- as_tibble(meta_df)
#   
#   # Create time intervals explicitly
#   new_meta <- meta_df %>%
#     mutate(time_interval = findInterval(pseudo, time_grid)) %>%
#     # Filter out any rows where time_interval is 0
#     filter(time_interval > 0) %>%
#     # Ensure time_interval exists before grouping
#     group_by(time_interval) %>%
#     summarise(
#       max_level1_anno = names(which.max(table(level1_anno))),
#       max_group = names(which.max(table(group))),
#       pseudo_time = mean(pseudo),
#       .groups = 'drop'  # Explicitly drop grouping
#     ) %>%
#     arrange(time_interval)
#   
#   # Create complete intervals
#   complete_intervals <- tibble(time_interval = 1:100)
#   
#   # Join and handle missing values
#   final_meta <- complete_intervals %>%
#     left_join(new_meta, by = "time_interval") %>%
#     mutate(
#       max_level1_anno = if_else(is.na(max_level1_anno), NA_character_, max_level1_anno),
#       max_group = if_else(is.na(max_group), NA_character_, max_group),
#       pseudo_time = if_else(is.na(pseudo_time), NA_real_, pseudo_time)
#     )
#   
#   return(final_meta)
# }
#new_meta3 <- process_metadata(meta3, time_grid)
new_meta2 <- process_metadata(meta2, time_grid)
new_meta1 <- process_metadata(meta1, time_grid)
new_meta <- rbind(new_meta1,new_meta2)

# Function to create a color mapping that includes gray for NA
create_color_map <- function(values, palette_name = "Set1") {
  unique_values <- unique(values[!is.na(values)])
  n_colors <- length(unique_values)
  
  if (n_colors == 0) {
    return(c("NA" = "gray"))
  } else if (n_colors == 1) {
    colors <- brewer.pal(3, palette_name)[1]
  } else if (n_colors == 2) {
    colors <- brewer.pal(3, palette_name)[1:2]
  } else if (n_colors <= 8) {
    colors <- brewer.pal(n_colors, palette_name)
  } else {
    colors <- colorRampPalette(brewer.pal(8, palette_name))(n_colors)
  }
  
  color_map <- setNames(colors, unique_values)
  if (any(is.na(values))) {
    color_map <- c(color_map, "NA" = "gray")
  }
  return(color_map)
}



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


#== some modification--------------
# sigGene2 <- cluster_assignments[cluster_assignments!=2]
# 
# new_meta1[54:100,c("max_level1_anno","max_group","pseudo_time")] <- NA
# new_meta <- rbind(new_meta1,new_meta2,new_meta3)
# mat[,c(54:100)] <- NA
# scaled_mat <- scale_rows_preserve_nan(mat)
# 
# scaled_mat <- scaled_mat[names(sigGene2),]

sigGene <- cluster_assignments
# Define the desired order and new names
cluster_order <- c(4,6,1,2,5,3,8,7)
new_cluster_names <- paste0("GC", 1:8)

# Create a named vector for renaming
rename_vector <- setNames(new_cluster_names, cluster_order)

# Reorder and rename the sigGene2 vector
sigGene2_ordered <- factor(sigGene, levels = cluster_order)
levels(sigGene2_ordered) <- rename_vector[levels(sigGene2_ordered)]

new_meta$max_level1_anno <- factor(new_meta$max_level1_anno, c("stem cells", "Transit","K14(-)", "K6+ cells"))
pdf("descriptive_results/rna_trajectory_gene/20250117_gene_traj.pdf",width = 6,height = 8)
hmAll2 <- Heatmap(scaled_mat,
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
hmAll2 <- draw(hmAll2)
dev.off()
# Export the table
# First, let's create a data frame with the gene clusters
gene_clusters <- data.frame(
  Gene = rownames(scaled_mat),
  Cluster = sigGene2_ordered
)
# Write the table to a CSV file
write.csv(gene_clusters, file = "descriptive_results/rna_trajectory_gene/20250117_gene_clusters.csv", row.names = FALSE)

write.csv(scaled_mat, file = "descriptive_results/rna_trajectory_gene/20250117_gene_scaled_mat.csv", row.names = FALSE)
#scaled_mat
saveRDS(hmAll2,"process/trajectory/20250117_gene_traj_heatmap.Rds")



