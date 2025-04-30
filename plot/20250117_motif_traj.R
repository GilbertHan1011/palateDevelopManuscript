motifObj <- readRDS("process/framework/motif/20250117_motifMat.Rds")

motifSeurat <- CreateSeuratObject(motifObj@assays@data$deviations)
count_exp <- exp(count)
pairCell <- read.csv("process/integration/20241222_FirR_pairing.csv",row.names = 1)
# Create a mapping dictionary
prefix_map <- c(
  "C2" = "E12_2",
  "C1" = "E12_1",
  "B2" = "E14_2",
  "B1" = "B1"
  # Add more mappings if needed
)

# Function to replace prefix
# Vectorized function to replace prefixes
# Modified function to remove the original prefix
replace_prefix <- function(x, prefix_map) {
  # Extract the prefix (C1, C2, etc.)
  old_prefixes <- sub("^([^_]+)_.*", "\\1", x)
  
  # Get the barcode part (everything after the underscore, before -1)
  barcodes <- sub("^[^_]+_(.*)-1", "\\1", x)
  
  # Get new prefixes
  new_prefixes <- prefix_map[old_prefixes]
  
  # Combine with barcode
  new_strings <- paste0(new_prefixes, "#", barcodes, "-1")
  return(new_strings)
}

# Apply the transformation
pairCell$ATAC <- replace_prefix(pairCell$ATAC, prefix_map)

intersect(pairCell$ATAC, colnames(motifSeurat)) %>% length()

atac_label <- read.csv("process/framework/cluster/20241211_knn_label.csv")

ident <- motifSeurat %>% colnames%>% strsplit("#") %>% lapply(., function(x) x[[1]]) %>% unlist
motifSeurat$ident <- ident


uniqueLabel <- atac_label$X %>% strsplit("_")
labels <- lapply(uniqueLabel, function(x) x[[1]]) %>% unlist() %>% unique

atac_label$name <- replace_prefix(atac_label$X, prefix_map)
intersect(atac_label$name, colnames(motifSeurat)) %>% length()
rownames(atac_label) <- atac_label$name
motifSeurat$label <- atac_label[colnames(motifSeurat),][2]

rnaPseudo <- seurat$pseudo

colnames(motifSeurat)
get_pseudotime <- function(atac_cells, pair_df, pseudotime_data) {
  # Convert pseudotime_data names to match RNA format in pair_df
  names(pseudotime_data) <- sub("_", "_", names(pseudotime_data), fixed=TRUE)
  
  # For each ATAC cell:
  # 1. Find its paired RNA cell from pair_df
  # 2. Look up that RNA cell's pseudotime
  pseudotimes <- sapply(atac_cells, function(atac_cell) {
    # Find the paired RNA cell
    rna_cell <- pair_df$RNA[pair_df$ATAC == atac_cell]
    
    # If we found a pair, get its pseudotime
    if (length(rna_cell) > 0) {
      return(pseudotime_data[rna_cell])
    } else {
      return(NA)
    }
  })
  
  return(pseudotimes)
}

pseudo <- get_pseudotime(colnames(motifSeurat), pairCell,rnaPseudo)
motifSeurat$pseudo <-as.numeric(pseudo)
boxplot(pseudo ~ label, data = motifSeurat@meta.data)

ggplot( motifSeurat@meta.data, aes(x = label, y = pseudo)) +
  geom_boxplot(fill = "lightblue", outlier.color = "red") +
  theme_minimal() +
  labs(title = "Your Title",
       x = "X Label",
       y = "Y Label") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels if needed

#motifSeurat$pseudotime <- pairCell$

# Define lineages and get cells for each
lineage1 <- c("stem cells", "Transit", "K6+ cells")
lineage2 <- c("K14(-)", "K6+ cells") 
cells1 <- colnames(motifSeurat)[motifSeurat$label%in%lineage1]
cells2 <- colnames(motifSeurat)[motifSeurat$label%in%lineage2]

# Create cell weight matrix indicating which cells belong to each lineage
cellsCol <- data.frame(colnames(motifSeurat))
colnames(cellsCol) <- "cells"
cellsCol$stem <- 0
cellsCol$stem[cellsCol$cells %in% cells1] <- 1

cellsCol$K14 <- 0 
cellsCol$K14[cellsCol$cells %in% cells2] <- 1
cellsCol <- cellsCol %>% column_to_rownames("cells")

# Create pseudotime matrix with values for each lineage
pseudotimeCol <- data.frame(colnames(motifSeurat))
colnames(pseudotimeCol) <- "cells"




pseudotimeCol <- data.frame(colnames(motifSeurat))
colnames(pseudotimeCol) <- "cells"
pseudotimeCol$stem <- pseudo_K6
pseudotimeCol$K14 <- pseudo_K14
pseudotimeCol <- pseudotimeCol %>% column_to_rownames("cells")

pseudotimeCol$stem <-pseudo_K6
pseudotimeCol$K14 <- pseudo_K14
pseudotimeCol <- pseudotimeCol %>% column_to_rownames("cells")
Pseudotime_Lineage_K6 <- pseudotime1$Pseudotime_Lineage_K6
names(Pseudotime_Lineage_K6) <- rownames(pseudotime1)
Pseudotime_Lineage_K14 <- pseudotime1$Pseudotime_Lineage_K14
names(Pseudotime_Lineage_K14) <- rownames(pseudotime1)
pseudo_K6 <- get_pseudotime(colnames(motifSeurat), pairCell,Pseudotime_Lineage_K6)
pseudo_K14 <- get_pseudotime(colnames(motifSeurat), pairCell,Pseudotime_Lineage_K14)

motifSeurat$pseudo_K6 <- as.numeric(pseudo_K6)
motifSeurat$pseudo_K14 <- as.numeric(pseudo_K14)


ggplot(motifSeurat@meta.data %>% filter(!is.na(pseudo_K6)), 
       aes(x = label, y = pseudo_K6)) +
  geom_boxplot(fill = "lightblue", outlier.color = "red") +
  theme_minimal() +
  labs(title = "Pseudotime Distribution by Label",
       x = "Cell Type",
       y = "Pseudotime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Adjust text alignment
        plot.title = element_text(hjust = 0.5),  # Center title
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95")) +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white") # Add mean points


# Set cells not in lineage to weight 0
cellsCol$stem[is.na(pseudo_K6)] <- 0
cellsCol$K14[is.na(pseudo_K14)] <- 0

# Set up parallel processing
BPPRARM <-BiocParallel::bpparam()
BPPRARM$workers<- 10

# Filter cells and prepare data for trajectory analysis
filterCell <- rowSums(cellsCol)!=0
count <- as.matrix(motifSeurat@assays$RNA@layers$counts)
count <- count[,filterCell]

pseudotimeCol

pseudotimeCol <- pseudotimeCol[filterCell,]
cellsCol <- cellsCol[filterCell,]

# Create batch correction matrix
batch <- factor(motifSeurat$ident)
U <- model.matrix(~batch)
U <- U[filterCell,]

pseudotimeCol_bk <- pseudotimeCol
pseudotimeCol_bk <- as.data.frame(pseudotimeCol_bk)
pseudotimeCol_bk <- pseudotimeCol_bk %>% column_to_rownames("cells")
pseudotimeCol$stem[is.na(pseudotimeCol$stem)] <- 0
pseudotimeCol$K14[is.na(pseudotimeCol$K14)] <- 0
colnames(count_exp) <- rownames(cellsCol)
rownames(count_exp) <- rownames(motifSeurat)
gam_motif <- fitGAM(counts = count_exp,
              pseudotime =pseudotimeCol,
              cellWeights = cellsCol,
              nknots = 8,
              U=U,
              parallel=TRUE,
              BPPARAM=BPPRARM,
              verbose = T)

rowData(gam_motif)$assocRes <- associationTest(gam_motif, lineages = TRUE, l2fc = 0)
assocRes <- rowData(gam_motif)$assocRes
View(assocRes)
line1Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.1)
]
line2Genes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.1)
]
# line3Genes <-  rownames(assocRes)[
#   which(p.adjust(assocRes$pvalue_3, "fdr") <= 0.01)
# ]

UpSetR::upset(UpSetR::fromList(list(line1 = line1Genes, line2 = line2Genes)))

yhatSmooth <- lapply(list(line1Genes,line2Genes),function(x) 
  predictSmooth(gam_motif, gene = x, nPoints = 100, tidy = FALSE))

geneUnion <- unique(c(line1Genes,line2Genes))
yhatSmoothAll <- predictSmooth(gam_motif, gene = geneUnion, nPoints = 100, tidy = FALSE)
yhatSmoothAllTidy <- predictSmooth(gam_motif, gene = geneUnion, nPoints = 500, tidy = TRUE)


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


# Check for infinite values
sum(is.infinite(mat))

# Check for NaN values
sum(is.nan(mat))


set.seed(1234)
hmAll<- Heatmap(scaled_mat,cluster_columns = F,show_column_names = F,
                show_row_names = T,
                col=colorRamp2(c(-2, 0, 4), c("blue", "white", "red")),na_col = "grey")
hmAll<- draw(hmAll)

