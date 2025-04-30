
pdf("descriptive_results/rna_trajectory_gene/20250117_trajectory_two.pdf",width = 8,height = 5)
mycolor<- c("#282828","#F2C9D5", "#B43E44", "#8B6D9C", "#496496", "#904869","#FADF92" )
names(mycolor) <- c( "K14(-)",  "K5(-)", "K6+ cells", "Mature K5(-) ", "Shh(+)","stem cells","Transit" )
emb <- seurat@reductions$scviumap@cell.embeddings
plot(seurat@reductions$scviumap@cell.embeddings, col = mycolor[as.character(seurat$level1_anno)], 
     pch=16, 
     asp = 1,
     cex = .5)
# Get the plot dimensions
plot_dims <- par("usr")
x_range <- plot_dims[2] - plot_dims[1]
y_range <- plot_dims[4] - plot_dims[3]

# Add legend for clusters (adjust the x and y values as needed)
legend(x = plot_dims[1],  # 5% from the left edge
       y = plot_dims[4],  # Top of the plot
       legend = names(mycolor), 
       fill = mycolor,
       bty = "n",  # No box around the legend
       xpd = TRUE) 
# Function to calculate cluster centers
get_cluster_center <- function(embedding, cluster) {
  cluster_cells <- embedding[seurat$level1_anno == cluster, ]
  return(colMeans(cluster_cells))
}

# Define the clusters you want to connect and their order
clusters_to_connect1 <- c("stem cells",  "K6+ cells")
clusters_to_connect2 <- c("K14(-)", "K6+ cells")
# clusters_to_connect3 <- c("K14(-)", "K6+ cells")

# Define colors for each path
path_colors <- RColorBrewer::brewer.pal(name = "Set3",n = 3)

# Function to draw arrows for a given path with offset
draw_arrows <- function(clusters, color, offset = c(0, 0)) {
  centers <- lapply(clusters, function(cluster) {
    center <- get_cluster_center(seurat@reductions$scviumap@cell.embeddings, cluster)
    return(center + offset)
  })
  for (i in 1:(length(centers) - 1)) {
    arrows(centers[[i]][1], centers[[i]][2],
           centers[[i+1]][1], centers[[i+1]][2],
           length = 0.15,  # Increased arrow head size
           col = color, 
           lwd = 4,  # Increased line width
           angle = 20)  # Adjusted arrow head angle
  }
}

# Calculate the direction of the offset
transit_center <- get_cluster_center(seurat@reductions$scviumap@cell.embeddings, "Transit")
k6_center <- get_cluster_center(seurat@reductions$scviumap@cell.embeddings, "K6+ cells")
direction <- k6_center - transit_center
perpendicular <- c(-direction[2], direction[1])  # Rotate 90 degrees
offset <- 0 * perpendicular / sqrt(sum(perpendicular^2))  # Scale to 10% of the plot range

# Draw arrows for each path
draw_arrows(clusters_to_connect1, path_colors[1], -offset/2)  # Offset in one direction
draw_arrows(clusters_to_connect2, path_colors[2], offset/2)   # Offset in the opposite direction
#draw_arrows(clusters_to_connect3, path_colors[3])

# Add legend for paths (adjust the x and y values as needed)
legend(x = plot_dims[2]-8,  # Right edge of the plot
       y = plot_dims[4],  # 10% above the top of the plot
       legend = c("stem->K6(+)",
                  "K14(-)->K6(+)"),
       col = path_colors,
       lty = 1,
       lwd = 4,
       bty = "n",  # No box around the legend
       xpd = TRUE,  # Allow plotting outside the plot region
       horiz = FALSE)  # Vertical legend
dev.off()


pseudotime1 <- read.csv("../2024.4_scATAC/processed_data/trajectory/5.1_rna_celloracle.pseudotime_run2.csv",row.names = 1)
seurat$pseudo <- pseudotime1$Pseudotime

pdf("descriptive_results/rna_trajectory_gene/20250117_trajectory_three_pseudo.pdf",width = 8,height = 5)
# Get the UMAP coordinates
umap_coords <- seurat@reductions$scviumap@cell.embeddings

# Create a color gradient
spectral_colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
color_palette <- colorRampPalette(spectral_colors)(100)

# Normalize pseudo values to 0-1 range for color mapping
pseudo_normalized <- (seurat$pseudo - min(seurat$pseudo)) / (max(seurat$pseudo) - min(seurat$pseudo))

# Create the plot
plot(umap_coords, 
     col = color_palette[ceiling(pseudo_normalized * 99) + 1],  # Map normalized values to color palette
     pch = 16, 
     asp = 1,
     cex = 0.5,
     xlab = "UMAP 1",
     ylab = "UMAP 2",
     main = "UMAP colored by Pseudotime")

# Add a color legend
color_legend <- as.raster(matrix(rev(color_palette), ncol = 1))
rasterImage(color_legend, 
            xleft = par("usr")[2], 
            xright = par("usr")[2] + (par("usr")[2] - par("usr")[1])*0.05, 
            ybottom = par("usr")[3], 
            ytop = par("usr")[4])

# Add text labels to the color legend
text(x = par("usr")[2] + (par("usr")[2] - par("usr")[1])*0.07, 
     y = c(par("usr")[3], mean(par("usr")[3:4]), par("usr")[4]),
     labels = round(c(min(seurat$pseudo), mean(range(seurat$pseudo)), max(seurat$pseudo)), 2),
     adj = 0)
# Define the clusters you want to connect and their order
clusters_to_connect1 <- c("stem cells", "K6+ cells")

clusters_to_connect2 <- c("K14(-)", "K6+ cells")

# Define colors for each path
path_colors <- RColorBrewer::brewer.pal(name = "Set1",n = 3)

# Function to draw arrows for a given path with offset
draw_arrows <- function(clusters, color, offset = c(0, 0)) {
  centers <- lapply(clusters, function(cluster) {
    center <- get_cluster_center(seurat@reductions$scviumap@cell.embeddings, cluster)
    return(center + offset)
  })
  for (i in 1:(length(centers) - 1)) {
    arrows(centers[[i]][1], centers[[i]][2],
           centers[[i+1]][1], centers[[i+1]][2],
           length = 0.15,  # Increased arrow head size
           col = color, 
           lwd = 4,  # Increased line width
           angle = 20)  # Adjusted arrow head angle
  }
}

# Calculate the direction of the offset
transit_center <- get_cluster_center(seurat@reductions$scviumap@cell.embeddings, "Transit")
k6_center <- get_cluster_center(seurat@reductions$scviumap@cell.embeddings, "K6+ cells")
direction <- k6_center - transit_center
perpendicular <- c(-direction[2], direction[1])  # Rotate 90 degrees
offset <- 0 * perpendicular / sqrt(sum(perpendicular^2))  # Scale to 10% of the plot range

# Draw arrows for each path
draw_arrows(clusters_to_connect1, path_colors[1], -offset/2)  # Offset in one direction
draw_arrows(clusters_to_connect2, path_colors[2], offset/2)   # Offset in the opposite direction
#draw_arrows(clusters_to_connect3, path_colors[3])

# Add legend for paths (adjust the x and y values as needed)
legend(x = plot_dims[2]-8,  # Right edge of the plot
       y = plot_dims[4],  # 10% above the top of the plot
       legend = c("stem->K6(+)",
                  "K14(-)->K6(+)"),
       col = path_colors,
       lty = 1,
       lwd = 4,
       bty = "n",  # No box around the legend
       xpd = TRUE,  # Allow plotting outside the plot region
       horiz = FALSE)  # Vertical legend
#dev.off()
dev.off()

# saveRDS(hmAll2,"descriptive_results/rna_trajectory_gene/20241021_trajHM.Rds")

