library(igraph)
library(RColorBrewer) # For better color palettes
setwd("../lh/2024_7.ATAC_E12/")
umap <- read.csv("../2024.4_scATAC/processed_data/embedding/RNA_scvi_umap.csv",row.names = 1)
label <- read.csv("../2024.4_scATAC/processed_data/cluster/rna_annotation_scvi_mod.csv",row.names = 1)


# Define the data as a multi-line string with consistent spacing
data <- read.csv("descriptive_results/rna_trajecotry/moscot_transition.csv", row.names = 1)
colnames(data) <- c("Transit", "K5(-)", "Mature K5(-) ", "Shh(+)", "K6+ cells", 
                    "stem cells", "K14(-)")
# Read the data into a matrix
mat <- as.matrix(data)

# Define the color palette based on provided colors
colorPalate3 <- c(
  "Transit" = "#FADF92",
  "K5(-)" = "#F2C9D5",
  "Mature K5(-) " = "#8B6D9C",
  "Shh(+)" = "#496496",
  "K6+ cells" = "#B43E44",
  "stem cells" = "#904869",
  "K14(-)" = "#282828"
)

# Create a directed graph from the adjacency matrix
G <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)

# Filter edges with weight < 0.05
G <- delete_edges(G, E(G)[weight < 0.05])

# Assign colors based on node names
V(G)$color <- colorPalate3[V(G)$name]

# Set node sizes based on both in and out degree
V(G)$size <- (degree(G, mode = "total") + 3) * 3

# Increase edge width based on weight
E(G)$width <- sqrt(E(G)$weight) * 8

# Add transparency to edges with higher contrast
edge_colors <- adjustcolor(ifelse(ends(G, E(G))[, 2] == "K6+ cells", 
                                  "firebrick", "steelblue"), alpha.f = 0.8)
E(G)$color <- edge_colors

# Set arrow parameters - bigger arrows
E(G)$arrow.size <- 0.6
E(G)$arrow.width <- 2.0
E(G)$curve <- 0.25  # Curved edges for better visibility

# Calculate mean UMAP coordinates for each cell type
# Note: This assumes umap and label dataframes are already in your R environment

# Combine UMAP and label data
combined_data <- data.frame(
  X0 = umap$X0,
  X1 = umap$X1,
  cell_type = label$level1_anno
)

# Calculate mean UMAP coordinates for each cell type
mean_coords <- aggregate(. ~ cell_type, data = combined_data[, c("cell_type", "X0", "X1")], FUN = mean)

# Create a layout matrix from the mean coordinates
manual_layout <- matrix(ncol = 2, nrow = length(V(G)$name))
rownames(manual_layout) <- V(G)$name

# Scale factor for adjusting the spread of nodes
scale_factor <- 0.05  # Adjust as needed

# Fill the layout matrix with mean UMAP coordinates
for (node in V(G)$name) {
  # Check for exact matches first
  if (node %in% mean_coords$cell_type) {
    idx <- which(mean_coords$cell_type == node)
    manual_layout[node, ] <- c(mean_coords$X0[idx], mean_coords$X1[idx])
  } else {
    # Try to match without whitespace differences or case sensitivity
    normalized_node <- gsub("\\s+", "", tolower(node))
    normalized_types <- gsub("\\s+", "", tolower(mean_coords$cell_type))
    matching_idx <- which(normalized_types == normalized_node)
    
    if (length(matching_idx) > 0) {
      manual_layout[node, ] <- c(mean_coords$X0[matching_idx[1]], mean_coords$X1[matching_idx[1]])
    } else {
      # If a node type doesn't exist in the data, use a fallback position
      warning(paste("Node type", node, "not found in UMAP data. Using fallback position."))
      manual_layout[node, ] <- c(0, 0)  # Fallback to origin
    }
  }
}

# Print the node positions for reference
print(mean_coords)
print(manual_layout)

# Options for final adjustments
flip_x <- FALSE  # Set to TRUE to flip the X axis
flip_y <- FALSE  # Set to TRUE to flip the Y axis
center_layout <- TRUE  # Center the layout around origin

# Apply flipping if needed
if (flip_x) {
  manual_layout[, 1] <- -manual_layout[, 1]
}
if (flip_y) {
  manual_layout[, 2] <- -manual_layout[, 2]
}

# Center the layout if requested
if (center_layout) {
  manual_layout[, 1] <- manual_layout[, 1] - mean(manual_layout[, 1])
  manual_layout[, 2] <- manual_layout[, 2] - mean(manual_layout[, 2])
}

# Manual adjustments for specific nodes (if needed)
# Uncomment and modify these lines to make fine adjustments to specific nodes
# manual_layout["Transit", ] <- manual_layout["Transit", ] + c(0.2, 0)  # Move Transit right by 0.2
# manual_layout["K6+ cells", ] <- manual_layout["K6+ cells", ] + c(0, -0.3)  # Move K6+ cells down by 0.3

# Apply scaling factor to control the spread
manual_layout <- manual_layout * scale_factor

# OPTIONAL: Save coordinates to reuse later
# write.csv(manual_layout, "descriptive_results/rna_trajecotry/node_positions.csv")

# Use the calculated layout
layout_choice <- manual_layout

# Create a legend for node types
legend_labels <- V(G)$name
legend_colors <- colorPalate3[legend_labels]

# Open PDF with higher resolution
pdf("descriptive_results/rna_trajecotry/plot/moscot_network2.pdf", width = 6, height = 6)

# Set up the plotting area with margins for the legend
par(mar = c(1, 1, 3, 8))

# Plot the network
plot(G, 
     layout = layout_choice,
     edge.color = E(G)$color, 
     edge.arrow.size = E(G)$arrow.size, 
     edge.arrow.width = E(G)$arrow.width,
     vertex.label.color = "black", 
     vertex.label.cex = 1.2,
     vertex.label.family = "sans",
     vertex.label.font = 2,  # Bold font
     vertex.frame.color = "black",
     vertex.frame.width = 2.0,
     main = "Cell Type Transition Network", 
     sub = "Edge thickness represents transition probabilities",
     margin = c(0, 0, 0, 0))

# Add a legend
legend("topright", 
       legend = legend_labels, 
       col = legend_colors,
       pch = 19, 
       pt.cex = 2.5, 
       bty = "n", 
       title = "Cell Types",
       text.col = "black",
       inset = c(-0.25, 0))

# Add a note about red edges
legend("bottomright", 
       legend = c("Transitions to K6
                  + cells", "Other transitions"), 
       col = c("firebrick", "steelblue"), 
       lwd = 3, 
       bty = "n",
       inset = c(0.05, 0.05))

dev.off()

