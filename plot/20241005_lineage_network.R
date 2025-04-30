library(igraph)

# Define the data as a multi-line string with consistent spacing
data <- read.csv("descriptive_results/rna_trajecotry/moscot_transition.csv",row.names = 1)
colnames(data) <- c("Transit", "K5(-)", "Mature K5(-) ", "Shh(+)", "K6+ cells", 
                    "stem cells", "K14(-)")
# Read the data into a matrix
mat <- as.matrix(data)



# Create a directed graph from the adjacency matrix
G <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)

# Filter edges with weight < 0.05
G <- delete_edges(G, E(G)[weight < 0.05])

# Set node colors and sizes
V(G)$color <- "lightblue"
V(G)$size <- (degree(G, mode = "out") + 1) * 10

# Set edge width based on weight
E(G)$width <- E(G)$weight * 10

# Set arrow size and edge curve
E(G)$arrow.size <- 1.5
E(G)$curve <- 0.2


# Choose a layout
layout_choice <- layout_with_fr(G)  # Try layout_with_kk, layout_in_circle, layout_with_graphopt

# Plot the network
plot(G, 
     layout = layout_choice,
     edge.color = "skyblue", 
     edge.arrow.size = 0.8, 
     edge.arrow.width = 0.8,
     vertex.label.color = "black", 
     vertex.label.cex = 1.2, 
     vertex.frame.color = "black",
     main = "Cell Type Network", 
     margin = c(0, 0, 0, 0))


# Highlight edges coming into "K6+ cells" with red
edge_colors <- ifelse(ends(G, E(G))[, 2] == "K6+ cells", "red", "skyblue")

# Choose a layout
layout_choice <- layout_with_fr(G)  # Try layout_with_kk, layout_in_circle, layout_with_graphopt

pdf("descriptive_results/rna_trajecotry/plot/moscot_network.pdf")
# Plot the network
plot(G, 
     layout = layout_choice,
     edge.color = edge_colors, 
     edge.arrow.size = 2, 
     edge.arrow.width = 0.8,
     vertex.label.color = "black", 
     vertex.label.cex = 1.2, 
     vertex.frame.color = "black",
     main = "Cell Type Network", 
     margin = c(0, 0, 0, 0))
dev.off()
