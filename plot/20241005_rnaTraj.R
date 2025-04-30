pseudotime1 <- read.csv("../2024.4_scATAC/processed_data/trajectory/5.1_rna_celloracle.pseudotime_run2.csv",row.names = 1)
dietPalate$pseudo <- pseudotime1$Pseudotime
FeaturePlot(dietPalate,"pseudo",label = T)+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/rna_trajecotry/plot/20241005_pseudotime.pdf",width = 6,height = 6)


pseudotime2 <- read.csv("../2024.4_scATAC/processed_data/trajectory/5.1_rna_celloracle.pseudotime_run2.csv",row.names = 1)