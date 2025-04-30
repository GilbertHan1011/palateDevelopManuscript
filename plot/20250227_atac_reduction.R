library(Seurat)

atac = readRDS("process/framework/obj/combined_E12_E14_signac.Rds")
cluster <- read.csv("process/framework/cluster/20241211_knn_label.csv",row.names = 1)
atac$knn_label <- cluster$knn_label
spectral_umap <- read.csv("process/framework/reduction/20241222_snapatac_umap.csv",row.names = 1)
spectral_drawgraph <- read.csv("process/framework/reduction/20250118_snap_harmony_drawgraph.csv",row.names = 1)

colnames(spectral_umap) <- c(1,2)
colnames(spectral_drawgraph) <- c(1,2)


prefix_map <- c(
  "E125_2:"="C2_",
  "E125_1:" = "C1_",
  "E145_2:"= "B2_",
  "E145_1:" = "B1_"
  # Add more mappings if needed
)

newnames <- rownames(spectral_umap)
for (new_prefix in names(prefix_map)) {
  old_prefix <- prefix_map[new_prefix]
  newnames <- gsub( new_prefix,paste0(old_prefix), newnames)
}

# Assign the new row names
rownames(spectral_umap) <- newnames

spectral_umap <- spectral_umap[colnames(atac),]
spectral_drawgraph <- spectral_drawgraph[colnames(atac),]
atac@reductions$spectral_umap <-  CreateDimReducObject(embeddings = as.matrix(spectral_umap),key = "umapSpec_", assay = "ATAC")
atac@reductions$spectral_drawgraph <-  CreateDimReducObject(embeddings = as.matrix(spectral_drawgraph),key = "faSpec_", assay = "ATAC")

rownames(atac@reductions$spectral_umap@cell.embeddings) <- rownames(atac@reductions$umap_after@cell.embeddings)

colorPalate3 <- c("#282828","#F2C9D5", "#B43E44", "#8B6D9C", "#496496", "#904869","#FADF92" )
#colorPalate4 <- colorRampPalette(colorPalate3)(7)
DimPlot(atac,reduction = "spectral_umap",group.by = "knn_label",cols = colorPalate3)
ggsave("descriptive_results/rna_atac_integration/20250227_snapatac_spectral.pdf",width = 8,height = 6)


DimPlot(atac,reduction = "spectral_drawgraph",group.by = "knn_label",cols = colorPalate3)
ggsave("descriptive_results/rna_atac_integration/20250227_snapatac_spectral_draw.pdf",width = 8,height = 6)



drawgraph_rna <- read.csv("../2024.4_scATAC/processed_data/embedding/RNA_scvi_20neibor_drawgraph.csv",row.names = 1)
colnames(drawgraph_rna) <- c(1,2)
rna@reductions$drawgraph <-  CreateDimReducObject(embeddings = as.matrix(drawgraph_rna),key = "fa_", assay = "ATAC")
DimPlot(rna,reduction = "drawgraph",group.by = "level1_anno",cols = colorPalate3)
ggsave("descriptive_results/rna_atac_integration/20250227_rna_draw.pdf",width = 8,height = 6)

DimPlot(rna,reduction = "scviumap",group.by = "level1_anno",cols = colorPalate3)
ggsave("descriptive_results/rna_atac_integration/20250227_rna_umap.pdf",width = 8,height = 6)

FeaturePlot(rna_full,"Krt5",reduction = "scvi_umap")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/rna_marker//20250225_Krt5_drawgraph.pdf",width = 6,height = 5)

FeaturePlot(rna_full,"Krt5",reduction = "scviumap")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("descriptive_results/rna_marker//20250225_Krt5_umap.pdf",width = 6,height = 5)

corrGene <- rna_full[c("Krt5","Krt14","Krt6a"),]@assays$MAGIC_SCT@data%>%t()%>%as.data.frame()
fullct <- rna_full$level1_anno
day <- rna_full$group
ggplot(corrGene,aes(x = Krt5, y =Krt14,color = fullct))+
  geom_point()+theme_bw()

# First, combine your data with annotation info
plot_data <- corrGene %>%
  cbind(fullct, day) %>%
  # Convert day to factor if it's not already
  mutate(day = factor(day, levels = unique(day)))

# Create the beautiful plot
ggplot(plot_data, aes(x = Krt5, y = Krt14, color = fullct)) +
  geom_point(alpha = 0.7, size = 1) +
  facet_wrap(~ day, scales = "free") +
  scale_color_manual(values = colorPalate3, name = "Cell Type") +
  labs(
    title = "Correlation between Krt5 and Krt14 expression",
    subtitle = "Split by developmental stage",
    x = expression(italic("Krt5") ~ "expression"),
    y = expression(italic("Krt14") ~ "expression")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "#F0F0F0"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10)
  )
ggsave("descriptive_results/rna_marker//20250227_scatter_Krt14_K5.pdf",width = 10,height = 4)


my_theme <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "#F0F0F0"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10)
  )
# Plot 1: Krt5 vs Krt14
p1 <- ggplot(plot_data, aes(x = Krt5, y = Krt14, color = fullct)) +
  geom_point(alpha = 0.7, size = 1) +
  facet_wrap(~ day, scales = "free") +
  scale_color_manual(values = colorPalate3, name = "Cell Type") +
  labs(
    title = expression(italic("Krt5") ~ "vs" ~ italic("Krt14") ~ "expression"),
    x = expression(italic("Krt5") ~ "expression"),
    y = expression(italic("Krt14") ~ "expression")
  ) +
  my_theme

# Plot 2: Krt6a vs Krt14
p2 <- ggplot(plot_data, aes(x = Krt6a, y = Krt14, color = fullct)) +
  geom_point(alpha = 0.7, size = 1) +
  facet_wrap(~ day, scales = "free") +
  scale_color_manual(values = colorPalate3, name = "Cell Type") +
  labs(
    title = expression(italic("Krt6a") ~ "vs" ~ italic("Krt14") ~ "expression"),
    x = expression(italic("Krt6a") ~ "expression"),
    y = expression(italic("Krt14") ~ "expression")
  ) +
  my_theme

# Plot 3: Krt6a vs Krt5
p3 <- ggplot(plot_data, aes(x = Krt6a, y = Krt5, color = fullct)) +
  geom_point(alpha = 0.7, size = 1) +
  facet_wrap(~ day, scales = "free") +
  scale_color_manual(values = colorPalate3, name = "Cell Type") +
  labs(
    title = expression(italic("Krt6a") ~ "vs" ~ italic("Krt5") ~ "expression"),
    x = expression(italic("Krt6a") ~ "expression"),
    y = expression(italic("Krt5") ~ "expression")
  ) +
  my_theme

# To display the plots one after another
p1

p2
ggsave("descriptive_results/rna_marker//20250227_scatter_Krt6_K14.pdf",width = 10,height = 4)

p3
ggsave("descriptive_results/rna_marker//20250227_scatter_Krt6_K5.pdf",width = 10,height = 4)

# Alternatively, arrange them in a grid (requires the patchwork package)
# install.packages("patchwork")
library(patchwork)
combined_plot <- (p1 / p2 / p3) + 
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'Relationships between Keratin gene expressions',
    subtitle = 'Faceted by developmental stage',
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

combined_plot
ggsave("descriptive_results/rna_marker//20250227_scatter_Krt6_K5_K14.pdf",width = 8,height = 15)



E12corrGene <- E12[c("Krt5","Krt14","Krt6a"),]@assays$MAGIC_SCT@data%>%t()%>%as.data.frame()
E12corrGene$ident = E12$celltype
E14corrGene <- E14[c("Krt5","Krt14","Krt6a"),]@assays$MAGIC_SCT@data%>%t()%>%as.data.frame()
E14corrGene$ident = E14$celltype
ggplot(E12corrGene,aes(x = Krt5, y =Krt14,color = ident))+
  geom_point()+theme_bw()
ggsave("result/expression/E12_K5K14_scatter.pdf",width = 6,height = 4)

ggplot(E14corrGene,aes(x = Krt5, y =Krt14,color = ident))+
  geom_point()+theme_bw()
ggsave("result/expression/E14_K5K14_scatter.pdf",width = 6,height = 4)


