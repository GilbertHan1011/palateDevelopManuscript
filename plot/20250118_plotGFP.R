grep("GFP",rownames(seurat),value = T)
FeaturePlot(seurat,"GFP",reduction = "scviumap",max.cutoff =0.1)
sum(seurat@assays$RNA@counts["GFP",]>0)
seurat$gfp <- "Non-GFP"
seurat$gfp[seurat@assays$RNA@counts["GFP",]>0] <- "GFP"
DimPlot(seurat,group.by = "gfp",reduction = "scviumap",cols = c("red","Deepskyblue3"),)
ggsave("descriptive_results/rna_marker/20250118_gfp_umap.pdf",width = 5,height = 4)

metadf <- seurat@meta.data[c("gfp","level1_anno")]
plot_data <- as.data.frame(table(metadf$gfp, metadf$level1_anno))
colnames(plot_data) <- c("gfp", "level1_anno", "count")


# Calculate proportions and order
# Calculate proportions and order
# Calculate correct proportions and order
plot_data_ordered <- plot_data %>%
  group_by(level1_anno) %>%
  mutate(
    group_total = sum(count),
    gfp_prop = count[gfp == "GFP"] / group_total
  ) %>%
  ungroup()




plot_data_ordered$level1_anno <- factor(plot_data_ordered$level1_anno,levels = c("stem cells","K14(-)","K6+ cells", 
                                                                                 "Transit", "K5(-)", "Shh(+)" ,"Mature K5(-) "
                                                                      ))
ggplot(plot_data_ordered, 
       aes(x = level1_anno, y = count, fill = gfp)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Cell Types by GFP Status",
       x = "Cell Type",
       y = "Proportion",
       fill = "GFP Status") +
  theme(
    # Make text bold and larger
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               face = "bold", 
                               size = 12),
    axis.text.y = element_text(face = "bold", 
                               size = 12),
    axis.title.x = element_text(face = "bold", 
                                size = 14, 
                                margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", 
                                size = 14, 
                                margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, 
                              face = "bold", 
                              size = 16,
                              margin = margin(b = 10)),
    legend.title = element_text(face = "bold", 
                                size = 12),
    legend.text = element_text(face = "bold", 
                               size = 11)
  ) +
  scale_fill_manual(values = c("Non-GFP" = "gray", "GFP" = "lightgreen")) +
  scale_y_continuous(labels = scales::percent)
ggsave("descriptive_results/rna_marker/20250118_gfp_prop_E12E14.pdf",width = 6,height = 6)



#==E12------------------


metadf <- seurat@meta.data[c("group","gfp","level1_anno")]


# Create summary with group
plot_data <- as.data.frame(table(metadf$gfp, metadf$level1_anno, metadf$group))
colnames(plot_data) <- c("gfp", "level1_anno", "group", "count")

# Calculate proportions
plot_data_ordered <- plot_data %>%
  group_by(level1_anno, group) %>%
  mutate(
    group_total = sum(count),
    gfp_prop = count[gfp == "GFP"] / group_total
  ) %>%
  ungroup()

# Set factor levels
plot_data_ordered$level1_anno <- factor(plot_data_ordered$level1_anno,
                                        levels = c("stem cells","K14(-)","K6+ cells", 
                                                   "Transit", "K5(-)", "Shh(+)" ,"Mature K5(-) "))

# Option 1: Facet by group
ggplot(plot_data_ordered, 
       aes(x = group, y = count, fill = gfp)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Cell Types by GFP Status and Group",
       x = "Cell Type",
       y = "Proportion",
       fill = "GFP Status") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 12)  # Format facet labels
  ) +
  scale_fill_manual(values = c("Non-GFP" = "gray", "GFP" = "lightgreen")) +
  scale_y_continuous(labels = scales::percent)

ggsave("descriptive_results/rna_marker/20250118_gfp_daygroup.pdf",width = 4,height = 6)

# Define invalid combinations
invalid_combinations <- data.frame(
  group = "E145",
  level1_anno = c("stem cells", "K14(-)")
)

# Filter out invalid combinations
plot_data_filtered <- plot_data_ordered %>%
  anti_join(invalid_combinations, by = c("group", "level1_anno"))
# Create the plot with filtered data
ggplot(plot_data_filtered, 
       aes(x = level1_anno, y = count, fill = gfp)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~group) +
  theme_minimal() +
  labs(title = "Proportion of Cell Types by GFP Status",
       x = "Cell Type",
       y = "Proportion",
       fill = "GFP Status") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_fill_manual(values = c("Non-GFP" = "gray", "GFP" = "lightgreen")) +
  scale_y_continuous(labels = scales::percent)
ggsave("descriptive_results/rna_marker/20250118_gfp_prop_E12E14_sep.pdf",width = 6,height = 6)

