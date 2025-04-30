# Plot Figure 2
library(ComplexHeatmap)
library(circlize)
csvfile <- list.files("process/regulation/20250110_celloracle/20250110_pertubation_score/",full.names = T, pattern = "*.csv")
name <- list.files("process/regulation/20250110_celloracle/20250110_pertubation_score/",full.names = F, pattern = "*.csv") %>% gsub(".csv","",.)

scoreData=lapply(csvfile,read.csv)
names(scoreData) <- name
combineData <- lapply(scoreData,function(x){
  x$PS
})
combineData <- do.call(rbind,combineData)
colnames(combineData) <- scoreData[[1]]$X
rownames(combineData) <- name
labeled_gene <- c("Arid3a","Grhl3","Grhl1","Irf6","Ahr","Klf7","Zbtb7a",
                  "Tfap2c","Cebpb","Klf5",'Sox11','Fosl2',"Klf4","Tgif1",
                  "Cebpb","Klf6",'Jund', 'Jun', 'Atf3')
write.csv(combineData,"descriptive_results/regulon_pertubation/pertubation_score/20250225_pertubation_data.csv")

# Find positions of labeled genes
position <- which(rownames(combineData) %in% labeled_gene)

# Improve column names by cleaning up the labels
colnames(combineData) <- c("Stem Cells", "Transit", "Shh(+)", "K5(-)", "Mature K5(-)", "K14(-)", "K6(+)")

# Create row annotation with marks - improved styling
har <- rowAnnotation(
    link = anno_mark(
        at = position,
        labels = labeled_gene,
        labels_gp = gpar(fontsize = 10, fontface = "bold", family = "sans"),
        link_width = unit(1, "cm"),
        link_gp = gpar(lwd = 1, col = "black"),
        padding = unit(0.5, "mm")
    )
)

#combineData <- combineData[,c("stem cells", "Transit", "Shh(+)","K5(-)", "Mature K5(-) ", 
#              "K14(-)", "K6+ cells")]

# Create and draw the heatmap with enhanced aesthetics
pdf("descriptive_results/regulon_pertubation/pertubation_score/20250225_pertubation_mat.pdf", 
    width = 8, height = 6)

Heatmap(combineData,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    
    # Enhanced color scheme
    col = colorRamp2(c(-0.003, -0.0015, 0, 0.0015, 0.003), 
                     c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")),
    
    # Improved column styling
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    
    # Enhanced titles
    column_title = "Perturbation Score Across Cell Types",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    name = "Perturbation\nScore",
    
    # Better spacing and borders
    heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        title_position = "topcenter",
        legend_height = unit(3, "cm")
    ),
    
    # Add border
    border = TRUE,
    
    # Add right annotation
    right_annotation = har
)

dev.off()
