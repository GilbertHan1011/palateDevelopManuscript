#rm(list=ls())
library(dplyr)
library(ggplot2)
##' input
dietPalate <- readRDS("process/framework/obj/RNA_E12E14_diet.Rds")
# dietPalate@reductions$umap <- dietPalate@reductions$scviumap
#
#
# markersAll <- FindAllMarkers(dietPalate,only.pos = T)
# FeaturePlot(dietPalate,"Pitx1")
# FeaturePlot(dietPalate,"Krt5")

#== rename ident--
Idents(dietPalate) <- dietPalate$level1_anno
newid <- c("stem cells", "Transit", "K6+ periderm", "K5(-)", "Shh(+)", "Mature K5(-) ",
           "K14(-)")
names(newid) <- c("stem cells", "Transit", "K6+ cells", "K5(-)", "Shh(+)", "Mature K5(-) ",
                  "K14(-)")
dietPalate <- RenameIdents(dietPalate,newid)
dietPalate$level1_anno <- Idents(dietPalate)

##' output


#== plot figure-----------------
# top_feature <- markersAll %>%  group_by(cluster) %>% arrange(p_val) %>% slice_head(n = 5)
# features <- unique(top_feature$gene)



# Combined list of all genes
all_genes <- c("Cxcl14", "Wnt3", "Wnt7b", "Irx3",
               "Krt5", "Krt14", "Gabrp", "Grhl3",
               "Krt6a", "Gabrp", "Cldn3", "Cldn23",
               "Casz1", "Sp6", "Fgf4", "Shh",
               "Pitx2", "Irx1", "Dsc3", "Fgf9",
               "Barx1", "Foxe1", "Sox2", "Col14a1",
               "Sim2", "Sox21", "Barx1", "Prox1")
all_genes <- unique(all_genes)
dp <- DotPlot(dietPalate, all_genes)
# DotPlot(dietPalate, features = features,cols = c("white","red"))+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),title = NULL)

dp_data <- dp$data
#dp_data_small <- dp_data[c(1:20),]
#devtools::install_github("Simon-Leonard/FlexDotPlot")
library(FlexDotPlot)
library(RColorBrewer)
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(17)

# dot_plot(dp_data[,c(3,4,1,2,5)], size_var = "pct.exp", col_var = "avg.exp.scaled",
#          size_legend = "Percent Expressed", col_legend = "Average Expression",
#          x.lab.pos = "bottom", display_max_sizes = FALSE)
# do_small <- dot_plot(dp_data_small[,c(3,4,1,2,5)], size_var = "pct.exp", col_var = "avg.exp.scaled", cols.use = rev(mycolor),
#          size_legend = "Percent Expressed", col_legend = "Average Expression",
#          x.lab.pos = "bottom", display_max_sizes = FALSE,do.return = T)
#do_small
#pdf("descriptive_results/rna_annotation_plot/20240930_marker_dotplot_wide.pdf",width = 10,height = 6)
dotplot=dot_plot(dp_data[,c(3,4,1,2,5)], shape_var = "pct.exp", col_var = "avg.exp.scaled", cols.use = rev(mycolor),
                 shape_legend = "Percent Expressed", col_legend = "Average Expression", x.lab.pos = "bottom",
                 dend_x_var = c("pct.exp","avg.exp.scaled"), dend_y_var = c("pct.exp","avg.exp.scaled"),
                 hclust_method = "ward.D2", do.return=T)

#dev.off()

pdf("descriptive_results/rna_annotation_plot/20250403_marker_dotplot_wide.pdf",width = 10,height = 6)
dotplot=dot_plot(dp_data[,c(3,4,1,2,5)], shape_var = "pct.exp",
                 col_var = "avg.exp.scaled", cols.use = rev(mycolor),
                 shape_legend = "Percent Expressed", col_legend = "Average Expression", x.lab.pos = "bottom",
                 dend_x_var = NULL,
                 dend_y_var = c("pct.exp","avg.exp.scaled"),
                 hclust_method = "ward.D2", do.return=T)
dev.off()
# pdf("descriptive_results/rna_annotation_plot/20240930_marker_dotplot_long.pdf",width = 6,height = 8)
# dotplot2=dot_plot(dp_data[,c(4,3,1,2,5)], shape_var = "pct.exp", col_var = "avg.exp.scaled",cols.use = rev(mycolor),
#                   shape_legend = "Percent Expressed", col_legend = "Average Expression", x.lab.pos = "bottom",
#                   dend_x_var = c("pct.exp","avg.exp.scaled"), dend_y_var = c("pct.exp","avg.exp.scaled"),
#                   hclust_method = "ward.D2", do.return=T,do.plot = T)
# dev.off()
#dotplot2+scale_color_viridis_c()



#== add heatmap--------------------------
gam <- readRDS("../2024.4_scATAC/processed_data/trajectory/20250107_gam_threelineage.Rds")
gene_clusters = read.csv("descriptive_results/rna_trajectory_gene/20250117_gene_clusters.csv")

scaled_mat= read.csv("descriptive_results/rna_trajectory_gene/20250117_gene_scaled_mat.csv")
#scaled_mat
hmAll2 <- readRDS("process/trajectory/20250117_gene_traj_heatmap.Rds")

hmAll2@annotation_legend_param


column_split <- rep(c("stem","K14"), each = 100)
column_split <- factor(column_split,levels = c("stem","K14"))

hmAll <- Heatmap(scaled_mat,
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 show_row_names = FALSE,
                 row_split = gene_clusters$Cluster,  # Assuming you have this from k-means
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
