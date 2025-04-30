library(ComplexHeatmap)
gene_clusters <- read.csv("descriptive_results/rna_trajectory_gene/10.21_gene_clusters.csv")
selected_gene <- intersect(gene_clusters$Gene, gene_clusters_score$Gene)
hmAll_score
cluster_order <- c(5,3,1,4,2)
sigGene2_ordered_gene <- factor(sigGene, levels = cluster_order)
levels(sigGene2_ordered_gene) <- rename_vector[levels(sigGene2_ordered_gene)]
gene_select <- c(
  "Sprr1a",
  "Krtdap",
  "Krt16",
  "Sprr1b",
  "Sprr4",
  "Klf2",
  "Snai2",
  "Vgll2",
  "Krt6a"
)
gene_selected <- intersect(gene_select, selected_gene)

position <- c()
for (i in gene_selected){
  position <- c(position, which(is.element(rownames(scaled_mat[selected_gene,]), i) == TRUE))
}

har <- rowAnnotation(link = anno_mark(at = position, labels = gene_selected,
                                      labels_gp = gpar(fontsize = 9, fontface = "bold"), link_width = unit(1.5, "cm")))

hmAll_score2 <- Heatmap(scaled_mat[selected_gene,],
                       cluster_columns = FALSE,
                       show_column_names = FALSE,
                       show_row_names = FALSE,
                       cluster_row_slices = F,
                       show_row_dend = F,
                       column_split = column_split,  # Splitting into 3 groups of 100
                       col = colorRamp2(c(-2, 0, 4), c("Deepskyblue3", "white", "red")),
                       border = TRUE,
                       na_col = "grey",
                       split = 2,
                       row_title_gp = gpar(fontsize = 10),
                       right_annotation = har,
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

hmAll_score2 <- draw(hmAll_score2)
pdf("descriptive_results/atac_trajectory/20250127_genescore_heatmap.pdf",width = 8,height = 6)
hmAll_score2
dev.off()



#== plot final plot of scATAC-------------

selectedMotif <- read.csv("process/trajectory/20250127_motifVarDev.csv",row.names = 1)
selectedMotif <- selectedMotif[selectedMotif$rank<=60,]$name
scaled_mat_motif <- read.csv("descriptive_results/rna_trajectory_gene/20250123_motif_scaled_mat.csv")
scaled_mat_motif <- hmAll2@ht_list$matrix_6@matrix

rownames(scaled_mat_motif) <- gsub("-","_",rownames(scaled_mat_motif))
selectedMotif <- intersect(selectedMotif,rownames(scaled_mat_motif))
hmAll_motif <- Heatmap(scaled_mat_motif[selectedMotif,],
                  cluster_columns = FALSE,
                  show_column_names = FALSE,
                  show_row_names = TRUE,
                  cluster_row_slices = F,
                  show_row_dend = F,
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

hmAll_motif <- draw(hmAll_motif)

mannual_select <- c("Arid3a_7",
                    "Cebpd_97",
                    "Cebpg_129",
                    "Cebpb_130",
                    "Mafb_136",
                    "Vax1_413",
                    "Trp53_692",
                    "Rel_697",
                    "Rela_698",
                    "Tcf7l2_732",
                    "Lef1_734",
                    "Tbx21_762",
                    "Tbx5_764",
                    "Mafa_792",
                    "Maf_793",
                    "Tead3_867",
                    "Tead4_868",
                    "Grhl1_390")

selectedMotif_final <- intersect(selectedMotif,mannual_select)


hmAll_motif2 <- Heatmap(scaled_mat_motif[mannual_select,],
                       cluster_columns = FALSE,
                       show_column_names = FALSE,
                       show_row_names = TRUE,
                       cluster_row_slices = F,
                       show_row_dend = F,
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
hmAll_motif2 <- draw(hmAll_motif2)

position <- c()
for (i in selectedMotif_final){
  position <- c(position, which(is.element(rownames(scaled_mat_motif[selectedMotif,]), i) == TRUE))
}

har_motif <- rowAnnotation(link = anno_mark(at = position, labels = selectedMotif_final,
                                      labels_gp = gpar(fontsize = 9, fontface = "bold"), link_width = unit(1.5, "cm")))




hmAll_motif3 <- Heatmap(scaled_mat_motif[selectedMotif,],
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        show_row_names = F,
                        cluster_row_slices = F,
                        show_row_dend = F,
                        column_split = column_split,  # Splitting into 3 groups of 100
                        col = colorRamp2(c(-2, 0, 4), c("Deepskyblue3", "white", "red")),
                        border = TRUE,
                        na_col = "grey",
                        row_title_gp = gpar(fontsize = 10),
                        row_title_rot = 0,  # Horizontal text
                        row_gap = unit(1, "mm"),  # Add some gap between row clusters
                        right_annotation = har_motif,
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

pdf("descriptive_results/atac_trajectory/20250127_motif_heatmap.pdf",width = 8,height = 6)
hmAll_motif3
dev.off()
#dev.off()