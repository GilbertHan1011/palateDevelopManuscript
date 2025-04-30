
```R
library(dplyr)

library(ggplot2)

##' input

dietPalate <- readRDS("process/framework/obj/RNA_E12E14_diet.Rds")

dietPalate@reductions$umap <- dietPalate@reductions$scviumap

  

Idents(dietPalate) <- dietPalate$level1_anno

markersAll <- FindAllMarkers(dietPalate,only.pos = T)

FeaturePlot(dietPalate,"Pitx1")

FeaturePlot(dietPalate,"Krt5")

  

#== rename ident--

newid <- c("stem cells", "Transit", "K6+ periderm", "K5(-)", "Shh(+)", "Mature K5(-) ",

"K14(-)")

names(newid) <- c("stem cells", "Transit", "K6+ cells", "K5(-)", "Shh(+)", "Mature K5(-) ",

"K14(-)")

dietPalate <- RenameIdents(dietPalate,newid)

dietPalate$level1_anno <- Idents(dietPalate)

  

##' output

write.csv(markersAll,"descriptive_results/rna_marker/20240930_level1_anno.csv")

  
  

#== plot figure-----------------

top_feature <- markersAll %>% group_by(cluster) %>% arrange(p_val) %>% slice_head(n = 5)

features <- unique(top_feature$gene)

dp <- DotPlot(dietPalate, features)

  

DotPlot(dietPalate, features = features,cols = c("white","red"))+

theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),title = NULL)

  

dp_data <- dp$data

#dp_data_small <- dp_data[c(1:20),]

#devtools::install_github("Simon-Leonard/FlexDotPlot")

library(FlexDotPlot)

library(RColorBrewer)

mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(17)

  

# dot_plot(dp_data[,c(3,4,1,2,5)], size_var = "pct.exp", col_var = "avg.exp.scaled",

# size_legend = "Percent Expressed", col_legend = "Average Expression",

# x.lab.pos = "bottom", display_max_sizes = FALSE)

# do_small <- dot_plot(dp_data_small[,c(3,4,1,2,5)], size_var = "pct.exp", col_var = "avg.exp.scaled", cols.use = rev(mycolor),

# size_legend = "Percent Expressed", col_legend = "Average Expression",

# x.lab.pos = "bottom", display_max_sizes = FALSE,do.return = T)

#do_small

pdf("descriptive_results/rna_annotation_plot/20240930_marker_dotplot_wide.pdf",width = 10,height = 6)

dotplot=dot_plot(dp_data[,c(3,4,1,2,5)], shape_var = "pct.exp", col_var = "avg.exp.scaled", cols.use = rev(mycolor),

shape_legend = "Percent Expressed", col_legend = "Average Expression", x.lab.pos = "bottom",

dend_x_var = c("pct.exp","avg.exp.scaled"), dend_y_var = c("pct.exp","avg.exp.scaled"),

hclust_method = "ward.D2", do.return=T)

  

dev.off()

  

pdf("descriptive_results/rna_annotation_plot/20240930_marker_dotplot_long.pdf",width = 6,height = 8)

dotplot2=dot_plot(dp_data[,c(4,3,1,2,5)], shape_var = "pct.exp", col_var = "avg.exp.scaled",cols.use = rev(mycolor),

shape_legend = "Percent Expressed", col_legend = "Average Expression", x.lab.pos = "bottom",

dend_x_var = c("pct.exp","avg.exp.scaled"), dend_y_var = c("pct.exp","avg.exp.scaled"),

hclust_method = "ward.D2", do.return=T,do.plot = T)

dev.off()

#dotplot2+scale_color_viridis_c()
```