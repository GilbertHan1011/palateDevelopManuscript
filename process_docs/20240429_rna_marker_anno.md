
```R
palate <- readRDS("data/palate_epi_MAGIC_periderm_MAGIC_harmony_renamed.rds")

  

FeaturePlot(palateE14,c("Krt5","Krt14"))

E12 <- palate[,palate$group=="E125"]

E14 <- palate[,palate$group=="E145"]


FeaturePlot(palate,c("Krt5","Krt14","Krt6a"))

p1 <- FeaturePlot(E12,c("Krt5","Krt14","Krt6a"),ncol = 3)&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

p2 <- FeaturePlot(E14,c("Krt5","Krt14","Krt6a"),ncol = 3)&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

  

ggsave(plot = p1/p2,filename = "result/expression/featureplot_E12E14.pdf",width = 10,height = 8)

  

FeaturePlot(palateE14,c("Krt5","Krt14","Krt6a"),ncol = 3)&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

  

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

  
  

ggplot(E14corrGene,aes(x = Krt5, y =Krt6a,color = ident))+

geom_point()+theme_bw()

ggsave("result/expression/E14_K5K6_scatter.pdf",width = 6,height = 4)

  
  

ggplot(E12corrGene,aes(x = Krt5, y =Krt6a,color = ident))+

geom_point()+theme_bw()

ggsave("result/expression/E12_K5K6_scatter.pdf",width = 6,height = 4)

  

ggplot(E12corrGene,aes(x = Krt14, y =Krt6a,color = ident))+

geom_point()+theme_bw()

ggsave("result/expression/E12_K14K6_scatter.pdf",width = 6,height = 4)

  

ggplot(E14corrGene,aes(x = Krt14, y =Krt6a,color = ident))+

geom_point()+theme_bw()

ggsave("result/expression/E14_K14K6_scatter.pdf",width = 6,height = 4)

  
  

allCorrGene <- palate[c("Krt5","Krt14","Krt6a"),]@assays$MAGIC_SCT@data%>%t()%>%as.data.frame()

allCorrGene$ident = palate$celltype

ggplot(allCorrGene,aes(x = Krt14, y =Krt6a,color = ident))+

geom_point()+theme_bw()

p1 <- ggplot(E12corrGene,aes(x = Krt5, y =Krt14,color = ident))+

geom_point()+theme_bw()+ggtitle("E12")

p2 <- ggplot(E14corrGene,aes(x = Krt5, y =Krt14,color = ident))+

geom_point()+theme_bw()+ggtitle("E14")

p3 <- ggplot(allCorrGene,aes(x = Krt5, y =Krt14,color = ident))+

geom_point()+theme_bw()+ggtitle("E12_E14")

p1| p2| p3

ggsave("result/expression/All_K5K14_scatter.pdf",width = 10,height = 3)

  
  

p1 <- ggplot(E12corrGene,aes(x = Krt5, y =Krt6a,color = ident))+

geom_point()+theme_bw()+ggtitle("E12")

p2 <- ggplot(E14corrGene,aes(x = Krt5, y =Krt6a,color = ident))+

geom_point()+theme_bw()+ggtitle("E14")

p3 <- ggplot(allCorrGene,aes(x = Krt5, y =Krt6a,color = ident))+

geom_point()+theme_bw()+ggtitle("E12_E14")

p1| p2| p3

ggsave("result/expression/All_K5K6_scatter.pdf",width = 10,height = 3)

  
  

p1 <- ggplot(E12corrGene,aes(x = Krt14, y =Krt6a,color = ident))+

geom_point()+theme_bw()+ggtitle("E12")

p2 <- ggplot(E14corrGene,aes(x = Krt14, y =Krt6a,color = ident))+

geom_point()+theme_bw()+ggtitle("E14")

p3 <- ggplot(allCorrGene,aes(x = Krt14, y =Krt6a,color = ident))+

geom_point()+theme_bw()+ggtitle("E12_E14")

p1| p2| p3

ggsave("result/expression/All_K14K6_scatter.pdf",width = 10,height = 3)

  

scviUmap <- read.csv("processed_data/embedding/RNA_scvi_umap.csv",row.names = 1)

  

colnames(scviUmap) <- c("umap_1","umap_2")

rownames(scviUmap) <- colnames(palate)

palate@reductions$df = CreateDimReducObject(embeddings = as.matrix(scviUmap), key = "umap_")

  

p1 <- FeaturePlot(E12,c("Krt5","Krt14","Krt6a"),ncol = 3,reduction = "df")&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

p2 <- FeaturePlot(E14,c("Krt5","Krt14","Krt6a"),ncol = 3,,reduction = "df")&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

p3 <- FeaturePlot(palate,c("Krt5","Krt14","Krt6a"),ncol = 3,,reduction = "df")&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

p1/p2/p3

#ggsave("result/expression/E14_K14K6_scatter.pdf",width = 6,height = 4)

  

FeaturePlot(palate,c("Shh","Tgfb2","Sp6","Cd44","Barx2","Sstr2"),ncol = 3,,reduction = "df")&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

  
  

FeaturePlot(palate,c("Shh","Tgfb2","Sp6","Cd44","Barx2","Sstr2"),ncol = 3,,reduction = "df")&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

  
  

FeaturePlot(palate,c("Sgk1","Ephb2","Pthlh","Krt6a","Tagln","Wnt11","Fgf21","Krt16","Gsg1l"),ncol = 3,,reduction = "df")&

scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))

  
  

leiden_rna <- read.csv("processed_data/cluster/rna_leiden_scvi.csv",row.names = 1)

  

palate@meta.data[colnames(leiden_rna)] <- leiden_rna

DimPlot(palate,group.by = c("leiden_0.5","group"),reduction = "df",label = T)

palate$level1_anno <- palate$leiden_0.5

Idents(palate) <- palate$leiden_0.5

newid <- c("K5- cells","stem cells","stem cells","transitional cells",

"K6+ cells","transitional cells","transitional cells","transitional cells","stem cells")

names(newid) <- c(0:8)%>%as.character()

palate <- RenameIdents(palate,newid)

palate$level1_anno <- Idents(palate)

  
  
  
  

DimPlot(palate)

  

DimPlot(palate,group.by = c("leiden_0.5","leiden","group"),reduction = "df",label = T)

DimPlot(palate,group.by = c("leiden_2"),reduction = "df",label = T)

palate$level2_anno <- palate$leiden_0.5

Idents(palate) <- palate$leiden_0.5

newid <- c("K5- cells","stem cells 1","stem cells 2","transitional cells 1",

"K6+ cells","transitional cells 2","transitional cells 3","transitional cells 4","stem cells 3")

names(newid) <- c(0:8)%>%as.character()

palate <- RenameIdents(palate,newid)

DimPlot(palate,reduction = "df")

palate$level2_anno <- Idents(palate)%>%as.character()

k14neg <- colnames(palate)[palate$leiden_2==16]

K6pos1 <- colnames(palate)[palate$leiden_2==31]

K6pos2 <- colnames(palate)[palate$leiden_2==15]

shhpos <- colnames(palate)[palate$leiden_2==24]

k5negMature <- colnames(palate)[palate$leiden_2==3]

  

palate$level2_anno[k14neg] <- "K14- cells"

palate$level2_anno[K6pos1] <- "K6+ cells 1"

palate$level2_anno[K6pos2] <- "K6+ cells 2"

palate$level2_anno[shhpos] <- "Shh+ cells"

palate$level2_anno[k5negMature] <- "Mature K5- cells "

DimPlot(palate,reduction = "df",group.by = "level2_anno")

  

palate$level1_anno <- as.character(palate$level1_anno)

palate$level1_anno[k14neg] <- "K14- cells"

palate$level1_anno[ colnames(palate)[palate$leiden==12]] <- "transitional cells"

  

palate$level1_anno[shhpos] <- "Shh+ cells"

palate$level1_anno[k5negMature] <- "Mature K5- cells "

DimPlot(palate,reduction = "df",group.by = c("level1_anno","level2_anno"))

  

DimPlot(palate,group.by = c("level1_anno","level2_anno"))

write.csv(palate@meta.data[c("level1_anno","level2_anno")],"processed_data/cluster/rna_annotation_scvi.csv")
```