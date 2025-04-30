```R
## GREAT analysis----

my_theme <- function()

{

theme(

plot.title = element_text(size = 20),

axis.line.x = element_line(color = "black"),

axis.line.y = element_line(color = "black"),

axis.title.x = element_text(size = 20),

axis.title.y = element_text(size = 20),

axis.text.x = element_text(size = 15, color = "black"),

axis.text.y = element_text(size = 15, color = "black"),

panel.background = element_blank(),

panel.grid.minor = element_blank(),

panel.grid.major = element_blank(),

legend.text = element_text(size = 15),

legend.key = element_rect(fill = "white"))

}

  

library(rGREAT)

res = great(gr, "GO:BP", "txdb:mm10")

res

grName <- table(gr@ranges@NAMES)%>%as.data.frame()

ggplot(grName,aes(x = Var1, y = Freq))+

geom_bar(stat = "identity")+theme_bw()+my_theme()

  

grNameTotal <- table(matchFile@rowRanges@ranges@NAMES)%>%as.data.frame()

grName$frac <- grName$Freq/grNameTotal$Freq

  

tb = getEnrichmentTable(res)

head(tb)

DotPlot(tb)

write.csv(tb,"processed_data/Arid3a_centric/archr_cisbp_enriched_rgreat_bp.csv")

plotVolcano(res)

pdf("result/arid3a_centric/5.5_tss_locate.pdf",width = 8,height = 4)

plotRegionGeneAssociations(res)

dev.off()

tss_distance <- getRegionGeneAssociations(res)

plotRegionGeneAssociations(res, term_id = "hemopoiesis")

shinyReport(res)

library(simplifyEnrichment)

sig_go_ids = tb$id[tb$p_adjust < 0.001]

cl = simplifyGO(sig_go_ids)
```