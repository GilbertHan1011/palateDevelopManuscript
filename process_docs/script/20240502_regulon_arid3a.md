```R
arid3aGene <- read.csv("processed_data/regulon/celloracle/5.1_arid3a_reg_gene.csv",row.names = 1)

arid3aGene$binaryp <- arid3aGene$p<0.001

arid3aGene$score <- arid3aGene$binaryp * arid3aGene$coef_mean

  
  

top_5_genes <- arid3aGene %>%

dplyr::filter(key %in% unique(c("K6+ cells","K14(-)","stem cells")))%>%

group_by(key) %>%

arrange(desc(score)) %>%

slice(1:10) %>%

ungroup()

plotGene <- arid3aGene%>%dplyr::filter(target %in% unique(top_5_genes$target))%>%

dplyr::filter(key %in% unique(c("stem cells","K6+ cells","K14(-)")))

  
  
  

ggplot(plotGene)+

geom_bar(aes(x = target, fill = factor(key), y =score), stat = "identity",position = "dodge")+

theme_bw()+theme(

axis.text.x = element_text(angle = 45, hjust = 1, size = 12),

axis.title.x = element_blank(),

legend.title = element_text(size = 14),

legend.text = element_text(size = 12)

)+ggtitle("Genes which Arid3a regulate")

ggsave("result/regulon/celloracle/5.2_K14_K6_compare.pdf",width = 10,height = 4)

  
  

top_5_genes <- arid3aGene %>%

group_by(key) %>%

arrange(desc(score)) %>%

slice(1:5) %>%

ungroup()

plotGene <- arid3aGene%>%dplyr::filter(target %in% unique(top_5_genes$target))

  

ggplot(plotGene)+

geom_bar(aes(x = target, fill = factor(key), y =score), stat = "identity",position = "dodge")+

theme_bw()+theme(

axis.text.x = element_text(angle = 45, hjust = 1, size = 12),

axis.title.x = element_blank(),

legend.title = element_text(size = 14),

legend.text = element_text(size = 12)

)+ggtitle("Genes which Arid3a regulate")

  

ggsave("result/regulon/celloracle/5.2_all_compare.pdf",width = 10,height = 4)

  
  
  

arid3aTarget <- read.csv("processed_data/regulon/celloracle/5.1_arid3a_target.csv",row.names = 1)

arid3aTarget$binaryp <- arid3aTarget$p<0.001

arid3aTarget$score <- arid3aTarget$binaryp * arid3aTarget$coef_mean

  
  

top_5_genes <- arid3aTarget %>%

group_by(key) %>%

arrange(desc(score)) %>%

slice(1:4) %>%

ungroup()

plotGene <- arid3aTarget%>%dplyr::filter(source %in% unique(top_5_genes$source))

  
  

ggplot(plotGene)+

geom_bar(aes(x = source, fill = factor(key), y =score), stat = "identity",position = "dodge")+

theme_bw()+theme(

axis.text.x = element_text(angle = 45, hjust = 1, size = 12),

axis.title.x = element_blank(),

legend.title = element_text(size = 14),

legend.text = element_text(size = 12)

)+ggtitle("Regulon which regulate Arid3a")

ggsave("result/regulon/celloracle/5.2_all_compare_target.pdf",width = 10,height = 4)

  
  

top_5_genes <- arid3aTarget %>%

group_by(key) %>%

dplyr::filter(key %in% unique(c("K6+ cells","K14(-)","stem cells")))%>%

arrange(desc(score)) %>%

slice(1:10) %>%

ungroup()

plotGene <- arid3aTarget%>%dplyr::filter(source %in% unique(top_5_genes$source))%>%

dplyr::filter(key %in% unique(c("K6+ cells","K14(-)","stem cells")))

  
  

ggplot(plotGene)+

geom_bar(aes(x = source, fill = factor(key), y =score), stat = "identity",position = "dodge")+

theme_bw()+theme(

axis.text.x = element_text(angle = 45, hjust = 1, size = 12),

axis.title.x = element_blank(),

legend.title = element_text(size = 14),

legend.text = element_text(size = 12)

)+ggtitle("Regulon which regulate Arid3a in stem cells, K14(-) and K6(+)")

ggsave("result/regulon/celloracle/5.2_K14_compare_target.pdf",width = 10,height = 4)
```