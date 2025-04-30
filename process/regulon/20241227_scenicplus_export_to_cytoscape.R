regulonDirect <- read.table("process//20241219_scenicplus_outs//eRegulon_direct.tsv",header = T)
regulonExtend <- read.table("process/20241219_scenicplus_outs/eRegulons_extended.tsv",header = T)
regulon <- rbind(regulonDirect,regulonExtend)

tf <- read.csv("process/20241219_scenicplus_outs/20241222_regulon/20241227_select_regulon.csv",row.names = 1)

E14_reg <- regulon %>% filter(Gene_signature_name %in% tf$K14) %>% 
  filter(triplet_rank < 8000)
write.csv(E14_reg,"process/20241219_scenicplus_outs/20241222_regulon/20241227_K14_reg1_to_cytoscape.csv",quote = F)
