regulon = read.csv("result/Arid3a/regulon/8.25_regulon_short.csv",row.names = NULL)
regulon_filter = regulon[regulon$X=="Arid3a",]
regulon_filter <- regulon_filter %>% filter(Enrichment.6!="[]")
strName = regulon_filter$Enrichment.6[3]
#library(jsonlite)


strFun <- function(strName){
  require(jsonlite)
  str_clean <- gsub("\\(", "[", strName)
  str_clean <- gsub("\\)", "]", str_clean)
  str_clean <- gsub("'", "\"", str_clean)
  
  # Convert the cleaned string into a list
  data_list <- fromJSON(paste0("[", str_clean, "]"))
  
  # Convert the list into a data frame
  df <- as.data.frame(data_list[1,,], stringsAsFactors = FALSE)
  df$V2 <- as.numeric(df$V2)
  return(df)
}

dfList <- lapply(regulon_filter$Enrichment.6,strFun)

names(dfList) <- paste0(regulon_filter$X.1,regulon_filter$Enrichment.5)

arid3a_1 <- dfList[[4]]
arid3a_1$source <- "Arid3a"
colnames(arid3a_1) <- c("target","importance","source")

write.csv(arid3a_1,"result/Arid3a/regulon/1_arid3a.csv",quote = F)



arid3a_5 <- dfList$`Arid3a_5_gluefrozenset({'top50', 'activating', 'input_glue.genes_vs_tracks.rankings'})`
processDf <- function(df){
  df$source <- "Arid3a"
  colnames(df) <- c("target","importance","source")
  return(df)
}
arid3a_5 <- processDf(arid3a_5)

arid3a_3 <- dfList$`Arid3a_3_gluefrozenset({'activating', 'top10perTarget', 'input_glue.genes_vs_tracks.rankings'})`

arid3a_3 <- processDf(arid3a_3)

arid3a_1_supp <- dfList$`Arid3a_1_suppfrozenset({'input_feather_test_short.genes_vs_tracks.rankings', 'top5perTarget', 'activating'})`
arid3a_1_supp <- processDf(arid3a_1_supp)

arid3a_3_supp <- dfList$`Arid3a_1_suppfrozenset({'input_feather_test_short.genes_vs_tracks.rankings', 'top5perTarget', 'activating'})`
arid3a_3_supp <- processDf(arid3a_3_supp)

arid3a_5_supp <- dfList$`Arid3a_5_suppfrozenset({'input_feather_test_short.genes_vs_tracks.rankings', 'top5perTarget', 'activating'})`
arid3a_5_supp <- processDf(arid3a_5_supp)


write.csv(arid3a_1,"result/Arid3a/regulon/1_arid3a.csv",quote = F)
write.csv(arid3a_3,"result/Arid3a/regulon/3_arid3a.csv",quote = F)
write.csv(arid3a_5,"result/Arid3a/regulon/5_arid3a.csv",quote = F)
write.csv(arid3a_1_supp,"result/Arid3a/regulon/1_supp_arid3a.csv",quote = F)
write.csv(arid3a_3_supp,"result/Arid3a/regulon/3_supp_arid3a.csv",quote = F)
write.csv(arid3a_5_supp,"result/Arid3a/regulon/5_supp_arid3a.csv",quote = F)




library(UpSetR)
UpSetR::upset(fromList(list(motif1 = arid3a_1$target, cuttag = arid3a_3$target, cuttag_combine = arid3a_5$target)),point.size = 5,set_size.numbers_size = 5,
              line.size = 2,set_size.show = 5,)

library(UpSetR)
UpSetR::upset(fromList(list(link = arid3a_1$target, flank = arid3a_1_supp$target)),point.size = 5,set_size.numbers_size = 5,
              line.size = 2,set_size.show = 5,)
