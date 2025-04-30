csvfile <- list.files("process/regulation/20250107_celloracle/20250107_pertubation_score/",full.names = T, pattern = "*.csv")
name <- list.files("process/regulation/20250107_celloracle/20250107_pertubation_score/",full.names = F, pattern = "*.csv") %>% gsub(".csv","",.)

scoreData=lapply(csvfile,read.csv)
names(scoreData) <- name
combineData <- lapply(scoreData,function(x){
  x$PS
  })
combineData <- do.call(rbind,combineData)
colnames(combineData) <- scoreData[[1]]$X
rownames(combineData) <- name
pdf("process/regulation/20250107_celloracle/purtubation_score.pdf",width = 6,height = 10)
ComplexHeatmap::Heatmap(combineData)
dev.off()
K6top5 <-combineData[,"K6+ cells"]
sort(K6top5,decreasing = T) %>% head(5) %>% names
K14top5 <-combineData[,"K14(-)"]
sort(K14top5,decreasing = T) %>% head(5) %>% names
