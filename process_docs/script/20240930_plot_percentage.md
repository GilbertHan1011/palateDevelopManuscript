```R
## Plot composition of RNA sample-------------------
plotBarData <- dietPalate@meta.data%>%

select(ident=level1_anno,time=orig.ident)

plotBarData$time <- factor(plotBarData$time,levels=c("E125_1","E125_2","E145_1","E145_2"))

mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(7)

ggplot(plotBarData, aes(x=time, fill=ident)) +

geom_bar(position = "fill")+theme_classic()+

theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1,size = 10,face = "bold"),title = NULL)+

scale_fill_manual(values=mycolor)

ggplot(plotBarData, aes(x=time, fill=ident)) +

geom_bar(position = "fill") +

theme_classic() +

theme(

axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10, face = "bold"), # Adjusted x-axis text angle

plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Added title formatting

legend.position = "right", # Moved legend to the right

legend.title = element_blank() # Removed legend title

) +

labs(

title = "Composition of RNA Sample Over Time", # Added plot title

x = "Time Point", # Added x-axis label

y = "Proportion" # Added y-axis label

) +

scale_fill_manual(values=mycolor)

  

unique_idents <- c("stem cells")

  

mycolor <- colorRampPalette(brewer.pal(8,'Set3'))(length(unique(plotBarData$ident)))

names(mycolor) <- unique(plotBarData$ident)

mycolor[unique_idents] <- "red" # Assign red color to unique idents

  

plotBarData <- dietPalate@meta.data %>%

select(ident = level1_anno, time = orig.ident)

plotBarData$time <- factor(plotBarData$time, levels = c("E125_1", "E125_2", "E145_1", "E145_2"))

  

ggplot(plotBarData, aes(x = time, fill = ident)) +

geom_bar(position = "fill") +

theme_classic() +

theme(

axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, face = "bold"), # Adjusted x-axis text angle

plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Added title formatting

legend.position = "right", # Moved legend to the right

legend.title = element_blank() # Removed legend title

) +

labs(

title = "Composition of Idents Over Time", # Added plot title

x = "Time Point", # Added x-axis label

y = "Proportion" # Added y-axis label

) +

scale_fill_manual(values = mycolor)

ggsave("descriptive_results/rna_annotation_plot/20240930_composition.pdf",width = 6,height = 6)

  
  

plotLineData <- plotBarData %>%

group_by(time, ident) %>%

summarise(count = n()) %>%

mutate(proportion = count / sum(count))

  

# Create a custom color palette

mycolor <- colorRampPalette(brewer.pal(8,'Set3'))(length(unique(plotBarData$ident)))

names(mycolor) <- unique(plotBarData$ident)

  

#== default-------

ggplot(plotLineData, aes(x = time, y = proportion, color = ident, group = ident)) +

geom_line(size = 1) +

geom_point(size = 3) +

theme_classic() +

theme(

axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, face = "bold"), # Adjusted x-axis text angle

plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Added title formatting

legend.position = "right", # Moved legend to the right

legend.title = element_blank() # Removed legend title

) +

labs(

title = "Trend of RNA Sample Composition", # Added plot title

x = "Time Point", # Added x-axis label

y = "Proportion" # Added y-axis label

) +

scale_color_manual(values = mycolor)

  
  

#== adjust point color-------

# Create a custom color palette for ident

ident_colors <- colorRampPalette(brewer.pal(8,'Set3'))(length(unique(plotBarData$ident)))

names(ident_colors) <- unique(plotBarData$ident)

# Create a color palette for time points

  

ggplot(plotLineData, aes(x = time, y = proportion, color = ident, group = ident)) +

geom_line(size = 1) +

geom_point(aes(color = time), size = 3) +

theme_classic() +

theme(

axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, face = "bold"), # Adjusted x-axis text angle

plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Added title formatting

legend.position = "right", # Moved legend to the right

legend.title = element_blank() # Removed legend title

) +

labs(

title = "Trend of RNA Sample Composition", # Added plot title

x = "Time Point", # Added x-axis label

y = "Proportion" # Added y-axis label

) +

scale_color_manual(values = c(ident_colors, time_colors))

  

#== add rect-------

ggplot(plotLineData, aes(x = time, y = proportion, color = ident, group = ident)) +

geom_line(size = 1) +

geom_point(size = 3) +

theme_classic() +

theme(

axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, face = "bold"), # Adjusted x-axis text angle

plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Added title formatting

legend.position = "right", # Moved legend to the right

legend.title = element_blank() # Removed legend title

) +

labs(

title = "Trend of RNA Sample Composition", # Added plot title

x = "Time Point", # Added x-axis label

y = "Proportion" # Added y-axis label

) +

scale_color_manual(values = mycolor) +

geom_rect(data = data.frame(xmin = c(0.5, 2.5), xmax = c(2.5, 4.5), ymin = -Inf, ymax = Inf),

aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),

fill = c("Deepskyblue", "green"), alpha = 0.1, inherit.aes = FALSE)

ggsave("descriptive_results/rna_annotation_plot/20240930_composition_with_rect.pdf",width = 6,height = 6)

  
  
  

plotBarData <- dietPalate@meta.data%>%

select(ident=level1_anno,time=orig.ident)

plotBarData$time <- factor(plotBarData$time,levels=c("E125_1","E125_2","E145_1","E145_2"))

  

# Specify the desired order for the idents

desired_order <- c("stem cells", "K14(-)", "Transit", "K6+ periderm",

"K5(-)", "Shh(+)", "Mature K5(-) ")

  

# Order the idents based on the specified order

plotBarData$ident <- factor(plotBarData$ident, levels = desired_order)

  

# Create a custom color palette for time

time_colors <- c("E125_1" = "deepskyblue3", "E125_2" = "deepskyblue2", "E145_1" = "lightcoral", "E145_2" = "red")

  

ggplot(plotBarData, aes(x = ident, fill = time)) +

geom_bar(position = "fill") +

theme_classic() +

theme(

axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, face = "bold"), # Adjusted x-axis text angle

plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Added title formatting

legend.position = "right", # Moved legend to the right

legend.title = element_blank() # Removed legend title

) +

labs(

title = "Composition of Time Over Ident", # Added plot title

x = "Propotion", # Updated x-axis label

y = "Ident" # Added y-axis label

) +

scale_fill_manual(values = time_colors)

  

ggsave("descriptive_results/rna_annotation_plot/20240930_composition_time.pdf",width = 6,height = 6)
```