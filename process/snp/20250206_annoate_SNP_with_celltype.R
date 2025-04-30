atac_peak <- read.csv("process/framework/peakSets/20241221_snapatac2_merged_peak.csv")
#atac_peak_hvg <- rtracklayer::import("process/framework/peakSets/8.25_hvg_peaks_17k.bed")
# Split the coordinates
peaks_df <- data.frame(coords = atac_peak$Peaks) %>%
  separate(coords, 
           into = c("chr", "range"), 
           sep = ":", 
           remove = FALSE) %>%
  separate(range, 
           into = c("start", "end"), 
           sep = "-")

peaks_df <- cbind(peaks_df,atac_peak)

# First, rename the columns in peaks_df to match df
peaks_df <- peaks_df %>%
  dplyr::rename(seqnames = chr,
         nearest_peak_start = start,
         nearest_peak_end = end)
peaks_df <- peaks_df %>%
  mutate(across(c(nearest_peak_start, nearest_peak_end), as.integer))
# Perform the left join
result_df_with_peak <- df %>%
  dplyr::left_join(peaks_df, 
            by = c("seqnames", 
                   "nearest_peak_start", 
                   "nearest_peak_end"))
result_df_with_peak <- result_df_with_peak[,c("seqnames", "start", "end", "width", "strand", "rsid", "peak_distance", 
                      "peak_overlap", "peak_position", "nearest_peak_start", "nearest_peak_end", 
                      "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", 
                      "geneStrand", "geneId", "transcriptId", "distanceToTSS", "ENSEMBL", 
                      "SYMBOL", "GENENAME", "coords",  "stem.cells", "Transit", 
                      "K5...", "Mature.K5...", "K6..cells", "K14...", "Shh...")]

write.csv(result_df_with_peak,"process/framework/snps/20250206_SNPs_annoated_peaks_genes_withcelltype.csv")

result_df_with_peak$stem.cells = result_df_with_peak$stem.cells=="true"
result_df_with_peak$Transit = result_df_with_peak$Transit=="true"
result_df_with_peak$K5... = result_df_with_peak$K5...=="true"
result_df_with_peak$Mature.K5... = result_df_with_peak$Mature.K5...=="true"
result_df_with_peak$K6..cells = result_df_with_peak$K6..cells=="true"
result_df_with_peak$K14... = result_df_with_peak$K14...=="true"
result_df_with_peak$Shh... = result_df_with_peak$Shh...=="true"

summary_data <- data.frame(
  Cell_Type = c("Stem Cells", "Transit", "K5+", "Mature K5+", "K6+ cells", "K14-", "Shh+"),
  Count = c(
    sum(result_df_with_peak$stem.cells, na.rm = TRUE),
    sum(result_df_with_peak$Transit, na.rm = TRUE),
    sum(result_df_with_peak$K5..., na.rm = TRUE),
    sum(result_df_with_peak$Mature.K5..., na.rm = TRUE),
    sum(result_df_with_peak$K6..cells, na.rm = TRUE),
    sum(result_df_with_peak$K14..., na.rm = TRUE),
    sum(result_df_with_peak$Shh..., na.rm = TRUE)
  )
)


# Set the factor levels to specify the order
summary_data$Cell_Type <- factor(summary_data$Cell_Type, 
                                 levels = c("Stem Cells", "Transit", "K5+", 
                                            "Mature K5+", "K6+ cells", "K14-", "Shh+"))

# Create the plot
ggplot(summary_data, aes(x = Cell_Type, y = Count, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Distribution of Cell Types",
       x = "Cell Type",
       y = "Count") +
  geom_text(aes(label = Count), vjust = -0.3) +
  scale_fill_brewer(palette = "Set3")


# Calculate proportions
total_rows <- nrow(result_df_with_peak)
summary_data$Proportion <- summary_data$Count / total_rows


# Create proportion plot
ggplot(summary_data, aes(x = Cell_Type, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "SNP_",
       x = "Cell Type",
       y = "Proportion") +
  geom_text(aes(label = sprintf("%.3f", Proportion)), vjust = -0.3) +
  scale_fill_brewer(palette = "Set3")


peaks_df <- peaks_df %>%
  mutate(across(c(stem.cells, Transit, K5..., Mature.K5..., 
                  K6..cells, K14..., Shh...),
                ~.x == "true"))



summary_data2 <- data.frame(
  Cell_Type = c("Stem Cells", "Transit", "K5+", "Mature K5+", "K6+ cells", "K14-", "Shh+"),
  Count = c(
    sum(peaks_df$stem.cells, na.rm = TRUE),
    sum(peaks_df$Transit, na.rm = TRUE),
    sum(peaks_df$K5..., na.rm = TRUE),
    sum(peaks_df$Mature.K5..., na.rm = TRUE),
    sum(peaks_df$K6..cells, na.rm = TRUE),
    sum(peaks_df$K14..., na.rm = TRUE),
    sum(peaks_df$Shh..., na.rm = TRUE)
  )
)


# Set the factor levels to specify the order
summary_data2$Cell_Type <- factor(summary_data2$Cell_Type, 
                                 levels = c("Stem Cells", "Transit", "K5+", 
                                            "Mature K5+", "K6+ cells", "K14-", "Shh+"))

# Create the plot
ggplot(summary_data2, aes(x = Cell_Type, y = Count, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Distribution of Cell Types",
       x = "Cell Type",
       y = "Count") +
  geom_text(aes(label = Count), vjust = -0.3) +
  scale_fill_brewer(palette = "Set3")
total_rows2 <- nrow(peaks_df)
summary_data2$Proportion <- summary_data2$Count / total_rows2

# Create proportion plot
ggplot(summary_data2, aes(x = Cell_Type, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Original propotion",
       x = "Cell Type",
       y = "Proportion") +
  geom_text(aes(label = sprintf("%.3f", Proportion)), vjust = -0.3) +
  scale_fill_brewer(palette = "Set3")


create_cell_type_plots <- function(df, 
                                   cell_type_order = c("Stem Cells", "Transit", "K5+", 
                                                       "Mature K5+", "K6+ cells", "K14-", "Shh+"),
                                   title_count = "Distribution of Cell Types",
                                   title_prop = "Original proportion") {
  
  # Create summary data
  summary_data <- data.frame(
    Cell_Type = c("Stem Cells", "Transit", "K5+", "Mature K5+", "K6+ cells", "K14-", "Shh+"),
    Count = c(
      sum(df$stem.cells, na.rm = TRUE),
      sum(df$Transit, na.rm = TRUE),
      sum(df$K5..., na.rm = TRUE),
      sum(df$Mature.K5..., na.rm = TRUE),
      sum(df$K6..cells, na.rm = TRUE),
      sum(df$K14..., na.rm = TRUE),
      sum(df$Shh..., na.rm = TRUE)
    )
  )
  
  # Set factor levels
  summary_data$Cell_Type <- factor(summary_data$Cell_Type, levels = cell_type_order)
  
  # Calculate proportion
  total_rows <- nrow(df)
  summary_data$Proportion <- summary_data$Count / total_rows
  
  # Create count plot
  count_plot <- ggplot(summary_data, aes(x = Cell_Type, y = Count, fill = Cell_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = title_count,
         x = "Cell Type",
         y = "Count") +
    geom_text(aes(label = Count), vjust = -0.3) +
    scale_fill_brewer(palette = "Set3")
  
  # Create proportion plot
  prop_plot <- ggplot(summary_data, aes(x = Cell_Type, y = Proportion, fill = Cell_Type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = title_prop,
         x = "Cell Type",
         y = "Proportion") +
    geom_text(aes(label = sprintf("%.3f", Proportion)), vjust = -0.3) +
    scale_fill_brewer(palette = "Set3")
  
  # Return both plots in a list
  return(list(count_plot = count_plot, proportion_plot = prop_plot))
}

#== snp_nearby_peak------------
result_df_with_peak_snp <- result_df_with_peak[result_df_with_peak$peak_overlap,]
plots1 <- create_cell_type_plots(peaks_df)
plots2 <- create_cell_type_plots(result_df_with_peak_snp,title_prop = "SNP_nearby")
grid.arrange(plots1$proportion_plot, plots2$proportion_plot, ncol = 2)

write.csv(result_df_with_peak_snp,"process/framework/snps/20250206_SNPs_annoated_peaks_genes_withcelltype_overlappedSNP.csv")
