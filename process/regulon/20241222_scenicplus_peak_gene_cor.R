clusterdf = read.csv("result/regulation/20241221_scenicplus/20241222_peak_cluster.csv")
unique(clusterdf$X) %>% length()
K14_high_peak <- clusterdf$X[clusterdf$Cluster%in%c(0,1)]
K14_low_peak <- clusterdf$X[clusterdf$Cluster%in%c(4)]
K14_high_gene <- c(gene_lists$`K14-Enriched`)
K14_low_gene <- c(gene_lists$`K14-Exclueded`)

peak_high_gene <- arid3a_reg1$V2[arid3a_reg1$V1 %in% K14_high_peak]
peak_low_gene <- arid3a_reg1$V2[arid3a_reg1$V1 %in% K14_low_peak]
intersect(peak_high_gene, K14_high_gene) %>% length()
intersect(peak_high_gene, K14_low_gene) %>% length()

intersect(peak_low_gene, K14_high_gene) %>% length()
intersect(peak_low_gene, K14_low_gene) %>% length()

bulk_peak <- read.csv("result/regulation/20241221_scenicplus/20241222_arid3a_peak_bulk.csv")
high_peak <- bulk_peak$X[bulk_peak$K14... > 0.5]
low_peak <- bulk_peak$X[bulk_peak$K14... < -0.5]

high_gene <- rownames(hmData_scale)[hmData_scale[,"K14(-)"] > 0.5]
low_gene <- rownames(hmData_scale)[hmData_scale[,"K14(-)"] < -0.5]

peak_high_gene2 <- arid3a_reg1$V2[arid3a_reg1$V1 %in% high_peak]
peak_low_gene2 <- arid3a_reg1$V2[arid3a_reg1$V1 %in% low_peak]

l1 <- intersect(high_gene, peak_high_gene2) %>% length()
l2 <- intersect(high_gene, peak_low_gene2) %>% length()


l3 <- intersect(low_gene, peak_high_gene2) %>% length()
l4 <- intersect(low_gene, peak_low_gene2) %>% length()


library(UpSetR)

# Create list of sets
sets <- list(
  "High Gene" = high_gene,
  "Peak High" = peak_high_gene2,
  "Low Gene" = low_gene,
  "Peak Low" = peak_low_gene2
)

upset(fromList(sets),
      nsets = 4,
      order.by = "freq",
      main.bar.color = "darkblue",
      sets.bar.color = "darkgreen")





# Fisher's Exact Test for pairwise comparisons
library(tidyverse)

# 1. For High Gene vs Peak High
high_contingency <- matrix(c(
  length(intersect(high_gene, peak_high_gene2)),  # overlap
  length(setdiff(high_gene, peak_high_gene2)),    # only in high_gene
  length(setdiff(peak_high_gene2, high_gene)),    # only in peak_high
  length(setdiff(union(high_gene, peak_high_gene2), 
                 intersect(high_gene, peak_high_gene2))) # in neither
), nrow = 2)

fisher_high <- fisher.test(high_contingency)

# 2. For Low Gene vs Peak Low
low_contingency <- matrix(c(
  length(intersect(low_gene, peak_low_gene2)),
  length(setdiff(low_gene, peak_low_gene2)),
  length(setdiff(peak_low_gene2, low_gene)),
  length(setdiff(union(low_gene, peak_low_gene2), 
                 intersect(low_gene, peak_low_gene2)))
), nrow = 2)

fisher_low <- fisher.test(low_contingency)

# Create summary table
enrichment_results <- data.frame(
  Comparison = c("High Gene vs Peak High", "Low Gene vs Peak Low"),
  Overlap = c(length(intersect(high_gene, peak_high_gene2)),
              length(intersect(low_gene, peak_low_gene2))),
  OddsRatio = c(fisher_high$estimate, fisher_low$estimate),
  PValue = c(fisher_high$p.value, fisher_low$p.value)
)

# Add FDR correction
enrichment_results$FDR <- p.adjust(enrichment_results$PValue, method = "BH")

# Print results
print(enrichment_results)

# Visualize results
library(ggplot2)

ggplot(enrichment_results, aes(x = Comparison, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("p = %.2e\nFDR = %.2e", PValue, FDR)), 
            vjust = -0.5) +
  theme_minimal() +
  labs(title = "Set Enrichment Analysis",
       y = "Odds Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
