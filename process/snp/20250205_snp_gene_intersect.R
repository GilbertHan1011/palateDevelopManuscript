library(stringr)
library(tidyr)
library(GenomicFeature)
snp <- read.table("data/SNP/SNP_summary/GWAS_catlog_with_thx_summary_mm10.bed")
snp_gr <- GRanges(
  seqnames = snp$V1,
  ranges = IRanges(start = snp$V2, end = snp$V3),
  rsid = snp$V4
)

atac_peak <- read.csv("process/framework/peakSets/20241221_snapatac2_merged_peak.csv")
atac_peak_hvg <- rtracklayer::import("process/framework/peakSets/8.25_hvg_peaks_17k.bed")
# Split the coordinates
peaks_df <- data.frame(coords = atac_peak$Peaks) %>%
  separate(coords, 
           into = c("chr", "range"), 
           sep = ":", 
           remove = FALSE) %>%
  separate(range, 
           into = c("start", "end"), 
           sep = "-")

# Convert to GRanges
peaks_gr <- GRanges(
  seqnames = peaks_df$chr,
  ranges = IRanges(
    start = as.numeric(peaks_df$start),
    end = as.numeric(peaks_df$end)
  )
)

library(ChIPseeker)
library(GenomicRanges)
peaks_overlap <- makeVennDiagram(
  list(Peaks1 = gr1, Peaks2 = gr2),
  NameOfPeaks = c("Set1", "Set2"),
  fill = c("blue", "red"),
  alpha = 0.5,
  counts.col = "black"
)

library(VennDiagram)

create_granges_venn <- function(gr1, gr2, 
                                name1 = "Set1", 
                                name2 = "Set2",
                                colors = c("blue", "red"),
                                alpha = 0.5) {
  
  # Input validation
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) stop("GenomicRanges required")
  if (!requireNamespace("VennDiagram", quietly = TRUE)) stop("VennDiagram required")
  if (!requireNamespace("grid", quietly = TRUE)) stop("grid required")
  
  # Ensure inputs are GRanges objects
  if (!is(gr1, "GRanges") || !is(gr2, "GRanges")) {
    stop("Both inputs must be GRanges objects")
  }
  
  # Calculate overlaps
  overlaps <- countOverlaps(gr1, gr2) > 0
  set1_unique <- sum(!overlaps)
  set2_unique <- length(gr2) - sum(overlaps > 0)
  shared <- sum(overlaps)
  
  # Create new plotting page
  grid::grid.newpage()
  
  # Create and return Venn diagram
  venn.plot <- VennDiagram::draw.pairwise.venn(
    area1 = length(gr1),
    area2 = length(gr2),
    cross.area = shared,
    category = c(name1, name2),
    fill = colors,
    alpha = alpha
  )
  
  # Return statistics as well
  stats <- list(
    set1_total = length(gr1),
    set2_total = length(gr2),
    set1_unique = set1_unique,
    set2_unique = set2_unique,
    shared = shared
  )
  
  return(list(plot = venn.plot, statistics = stats))
}

pdf("result/SNP/20250205_snp_atac_peak/snp_peak_venn.pdf")
result1 <- create_granges_venn(snp_gr, peaks_gr, 
                              name1 = "SNPs", 
                              name2 = "scATAC")
dev.off()

gr_extend_500 <- extendGR(snp_gr,upstream = 500,downstream = 500)
pdf("result/SNP/20250205_snp_atac_peak/snp_peak_venn_extended.pdf")
result2 <- create_granges_venn(gr_extend_500, peaks_gr, 
                               name1 = "SNPs", 
                               name2 = "scATAC")
dev.off()

pdf("result/SNP/20250205_snp_atac_peak/snp_peak_venn_hvg.pdf")
result3 <- create_granges_venn(snp_gr, atac_peak_hvg, 
                               name1 = "SNPs", 
                               name2 = "scATAC")
dev.off()

pdf("result/SNP/20250205_snp_atac_peak/snp_peak_venn_hvg_extend.pdf")
result3 <- create_granges_venn(gr_extend_500, atac_peak_hvg, 
                               name1 = "SNPs", 
                               name2 = "scATAC")
dev.off()


library(GenomicRanges)
annotate_snps_detailed <- function(snp_gr, peak_gr, max_distance = 5000) {
  # Find nearest peaks and distances
  nearest_peaks <- nearest(snp_gr, peak_gr)
  distances <- distanceToNearest(snp_gr, peak_gr)
  
  # Get peak coordinates
  peak_coords <- peak_gr[nearest_peaks]
  
  # Create new metadata columns
  snp_gr$peak_distance <- mcols(distances)$distance
  snp_gr$peak_overlap <- overlapsAny(snp_gr, peak_gr)
  snp_gr$peak_position <- ifelse(
    start(snp_gr) < start(peak_coords),
    "upstream",
    ifelse(end(snp_gr) > end(peak_coords),
           "downstream",
           "within")
  )
  
  # Add peak coordinates
  snp_gr$nearest_peak_start <- start(peak_coords)
  snp_gr$nearest_peak_end <- end(peak_coords)
  
  # Filter by distance
  if (!is.null(max_distance)) {
    far_peaks <- snp_gr$peak_distance > max_distance
    snp_gr$peak_distance[far_peaks] <- NA
    snp_gr$peak_position[far_peaks] <- "distant"
  }
  
  return(snp_gr)
}

# Use the function
annotated_snps <- annotate_snps_detailed(snp_gr, peaks_gr)


library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)  # or your preferred TxDb
library(org.Mm.eg.db)

# Get TxDb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate SNPs
annotated_snps_gene <- annotatePeak(annotated_snps,
                          tssRegion = c(-3000, 3000),
                          TxDb = txdb,
                          annoDb = "org.Mm.eg.db")

# Convert to GRanges with annotations
annotated_snps_gene_gr <- as.GRanges(annotated_snps_gene)
df <- as.data.frame(annotated_snps_gene_gr)
write.csv(df,"process/framework/snps/20250205_SNPs_annoated_peaks_genes.csv",quote = F)
rtracklayer::export(annotated_snps_gene_gr, "process/framework/snps/20250205_SNPs_annoated_peaks_genes.bed")

