
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)

# Get the TxDb object for mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Extract gene coordinates
genes <- genes(txdb)

# Add gene symbols
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = genes$gene_id,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
# Create GRanges with upstream 5000bp and gene body
genes_with_upstream <- GRanges(
  seqnames = seqnames(genes),
  ranges = IRanges(
    start = ifelse(strand(genes) == "+",
                   start(genes) - 5000,
                   start(genes)),
    end = ifelse(strand(genes) == "+",
                 end(genes),
                 end(genes) + 5000)
  ),
  strand = strand(genes),
  gene_id = genes$gene_id,
  symbol = gene_symbols
)
# Remove any negative coordinates
genes_with_upstream <- genes_with_upstream[start(genes_with_upstream) > 0]


# To get separate regions for upstream and gene body
create_separate_regions <- function(genes) {
  # Upstream regions
  upstream <- GRanges(
    seqnames = seqnames(genes),
    ranges = IRanges(
      start = ifelse(strand(genes) == "+",
                     start(genes) - 5000,
                     end(genes)),
      end = ifelse(strand(genes) == "+",
                   start(genes),
                   end(genes) + 5000)
    ),
    strand = strand(genes),
    gene_id = genes$gene_id,
    symbol = gene_symbols,
    region_type = "upstream"
  )
  
  # Gene body regions
  gene_body <- GRanges(
    seqnames = seqnames(genes),
    ranges = IRanges(
      start = start(genes),
      end = end(genes)
    ),
    strand = strand(genes),
    gene_id = genes$gene_id,
    symbol = gene_symbols,
    region_type = "gene_body"
  )
  
  # Combine both
  all_regions <- c(upstream, gene_body)
  return(all_regions[start(all_regions) > 0])
}

# Get separate regions
separate_regions <- create_separate_regions(genes)

# Create simple intersection function
annotate_gene_regions <- function(gene_regions, snp_gr) {
  # Find overlaps
  overlaps <- findOverlaps(gene_regions, snp_gr)
  
  # Create result vector (default 0)
  result <- rep(0, length(gene_regions))
  
  # Set 1 for gene regions that overlap with SNPs
  result[queryHits(overlaps)] <- 1
  
  # Add annotation to gene regions
  gene_regions$has_snp <- result
  
  return(gene_regions)
}

# Use the function
annotated_genes <- annotate_gene_regions(separate_regions, snp_gr)

annotated_genes_df <- as.data.frame(annotated_genes)
write.csv(annotated_genes_df,"process/SNP/SNP2ell/20250206_prepare_gene_snp.csv")
