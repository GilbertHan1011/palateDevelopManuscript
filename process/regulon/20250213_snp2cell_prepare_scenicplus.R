atac_peak_scenic <- read.csv("process/framework/peakSets/20250213_scenicplus_peak.csv",row.names = 1)
# Split the coordinates
peaks_df <- data.frame(coords = atac_peak_scenic$RegionIDs) %>%
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
# First, add a new metadata column to peaks_gr for SNP overlap
peaks_gr$has_snp <- FALSE  # Initialize all to FALSE

# Find overlaps between peaks and SNPs
overlaps <- findOverlaps(peaks_gr, snp_gr)

# Mark peaks that have overlapping SNPs as TRUE
peaks_gr$has_snp[queryHits(overlaps)] <- TRUE

# If you also want to know which SNP IDs overlap with each peak
# Create a new metadata column for SNP IDs
peaks_gr$overlapping_snps <- ""

# For each peak that has overlaps, collect the SNP IDs
peak_snps <- split(snp_gr$rsid[subjectHits(overlaps)], queryHits(overlaps))
peaks_gr$overlapping_snps[as.numeric(names(peak_snps))] <- sapply(peak_snps, paste, collapse=",")

# View results
head(peaks_gr)
scenicdf <- as.data.frame(peaks_gr)
scenicdf <- scenicdf[,c(1:6)]
write.csv(scenicdf,"process/framework/snps/20250213_SNPs_annoated_peaks_scenic.csv",quote = F)
