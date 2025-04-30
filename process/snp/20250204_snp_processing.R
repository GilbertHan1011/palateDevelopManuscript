library(tidyverse)
library(data.table)
fileList <- list.files("data/SNP/", pattern = "*.tsv",full.names = T)
fileShort <- list.files("data/SNP/", pattern = "*.tsv") %>% gsub(".tsv","",.)
files <- lapply(fileList,data.table::fread)
names(files) <- fileShort
fileBind <- do.call(rbind,files)
fileBed <- fileBind[,c("CHR_ID","CHR_POS","SNPS")]
length(unique(fileBed$SNPS))
bedSum <- data.frame(
  chr = paste0("chr", fileBed$CHR_ID),  # Add 'chr' prefix
  start = fileBed$CHR_POS - 1,          # BED is 0-based, so subtract 1
  end = fileBed$CHR_POS,                # Original position as end
  name = fileBed$SNPS                  # SNP ID as name
) %>% unique
write.table(bedSum,"data/SNP/GWAS_catlog_summary.bed",sep = "\t",col.names = F, row.names = F, quote = F)

##
# Install and load rsnps
if (!require("rsnps")) {
  install.packages("rsnps")
}
library(rsnps)

# Get SNP information
#snps <- ncbi_snp_query(bedSum$name[1:10])

# find where snp info missing
missing_bed <- bedSum[is.na(bedSum$start),]
unmissing_bed <- bedSum[!is.na(bedSum$start),]
rs_bed <- missing_bed[grep("rs",missing_bed$name),]
non_rs_bed <-  missing_bed[grep("rs",missing_bed$name,invert = T),]
snps <- ncbi_snp_query(rs_bed$name)
snps_bed <- data.frame(chr = paste0("chr", snps$chromosome),
                       start = snps$bp - 1,
                       end = snps$bp,
                       name = snps$query)
# Create a function to process the positions
process_positions <- function(pos_str) {
  # Split chromosome and position
  parts <- strsplit(pos_str, ":")[[1]]
  
  # Handle cases where 'chr' prefix is missing
  chr <- parts[1]
  if (!grepl("^chr", chr, ignore.case = TRUE)) {
    chr <- paste0("chr", chr)
  }
  # Standardize chr format (lowercase 'chr')
  chr <- tolower(chr)
  
  # Get position
  pos <- as.numeric(parts[2])
  
  # Return processed values
  return(c(chr, pos))
}

# Process all positions
result <- do.call(rbind, lapply(non_rs_bed$name, process_positions))

# Create data frame
bed_df <- data.frame(
  chr = result[,1],
  start = as.numeric(result[,2]) - 1,  # BED format is 0-based
  end = as.numeric(result[,2]),
  name = paste0("region_", 1:nrow(result))  # Create names if not provided
)

# Sort by chromosome and position
bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]
bed_df <- drop_na(bed_df)
final_bed <- do.call( rbind,list(unmissing_bed,snps_bed, bed_df))
write.table(final_bed,"data/SNP/GWAS_catlog_summary.bed",sep = "\t",col.names = F, row.names = F, quote = F)
write.table(fileBind,"data/SNP/GWAS_catlog_summary_full.tsv",sep = "\t",col.names = F, row.names = F, quote = F)


#== transform snp from thx------
hg19_thx <- read.table("data/SNP/NSCPO_hg19.bed",sep = "\t")
hg19_thx$name <- strsplit(hg19_thx$V4,"_") %>% lapply(`[`, 2) %>% unlist()
hg19_thx <- hg19_thx[,c(1,2,3,5)]
colnames(hg19_thx) <- colnames(final_bed)
names(hg19_thx)
summary_bed <- rbind(final_bed,hg19_thx)
summary_bed <- summary_bed %>% unique
#write.table(fileBind,"data/SNP/GWAS_catlog_summary_full.tsv",sep = "\t",col.names = F, row.names = F, quote = F)


duplicate_row <- summary_bed[summary_bed$name%in%summary_bed$name[duplicated(summary_bed$name)],]
duplicate_row <- duplicate_row[order(duplicate_row$chr, duplicate_row$start), ]

summary_bed$start[summary_bed$start == summary_bed$end] = summary_bed$start[summary_bed$start == summary_bed$end]-1
summary_bed <- summary_bed[!duplicated(summary_bed$name),]
write.table(summary_bed,"data/SNP/GWAS_catlog_with_thx_summary.bed",sep = "\t",col.names = F, row.names = F, quote = F)
