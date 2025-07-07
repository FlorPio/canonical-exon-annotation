library(GenomicRanges)
# Download MANE Select data for human (GRCh38)
mane_url <- "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz"
temp_file <- tempfile(fileext = ".gz")
download.file(mane_url, temp_file)

#Read the gzipped file
con <- gzfile(temp_file, "r")
lines <- readLines(con)
close(con)
# The first row is the header with columns names (with #)
header <- unlist(strsplit(lines[1], "\t"))
# Quit the # from the first name
header[1] <- gsub("^#", "", header[1])
# Filtering rows which are not comments or header
data_lines <- lines[!grepl("^#", lines)]
# Convert the data lines into a data frame
mane_data <- data.frame(matrix(ncol = length(header), nrow = length(data_lines)))
colnames(mane_data) <- header
#Fill the data frame
for (i in 1:length(data_lines)) {
  row_data <- unlist(strsplit(data_lines[i], "\t"))
  mane_data[i, ] <- row_data
}


#quit the version number from the RefSeq_nuc column
mane_data$RefSeq_nuc <- gsub("\\..*", "", mane_data$RefSeq_nuc)

#Exctract the MANE IDs (without version)
mane_ids <- mane_data$RefSeq_nuc

#Ad chr prefix to the chromosome names
mane_data$GRCh38_chr <- paste0("chr", mane_data$GRCh38_chr) 

#Convert to GRanges 
mane_data_gr <- makeGRangesFromDataFrame(mane_data, 
                                         seqnames.field = "GRCh38_chr", 
                                         start.field = "chr_start", 
                                         end.field = "chr_end", 
                                         strand.field = "chr_strand",
                                         keep.extra.columns = TRUE)

#Load necessary libraries
#Load TxDb for the human genome (UCSC hg38) to get exons of transcripts
library(TxDb.Hsapiens.UCSC.hg38.refGene)
Tx_38 <- TxDb.Hsapiens.UCSC.hg38.refGene
columns(Tx_38)

exonsBytranscript <- exonsBy(Tx_38, by="tx", use.names = TRUE)
exonsBytranscript[1:2]
# Filter exonsBytranscript to keep only those transcripts that are in the MANE dataset
exons_mane <- exonsBytranscript[names(exonsBytranscript) %in% mane_ids]
exons_mane[1:2]
# Convert the exons to a GRangesList object in which each transcript's exons are grouped together
exons_mane_gr <- unlist(exons_mane)

# Add metadata: gene and transcript to the Tx GRanges object
mcols(exons_mane_gr)$gene_symbol <- 
  mane_data$symbol[match(names(exons_mane_gr), mane_data$RefSeq_nuc)]

mcols(exons_mane_gr)$transcript_id <- names(exons_mane_gr)

#Create a data frame 
exons_df <- data.frame(
  chr = seqnames(exons_mane_gr),
  start = start(exons_mane_gr),
  end = end(exons_mane_gr),
  gene_symbol = mcols(exons_mane_gr)$gene_symbol,
  exon = mcols(exons_mane_gr)$exon_rank,
  transcript_id = mcols(exons_mane_gr)$transcript_id
  
)

missing_transcripts <- setdiff(mane_data$RefSeq_nuc, names(exons_mane))
# Show result
head(exons_df)

# Save the data frame as CSV
write.csv(exons_df, "exons_mane.csv", row.names = FALSE)

# Save the data frame as a txt object 
write.table(exons_df, "exons_mane.txt", row.names = FALSE)
