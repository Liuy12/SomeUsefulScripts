## ----- generate_filter_gtf
## <<generate_filter_gtf.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))

input_gtf=args[1]
input_truth=args[2]
min_abund=args[3]
output_gtf=args[4]
min_abund<- as.numeric(min_abund)

##### usage
if (!file.exists(input_gtf)){
  writeLines ("#################################################################################
              \n\nUsage:\nRscript generate_filter_gtf.R input_gtf input_truth min_abund output_gtf
              
              input_gtf: path to input gtf 
              input_truth: path to table containing isoform level abundance
              min_abund: percentage of minnimum abundance
              output_gtf: path to output gtf
              ##################################################################################")
  quit()
}

print(input_gtf)
print(input_truth)
print(min_abund)
print(output_gtf)

source('http://bioconductor.org/biocLite.R')
requiredPacks <- 'rtracklayer'
for(i in requiredPacks){
  if(!require(i, character.only = T)) biocLite(i)
  library(i, character.only = T)
}

## Read gtf file and truth
gtf_file <- import(input_gtf)
truth <- read.delim(input_truth, header = TRUE, as.is = TRUE)

## Determine which transcripts to exclude
truth$avg <- apply(truth[,-1], 1, mean)
tx_remove <- truth$transcript_id[truth$avg <= quantile(truth$avg, min_abund)]

## Subset gtf file
new_gtf <- gtf_file[-which(gtf_file$transcript_id %in% tx_remove), ]
new_gtf$exon_number <- as.character(new_gtf$exon_number)

## Write to files
export(new_gtf, output_gtf)