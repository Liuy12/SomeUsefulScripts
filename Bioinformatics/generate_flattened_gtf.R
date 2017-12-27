## ----- generate_flattened_gtf
## <<generate_flattened_gtf.R>>

## Generate manually flattened gtf files to use with featureCounts
args <- (commandArgs(trailingOnly = TRUE))
input_gtf=args[1]
output_gtf=args[2]
ignore_strand=args[3]

##### usage
if (!file.exists(input_gtf)){
  writeLines ("#################################################################################
              \n\nUsage:\nRscript generate_flattened_gtf.R input_gtf output_gtf ignore_strand
              
              input_gtf: path to input gtf 
              output_gtf: path to output gtf
              ignore_strand: ignore strand or not. TRUE/FALSE
              ##################################################################################")
  quit()
}

if(file.exists(output_gtf))
  stop('Output gtf exsits and will be replaced')

print(input_gtf)
print(output_gtf)
print(ignore_strand)

t1 <- Sys.time()
#### loading required libraries
source('http://bioconductor.org/biocLite.R')
requiredPacks <- c('GenomicRanges', 'rtracklayer')
for(i in requiredPacks){
  if(!require(i, character.only = T)) biocLite(i)
  library(i, character.only = T)
}
## Import gtf file
gtf0 <- import(input_gtf)

## Keep only exons
idx <- mcols(gtf0)$type == "exon" 
gtf0 <- gtf0[idx]

# Split by gene
gtf.g <- split(gtf0, mcols(gtf0)$gene_id)

# Flatten per gene
gtf.g.flat <- disjoin(gtf.g)

gtf.g.flat <- unlist(gtf.g.flat)
mcols(gtf.g.flat)$gene_id <- names(gtf.g.flat)

# Flatten overall
gtf.o.flat <- disjoin(gtf.g.flat, ignore.strand = as.logical(ignore_strand))

# Remove overlapping parts
co <- countOverlaps(gtf.o.flat, gtf.g.flat)
gtf.o.flat <- gtf.o.flat[co == 1]

# Find the gene names
fo <- findOverlaps(gtf.o.flat, gtf.g.flat)
mcols(gtf.o.flat)$gene_id <- mcols(gtf.g.flat)$gene_id[subjectHits(fo)]
mcols(gtf.o.flat)$transcript_id <- mcols(gtf.o.flat)$gene_id
mcols(gtf.o.flat)$type <- "exon"

gtf.o.flat.sort <- gtf.o.flat[order(mcols(gtf.o.flat)$gene_id, decreasing = FALSE)]

exon.id <- split(mcols(gtf.o.flat.sort)$gene_id, mcols(gtf.o.flat.sort)$gene_id)

exon.id.new <- unlist(lapply(exon.id, function(g){ 
  seq(1, length(g)) 
}))

mcols(gtf.o.flat.sort)$exon_number <- exon.id.new
mcols(gtf.o.flat.sort)$exon_id <- paste0(mcols(gtf.o.flat.sort)$gene_id,
                                         ":", sprintf( "%03d", exon.id.new))
export(gtf.o.flat.sort, output_gtf, format = "gtf")

## Fix gff file for DEXSeq. Note! Still need to change the strand to . instead of *
# gtf.o.flat.sort$type[gtf.o.flat.sort$type == "exon"] <- "exonic_part"
# mcols(gtf.o.flat.sort)$exonic_part_number <- sprintf("%03d", gtf.o.flat.sort$exon_number)
# tmp <- gtf.o.flat.sort
# mcols(tmp)$group <- paste0("transcripts \"", tmp$transcript_id, 
#                            "\"; exonic_part_number \"", tmp$exonic_part_number, 
#                            "\"; gene_id \"", tmp$gene_id, "\"")
# 
# export(tmp, gsub("\\.gtf$", ".gff", output_gtf), format = "gff")
