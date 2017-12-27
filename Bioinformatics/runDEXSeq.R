## ----- DEXSeq workflow
## <<runDEXSeq.R>>

args <- (commandArgs(trailingOnly = TRUE))
countfile=args[1]
sampleTable=args[2]
flattenedfile=args[3]
output=args[4]


##### usage
if (!file.exists(countfile)){
  writeLines ("#################################################################################
              \n\nUsage:\nRscript runDEXSeq.R countfile sampleTable flattenedfile output
              
              countfile: path to count table from featurecount
              sampleTable: path to sample table indicating grouping information
              flattenedfile: path to flattened gtf annotation
              output: path to output folder
              ##################################################################################")
  quit()
}


t1 <- Sys.time()
#### loading required libraries
source('http://bioconductor.org/biocLite.R')
requiredPacks <- c('BiocParallel', 'DEXSeq', 'dplyr')
for(i in requiredPacks){
  if(!require(i, character.only = T)) biocLite(i)
  library(i, character.only = T)
}

## function
DEXSeqDataSetFromFeatureCounts <- function (countfile, sampleData, 
                                            design = ~sample + exon + condition:exon, flattenedfile = NULL) 
  
{
  # Take a fcount file and convert it to dcounts for dexseq
  message("Reading and adding Exon IDs for DEXSeq")
  read.table(countfile,skip = 2, stringsAsFactors = F) %>% dplyr::arrange(V1,V3,V4) %>% dplyr::select(-(V2:V6)) -> dcounts
  sampleData <- read.table(sampleData, row.names = 1, sep = '\t', stringsAsFactors = F)
  colnames(dcounts) <- c("GeneID", rownames(sampleData) )
  id <- as.character(dcounts[,1])
  #n <- id
  #split(n,id) <- lapply(split(n ,id), seq_along )
  id <- strsplit(id, ':')
  id1 <- sapply(id, function(i) i[1])
  n <- sapply(id, function(i) i[2])
  rownames(dcounts) <- sprintf("%s%s%03.f",id1,":E",as.numeric(n))
  dcounts <- dcounts[,2:ncol(dcounts)]
  
  dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name 
  
  ## get genes and exon names out
  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, "[[", 2)
  genesrle <- sapply(splitted, "[[", 1)
  
  ## parse the flattened file
  if (!is.null(flattenedfile)) {
    aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, 
                             header = FALSE)
    colnames(aggregates) <- c("chr", "source", "class", "start", 
                              "end", "ex", "strand", "ex2", "attr")
    aggregates$strand <- gsub("\\.", "*", aggregates$strand)
    aggregates <- aggregates[which(aggregates$class == "exon"), # exonic_part
                             ]
    aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
    aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", 
                              aggregates$attr)
    # trim the gene_ids to 255 chars in order to match with featurecounts
    longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
    warning(paste0(longIDs, 
                   " aggregate geneIDs were found truncated in featureCounts output"), 
            call. = FALSE)
    aggregates$gene_id <- substr(aggregates$gene_id,1,255)
    
    transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", 
                        aggregates$attr)
    transcripts <- strsplit(transcripts, "\\+")
    exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1", # exonic_part_number
                    aggregates$attr)
    exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, 
                                                              end = aggregates$end), strand = aggregates$strand)
    names(exoninfo) <- sprintf("%s%s%03.f",aggregates$gene_id,":E",as.numeric(exonids))
    names(transcripts) <- names(exoninfo) 
    # if (!all(rownames(dcounts) %in% names(exoninfo))) {
    #   stop("Count files do not correspond to the flattened annotation file")
    # }
    matching <- match(rownames(dcounts), names(exoninfo))
    cat('number of matched items: ', length(matching), 'out of ', nrow(dcounts), '\n')
    stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
    stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
    dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, 
                         genesrle, exoninfo[matching], transcripts[matching])
    return(dxd)
  }
  else {
    dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons, 
                         genesrle)
    return(dxd)
  }
  
}


#################################### 
##### main script

dxd <- DEXSeqDataSetFromFeatureCounts(
  countfile, 
  sampleTable, 
  design = ~sample + exon + condition:exon, 
  flattenedfile=flattenedfile)

BPPARAM = MulticoreParam(workers=detectCores()/4)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
dxd = testForDEU( dxd, BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )
#table ( dxr1$padj < 0.1 )
#table ( tapply( dxr1$padj < 0.01, dxr1$groupID, any ) )
#plotMA( dxr1, cex=0.8 )
#plotDEXSeq( dxr1, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE,
#            cex.axis=1.2, cex=1.3, lwd=2 )
DEXSeqHTML( dxr1, FDR=0.01, path = paste0(output, "/DEXSeqReport"), color=c("#FF000080", "#0000FF80") )
## control fdr at gene level
perGene <- perGeneQValue(dxr1)

### saving results
write.table(dxr1, paste0(output, "/deutestStats.txt"), sep = '\t', quote =F, row.names = T, col.names = T)
write.table(perGene, paste0(output, '/deutestStats_pergene.txt'), sep = '\t', quote =F, row.names = T, col.names = F)
save.image(paste0(output, '/data.rdata'))

t2 <- Sys.time()
cat('Total time spent: ', as.double(t2-t1), units(t2-t1) , '\n')