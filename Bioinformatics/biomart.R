library(biomaRt)

## listMarts will display all available BioMart web services
listMarts()
## use archive verison of biomart
listEnsemblArchives()
## use a previous version of ensemble, e.g. v91
ensembl <- useMart("ensembl",host = "https://dec2017.archive.ensembl.org")
## otherwise use default
ensembl <- useMart("ensembl")

## we look at which datasets are available
## in the selected BioMart by using the function listDatasets
listDatasets(ensembl)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

## list attributes you can retrive from
listAttributes(ensembl)

## retrive all gene symbols
geneSymbols <- unique( getBM(attributes = "hgnc_symbol",
                                values = "*", 
                                mart = ensembl) )
### get 3'utr sequence of a selected set of genes
seq <- getSequence(id=geneSymbols$hgnc_symbol, type="hgnc_symbol", seqType="3utr", mart = ensembl)
### get peptide sequence for BRCA1
seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="peptide", mart = ensembl)

## gene id conversion between ref-seq and gene symbol
results <- getBM(attributes = c('ref_mrna',"hgnc_symbol"),
                 filters = "refseq_mrna", values = 'NM110332',
                 mart = ensembl)


#########################

### retriving gene annotation information for mouse samples
ensembl <- useMart('ensembl',dataset='mmusculus_gene_ensembl')
results <- getBM(attributes = c('mgi_symbol','ensembl_gene_id','gene_biotype', 'mgi_description'), filters = 'mgi_symbol', values = dataMat$GeneID, mart = ensembl)

### retriving gene annotation information for human samples
ensembl <- useMart('ensembl',dataset='hsapiens_gene_ensembl')
results <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id','gene_biotype', 'description'), filters = 'hgnc_symbol', values = dataMat$GeneID, mart = ensembl)

####### retrive homology information for mouse and human
####### match between mouse gene symbol, human gene symbol
musGenes <- c("Hmmr", "Tlx3", "Cpeb4")

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #### retrive information between two linked dataset
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- genesV2
  return(humanx)
}

