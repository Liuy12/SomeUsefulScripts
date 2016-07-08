library(biomaRt)

## listMarts will display all available BioMart web services
listMarts()

ensembl <- useMart("ensembl")

## we look at which datasets are available
## in the selected BioMart by using the function listDatasets
listDatasets(ensembl)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

## list attributes you can retrive from
listAttributes(ensembl)


geneSymbols <- unique( getBM(attributes = "hgnc_symbol",
                                values = "*", 
                                mart = ensembl) )

seq <- getSequence(id=geneSymbols$hgnc_symbol, type="hgnc_symbol", seqType="3utr", mart = ensembl)

seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="peptide", mart = ensembl)


