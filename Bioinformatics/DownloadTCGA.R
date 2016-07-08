library(TCGAbiolinks)

LIHC_RNASeqV2_version <- TCGAquery_Version(tumor = "lihc", 
                                           platform = "illuminahiseq_rnaseqv2")

query <- TCGAquery(tumor = c("LIHC"), level = 3,
                   platform = c("IlluminaHiSeq_RNASeqV2")
)

TCGAdownload(query, path = "Downloads/", type = "rsem.genes.results")

