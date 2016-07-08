library(data.table)

inforTable <- fread('Downloads/unc.edu_LIHC.IlluminaHiSeq_RNASeqV2.1.15.0.sdrf.txt', data.table = F)
filenames <- list.files(path = "./",full.names=T)
dataMat <- do.call("cbind", lapply(filenames, fread, header = TRUE,data.table =F))
dataMat <- dataMat[,seq(2,1696,4)]

filenames1 <- sapply(filenames, function(i) {
  index <- grep(substr(i, 12, 47), inforTable[,1])
  if(length(index))
    inforTable[index[1],2]
})

colnames(dataMat) <- filenames1

dataMat <- constructPairedMatrix(dataMat)

constructPairedMatrix <- function(dataMat){
  dataMat <- dataMat[,!duplicated(substr(colnames(dataMat), 9, 15))]
  pID <- substr(colnames(dataMat), 9, 12)
  dupInd <- pID %in% pID[duplicated(pID)]
  dataMat <- dataMat[,dupInd]
  pID <- pID[dupInd]
  TNind <- as.numeric(substr(colnames(dataMat), 14, 15))
  PairedSamples <- data.frame(row.names = rownames(dataMat))
  for(i in unique(pID)){
    index <- grep(i, pID)
    if(any(TNind[index] > 0 & TNind[index] < 10) & any(TNind[index] > 9 & TNind[index] < 20)){
      indexN <- index[which(TNind[index] > 0 & TNind[index] < 10)]
      indexT <- index[which(TNind[index] > 9 & TNind[index] < 20)]
      PairedSamples <- cbind(PairedSamples, dataMat[,c(indexN[1], indexT[1])])
    }
  }
  return(PairedSamples)
}

