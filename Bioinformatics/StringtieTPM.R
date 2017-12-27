library(plyr)
library(dplyr)
library(doMC)

getTPMall <- function(path){
  gtffiles <- list.files(path, '*mergedTrans.gtf', recursive = T)
  #filenames <- sapply(strsplit(gtffiles, '/'), function(i) i[1])
  registerDoMC(detectCores()/4)
  tpmMat <- foreach(i=1:length(gtffiles)) %dopar% getTPM(gtffiles[i])
  # for(i in 1:length(gtffiles)){
  #   assign(filenames[i], getTPM(gtffiles[i]))
  # }
  #filenames <- lapply(filenames, function(i) get(i))
  #mergedata <- join_all(filenames, type = 'left', by = 'id')
  mergedata <- join_all(tpmMat, type = 'left', by = 'id')
  return(mergedata)
}

getTPM <- function(dataMat){
  cat(dataMat, '\n')
  cat('loading gtf file', '\n')
  dataMat <- read.table(dataMat, sep = '\t', skip=2, stringsAsFactors=F)
  dataMat <- dataMat %>% filter(V3 == 'transcript')
  tpmMat <- as.data.frame(matrix('NA', nrow = nrow(dataMat), ncol=2), stringsAsFactors=F)
  colnames(tpmMat) <- c('id', 'tpm')
  for(i in 1:nrow(dataMat)){
    cat(i, '\n')
    idx <- grep('transcript_id', dataMat$V9[i])
    if(length(idx)){
      str1 <- strsplit(dataMat$V9[i], ';')
      idx <- grep('transcript_id', str1[[1]])
      str2 <- gsub('\\s*transcript_id\\s*','', str1[[1]][idx])
      idx1 <- grep('TPM', str1[[1]])
      str3 <- as.numeric(gsub('\\s*TPM\\s*','', str1[[1]][idx1]))
      data.table::set(tpmMat, i, 1L:2L, as.list(c(str2, str3)))
    }
  }
  tpmMat <- tpmMat[tpmMat$id != 'NA',]
  tpmMat$tpm <- as.numeric(tpmMat$tpm)
  return(tpmMat)
}


t1 <- Sys.time()
tmp <- getTPM(gtffiles[1])
t2 <- Sys.time()
t2 -t1

t1 <- Sys.time()
tmp <- getTPM1(gtffiles[1])
t2 <- Sys.time()
t2 -t1

t1 <- Sys.time()
tmp <- getTPMall(path)
t2 <- Sys.time()
t2 -t1





