PCAcovariate <- function(normdata, sampleinfo, output, npcs){
  library(ggplot2)
  library(heplots)
  library(gridExtra)
  pcares <- Principalstats(normdata, method = 'pca', cvCutoff = 0.1, npcs=npcs)
  gps <- vector(mode = "list", length = npcs*ncol(sampleinfo))
  for(i in 1:npcs){
    for(j in 1:ncol(sampleinfo)){
      if(is.numeric(sampleinfo[[j]])) cor1 <- cor.test(pcares[[1]][,i],sampleinfo[[j]],method = 's') else {
        df <- data.frame(y=pcares[[1]][,i],x=sampleinfo[[j]])
        mod1 <- aov(y~x,data=df)
        cor1 <- data.frame(p.value= summary(mod1)[[1]][1,5],estimate=sqrt(etasq(mod1,partial = F)[1,1]))
      }
      if(round(cor1$p.value,3) <=0.01) labcol <- "red" else labcol <- 'black'
      df <- data.frame(x=sampleinfo[[j]], y = pcares[[1]][,i],stringsAsFactors = F)
      gp <- ggplot() + geom_point(aes(x=x,y=y),data=df) +
        theme_classic() + labs(x=colnames(sampleinfo)[j],y=paste0('PC',i,' (',pcares[[2]][i],'%)')) +
        theme(axis.title = element_text(face = 'bold',size=15),axis.text = element_text(size=10)) +
        annotate("text",  x=Inf, y = Inf,color=labcol,size=10,label = paste0('cor: ',formatC(cor1$estimate,digits = 2, format = 'e'), '; pval: ',formatC(cor1$p.value,digits = 2, format = 'e')), vjust=1, hjust=1)
      gps[[(ncol(sampleinfo)*(i-1) + j)]] <- gp
    }
  }
  pdf(paste0(output, "/PCAcovariate.pdf"),useDingbats = F,height = npcs*6,width = ncol(sampleinfo)*6)
  grid.arrange(grobs=gps,nrow=npcs,ncol=ncol(sampleinfo))
  dev.off()
  png(paste0(output, "/PCAcovariate.png"),height = npcs*6,width = ncol(sampleinfo)*6,res = 300, units = 'in')
  grid.arrange(grobs=gps,nrow=npcs,ncol=ncol(sampleinfo))
  dev.off()
}


Principalstats <- function(dataMat, method, cvCutoff, npcs){
  dataMat <- log2(dataMat)
  cv.gene <- apply(dataMat, 1, function(x)
    sd(x) / mean(x))
  dataMat <- dataMat[which(cv.gene > cvCutoff),]
  dataMat <- scale(dataMat)
  if (method == 'mds') {
    dd <- dist(t(dataMat))
    mds.result <- cmdscale(dd, k = npcs, eig = TRUE)
    ppoints <- mds.result$points
    eig <- mds.result$eig
    percent <- round(eig/sum(eig) * 100, 1)
  }
  else{
    pca.result <- prcomp(t(dataMat))
    ppoints <- pca.result$x[,1:npcs]
    percent<-round((pca.result$sdev^2/sum(pca.result$sdev^2))*100,1)
  }
  ppoints <- as.data.frame(ppoints)
  return(list(ppoints, percent))
}

