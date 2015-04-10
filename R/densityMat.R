densityMat <- function(dataMatrix){
  x<-ggplot()
  for(i in 1:ncol(dataMatrix)){
    a<-as.data.frame(dataMatrix[,i])
    names(a)<-'x'
    x<-x+geom_density(aes(x=x),data=a,colour=i)
  }
}
