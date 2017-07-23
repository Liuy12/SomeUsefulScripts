densityMat <- function(dataMatrix){
  x<-ggplot()
  for(i in 1:ncol(dataMatrix)){
    a<-as.data.frame(dataMatrix[,i])
    names(a)<-'value'
    x<-x+geom_density(aes(x=value),data=a,colour=colors()[i])
  }
  return(x)
}
