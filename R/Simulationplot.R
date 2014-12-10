Simulationplot<-function(data1,data2,plottp,frame=c(1,1),misseddata=NA,NumberofReads=NA,Genelength=NA)
{
  par(mfrow=frame)
  for(i in 1:sum(frame))
  {
    Simulationdata<-data1[,i]
    MissInIntron<-misseddata
    Simulationdata_intron<-data2[,i]
    #RPKM norm
    #Log2 transform after
    Simulationdata_RPKM<-Simulationdata/NumberofReads[2*i,7]
    Simulationdata_RPKM<-Simulationdata_RPKM/as.numeric(Genelength[,2])
    Simulationdata_RPKM<-10^9*Simulationdata_RPKM
    Simulationdata_intron_RPKM<-Simulationdata_intron/NumberofReads[2*i,7]
    Simulationdata_intron_RPKM<-Simulationdata_intron_RPKM/as.numeric(Genelength[-MissInIntron,2])
    Simulationdata_intron_RPKM<-10^9*Simulationdata_intron_RPKM
    if(plottp=='hist')
    {
      hist(log2(Simulationdata_RPKM),breaks=400,xlab=NA,main=NA,col='blue4',xlim=c(-10,15))
      hist(log2(Simulationdata_intron_RPKM),breaks=400,xlab=NA,main=NA,col=rgb(red=0,green=1,blue=0,alpha=0.4),add=T)
    }
    if(plottp=='density')
    {
      plot(density(log2(Simulationdata_RPKM)),col='blue',xlab=NA,main=NA)
      lines(density(log2(Simulationdata_intron_RPKM)),col='green',xlab=NA,main=NA)
      mean.intron<-mean(Simulationdata_intron_RPKM)
      variance.intron<-var(Simulationdata_intron_RPKM)
      lines(function(x)dnorm(x,mean=mean.intron,sd=sqrt(variance.intron),col='red',xlan=NA,main=NA)
    }
    if(plottp=='pairwise')
      
    {
      par(pch=16)
      plot(log2(Simulationdata_RPKM[-MissInIntron]),log2(Simulationdata_intron_RPKM),type='p',xlab=NA,main=NA,ylab=NA)
    }
  }
}