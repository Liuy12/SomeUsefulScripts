##Generate MAplot based on two dataset or prodived M and A values 
#This can be used after the DESeq procedure 
MyMAplot<-function(adjust.col=F,alpha=1,plot=T,data1=NA,data2=NA,Mr=NA,Ar=NA,conint.plot=T,bar.plot=F,col.point=1,pch=19,col.bar=NA,length.arrow=0.25,conint=c(0.99,0.95),mad=F,cex=0.45,pval=NA,pcoff=0.05)
{
  #adjust transparency of colors of the bar 
  if(adjust.col==T)
  {
    mycols<-adjustcolor(palette(),alpha.f=alpha)
    palette(mycols)
  } 
  #show error imformation if input data is not appropriate 
  if(missing(data1)&missing(data2))
  {
    if(missing(Mr)|missing(Ar))
      stop('at least two datasets should be provided')
    if(length(Mr)!=length(Ar))
      stop('the length of the input must be equal')
    #log2 transform the M and A values provided based on non-zero data
    M<-log2(Mr[Ar!=0&Mr!=0&is.finite(Mr)])
    A<-log2(Ar[Ar!=0&Mr!=0&is.finite(Mr)])
    if(!missing(pval))
      MA<-cbind(A,M,pval[Ar!=0&Mr!=0&is.finite(Mr)])
  }
  else{
    if(length(data1)!=length(data2))
      stop('the length of the input must be equal')
    #calcualte M and A values based on data provided
    M<-log2(data1[data1!=0&data2!=0])-log2(data2[data1!=0&data2!=0])
    A<-log2((data1[data1!=0&data2!=0]+data2[data1!=0&data2!=0])/2)
    if(!missing(pval))
      MA<-cbind(A,M,pval[data1!=0&data2!=0])
    else
      MA<-cbind(A,M)
  }
  #make a standard MAplot 
  if(!conint.plot&plot)
  {
    if(missing(pval))
      plot(A,M,pch=pch,xlab=NA,ylab=NA,col=col.point)
    else
      plot(A,M,pch=pch,xlab=NA,ylab=NA,col=ifelse(MA[,3]>pcoff,col.point,rgb(1,0,0,alpha=alpha)))
    abline(h=0,col=rgb(1,0,0,alpha=0.5),lwd=4)
  }
  #caculate mean and variance on a 0.1 intervel to make a confidence intervel plot 
  else{
    Simulation.mean<-c()
    Simulation.conint<-c()
    for(i in seq(round(min(A),1),round(max(A),1),0.1))
    {
      temp1<-mean(MA[MA[,1]>i&MA[,1]<i+0.1,2])
      Simulation.mean<-c(Simulation.mean,temp1)
      if(conint==0.99)
        #whether to use a mad method,mad=median(|(x-median(x))|) 
        if(mad)
          temp2<-mad(MA[MA[,1]>i&MA[,1]<i+0.1,2])*3
      else
        temp2<-sd(MA[MA[,1]>i&MA[,1]<i+0.1,2])*3
      if(conint==0.95)
        if(mad)
          temp2<-mad(MA[MA[,1]>i&MA[,1]<i+0.1,2])*2
      else
        temp2<-sd(MA[MA[,1]>i&MA[,1]<i+0.1,2])*2
      Simulation.conint<-c(Simulation.conint,temp2)
    }
    #replace na with 0 in interval gap
    index<-which(is.na(Simulation.conint))
    Simulation.mean<-replace(Simulation.mean,index,0)
    Simulation.conint<-replace(Simulation.conint,index,0)
    if(plot)
    {
      if(missing(pval))
        plot(A,M,pch=pch,xlab=NA,ylab=NA,col=col.point,cex=cex,ylim=range(min(M),max(M),max(Simulation.conint),min(-Simulation.conint)))
      else
        plot(A,M,pch=pch,xlab=NA,ylab=NA,col=ifelse(MA[,3]>pcoff,col.point,rgb(1,0,0,alpha=alpha)),cex=cex,ylim=range(min(M),max(M),max(Simulation.conint),min(-Simulation.conint)))
      abline(h=0,col=rgb(1,0,0,alpha=0.5),lwd=4)
    } 
    #generate barplot break points
    barp<-seq(min(round(A,1)),max(round(A,1)),0.1)
    if(!bar.plot&plot)
    {
      temp<-Simulation.conint
      temp1<-temp[-(1:(which(temp==max(temp))-1))]
      temp2<-temp1[temp1!=0]
      barp<-barp[-(1:(which(temp==max(temp))-1))]
      barp<-barp[temp1!=0]
      lo1<-loess(temp2~barp)
      lines(barp,predict(lo1),xlab=NA,ylab=NA,main=NA,col='green',lwd=2)
      temp<- -Simulation.conint
      temp1<-temp[-(1:(which(temp==min(temp))-1))]
      temp2<-temp1[temp1!=0]
      lo2<-loess(temp2~barp)
      lines(barp,predict(lo2),xlab=NA,ylab=NA,main=NA,col='green',lwd=2)
    }
    if(bar.plot&plot)
    {
      rect(barp,0,barp+0.08,Simulation.mean,col=col.bar)
      arrows(barp+0.04,0,barp+0.04,Simulation.conint,angle=90,code=2,length=length.arrow)
      arrows(barp+0.04,0,barp+0.04,-Simulation.conint,angle=90,code=2,length=length.arrow)
      temp<-Simulation.conint
      temp1<-temp[-(1:(which(temp==max(temp))-1))]
      temp2<-temp1[temp1!=0]
      barp<-barp[-(1:(which(temp==max(temp))-1))]
      barp<-barp[temp1!=0]
      lo1<-loess(temp2~barp)
      lines(barp,predict(lo1),xlab=NA,ylab=NA,main=NA,col='green',lwd=2)
      temp<- -Simulation.conint
      temp1<-temp[-(1:(which(temp==min(temp))-1))]
      temp2<-temp1[temp1!=0]
      lo2<-loess(temp2~barp)
      lines(barp,predict(lo2),xlab=NA,ylab=NA,main=NA,col='green',lwd=2)
    }
  }
}