###########################################################################################
################ plot of cluster centers and centered gene expression #####################
################ Hierarchical clustering on the centers of kmeans result with p value #####
###########################################################################################
#cluster plot function
#xbrk is the break points of x axis
clusterplot<-function(mfrow,xbrk,clusterres,data,cent.col='red',xlables=xbrk)
{
  par(mfrow=mfrow)
  for(i in 1:prod(mfrow))
  {
    plot(xbrk,clusterres$centers[i,],ylim=range(data[clusterres$cluster==i,]),xlab='Time in hours',ylab='Expression',main=paste('Cluster',i,':',clusterres$size[i]),lwd=1,type='l',xaxt='n')
    for(j in 1:nrow(data[clusterres$cluster==i,]))
    {
      lines(xbrk,data[clusterres$cluster==i,][j,])
    }
    lines(xbrk,clusterres$centers[i,],col=cent.col,lwd=1)
    axis(1,at=xbrk,labels=xlables)
  }
}