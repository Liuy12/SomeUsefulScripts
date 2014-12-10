## Correlation test for miRNA and its target genes
# miRNA with enrichment score from Toppgene
miRNAcandi<-function(clustercenter,miRNAenrich,miRNAexpr){
  temp1<-c()
  for(i in 1:nrow(miRNAenrich))
  {
    temp<-grep(miRNAenrich[i,3],rownames(miRNAexpr))
    if(length(temp)==1)
    {
      temp2<-miRNAcor(rowindex=temp,a=miRNAexpr,b=clustercenter,c=miRNAenrich,i=i)
      temp1<-rbind(temp1,temp2)
    }
    if(length(temp)==2)
    {
      temp2<-miRNAcor(rowindex=temp[1],a=miRNAexpr,b=clustercenter,c=miRNAenrich,i=i)
      temp1<-rbind(temp1,temp2)
      temp2<-miRNAcor(rowindex=temp[2],a=miRNAexpr,b=clustercenter,c=miRNAenrich,i=i)
      temp1<-rbind(temp1,temp2)
    }
  }
  return(temp1)
}
miRNAcor<-function(rowindex,a,b,c,i)
{
  cortest<-cor.test(as.numeric(a[rowindex,]),b)
  if(!is.na(cortest$estimate)&cortest$estimate< -0.8& cortest$p.value<0.05)
  {
    temp2<-matrix(ncol=4)
    rownames(temp2)<-c[i,3]
    temp2[,1]<-c[i,5]
    temp2[,2]<-cortest$estimate
    temp2[,3]<-cortest$p.value
    temp2[,4]<-as.character(c[i,4])
    return(temp2)
  }
}