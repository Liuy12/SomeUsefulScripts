heatmap.my <- function(Exprs, sel=F, thres_mean, thres_var, numbreaks=100, col = c("blue","white","red"), 
                       breakratio = c(2,1,2)){
  suppressPackageStartupMessages(invisible(require('gplots', quietly=TRUE)))
  options(warn=-1)
  if(sel){
    gene_mean <- apply(Exprs,1,mean)
    gene_var <- apply(Exprs,1,var)
    Exprs <- Exprs[gene_mean>thres_mean & gene_var>thres_var,]
  }
  Exprs_scale <- t(scale(t(Exprs)))
  # lmat is a matrix describing how the screen is to be broken up. By default, heatmap.2 divides the screen into a four element grid, so lmat is a 2x2 matrix. The number in each element of the matrix describes what order to plot the next four plots in. Heatmap.2 plots its elements in the following order:
  # 1 Heatmap,
  # 2 Row dendrogram,
  # 3 Column dendrogram,
  # 4 Key
  lmat <- rbind(c(0,4),c(2,1),c(0,3))
  lwid <- c(1,4)
  lhei <- c(1,4,0.1)
  rg <- range(Exprs_scale)
  bp <- c(rg[1]+(breakratio[1]/sum(breakratio))*diff(rg), rg[2] - (breakratio[3]/sum(breakratio))*diff(rg))
  bk <- c(seq(rg[1],bp[1],length=numbreaks), seq(bp[1],bp[2],length=numbreaks),seq(bp[2],rg[2],length=numbreaks))
  hmcols<- colorRampPalette(col)(length(bk)-1)
  heatmap.2(Exprs, Colv=F,dendrogram = 'row',trace='none',scale='row',density.info='none',lmat=lmat,lwid=lwid,lhei=lhei,labRow=F,col=hmcols,breaks=bk) 
  options(warn=0)
}
