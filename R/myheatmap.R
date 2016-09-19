heatmap.my <- function(Exprs, sel=F, thres_mean, thres_var, numbreaks=100, col = c("blue","white","red"), 
                       breakratio = c(2,1,2), colsidebar, Colv=F, Rowv=T, scale= 'row', labRow=F, 
                       labCol=F, dendrogram = 'row'){
  suppressPackageStartupMessages(invisible(require('gplots', quietly=TRUE)))
  if(labRow)
    labRow <- rownames(Exprs)
  if(labCol)
    labCol <- colnames(Exprs)
  if(sel){
    gene_mean <- apply(Exprs,1,mean)
    gene_var <- apply(Exprs,1,var)
    Exprs <- Exprs[gene_mean>thres_mean & gene_var>thres_var,]
  }
  if(scale == 'row')
    Exprs_scale <- t(scale(t(Exprs)))
  else
    Exprs_scale <- Exprs
  Exprs_scale[is.na(Exprs_scale)] <- min(Exprs_scale, na.rm = TRUE)
  # lmat is a matrix describing how the screen is to be broken up. By default, heatmap.2 divides the screen into a four element grid, so lmat is a 2x2 matrix. 
  # The number in each element of the matrix describes what order to plot the next four plots in. Heatmap.2 plots its elements in the following order:
  # 1 Heatmap,
  # 2 Row dendrogram,
  # 3 Column dendrogram,
  # 4 Key
  if(missing(colsidebar)){
  lmat <- rbind(c(0,4), c(0,3), c(2,1))
  lwid <- c(1,4)
  lhei <- c(1,0.1,4)
  if(class(Colv) == 'dendrogram'){
    lhei <- c(1,1,4)
    dendrogram <- 'both'
  }
  }
  else{
    if(class(Colv) == 'dendrogram'){
      # 4 is column dendrogram, 5 is key, 1 is colcolorkey
      lmat <- rbind(c(0,5),c(0,4), c(3,2),c(0,1))
      lwid <- c(1,4)
      lhei <- c(1,1, 4,0.25)
      dendrogram <- 'both'
    }
    else{
      if(Colv){
        lmat <- rbind(c(0,5),c(0,4), c(3,2),c(0,1))
        lwid <- c(1,4)
        lhei <- c(1,1, 4,0.25)
        dendrogram <- 'both'
      }
      lmat <- rbind(c(0,5),c(0, 1), c(3,2),c(0,4))
      lwid <- c(1,4)
      lhei <- c(1,0.25,4,0.1)
    }

  }
  rg <- quantile(Exprs_scale,na.rm=T)
  rg_diff <- rg[4]-rg[2]
  rg_max <- max(abs(rg))
  Exprs_sd <- sd(Exprs_scale)
  Exprs_mean <- mean(Exprs_scale)
  if(rg_max > max(abs(c(Exprs_mean + 3*Exprs_sd, Exprs_mean - 3*Exprs_sd)))){
    rg_iqr <- max(abs(c(rg[2], rg[4])))
    bp <- c((breakratio[1]/sum(breakratio))*rg_diff - rg_iqr, rg_iqr - (breakratio[3]/sum(breakratio))*rg_diff)
    bk <- unique(c(seq(-rg_max, -rg_iqr, length= numbreaks), seq(-rg_iqr,bp[1],length = numbreaks), seq(bp[1],bp[2],length=numbreaks),seq(bp[2],rg_iqr,length=numbreaks), 
            seq(rg_iqr, rg_max, length = numbreaks)))
    hmcols<- colorRampPalette(col)(length(bk)-1)
  }
  else{
    rg <- range(Exprs_scale, na.rm=T)
    bp <- c((breakratio[1]/sum(breakratio))*diff(rg) - rg_max, rg_max - (breakratio[3]/sum(breakratio))*diff(rg))
    bk <- c(seq(-rg_max,bp[1],length=numbreaks), seq(bp[1],bp[2],length=numbreaks),seq(bp[2],rg_max,length=numbreaks))
    bk <- bk[!duplicated(bk)]
    hmcols<- colorRampPalette(col)(length(bk)-1)
  }
  heatmap.2(Exprs, Colv=Colv,Rowv=Rowv, dendrogram = dendrogram,trace='none',scale=scale ,density.info='none',
            lmat=lmat,lwid=lwid,lhei=lhei,labRow=labRow,labCol=labCol,col=hmcols,breaks=bk, 
            ColSideColors=colsidebar) 
}
