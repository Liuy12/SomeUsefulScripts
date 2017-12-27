args =  commandArgs(TRUE)
filePath=args[1]
sampleIn=args[2]
control=args[3]
case=args[4]
exprCutoff=args[5]
fc=args[6]
pVal=args[7]
pFlag=args[8]
paired=args[9]
outPath=args[10]
exprCutoff=as.numeric(exprCutoff)
fc=as.numeric(fc)
pVal=as.numeric(pVal)
pFlag=as.numeric(pFlag)
paired=as.numeric(paired)


##### usage
if (file.exists(filePath) == FALSE){
  writeLines ("#################################################################################
              \n\nUsage:\nRscript  runEdgeR.R GeneCount_raw.tsv sampleIn controlName caseName  
              exprCutoff foldChangeCutoff pValCutoff pFlag pairFlag outputDir
              
              GeneCount_raw.tsv: gene count matrix output by MAPRSeq
              sampleIn: txt file containing sample information
              controlName: character representing control samples
              caseName: character representing case samples
              exprCutoff: expression cutoff value for DEGs
              foldChangeCutoff: fold change cutoff for DGEs
              pValCufoff: p value cutoff for DEGs
              pFlag: whether to use padjust value or raw p value
              pairFlag: whether this is a paired comparison
              outputDir: output directory
              #################################################################################\n\nFor sampleIn file - Please make file in this format:\nsamples\tdiagnosis\nsample1_name\tcontrol1\nsample2_name\tcontrol1\nsample3_name\tcase1\nsample4_name\tcase1\nsample5_name\tcontrol2\nsample6_name\tcontrol2\nsample7_name\tcase2\nsample8_name\tcase2\n\n\n#################################################################################\n\n1. Please make sure controlName and CaseName match the diagnosis column in sampleIn file\n2. Please make sure number of rows in sampleIn file should equal to the total number of samples. If you have more than two conditions, this should include all conditions rather than just the conditions specified as case and control\n3. You can compare any two sets of control and case using this script\n\n");
  quit()
}

t1 <- Sys.time()
#### loading required libraries
source('http://bioconductor.org/biocLite.R')
##'threejs', 'htmlwidgets', 'metricsgraphics', 'rCharts'
requiredPacks <- c('edgeR', 'data.table', 'gplots', 'ggplot2', 'threejs', 'htmlwidgets', 'metricsgraphics', 'rCharts')
for(i in requiredPacks){
  if(!require(i, character.only = T)) biocLite(i)
  library(i, character.only = T)
}

#############################################################
#### Functions
heatmap.my <- function(Exprs, sel=F, thres_mean, thres_var, numbreaks=100, col = c("blue","white","red"), 
                       breakratio = c(2,1,2), colsidebar, Colv=F, Rowv=T, scale= 'row', labRow=F, 
                       labCol=F, dendrogram = 'row'){
  suppressPackageStartupMessages(invisible(require('gplots', quietly=TRUE)))
  if(labRow)
    labRow <- rownames(Exprs)
  if(labCol)
    labCol <- paste0('sample', 1:ncol(Exprs))
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


Principalstats <- function(dataMat, design, method, cvCutoff){
  dataMat <- log2(dataMat)
  cv.gene <- apply(dataMat, 1, function(x)
    sd(x) / mean(x))
  dataMat <- dataMat[which(cv.gene > cvCutoff),]
  dataMat <- scale(dataMat)
  if (method == 'mds') {
    dd <- dist(t(dataMat))
    mds.result <- cmdscale(dd, k = 3, eig = TRUE)
    ppoints <- mds.result$points
    eig <- mds.result$eig
    percent <- round(eig/sum(eig) * 100, 1)
  }
  else{
    pca.result <- prcomp(t(dataMat))
    ppoints <- pca.result$x[,1:3]
    percent<-round((pca.result$sdev^2/sum(pca.result$sdev^2))*100,1)
  }
  ppoints <-
    cbind(as.data.frame(ppoints), design, paste('S', 1:ncol(dataMat), sep = ''))
  colnames(ppoints) <- c('PC1', 'PC2', 'PC3', 'Design', 'Samplename')
  ppoints <- as.data.frame(ppoints)
  return(list(ppoints, percent))
} 


MAplot <- function(stats, ylim, padj=1, pcuff=0.1, lfccuff=1, exprCutoff,linecol='red3',
                   xlab='mean of normalized counts', ylab=expression(log[2]~fold~change), shape, interactive=T, geneName, outPath, case, control)
{
  if(padj==1)
    col = ifelse(stats$padj<=pcuff&abs(stats$log2FoldChange)>=lfccuff&(stats$log2baseMeanControl>=log2(exprCutoff)|stats$log2baseMeanCase>=log2(exprCutoff)), "red", "gray32")
  else if(padj==2)
    col = ifelse(stats$padj1<=pcuff&abs(stats$log2FoldChange)>=lfccuff&(stats$log2baseMeanControl>=log2(exprCutoff)|stats$log2baseMeanCase>=log2(exprCutoff)), "red", "gray32")
  else
    col = ifelse(stats$pval<=pcuff&abs(stats$log2FoldChange)>=lfccuff&(stats$log2baseMeanControl>=log2(exprCutoff)|stats$log2baseMeanCase>=log2(exprCutoff)), "red", "gray32")
  y = stats$log2FoldChange
  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(y[is.finite(y)]), probs=0.99) * 1.1
  if (missing(shape))
    shape = ifelse(y<ylim[1], 6, ifelse(y>ylim[2], 2, 16) )
  stats$log2FoldChange1 = pmax( ylim[1], pmin(ylim[2], y) )
  stats$geneName <- geneName
  gp <- ggplot() + geom_point( data=stats,aes( x=log2baseMean, y=log2FoldChange1 ), color=col, shape=shape ) +
    ylim(ylim) + geom_hline(yintercept=0,colour=linecol,size=1)  +
    labs( x=xlab, y=ylab )
  if(length(unique(col)) == 1)
    color_rg <- 'gray32'
  else
    color_rg <- c('red', 'gray32')
  col[col == 'red'] <- 'DE'
  col[col == 'gray32'] <- 'Not DE'
  stats$col <- col
  mp <-
    metricsgraphics::mjs_plot(stats, log2baseMean, log2FoldChange, decimals = 2) %>%
    metricsgraphics::mjs_point(
      color_accessor = col, color_range = color_rg, color_type = "category", x_rug =
        FALSE, y_rug = FALSE
    ) %>%
    metricsgraphics::mjs_add_baseline(y_value = 0, label = 'baseline') %>%
    metricsgraphics::mjs_labs(x_label = 'log2 mean of normlized counts', y_label = "Log2 fold change") %>%
    metricsgraphics::mjs_add_mouseover("function(d) {
                                       $('{{ID}} svg .mg-active-datapoint')
                                       .text('Gene Name: ' +  d.point.geneName + ',' + ' Log2 intensity: ' + d.point.log2baseMean + ',' + ' Log2 fold change: ' + d.point.log2FoldChange);
}")
  ggsave(paste0(outPath, case, '_vs_', control, '_MAplot.pdf'), gp)
  ggsave(paste0(outPath, case, '_vs_', control, '_MAplot.png'), gp)
  #### MAplot interactive.
  htmlwidgets::saveWidget(mp, file = paste0(outPath, case, '_vs_', control, '_MAplot.html'), background = 'none')
}
#### end of functions



################################# main script
#### load datasets
message(paste0('Loading datasets!', '\n'))
dataMat <- fread(filePath, data.table = F, stringsAsFactors = F, check.names = F)
sampleInfo <- fread(sampleIn, data.table = F, stringsAsFactors = F, check.names = F)
s1 <- sapply(sampleInfo$samples, function(i) grep(paste0('^', i, '$'), colnames(dataMat)))
inforMat <- dataMat[,1:(ncol(dataMat)-length(s1))]
dataMat <- dataMat[,(ncol(dataMat)-length(s1)+1):ncol(dataMat)]
if(length(intersect(sampleInfo[,1], colnames(dataMat))) != ncol(dataMat)){
  stop('Check sample info file!')
}
  
#### filter genes based on abundance
message(paste0('Filtering genes based on abundance!', '\n'))
cat('number of genes before filtering: ', nrow(dataMat), '\n')

idx1 <- grep(paste0('^', control, '$'), sampleInfo$diagnosis)
idx1 <- sapply(sampleInfo$samples[idx1], function(i) grep(paste0('^', i, '$'), colnames(dataMat)))
idx2 <- grep(paste0('^', case, '$'), sampleInfo$diagnosis)
idx2 <- sapply(sampleInfo$samples[idx2], function(i) grep(paste0('^', i, '$'), colnames(dataMat)))
isExpr <- apply(cpm(dataMat[,idx1]), 1, function(x) mean(x) > 2) | apply(cpm(dataMat[,idx2]), 1, function(x) mean(x) > 2)
dataMat <- dataMat[isExpr,]
inforMat <- inforMat[isExpr,]
cat('number of genes after filtering: ', nrow(dataMat), '\n')

#### run edgeR
message(paste0('Carrying out differential expression analysis using edgeR!', '\n'))
cds <- DGEList(counts = dataMat, group = sampleInfo$diagnosis)
cds <- calcNormFactors(cds)
if(paired==1){
  design <- model.matrix(~pairinfo+diagnosis, data = sampleInfo)
  cds <- estimateGLMRobustDisp(cds,design)
  fit <- glmFit(cds, design)
  lrt <- glmLRT(fit)
}else{
  design <- model.matrix(~0+diagnosis, data = sampleInfo)
  
  colnames(design) <- levels(cds$samples$group)
  cds <- estimateGLMRobustDisp(cds,design)
  fit <- glmFit(cds, design)
  lrt<- eval(parse(text = paste0("glmLRT(fit,contrast=makeContrasts(",case,"-",control,",levels=design))")))
}
#### print top DEGs
cat(paste0('Top DEGs:', '\n'))
topTags(lrt)

#### box plot 
message(paste0('Generating figures!', '\n'))
normCounts <- cpm(cds)
tmp <- stack(as.data.frame(log2(normCounts)))
gp <- ggplot(tmp) +
  geom_boxplot(aes(x = ind, y = values)) +
  labs(x='', y='Log2 normalized CPM (counts per million)') +
  scale_x_discrete(labels=c(paste0('sample', 1:ncol(dataMat))))
ggsave(paste0(outPath, case, '_vs_', control, '_boxplot.pdf'), gp)
ggsave(paste0(outPath, case, '_vs_', control, '_boxplot.png'), gp)

#### correlation of samples
tmp2 <- cor(log2(normCounts+0.1))
pdf(paste0(outPath, case, '_vs_', control, '_corrplot.pdf'))
heatmap.2(tmp2, Rowv = F, Colv = F, labCol = paste('sample', 1:ncol(tmp2), sep = ''),labRow = paste('sample', 1:ncol(tmp2), sep = ''), dendrogram = 'none', scale = 'none', trace = 'none', density.info = 'none', col=colorRampPalette(c("blue","white","red")))
dev.off()
png(paste0(outPath, case, '_vs_', control, '_corrplot.png'))
heatmap.2(tmp2, Rowv = F, Colv = F, labCol = paste('sample', 1:ncol(tmp2), sep = ''),labRow = paste('sample', 1:ncol(tmp2), sep = ''), dendrogram = 'none', scale = 'none', trace = 'none', density.info = 'none', col=colorRampPalette(c("blue","white","red")))
dev.off()

#### heatmap of high CV genes
highcv <- order(apply(log2(normCounts), 1, function(x) mean(x)/sd(x)), decreasing = T)[1:1000]
pdf(paste0(outPath, case, '_vs_', control, '_heatmapHighcv.pdf'))
heatmap.my(Exprs = as.matrix(log2(normCounts[highcv,])), labCol = T, Colv = F, dendrogram = 'row')
dev.off()
png(paste0(outPath, case, '_vs_', control, '_heatmapHighcv.png'))
heatmap.my(Exprs = as.matrix(log2(normCounts[highcv,] + 0.1)), labCol = T, Colv = F, dendrogram = 'row')
dev.off()

#### hierarchy clustering
hc <- hclust(as.dist(1-cor(log2(normCounts[highcv,] + 0.1))), method ='average')
hc$call <- NULL
pdf(paste0(outPath, case, '_vs_', control, '_hclustHighcv.pdf'))
plot(hc)
dev.off()
png(paste0(outPath, case, '_vs_', control, '_hclustHighcv.png'))
plot(hc)
dev.off()

#### PCA plot
pcaStats <- Principalstats(normCounts, sampleInfo$diagnosis, 'pca', 0.1)
gp <- ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = Design), data = as.data.frame(pcaStats[[1]])) +
  labs(x = paste0("Principal component 1 (", pcaStats[[2]][1], '%)'), y = paste0("Principal component 2 (", pcaStats[[2]][2], '%)')) +
  annotate("text", x = pcaStats[[1]][,1], y = pcaStats[[1]][,2], label = paste0('Sample', 1:nrow(as.data.frame(pcaStats[[1]]))))
ggsave(paste0(outPath, case, '_vs_', control, '_pcaplot.pdf'), gp)
ggsave(paste0(outPath, case, '_vs_', control, '_pcaplot.png'), gp)

### 3d pca plot interactive
xlab <- paste('PC1 (', pcaStats[[2]][1],'%)', sep = '')
ylab <- paste('PC2 (', pcaStats[[2]][2], '%)', sep = '')
zlab <- paste('PC3 (', pcaStats[[2]][3], '%)', sep = '')
colors <-
  c(
    '#00FFFF', '#FFE4C4', '#D2691E', '#6495ED', '#9932CC', '#8B0000',
    '#FF00FF', '#FFD700'
  )
if(length(unique(pcaStats[[1]]$Design)) > 8)
  colors <- rainbow_hcl(length(unique(pcaStats[[1]]$Design)))
col <- colors[1:length(unique(pcaStats[[1]]$Design))]
col <- rep(col, c(as.numeric(table(pcaStats[[1]]$Design))))
colnames(pcaStats[[1]])[1:3] <- c(xlab, ylab, zlab)
scatter3d <- threejs::scatterplot3js(as.matrix(pcaStats[[1]][,1:3]),
                                     labels = colnames(dataMat),
                                     color = col,
                                     renderer = 'auto') # bg = '#ecf0f5'
htmlwidgets::saveWidget(scatter3d, paste0(outPath, case, '_vs_', control, '_pcaplot3d.html'), background = 'none')

#### MA plot
idx1 <- grep(paste0('^', control, '$'), sampleInfo$diagnosis)
idx1 <- sapply(sampleInfo$samples[idx1], function(i) grep(paste0('^', i, '$'), colnames(normCounts)))
idx2 <- grep(paste0('^', case, '$'), sampleInfo$diagnosis)
idx2 <- sapply(sampleInfo$samples[idx2], function(i) grep(paste0('^', i, '$'), colnames(normCounts)))

testStats <- data.frame(log2baseMean = apply(log2(normCounts[,c(idx1, idx2)]+0.1), 1, mean),
                        log2baseMeanControl = apply(log2(normCounts[,idx1]+0.1), 1, mean),
                        log2baseSDControl = apply(log2(normCounts[,idx1]+0.1), 1, sd),
                        log2baseMeanCase = apply(log2(normCounts[,idx2]+0.1), 1, mean),
                        log2baseSDControl = apply(log2(normCounts[,idx2]+0.1), 1, sd),
                        log2FoldChange = lrt$table$logFC,
                        Foldchange = 2^(lrt$table$logFC),
                        pval = lrt$table$PValue,
                        padj = p.adjust(lrt$table$PValue, method = 'BH')
                        )
tmp <- rep(NA, nrow(testStats))
tmp[testStats$log2baseMean >= log2(exprCutoff)] <- p.adjust(testStats$pval[(testStats$log2baseMeanControl>=log2(exprCutoff)|testStats$log2baseMeanCase>=log2(exprCutoff))], method = 'BH')
testStats$padj1 <- tmp
testStats$t.pval <- sapply(1:nrow(normCounts), function(i) {tmp <- try(t.test(log2(normCounts[i,idx2]+0.1), log2(normCounts[i,idx1]+0.1))); if(class(tmp) =='htest') tmp$p.value else NA})
testStats[,1:7] <- round(testStats[,1:7], digits = 3)
MAplot(testStats, padj=pFlag, pcuff=pVal, lfccuff=log2(fc), exprCutoff = exprCutoff, linecol='red3',
                      xlab='mean of normalized CPM (counts per million)', ylab=expression(log[2]~fold~change), geneName=inforMat$GeneName, outPath = outPath, case = case, control = control)

#### Volcano plot
if(pFlag==1) idx <- which(abs(testStats$log2FoldChange)>=log2(fc)&testStats$padj<=pVal&testStats$log2baseMean>=log2(exprCutoff)) else if(pFlag==2) 
  idx <- which(abs(testStats$log2FoldChange)>=log2(fc)&testStats$padj1<=pVal&(testStats$log2baseMeanControl>=log2(exprCutoff)|testStats$log2baseMeanCase>=log2(exprCutoff))) else 
    idx <- which(abs(testStats$log2FoldChange) >= log2(fc) & testStats$pval <= pVal & (testStats$log2baseMeanControl>=log2(exprCutoff)|testStats$log2baseMeanCase>=log2(exprCutoff))) 
    
DE <- rep('None', nrow(testStats))
DE[idx[testStats$log2FoldChange[idx] >= log2(fc)]] <- 'Up'
DE[idx[testStats$log2FoldChange[idx] <= log2(fc)]] <- 'Down'
DE <- as.factor(DE)
testStats['DE'] <- DE
if(pFlag==1) testStats['logpval'] <- -1*log10(testStats$padj) else if(pFlag==2)
  testStats['logpval'] <- -1*log10(testStats$padj1) else
    testStats['logpval'] <- -1*log10(testStats$pval)

if(length(unique(DE)) ==1) 
  color <- 'black' else if('Up' %in% DE & 'Down' %in% DE) 
    color <- c('green', 'black', 'red') else if('Up' %in% DE & !('Down' %in% DE)) 
      color <- c('black', 'red') else
        color <- c('green', 'black')
  
gp <- ggplot() +
  geom_point(aes(x = log2FoldChange, y = logpval, col = DE), data = testStats) +
  scale_color_manual(values = color) +
  labs(title='Volcano plot', x='Log2 fold change', y='-Log10 p vlaue')
ggsave(paste0(outPath, case, '_vs_', control, '_Volcanoplot.pdf'))
ggsave(paste0(outPath, case, '_vs_', control, '_Volcanoplot.png'))

testStats$GeneName <- inforMat$GeneName
# if(length(unique(DE)) == 1)
#   color_rg <- 'grey32'else color_rg <- c('blue', 'grey32','red')
mp <-
  metricsgraphics::mjs_plot(testStats, log2FoldChange, logpval, decimals = 2) %>%
  metricsgraphics::mjs_point(
    color_accessor = DE, color_range = color, color_type = "category", x_rug =
      FALSE, y_rug = FALSE
  ) %>%
  metricsgraphics::mjs_labs(x_label = 'log2 fold change', y_label = "-log10 p value") %>%
  metricsgraphics::mjs_add_mouseover("function(d) {
                                     $('{{ID}} svg .mg-active-datapoint')
                                     .text('Gene Name: ' +  d.point.GeneName + ',' + ' Log2 fold change: ' + d.point.log2FoldChange + ',' + ' Log10 p value: ' + d.point.logpval);
                                     }")
  htmlwidgets::saveWidget(mp, file = paste0(outPath, case, '_vs_', control, '_Volcanoplot.html'), background = 'none')


## scatter plot
gp <- ggplot() +
  geom_point(aes(x = log2baseMeanControl, y = log2baseMeanCase, col = DE), data = testStats) +
  scale_color_manual(values = color) +
  labs(title='scatter plot', x=paste('log2 ',control, ' CPM', sep =' '), y=paste('log2', case, ' CPM', sep = ' '))
ggsave(paste0(outPath, case, '_vs_', control, '_scatterplot.pdf'))
ggsave(paste0(outPath, case, '_vs_', control, '_scatterplot.png'))

# if(length(unique(tmp3$DE)) == 1)
#   color_rg <- 'grey32'else color_rg <- c('blue', 'grey32','red')
mp <-
  metricsgraphics::mjs_plot(testStats, log2baseMeanControl, log2baseMeanCase, decimals = 2) %>%
  metricsgraphics::mjs_point(
    color_accessor = DE, color_range = color, color_type = "category", x_rug =
      FALSE, y_rug = FALSE
  ) %>%
  metricsgraphics::mjs_labs(x_label = paste('log2 ',control, ' CPM', sep = ' '), y_label = paste('log2', case, ' CPM', sep = ' ')) %>%
  metricsgraphics::mjs_add_mouseover("function(d) {
                                     $('{{ID}} svg .mg-active-datapoint')
                                     .text('Gene Name: ' +  d.point.GeneName + ',' + ' type: ' + d.point.DE);
                                     }")
  htmlwidgets::saveWidget(mp, file = paste0(outPath, case, '_vs_', control, '_scatterplot.html'), background = 'none')
  
#### cutoff table
pLevel <- c(1, 0.1, 0.05, 0.01)
fdrLevel <- c(0.1, 0.05, 0.01, 0.001)
fcLevel <- c(0, log2(1.5), log2(2), log2(3), log2(4))
out <- matrix(data = NA, nrow = 8, ncol = 5)
for (i in 1:4) {
  for (j in 1:5) {
    tdx <- which(abs(testStats$log2FoldChange) >= fcLevel[j] & testStats$pval < pLevel[i] & testStats$log2baseMean >= log2(exprCutoff))
    out[i, j] <- length(tdx)
  }
}
if(pFlag ==1){
  for (i in 1:4) {
    for (j in 1:5) {
      tdx <- which(abs(testStats$log2FoldChange) >= fcLevel[j] & testStats$padj < fdrLevel[i] & (testStats$log2baseMeanControl>=log2(exprCutoff)|testStats$log2baseMeanCase>=log2(exprCutoff)))
      out[i+4, j] <- length(tdx)
    }
  }
} else if(pFlag==2){
  for (i in 1:4) {
    for (j in 1:5) {
      tdx <- which(abs(testStats$log2FoldChange) >= fcLevel[j] & testStats$padj1 < fdrLevel[i] & (testStats$log2baseMeanControl>=log2(exprCutoff)|testStats$log2baseMeanCase>=log2(exprCutoff)))
      out[i+4, j] <- length(tdx)
    }
  }
}


cutoff <- c("pval<1", "pval<0.1", "pval<0.05", "pval<0.01",
            "padj<0.1", "padj<0.05", "padj<0.01", "padj<0.001")
out <- cbind(cutoff,out)
colnames(out) <- c("cutoff", "FC>1", "FC>1.5", "FC>2", "FC>3", "FC>4")
out <- as.data.frame(out, stringsAsFactors = F)
for(i in 2:6)
  storage.mode(out[,i]) <- 'integer'

normCounts <- cbind(inforMat, normCounts)
res <- cbind(inforMat, testStats[,1:11])
sig <- cbind(inforMat[idx,], testStats[idx,1:11])
####

message(paste0('Saving results!', '\n'))
write.table(normCounts, file = paste0(outPath, case, '_vs_', control, '_normCPM.txt'), append = FALSE, quote = FALSE, 
            sep = "\t", row.names = F)
write.table(res, file = paste0(outPath, case, '_vs_', control, '_DEstatAll.txt'), append = FALSE, quote = FALSE, 
            sep = "\t", row.names = F)
write.table(sig, file = paste0(outPath, case, '_vs_', control, '_DEstat.txt'), append = FALSE, quote = FALSE, 
            sep = "\t", row.names = F)
write.table(out, file = paste0(outPath, case, '_vs_', control, '_DEtab.txt'), append = FALSE, quote = FALSE, 
            sep = "\t", row.names = F)
write.table(sampleInfo, file = paste0(outPath, case, '_vs_', control, '_sampleInfo.txt'), append = FALSE, quote = FALSE, 
            sep = "\t", row.names = F)
fileConn<-file(paste0(outPath,case , '_vs_',control,'_parameters.txt'))
writeLines(paste0('filePath: ', filePath, '\n',
                  'sampleInfo: ', sampleIn, '\n',
                  'control: ', control, '\n',
                  'case: ', case, '\n',
                  'expression cutoff: ', exprCutoff, '\n',
                  'fold change cutoff: ', fc, '\n',
                  'p value cutoff:', pVal, '\n',
                  'p Flag: ', pFlag, '\n',
                  'paired: ', paired, '\n',
                  'output path:', outPath, '\n'), fileConn)
close(fileConn)
unlink(paste0(outPath, '*_files'), recursive = T)
t2 <- Sys.time()
cat('Total time spent: ', as.double(t2-t1), units(t2-t1) , '\n')