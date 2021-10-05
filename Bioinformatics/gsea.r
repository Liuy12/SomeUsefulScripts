
#######################################################################################################
#########################plot summerizing figure including top pos/neg enriched gene sets
gseaSumPlot <- function(pos_path, neg_path, pos_sel, neg_sel, outpath){
  require(readxl)
  require(ggplot2) 

  gseaResUp <- read.delim(pos_path, stringsAsFactors = F)
  gseaResDown <- read.delim(neg_path, stringsAsFactors = F)
  
  idx_up <- sapply(pos_sel, function(i) {
    grep(paste0('^', i, '$'), gseaResUp$NAME)
  })
  print(idx_up)
  cat('retrived ', length(idx_up), ' number of up-regulated gene-sets', '\n')
  idx_down <- sapply(neg_sel, function(i) {
    grep(paste0('^', i, '$'), gseaResDown$NAME)
  })
  print(idx_down)
  cat('retrived ', length(idx_down), ' number of down-regulated gene-sets')
  if(length(idx_up) != length(pos_sel) | length(idx_down) != length(neg_sel))
    stop('Please check spelling of selected gene sets')
  gseaResUp <- gseaResUp[idx_up,]
  gseaResDown <- gseaResDown[idx_down,]
  #gseaResUp <- gseaResUp[order(gseaResUp$NES, decreasing = T),]
  #gseaResDown <- gseaResDown[order(gseaResDown$NES, decreasing = F),]
  tmp <- rbind(gseaResUp,gseaResDown)
  tmp$NAME <- sapply(tmp$NAME, function(i) {
    sp <- strsplit(i, '_')[[1]]
    sp <- sp[-1]
    sp <- tolower(paste(sp, collapse = ' '))
    sp <- paste(toupper(substr(sp, 1,1)), substr(sp,2,nchar(sp)), sep = '')
  })
  # Plotting
  p <- ggplot(tmp, aes(NES, NAME)) +
    geom_point(aes(colour=FDR.q.val, size=SIZE)) +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1)) +
    geom_vline(xintercept=0, size=0.5, colour="gray50") +
    theme(panel.background=element_rect(fill="gray95", colour="gray95"),
          panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
          axis.title.y=element_blank()) +
    #  expand_limits(x=c(-3,3)) +
    #  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
    scale_y_discrete(limits=rev(tmp$NAME))
  ggsave(paste0(outpath, '/gsea_sel', '.pdf'), p, width = 15, height = 10, useDingbats=FALSE)
}


#######################################################################################################
######################### Generate heatmap including one selected gene sets for pos and neg each
gseaHeatmap <- function(dataMatpath, testStatPath, msigdbpath, numgenes, upGSname,downGSname,colAnn, color=c('blue','white','red'), nAnno =6, rankCol=TRUE, outpath){
  require(GSA)
  require(NMF)
  dataMat <- data.table::fread(dataMatpath, data.table = F)
  testStats <- data.table::fread(testStatPath, data.table = F)
  ## load msigdb gene sets annotation
  msigdb <- GSA.read.gmt(msigdbpath)
  or <- order(testStats$Foldchange, decreasing = T)
  dataMat_sel <- dataMat[c(or[1:numgenes], tail(or, numgenes)),]
  idx <- grep(paste0('^', upGSname, '$'), msigdb$geneset.names)
  ann_gs1 <- sapply(dataMat_sel$GeneName, function(i) {
    idx1 <- grep(paste0('^', i, '$'), msigdb$genesets[[idx]])
    if(length(idx1)) 'hit' else 'off'
  })
  idx <- grep(paste0('^', downGSname, '$'), msigdb$geneset.names)
  ann_gs2 <- sapply(dataMat_sel$GeneName, function(i) {
    idx1 <- grep(paste0('^', i, '$'), msigdb$genesets[[idx]])
    if(length(idx1)) 'hit' else 'off'
  })
  
  rowAnn <- data.frame(upGSname=ann_gs1,
                       downGSname=ann_gs2)
  colAnn <- data.frame(Type=colAnn)
  annColors <- list(upGSname = c('black', 'white'),
                    downGSname = c('green','white'),
                    Type = c('yellow', 'yellow4','greenyellow', 'springgreen4', 'springgreen')[1:length(unique(colAnn$Type))])
  # examine whether there are extreme values in the dataset
  # if extreme values are detected, then modified the color ratio 
  # to make heatmap appear more contract
  Exprs_scale <- t(scale(t(dataMat_sel[,-c(1:nAnno)]+0.1)))
  Exprs_scale[is.na(Exprs_scale)] <- min(Exprs_scale, na.rm = TRUE)
  rg <- quantile(Exprs_scale,na.rm=T)
  rg_diff <- rg[4]-rg[2]
  rg_max <- max(abs(rg))
  Exprs_sd <- sd(Exprs_scale)
  Exprs_mean <- mean(Exprs_scale)
  breakratio = c(1,1,1)
  numbreaks=100
  col = color
  if(rg_max > max(abs(c(Exprs_mean + 3*Exprs_sd, Exprs_mean - 3*Exprs_sd)))){
    rg_iqr <- max(abs(c(rg[2], rg[4])))
    bp <- c((breakratio[1]/sum(breakratio))*rg_diff - rg_iqr, rg_iqr - (breakratio[3]/sum(breakratio))*rg_diff)
    bk <- unique(c(seq(-rg_max, -rg_iqr, length= numbreaks), seq(-rg_iqr,bp[1],length = numbreaks), seq(bp[1],bp[2],length=numbreaks),seq(bp[2],rg_iqr,length=numbreaks), 
                   seq(rg_iqr, rg_max, length = numbreaks)))
    hmcols<- colorRampPalette(col)(length(bk)-1)
  } else{
    rg <- range(Exprs_scale, na.rm=T)
    bp <- c((breakratio[1]/sum(breakratio))*diff(rg) - rg_max, rg_max - (breakratio[3]/sum(breakratio))*diff(rg))
    bk <- c(seq(-rg_max,bp[1],length=numbreaks), seq(bp[1],bp[2],length=numbreaks),seq(bp[2],rg_max,length=numbreaks))
    bk <- bk[!duplicated(bk)]
    hmcols<- colorRampPalette(col)(length(bk)-1)
  }
  
  tmp <- log2(dataMat_sel[,-c(1:nAnno)]+0.1)
  nmf.options(grid.patch=TRUE)
  #pdf(paste0(outpath, '/gseaheatmap_', numgenes, '_genes.pdf'), height =4, width = 9, useDingbats = F)
  pdf(paste0(outpath, '/gseaheatmap_', numgenes, '_genes.pdf'), useDingbats = F, height = 4, width = 5)
  if(rankCol){
    colrk <- order(sapply(1:(ncol(tmp)-nAnno), function(i) cor(tmp[,i], testStats$Foldchange[c(or[1:numgenes], tail(or, numgenes))])), decreasing = F)
    aheatmap(log2(dataMat_sel[,-c(1:nAnno)][,colrk]+0.1), annRow = rowAnn, annCol = colAnn ,scale = 'row', labCol = NA, annColors = annColors, Rowv = NA, Colv = NA, color = hmcols, breaks = bk)
  }
  else
    aheatmap(log2(dataMat_sel[,-c(1:nAnno)]+0.1), annRow = rowAnn, annCol = colAnn ,scale = 'row', labCol = NA, annColors = annColors, Rowv = NA, Colv = NA, color = hmcols, breaks = bk)
  dev.off()
}

## Script by Thomas Kuilman
## path argument: path to output folder of analysis (e.g. PATH/my_analysis.GseaPreranked.1470948568349)
## gene.set argument: name of the gene set (e.g. V$AP1_Q2).
## It is used in a grep command, so multiple matching is possible.
## Also, R regular expressions can be handled, e.g. "IL2[0-9]$"
## Leading "V$" from gene set names are stripped to allow using the grep command.
## In case of multiple grep matches a warning is given and the first option is plotted.
## class.name: the name of the class / variable to which genes have been correlated (e.g. drug-treatment)
## metric.range: the range of the metric; defaults to [min(DEFINED RANGE), max(DEFINED RANGE)]
## replotGSEA('my_analysis_2.GseaPreranked.1512668450791/', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', class.name = 'HALLMARK', outpath)


replotGSEA <- function(path, gene.set, class.name, metric.range, outpath=NULL) {
  
  if (missing(path)) {
    stop("Path argument is required")
  }
  # if (!file.exists(path)) {
  #   stop("The path folder could not be found. Please change the path")
  # }
  if (missing(gene.set)) {
    stop("Gene set argument is required")
  }
  
  ## Load .rnk data
  path.rnk <- list.files(path = file.path(path, "edb"),
                         pattern = ".rnk$", full.names = TRUE)
  gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
  colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
  if (missing(metric.range)) {
    metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
  }  
  
  ## Load .edb data
  path.edb <- list.files(path = file.path(path, "edb"),
                         pattern = ".edb$", full.names = TRUE)
  gsea.edb <- read.delim(file = path.edb,
                         header = FALSE, stringsAsFactors = FALSE)
  gsea.edb <- unlist(gsea.edb)
  gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
  gsea.metric <- unlist(strsplit(gsea.metric, " "))
  gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
  gsea.metric <- gsub("METRIC=", "", gsea.metric)
  gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]
  
  # Select the right gene set
  if (length(gsea.edb) == 0) {
    stop(paste("The gene set name was not found, please provide",
               "a correct name"))
  }
  if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)) > 1) {
    warning(paste("More than 1 gene set matched the gene.set",
                  "argument; the first match is plotted"))
  }
  gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)[1]]
  
  # Get template name
  gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
  gsea.edb <- unlist(strsplit(gsea.edb, " "))
  gsea.template <- gsea.edb[1]
  
  # Get gene set name
  gsea.gene.set <- gsea.edb[2]
  gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
  
  # Get enrichment score
  gsea.enrichment.score <- gsea.edb[3]
  gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
  
  # Get gene set name
  gsea.normalized.enrichment.score <- gsea.edb[4]
  gsea.normalized.enrichment.score <- gsub("NES=", "",
                                           gsea.normalized.enrichment.score)
  
  # Get nominal p-value
  gsea.p.value <- gsea.edb[5]
  gsea.p.value <- gsub("NP=", "", gsea.p.value)
  gsea.p.value <- as.numeric(gsea.p.value)
  
  # Get FDR
  gsea.fdr <- gsea.edb[6]
  gsea.fdr <- gsub("FDR=", "", gsea.fdr)
  gsea.fdr <- as.numeric(gsea.fdr)
  
  # Get hit indices
  gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
  gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
  gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
  gsea.hit.indices <- as.integer(gsea.hit.indices)
  
  # Get ES profile
  gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
  gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
  gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
  gsea.es.profile <- as.numeric(gsea.es.profile)
  
  if(!is.null(outpath))
    pdf(paste0(outpath, '/gsea_', gene.set, '.pdf'),width = 6.9, height = 6.9, useDingbats = F)
  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
  #layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
       c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       main = list(gsea.gene.set, font = 1, cex = 1),
       panel.first = {
         abline(h = seq(round(min(gsea.es.profile), digits = 1),
                        max(gsea.es.profile), 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score), adj = c(0, 0))
  } else {
    text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
  }
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk$metric)))
  abline(v = gsea.hit.indices, lwd = 0.75)
  
  par(mar = c(0, 5, 0, 2))
  rank.colors <- gsea.rnk$metric - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(gsea.rnk$metric)))
  box()
  text(length(gsea.rnk$metric) / 2, 0.7,
       labels = ifelse(!missing(class.name), class.name, gsea.template))
  text(length(gsea.rnk$metric) * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(length(gsea.rnk$metric) * 0.99, 0.7, "Negative", adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(round(gsea.rnk$metric, digits = 2))
  plot(gsea.rnk$metric, type = "n", xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
       ylim = metric.range, yaxs = "i",
       ylab = if(gsea.metric == "None") {"Ranking metric"} else {gsea.metric},
       panel.first = abline(h = seq(metric.range[1] / 2,
                                    metric.range[2] - metric.range[1] / 4,
                                    metric.range[2] / 2), col = "gray95", lty = 2))
  
  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
          xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
          ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
          ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric), space = 0, add = TRUE)
  box()
  if(!is.null(outpath))
    dev.off()
  # Reset to default
  par(def.par)
}


