RunDEseq_v2 <- function(data, conds_all, t1, t2, FC, pVal, pFlag) {
  # Differential gene expression analysis by using DESeq package. 
  # sig <- runDEseq(data, conds_all, t1, t2, 2, 0.05, 0)
  # Args:
  #  data: input read counts in matrix format, column is sampples, and row is genes.
  #  conds_all: The experimental group of each sample
  #  t1: The control group name
  #  t2: The test group name
  #  FC: The fold change cutoff
  #  pVal: The raw pVal or adjusted pVal cutoff
  #  pFlag: 1: use adjusted pVal, 0: use raw pVal
  # 
  # Returns:
  #  sig: The signature gene list using criterions above
  #  
  # Output file generated:
  #  t1_vs_t2_result_DEseq_all.txt: DESeq result of all genes
  #  t1_vs_t2_result_DEseq_filtered.txt: DESeq result of DE genes
  #  t1_vs_t2_volcano.pdf: volcano plot, red: up-regulated genes, green: down-regulated genes
  #  t1_vs_t2_heatmap_all.pdf: heatmap of genes (in log2 intensity) that standard deviation > 0.5 
  #  t1_vs_t2_heatmap_filtered.pdf: heatmap of DE genes
  #  t1_vs_t2_MAplot.pdf: maplot of all genes, in red are the DE genes
  #  t1_vs_t2_cutoffNum.txt: number of DE genes using different cutoff criterions
  
  # Hung-I Harry Chen, April 2015, GCCRI/UTHSCSA
  # 
  library(DESeq2)
  library(calibrate)
  library(gplots)
  rownames(data) <- data[, 1]
  data <- data[, -1]
  if(tail(rownames(data), 1) == "alignment_not_unique") {
    data <- head(data, -5)
  }
  rpkm <- read.delim('Norm_logRPKM_allSamples_Genes_ReadCount.txt', header=TRUE, stringsAsFactors=TRUE)
  rownames(rpkm) <- rpkm[, 1]
  rpkm <- rpkm[, -1]
  ifelse(length(which(conds_all == t1)) == 1, 
         rpkm_avg1 <- rpkm[, which(conds_all == t1)],
         rpkm_avg1 <- apply(rpkm[, which(conds_all == t1)], 1, mean))    
  ifelse(length(which(conds_all == t2)) == 1, 
         rpkm_avg2 <- rpkm[, which(conds_all == t2)],
         rpkm_avg2 <- apply(rpkm[, which(conds_all == t2)], 1, mean))    
    

  dx <- c(which(conds_all == t1), which(conds_all == t2))
  countsTable <- data[, dx]
  conds <- conds_all[dx]
  colData <- data.frame(condition = factor(conds))
  cds_in <- DESeqDataSetFromMatrix(countsTable, colData, formula(~ condition))
  cds <- DESeq(cds_in)
  res <- results(cds)
  res$id <- rownames(res)
  res$baseMeanA <- rowMeans(counts(cds,normalized=TRUE)[,cds$condition == t1], 
                            na.rm = TRUE, dims = 1)
  res$baseMeanB <- rowMeans(counts(cds,normalized=TRUE)[,cds$condition == t2], 
                            na.rm = TRUE, dims = 1)
  res <- res[c(7,1,8,9,2:6)]
  ## export normalized data
  normData <- counts(cds, normalized = TRUE)
  ProbeIDs <- rownames(normData)
  colnames(normData) <- colnames(countsTable)
  normData <- cbind(ProbeIDs, normData)
  write.table(normData, file = "allSamples_normData.txt", append = FALSE, quote = FALSE, 
              sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"))
  
  ## create sub folder
  type <- paste(t1, 'vs', t2, 'DEseq2', sep="_")
  dir.create(type, showWarning = FALSE)
  
  ## plot heatmap of all genes
  if (length(dx) >= 3) {
    vsd <- getVarianceStabilizedData(cds)
    colnames(vsd) <- colnames(countsTable) 
    idx <- which(apply(log2(countsTable + 1), 1, sd) > 0.5)
    tmp <- vsd[idx, ]
    hr <- hclust(as.dist(1 - cor(t(tmp), method = "pearson")), method = "average")
    pdf_filename <- paste(t1, 'vs', t2, length(idx), 'heatmap_all.pdf', sep="_")
    pdf(paste(type, pdf_filename, sep = "/"))
    mycol <- colorpanel(n=99, low="green", mid="black", high="red")
    heatmap.2(tmp, Rowv = as.dendrogram(hr), Colv = FALSE, col = mycol,
              scale = "row", key = TRUE, dendrogram = "row", keysize = 1, symkey = FALSE,
              density.info = "none", trace = "none", mar = c(8,5))
    dev.off()
    ##Sample clustering
    pdf_filename <- paste(t1, 'vs', t2, 'sampleCluster.pdf', sep="_")
    pdf(paste(type, pdf_filename, sep = "/"))
    dists <- dist(t(vsd))
    heatmap(as.matrix(dists), symm = TRUE, cexCol = 0.9, cexRow = 0.9, mar = c(10,5),
             main="Sample Cluster on log2 normalized counts")
    dev.off()
  }
  ## Sample pair scatter plot
  if (length(dx) >= 3) {
    panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor) 
    {
      usr <- par("usr"); on.exit(par(usr)) 
      par(usr = c(0, 1, 0, 1)) 
      r <- abs(cor(x, y)) 
      txt <- format(c(r, 0.123456789), digits = digits)[1] 
      txt <- paste(prefix, txt, sep="") 
      if(missing(cex.cor)) cex <- 0.8 / strwidth(txt) 
      test <- cor.test(x ,y) 
      Signif <- symnum(test$estimate, corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.95, 0.96, 0.97, 0.98, 0.99, 1),
                       symbols = c(" ","*", "**", "***", "****", "*****"))
      text(0.5, 0.5, txt, cex = cex * r) 
      text(.8, .8, Signif, cex = cex * 0.5, col = 2) 
    }
    panel.scatter <-function (x, y, col = "black", bg = NA, pch = 20, 
                              cex = 0.4, col.smooth = "red", span = 2/3, iter = 3, ...) {
      points(x, y, pch = pch, col = col, bg = bg, cex = cex)
      ok <- is.finite(x) & is.finite(y)
      if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
    }
    png_filename <- paste(t1, 'vs', t2, 'ScatterPairs.png', sep = "_")
    png(filename = paste(type, png_filename, sep = "/"), width = 800, height = 800)
    pairs(vsd, lower.panel = panel.scatter, upper.panel=panel.cor,
          main="Normalized gene read counts pairwise scatter")
    dev.off()
  }
  ## run DESeq DE test

  if (pFlag == 1) {
    idx <- which(abs(res$log2FoldChange) >= log2(FC) & res$padj < pVal & 
                 pmax(rpkm_avg1,rpkm_avg2) >= 0)
#     idx_all <- which(abs(all$log2FoldChange) >= log2(FC) & all$padj < pVal & 
#                    pmax(rpkm_avg1,rpkm_avg2) >= 1)
  }
  else {
    idx <- which(abs(res$log2FoldChange) >= log2(FC) & res$pvalue < pVal & 
                   pmax(rpkm_avg1,rpkm_avg2) >= 0)
#     idx_all <- which(abs(all$log2FoldChange) >= log2(FC) & all$pval < pVal & 
#                        pmax(rpkm_avg1,rpkm_avg2) >= 1)
  }
  sig <- res[idx, ]
  sig <- sig[is.na(sig$pval) != "TRUE", ]
  ## heatmap filtered genes
  if (length(dx) >= 3) {
    tmp <- vsd[idx, ]
    hr <- hclust(as.dist(1 - cor(t(tmp), method = "pearson")), method = "average")
    pdf_filename <- paste(t1, 'vs', t2, 'heatmap_filtered.pdf', sep = "_")
    pdf(paste(type, pdf_filename, sep = "/"))
    mycol <- colorpanel(n = 99, low = "green", mid = "black", high = "red")
    heatmap.2(tmp, Rowv = as.dendrogram(hr), col = mycol,
              scale = "row", key = TRUE, dendrogram = "both", keysize = 1, symkey = FALSE,
              density.info = "none", trace = "none", mar = c(8,5))
    dev.off()
  }
  ## MA plot
  pdf_filename <- paste(t1, 'vs', t2, 'MAplot.pdf', sep = "_")
  pdf(paste(type, pdf_filename, sep = "/"))
  if (pFlag == 1) {
    plot(res$baseMean, res$log2FoldChange, log = "x", pch = 20, cex = .4, 
         col = ifelse((res$padj < pVal) & (abs(res$log2FoldChange) > log2(FC)), "red", "black"))
  }
  else {
    plot(res$baseMean, res$log2FoldChange, log = "x", pch = 20, cex = .4, 
         col = ifelse((res$pval < pVal) & (abs(res$log2FoldChange) > log2(FC)), "red", "black"))  
  }
  dev.off()
  ## Volcano plot

  log2FC <- res$log2FoldChange
  is.na(log2FC) <- sapply(log2FC, is.infinite)
  pdf_filename <- paste(t1,'vs',t2,'volcano.pdf',sep = "_")
  pdf(paste(type, pdf_filename, sep = "/"))
  with(res, plot(log2FoldChange, -log10(pvalue), 
                 xlim = c(min(log2FC, na.rm = T) - 3, max(log2FC, na.rm = T) + 3), 
                 pch = 20, cex = 0.4, main = "Volcano plot"))
  with(subset(res[idx, ],  log2FoldChange > log2(FC)), points(log2FoldChange, -log10(pvalue), 
       pch = 20, cex = 0.4, col = "red"))
  with(subset(res[idx, ],  log2FoldChange < -1*log2(FC)), points(log2FoldChange, -log10(pvalue), 
       pch = 20, cex = 0.4, col = "green"))
  with(sig[order(sig$pvalue)[1:15], ], textxy(log2FoldChange, -log10(pvalue), 
       labs = id, cex = .6))
  dev.off()
  
  ## scatter plot
  pdf_filename <- paste(t1, 'vs', t2, 'Scatter.pdf', sep = "_")
  pdf(paste(type, pdf_filename, sep = "/"))
  with(res, plot(log2(baseMeanB), log2(baseMeanA),  
                 pch = 20, cex = 0.4, main = "Scatter plot", 
                 xlab = paste("log2", t1), ylab = paste("log2", t2)))
  with(subset(res[idx, ],  log2FoldChange > log2(FC)), points(log2(baseMeanB), log2(baseMeanA), 
       pch = 20, cex = 0.4, col = "red"))
  with(subset(res[idx, ],  log2FoldChange < -1*log2(FC)), points(log2(baseMeanB), log2(baseMeanA), 
       pch = 20, cex = 0.4, col = "green"))
  with(sig[order(sig$pval)[1:15], ], textxy(log2(baseMeanB), 
       log2(baseMeanA), labs = id, cex = .6))
  abline(1, 1, col = "blue", lty = 2)
  abline(-1, 1, col = "blue", lty = 2)
  dev.off()
  
  ## cutoff table
  pLevel <- c(1, 0.1, 0.05, 0.01)
  fdrLevel <- c(0.1, 0.05, 0.01, 0.001)
  fcLevel <- c(0, log2(1.5), log2(2), log2(3), log2(4))
  out <- matrix(data = NA, nrow = 8, ncol = 5)
  for (i in 1:4) {
    for (j in 1:5) {
      tdx <- which(abs(res$log2FoldChange) >= fcLevel[j] & res$pvalue < pLevel[i] & 
                     pmax(rpkm_avg1,rpkm_avg2) >= 0)
      out[i, j] <- length(tdx)
    }
  }
  for (i in 1:4) {
    for (j in 1:5) {
      tdx <- which(abs(res$log2FoldChange) >= fcLevel[j] & res$padj < fdrLevel[i] & 
                     pmax(rpkm_avg1,rpkm_avg2) >= 0)
      out[i+4, j] <- length(tdx)
    }
  }
  cutoff <- c("pval < 1", "pval < 0.1", "pval < 0.05", "pval < 0.01",
              "adj.pVal < 0.1", "adj.pVal < 0.05", "adj.pVal < 0.01", "adj.pVal < 0.001")
  out <- cbind(cutoff,out)
  colnames(out) <- c("cutoff", "FC > 1", "FC > 1.5", "FC > 2", "FC > 3", "FC > 4")
  
  ## export DE result
  colnames(res)[3] <- paste("baseMean1_", t1, sep = "")
  colnames(res)[4] <- paste("baseMean2_", t2, sep = "")
  colnames(res)[5] <- paste("log2FC", "2over1", sep = "_")
  all <- res
  all <- cbind(all, counts(cds, normalized = TRUE))
  sig <- res[idx, ]
  sig <- sig[is.na(sig$pvalue) != "TRUE", ]
  all_sig <- all[idx, ]
  fileOut <- paste(type, "_result_DEseq_all.txt", sep = "")
  fileOut2 <- paste(type, "_result_DEseq_filtered.txt", sep = "")
  fileOut3 <- paste(type, "_cutoffNum.txt", sep="")
  fileOut4 <- paste(type, "_result_DEseq_all_withNormData.txt", sep = "")
  fileOut5 <- paste(type, "_result_DEseq_filtered_withNormData.txt", sep = "")
write.table(res[order(res$pvalue), ], file = paste(type, fileOut, sep = "/"), append = FALSE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))
write.table(sig[order(sig$pvalue), ], file = paste(type, fileOut2, sep = "/"), append = FALSE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))
write.table(out, file = paste(type, fileOut3, sep = "/"), append = FALSE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))
write.table(all, file = paste(type, fileOut4, sep = "/"), append = FALSE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))
write.table(all_sig, file = paste(type, fileOut5, sep = "/"), append = FALSE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))
  

  
  return(sig)  
}