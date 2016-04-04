RunDEseq <- function(data, conds_all, t1, t2, FC, pVal, pFlag) {
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
  library(DESeq)
  library(calibrate)
  library(gplots)
  dx <- c(which(conds_all == t1), which(conds_all == t2))
  countsTable <- data[, dx]
  conds <- conds_all[dx]
  cds <- newCountDataSet(countsTable, conds)
  cds <- estimateSizeFactors(cds)
  # sizeFactors(cds)
  if (length(dx) < 3) {
    cds <- estimateDispersions(cds, method="blind", sharingMode = "fit-only", fitType = "local")
  }
  else if (length(dx) < 6) {
    cds <- estimateDispersions(cds, method="pooled", sharingMode = "fit-only", fitType = "local")
  }
  else {
    cds <- estimateDispersions(cds)
  }
  
  ## plot heatmap of all genes
  if (length(dx) >= 3) {
    vsd <- getVarianceStabilizedData(cds)
    idx <- which(apply(log2(countsTable + 1), 1, sd) > 0.5)
    tmp <- vsd[idx, ]
    hr <- hclust(as.dist(1 - cor(t(tmp), method = "pearson")), method = "average")
    pdf_filename <- paste(t1, 'vs', t2, 'heatmap_all.pdf', sep="_")
    pdf(pdf_filename)
    mycol <- colorpanel(n=99, low="green", mid="black", high="red")
    heatmap.2(tmp, Rowv = as.dendrogram(hr), Colv = FALSE, col = mycol,
              scale = "row", key = TRUE, dendrogram = "row", keysize = 1, symkey = FALSE,
              density.info = "none", trace = "none", mar = c(8,5))
    dev.off()
  }
  
  res <- nbinomTest(cds, t1, t2)
  rownames(res) <- rownames(countsTable)
  if (length(dx) >= 3) {
    vsd <- vsd[order(res$pval), ]
  }
  res <- res[order(res$pval), ]
  
  if (pFlag == 1) {
    idx <- which(abs(res$log2FoldChange) >= log2(FC) & res$padj < pVal & 
                   pmax(res$baseMeanA,res$baseMeanB) >= 50)
  }
  else {
    idx <- which(abs(res$log2FoldChange) >= log2(FC) & res$pval < pVal & 
                   pmax(res$baseMeanA,res$baseMeanB) >= 50)
  }
  ## heatmap filtered genes
  if (length(dx) >= 3) {
    tmp <- vsd[idx, ]
    hr <- hclust(as.dist(1 - cor(t(tmp), method = "pearson")), method = "average")
    pdf_filename <- paste(t1, 'vs', t2, 'heatmap_filtered.pdf', sep = "_")
    pdf(pdf_filename)
    mycol <- colorpanel(n = 99, low = "green", mid = "black", high = "red")
    heatmap.2(tmp, Rowv = as.dendrogram(hr), Colv = FALSE, col = mycol,
              scale = "row", key = TRUE, dendrogram = "row", keysize = 1, symkey = FALSE,
              density.info = "none", trace = "none", mar = c(8,5))
    dev.off()
  }
  ## MA plot
  pdf_filename <- paste(t1, 'vs', t2, 'MAplot.pdf', sep = "_")
  pdf(pdf_filename)
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
  pdf(pdf_filename)
  with(res, plot(log2FoldChange, -log10(pval), 
                 xlim = c(min(log2FC, na.rm = T) - 3, max(log2FC, na.rm = T) + 3), 
                 pch = 20, cex = 0.4, main = "Volcano plot"))
  with(subset(res[idx, ],  log2FoldChange > log2(FC)), points(log2FoldChange, -log10(pval), 
                                                              pch = 20, cex = 0.4, col = "red"))
  with(subset(res[idx, ],  log2FoldChange < -1*log2(FC)), points(log2FoldChange, -log10(pval), 
                                                                 pch = 20, cex = 0.4, col = "green"))
  with(res[idx[1:10], ], textxy(log2FoldChange, -log10(pval), labs = id, cex = .6))
  dev.off()
  
  ## scatter plot
  pdf_filename <- paste(t1, 'vs', t2, 'Scatter.pdf', sep = "_")
  pdf(pdf_filename)
  with(res, plot(log2(baseMeanB), log2(baseMeanA),  
                 pch = 20, cex = 0.4, main = "Scatter plot", 
                 xlab = paste("log2", t2), ylab = paste("log2", t2)))
  with(subset(res[idx, ],  log2FoldChange > log2(FC)), points(log2(baseMeanB), log2(baseMeanA), 
                                                              pch = 20, cex = 0.4, col = "red"))
  with(subset(res[idx, ],  log2FoldChange < -1*log2(FC)), points(log2(baseMeanB), log2(baseMeanA), 
                                                                 pch = 20, cex = 0.4, col = "green"))
  with(res[idx[1:10], ], textxy(log2(baseMeanB), log2(baseMeanA), labs = id, cex = .6))
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
      tdx <- which(abs(res$log2FoldChange) >= fcLevel[j] & res$pval < pLevel[i] & 
                     pmax(res$baseMeanA, res$baseMeanB) >= 50)
      out[i, j] <- length(tdx)
    }
  }
  for (i in 1:4) {
    for (j in 1:5) {
      tdx <- which(abs(res$log2FoldChange) >= fcLevel[j] & res$padj < fdrLevel[i] & 
                     pmax(res$baseMeanA, res$baseMeanB) >= 50)
      out[i+4, j] <- length(tdx)
    }
  }
  cutoff <- c("pval<1", "pval<0.1", "pval<0.05", "pval<0.01",
              "padj<0.1", "padj<0.05", "padj<0.01", "padj<0.001")
  out <- cbind(cutoff,out)
  colnames(out) <- c("cutoff", "FC>1", "FC>1.5", "FC>2", "FC>3", "FC>4")
  
  ## export DE result
  colnames(res)[3] <- paste("baseMean1_", t1, sep = "")
  colnames(res)[4] <- paste("baseMean2_", t2, sep = "")
  colnames(res)[5] <- paste("FC", "2over1", sep = "_")
  colnames(res)[6] <- paste("log2FC", "2over1", sep = "_")
  sig <- res[idx, ]
  sig <- sig[is.na(sig$pval) != "TRUE", ]
  
  type <- paste(t1, 'vs', t2, sep="_")
  fileOut <- paste(type, "_result_DEseq_all.txt", sep = "")
  fileOut2 <- paste(type, "_result_DEseq_filtered.txt", sep = "")
  fileOut3 <- paste(type, "_cutoffNum.txt", sep="")
  write.table(res, file = fileOut, append = FALSE, quote = FALSE, 
              sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"))
  write.table(sig, file = fileOut2, append = FALSE, quote = FALSE, 
              sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"))
  write.table(out, file = fileOut3, append = FALSE, quote = FALSE, 
              sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"))
  return(sig)  
}