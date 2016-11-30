library(edgeR)
library(calibrate)
library('gplots')

################# Run edgeR with raw count
RawCount <- read.table('TCGA_LIHC_paired_N-T.txt', header = TRUE, row.names = 1 , sep = '\t')
y <- DGEList(counts=RawCount)
y <- calcNormFactors(y)
Patient_design <- paste('patient', rep(1:50, each =2), sep = '')
Tumor_design <- rep(c('N', 'T'), length = 100)
design <- model.matrix(~Patient_design+Tumor_design)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)

TestStats_edgeR <- as.data.frame(topTags(lrt, n = nrow(RawCount), sort.by = 'none', p.value = 1))

DE_index <- with(TestStats_edgeR, which(abs(logFC) > 1 & FDR < 0.05 & logCPM > 0))
NormCount_edgeR <- cpm(y)

summary(de <- decideTestsDGE(lrt))

# MAplot
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

## Volcano plot
idx <- which(abs(TestStats_edgeR$logFC) >= 1 & TestStats_edgeR$PValue < 0.05 & 
               TestStats_edgeR$logCPM >= 0)
log2FC <- TestStats_edgeR$logFC
is.na(log2FC) <- sapply(log2FC, is.infinite)
pdf_filename <- paste('Tumor','vs','Normal','volcano TCGA 6.10.16.pdf',sep = "_")
pdf(pdf_filename)
with(TestStats_edgeR, plot(logFC, -log10(PValue), 
                           xlim = c(min(logFC, na.rm = T) - 3, max(logFC, na.rm = T) + 3), 
                           pch = 20, cex = 0.4, main = "Volcano plot"))
with(subset(TestStats_edgeR[idx, ],  logFC > 1), points(logFC, -log10(PValue), 
                                                        pch = 20, cex = 0.4, col = "red"))
with(subset(TestStats_edgeR[idx, ],  logFC < -1), points(logFC, -log10(PValue), 
                                                         pch = 20, cex = 0.4, col = "green"))
with(TestStats_edgeR[idx[1:10], ], textxy(logFC, -log10(PValue), labs = rownames(TestStats_edgeR[idx[1:10],]), cex = .6))
dev.off()


##Heatmap2
hr <- hclust(as.dist(1 - cor(t(NormCount_edgeR[DE_index,c(seq(1,49,2), seq(2,50,2))]), method = "pearson")), method = "average")
mycol <- colorpanel(n = 99, low = "green", mid = "black", high = "red")
heatmap.2(NormCount_edgeR[DE_index,c(seq(1,99,2), seq(2,100,2))], Rowv = as.dendrogram(hr), Colv = FALSE, col = mycol,
          scale = "row", key = TRUE, dendrogram = "row", keysize = 1, symkey = FALSE,
          density.info = "none", trace = "none", mar = c(8,5))

##Export Results
write.table(TestStats_edgeR[DE_index,], 'DE_EdgeR TCGA 6.10.16.txt', quote = F, sep = '\t', row.names = T, col.names = T)
write.table(TestStats_edgeR, 'TestStat_EdgeR TCGA 6.10.16.txt', quote = F, sep = '\t', row.names = T, col.names = T)