#### load library
library(data.table)
library(scde)
library(ggplot2)

#### set working directory
setwd('../Dropbox/CoProject/Liu Yang/')

#### load dataset
dataMat <- data.frame(fread('DigitalGeneExpression_matrix_single_cell_pools.txt', data.table = F, stringsAsFactors = F), row.names = 1)
design <- read.table('design_copy.txt', sep = '\t', header = T, stringsAsFactors = F)

#### remove cells with low overall read coverage and genes with zero expression
dataMat_sel <- clean.counts(dataMat, min.lib.size=1000, min.reads = 1, min.detected = 1)
index_rm <- which(!colnames(dataMat) %in% colnames(dataMat_sel))
design_sel <- design$Design[-index_rm]


#### building error models
t1 <- Sys.time()
err_mod <- scde.error.models(counts = dataMat_sel, groups = as.factor(design_sel), n.cores = 2,
                             threshold.segmentation=T, save.crossfit.plots=F, 
                             save.model.plots=F,verbose=1, min.size.entries = 100)
t2 <- Sys.time()
t2 - t1

##### exclude cells showing abnormal fits
valid_cells <- err_mod$corr.a >0
err_mod <- err_mod[valid_cells, ]
dataMat_sel <- dataMat_sel[, valid_cells]

##### retrive normalized read counts
norm_counts <- as.matrix(2^(scde.expression.magnitude(models = err_mod, counts = dataMat_sel)))

exprs_prior <- scde.expression.prior(models = err_mod,
                                     counts = dataMat_sel,
                                     length.out=400,
                                     show.plot=F)
names(design_sel) <- row.names(err_mod)

###### test for differential expression
ediff <- scde.expression.difference(models = err_mod, counts = dataMat_sel, 
                                    prior =  exprs_prior, groups = as.factor(design_sel),
                                    n.randomizations=100, n.cores = 2,verbose=1, 
                                    return.posteriors = F)


TestStat <- data.frame(
  AveExpr = apply(norm_counts, 1, mean),
  logfc = ediff$mle,
  pval = 2*pnorm(-abs(ediff$Z))
)
TestStat$padj <- p.adjust(TestStat$pval, method = 'BH')
DE_index <- which(with(TestStat, AveExpr > quantile(AveExpr)[3] &
                         logfc > 1 &
                         padj > 0.01))


write.csv(TestStat,  'TestStat.csv', quote = F)
write.csv(TestStat[DE_index,],  'DEstat.csv', quote = F)

##### pca analysis
dataMat <- log2(norm_counts + 0.01)
cv.gene <- apply(dataMat, 1, function(x)
  sd(x) / mean(x))
dataMat <- dataMat[which(cv.gene > 0.1),]
dataMat <- scale(dataMat)
pca.result <- prcomp(t(dataMat))
ppoints <- pca.result$x[,1:3]
percent<-round((pca.result$sdev^2/sum(pca.result$sdev^2))*100,1)
ppoints <-
  cbind(as.data.frame(ppoints), design_sel, paste('S', 1:ncol(dataMat), sep = ''))
colnames(ppoints) <- c('PC1', 'PC2', 'PC3', 'Design', 'Samplename')

gp <- ggplot() + geom_point(aes(x = PC1, y = PC2, color = Design), data = as.data.frame(ppoints))
ggsave('pcaPlot2d.pdf', gp)





