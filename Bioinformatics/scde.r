#### load library
library(data.table)
library(scde)

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
err_mod <- scde.error.models(counts = dataMat_sel, groups = as.factor(design_sel), n.cores = 2,
                             threshold.segmentation=T, save.crossfit.plots=F, 
                             save.model.plots=F,verbose=1, min.size.entries = 2000)
valid_cells <- err_mod$corr.a >0
err_mod <- err_mod[valid_cells, ]
counts <- counts[, valid_cells]
norm_counts <- as.matrix(2^(scde.expression.magnitude(models = err_mod, counts = counts)))

