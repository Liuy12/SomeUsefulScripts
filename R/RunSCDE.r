##################################################################################################
################ scde wrapper function
##### RunSCDE(exprs_path, design_path, fc, pvalue, pflag, output_path, cores)
#####
################ prerequsite
##### scde 
##### data.table
#####
################ input:
##### exprs_path: path of expression matrix, rows indicating genes, columns indicating samples
##### design_path: path of design matrix, a file containing single column indicating group design 
#####              information (number of rows equals to number of columns in exprs_path file). 
#####              For instance, 'C1', 'C2'
##### fc: fold change cutoff for DEGs, default 2
##### pvalue: p value cutoff for DEGs, default 0.01
##### pflag: whether to use adjusted p value or pvalue, 1 indicates adjusted pvalue, 0 indicats 
#####        raw p value. Default 1
##### output_path: path to save output files
##### cores: number of cores to use for scde, default 1
#####
################ output: 
##### TestStat.csv: csv file containing four columns. 1, average expression across all cells for each gene; 
#####               2, log2 fold change; 3, p value; 4 adjusted p value. 
##### DEstat.csv: same columns as TestStat.csv but only include DE genes that pass pre-specified criteria.
##### Normcounts.csv: csv file containing normalized read counts

RunSCDE <- function(exprs_path, design_path, fc = 2, pvalue = 0.01, pflag = 1, output_path, cores = 1){
  library('scde')
  library('data.table')
  dataMat <- data.frame(fread(exprs_path, data.table = F, stringsAsFactors = F), row.names = 1)
  design <- read.table(design_path, sep = '\t', stringsAsFactors = F)
  
  message('remove cells with low overall read coverage and genes with zero expression')
  dataMat_sel <- clean.counts(dataMat, min.lib.size=1000, min.reads = 1, min.detected = 1)
  index_rm <- which(!colnames(dataMat) %in% colnames(dataMat_sel))
  if(length(index_rm))
    design_sel <- design[[1]][-index_rm]
  else
    design_sel <- design[[1]]
  cat(length(index_rm), ' number of cells dropped during this step', '\n', sep = '')
  
  message('building error models')
  t1 <- Sys.time()
  err_mod <- scde.error.models(counts = dataMat_sel, groups = as.factor(design_sel), n.cores = 2,
                               threshold.segmentation=T, save.crossfit.plots=F, 
                               save.model.plots=F,verbose=0, min.size.entries = 100)
  t2 <- Sys.time()
  cat('Total time consumes during this step: ', t2-t1, '\n' ,sep = '')
  
  message('Exclude cells showing abnormal fits')
  valid_cells <- err_mod$corr.a >0
  err_mod <- err_mod[valid_cells, ]
  dataMat_sel <- dataMat_sel[, valid_cells]
  design_sel <- design[[1]][valid_cells]
  cat(sum(!valid_cells), ' number of cells dropped during this step', '\n', sep = '')
  
  ##### retrive normalized read counts
  norm_counts <- as.matrix(2^(scde.expression.magnitude(models = err_mod, counts = dataMat_sel)))
  
  message('Estimate prior distribution for gene expression manitudes')
  t1 <- Sys.time()
  exprs_prior <- scde.expression.prior(models = err_mod,
                                       counts = dataMat_sel,
                                       length.out=400,
                                       show.plot=F)
  names(design_sel) <- row.names(err_mod)
  t2 <- Sys.time()
  cat('Total time consumes during this step: ', t2-t1, '\n' ,sep = '')
  
  message('Test for differential expression')
  t1 <- Sys.time()
  ediff <- scde.expression.difference(models = err_mod, counts = dataMat_sel, 
                                      prior =  exprs_prior, groups = as.factor(design_sel),
                                      n.randomizations=100, n.cores = cores,verbose=1, 
                                      return.posteriors = F)
  t2 <- Sys.time()
  cat('Total time consumes during this step: ', t2-t1, '\n' ,sep = '')
  
  TestStat <- data.frame(
    AveExpr = apply(norm_counts, 1, mean),
    logfc = ediff$mle,
    pval = 2*pnorm(-abs(ediff$Z))
  )
  TestStat$padj <- p.adjust(TestStat$pval, method = 'BH')
  if(pflag == 1)
    DE_index <- which(with(TestStat, AveExpr > quantile(AveExpr)[3] &
                           abs(logfc) > log2(fc) &
                           padj < pvalue))
  else
    DE_index <- which(with(TestStat, AveExpr > quantile(AveExpr)[3] &
                             abs(logfc) > log2(fc) &
                             pval < pvalue))
  
  
  write.csv(TestStat,  paste(output_path, 'TestStat.csv', sep = '/'), quote = F)
  write.csv(TestStat[DE_index,],  paste(output_path, 'DEstat.csv', sep = '/'), quote = F)
  write.csv(norm_counts, paste(output_path, 'Normcounts.csv', sep = '/'), quote = F)
  message('Done!')
}

