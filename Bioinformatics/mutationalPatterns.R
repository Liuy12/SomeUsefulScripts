args <- commandArgs(TRUE)
vcf_files <- args[1]
sample_names <- args[2]
genome_build <- args[3]
nsig <- as.numeric(args[4])
output <- args[5]

vcf_files <- read.delim(vcf_files, header = F, stringsAsFactors = F)
sample_names <- read.delim(sample_names, header = F, stringsAsFactors = F)
library(MutationalPatterns)
library(ggplot2)
## load reference genome
if(genome_build == 'hg19'){
  library('BSgenome.Hsapiens.UCSC.hg19')
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
} else {
  library('BSgenome.Hsapiens.UCSC.hg38')
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
}
message(paste0('Loading vcf files as granges object', '/n'))
vcfs <- read_vcfs_as_granges(vcf_files$V1, sample_names$V1, ref_genome)
## count mutation type occurrence
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
## overall mutation spectrum
message(paste0('saving overall mutation spectrum figure', '/n'))
gp <- plot_spectrum(type_occurrences, CT = TRUE)
ggsave(paste0(output, '/overall_mutation_spectrum.pdf'), gp, useDingbats=F)
ggsave(paste0(output, '/overall_mutation_spectrum.png'), gp)
## per-sample mutation spectrum
message(paste0('saving per-sample mutation spectrum figure', '/n'))
n <- min(20,length(vcfs))
gp <- plot_spectrum(type_occurrences[1:n,], by = 1:n, CT = TRUE, legend = F)
ggsave(paste0(output, '/persample_mutation_spectrum.pdf'), gp, useDingbats=F)
ggsave(paste0(output, '/persample_mutation_spectrum.png'), gp)
## 96 mutational profiles
message(paste0('saving mutation signature figure', '/n'))
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
gp <-  plot_96_profile(mut_mat[,1:n], condensed = TRUE)
ggsave(paste0(output,'/mutation_signature.pdf'), gp, useDingbats=F)
ggsave(paste0(output,'/mutation_signature.png'), gp)

## de novo mutation signature extracting using NMF
mut_mat <- mut_mat + 0.0001
# library("NMF")
# estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
# pdf('nmf_eval.png');plot(estimate);dev.off()
## extract de novo signatures
message(paste0('saving de novo signatures', '\n'))
nmf_res <- extract_signatures(mut_mat, rank = nsig, nrun = 50)
colnames(nmf_res$signatures) <- paste0('Signature ', 1:nsig)
rownames(nmf_res$contribution) <- paste0('Signature ', 1:nsig)
## plot extracted signatures
gp <- plot_96_profile(nmf_res$signatures, condensed = TRUE)
ggsave(paste0(output,'/denovo_mutation_signatures.pdf'),gp,useDingbats=F)
ggsave(paste0(output,'/denovo_mutation_signatuers.png'),gp)
## plot contribution heatmap
gp <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=TRUE)
ggsave(paste0(output,'/relative_contribution_heatmap.pdf'),gp, useDingbats=F)
ggsave(paste0(output,'/relative_contribution_heatmap.png'),gp)

## cosmic signatures
message(paste0('cosmic signatures', '\n'))
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/",
                "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
# Plot heatmap with specified signature order
gp <- plot_cosine_heatmap(cos_sim_samples_signatures,
                          col_order = cosmic_order,
                          cluster_rows = TRUE)
ggsave(paste0(output, '/cos_sim_cosmic_signature.pdf'),gp, useDingbats=F)
ggsave(paste0(output, '/cos_sim_cosmic_signature.png'),gp)

cos_sim_samples_signatures = cos_sim_matrix(nmf_res$signatures, cancer_signatures)
# Plot heatmap with specified signature order
gp <- plot_cosine_heatmap(cos_sim_samples_signatures,
                          col_order = cosmic_order,
                          cluster_rows = TRUE)
ggsave(paste0(output, '/cos_sim_cosmic_signature_denovo_sig.pdf'),gp, useDingbats=F)
ggsave(paste0(output, '/cos_sim_cosmic_signature_denovo_sig.png'),gp)





