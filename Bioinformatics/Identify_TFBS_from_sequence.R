### Firstly, install Biostrings package from Bioconductor
### remember you only need to do this once
source('http://www.bioconductor.org/biocLite.R')
biocLite('Biostrings')

#### load package into your working space
library(Biostrings)

#### Position frequency matrix of TEAD1; You can retrive
#### this information form JASPAR database
PWM_TEAD1 <- matrix(as.integer(c( 1, 9, 0, 12,  0,  0,  0,  0,  5,  1,  2,  0,
                       6,  0, 12,  0,  0,  0, 12, 11,  0,  7,  4,  2,
                       1,  3,  0,  0,  0,  0,  0,  0,  0,  4,  3,  8,
                       4,  0,  0,  0, 12, 12,  0,  1,  7,  0,  3,  2)), nrow = 4, ncol = 12, byrow = T)
rownames(PWM_TEAD1) <- c('A', 'C', 'G', 'T')

#### Calculate Position weight matrix that will be used for pattern matching
PWM_TEAD1 <- PWM(PWM_TEAD1)
# PWM_TEAD1 <- PWM_TEAD1 + 1 
# PWM_TEAD1 <- PWM_TEAD1/apply(PWM_TEAD1, 2, sum)
# PWM_TEAD1 <- log(PWM_TEAD1/0.25)

#### Load sequence information for CCL2 or any other genes
CCL2_fasta_human <- readDNAStringSet('Desktop/TFBScode/CCL2_human.fasta', use.names = T)
CCL2_fasta_mouse <- readDNAStringSet('Desktop/TFBScode/CCL2.fasta', use.names = T)

#### Apply matchPWM to position weight matrix to identify matched pattern
hits_human <- matchPWM(PWM_TEAD1, CCL2_fasta_human[[1]], with.score = TRUE, min.score = '85%')
hits_mouse <- matchPWM(PWM_TEAD1, CCL2_fasta_mouse[[1]], with.score = TRUE, min.score = '85%')
