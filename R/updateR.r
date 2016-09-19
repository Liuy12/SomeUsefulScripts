### retrive package information from old version of R
oldPackages <- row.names(installed.packages())
if(!'rsconnect' %in% oldPackages)
  install.packages('rsconnect')
library(rsconnect)
dir.create('~/tmpDir/')
setwd('~/tmpDir/')
write.table(paste('library(', oldPackages, ')', sep = ''), './tmpfile.r', sep = '\n', quote = F, row.names = F, col.names = F)
packinfo <- appDependencies(appDir = './')
write.table(packinfo, './packinfo.txt', quote = F, row.names = F, sep = '\t')

### after install the new version of R
### restart rstudio/r session to load the new version of R
### then install packages for new version of R
setwd('~/tmpDir/')
packinfo <- read.table('./packinfo.txt', sep = '\t', header = T)
unlink('~/tmpDir/')
diffPackages <- packinfo[!packinfo$package %in% rownames(installed.packages()),]

### install cran packages 
install.packages(diffPackages[diffPackages$source == 'CRAN',]$package)

### install bioconductor packages
source('http://www.bioconductor.org/biocLite.R')
biocLite(diffPackages[diffPackages$source == 'Bioconductor',]$package)

### install github packages 
### this might not be very accurate 
### but is the best shot I have so far
if(!'githubinstall' %in% diffPackages$package)
  install.packages('githubinstall')
library(githubinstall)
githubinstall(diffPackages[diffPackages$source == 'github',]$package, ask = F)

### for packages that can not be installed by abovementioned methods, 
### you might need to do it manually
devtools::install_github('nghiavtr/BPSC')
devtools::install_github('RGLab/MAST')
devtools::install_github('Liuy12/MBDDiff')
devtools::install_github('ramnathv/rCharts')
devtools::install_github('skardhamar/rga')
devtools::install_github('ramnathv/slidify')
devtools::install_github('ramnathv/slidifyLibraries')
install.packages('threejs')
biocLite('scde')
biocLite('ROTS')



