########## dotplot and barplot of FDR and Z score for top upstream regulators as predicted by IPA
#### input: input .csv file from IPA upregulator analysis
#### output: path to output directory
#### nup: number of activated regulators
#### ndown: number of inhibited regulators
twoplots <- function(input, output, nup, ndown){
  library(ggplot2)
  dataMat <- read.csv(input,check.names = F,header = T,stringsAsFactors = F)
  dataMat <- dataMat[order(dataMat$`p-value of overlap`,decreasing = F),]
  dataMat$targets <- sapply(dataMat$`Target molecules in dataset`,function(i) length(strsplit(i,',')[[1]]))
  dataMat$targets_cate <- sapply(1:nrow(dataMat), function(i){
    if(dataMat$targets[i] <= 40) '<=40' else 
      if(dataMat$targets[i]>40 &dataMat$targets[i]<=80) '41-80' else 
        if(dataMat$targets[i]>80 & dataMat$targets[i]<=120) '81-120' else '>120'
  })
  dataMat$targets_cate <- factor(dataMat$targets_cate,levels = c('<=40','41-80','81-120','>120'))
  dataMat <- dataMat[c(which(dataMat$`Activation z-score`>0)[1:nup],which(dataMat$`Activation z-score`<0)[1:ndown]) ,]
  dataMat <- dataMat[order(log10(dataMat$`p-value of overlap`),decreasing = T),]
  gp <- ggplot() + geom_point(aes(x=1:(nup+ndown),y=-log10(dataMat$`p-value of overlap`),size=dataMat$targets_cate),color='red') + 
    theme_linedraw() + scale_x_continuous(breaks = 1:(nup+ndown),labels = dataMat$`Upstream Regulator`) +
    labs(x='Upstream Regulator',y='-Log10(FDR)') + guides(size=guide_legend(title = 'Number of target genes')) +
    theme(axis.text.x = element_text(angle = 90,size=10),legend.title  = element_text(size=10,face='bold'),axis.title = element_text(size=15,face = 'bold')) + coord_flip()
  ggsave(paste0(output,'/dotplot.png'),gp,dpi = 600,width = 12, height = 8)
  dataMat <- dataMat[order(dataMat$`Activation z-score`,decreasing = F),]
  dataMat$activate <- ifelse(dataMat$`Activation z-score`>0,'Activated','Inhibited')
  gp1 <- ggplot() + geom_col(aes(x=1:(nup+ndown),y=dataMat$`Activation z-score`,fill=dataMat$activate)) + 
    theme_linedraw() + scale_x_continuous(breaks = 1:(nup+ndown),labels = dataMat$`Upstream Regulator`) +
    labs(x='Upstream Regulator',y='Activation Z-score') + guides(fill=guide_legend(title = 'Predicted State')) +
    theme(axis.text.x = element_text(angle = 90,size=10),legend.title  = element_text(size=10,face='bold'),axis.title = element_text(size=15,face = 'bold')) + coord_flip()
  ggsave(paste0(output,'/barplot.png'),gp1,dpi=600,width = 12, height = 8)
  }