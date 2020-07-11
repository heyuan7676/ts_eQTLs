setwd('/Users/Yuan/Desktop/plots/')
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gplots)
library(RColorBrewer) 
library(distances)
library(scales)


lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

tfs = read.table('data/enriched_TFs_tpm/interesting_TFs.txt', sep='\t', header=F, stringsAsFactors = F)
tfs = tfs$V1
tfs = c("ENSG00000101076.16")

tfs_tpm_across_tis = list()
for(tfi in tfs){
  tfs_tpm_across_tis[[tfi]] = list()
}

for(tis in gsub("-", "", tissues)){
  fn=paste0('data/enriched_TFs_tpm/interesting_TFs_in_',tis,'.txt')
  tf_tis = read.table(fn, sep='\t', stringsAsFactors = F, row.names = 1)
  for(tfi in tfs){
    tfs_tpm_across_tis[[tfi]] = lappend(tfs_tpm_across_tis[[tfi]], as.numeric(tf_tis[tfi,seq(2, dim(tf_tis)[2]) ]))
  }
}


library(gridExtra)
library(ggplot2)

p <- list()
for(i in seq(1, 1)){
  tissue_vec = c()
  tpm_vec = c()
  
  tfi = tfs[i]
  for(tis in seq(1,length(tissues))){
    tpm_vec = c(tpm_vec, as.numeric(tfs_tpm_across_tis[[tfi]][[tis]]))
    tissue_vec = c(tissue_vec, rep(tissues[tis], length(tfs_tpm_across_tis[[tfi]][[tis]])))
  }
  
  df_tp = data.frame("TPM" = tpm_vec, "tissue" = tissue_vec)
  df_tp$logTPM = log10(df_tp$TPM + 1)
  
  g = qplot(tissue, logTPM, data=df_tp, geom=c("boxplot"), fill = I(gtex_col$tissue_color_hex), xlab = "", outlier.size=0.01) +
    ggtitle(tfi) + 
    theme(plot.title = element_text(size=10), axis.title = element_text(size=8), axis.text = element_text(size=8))+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2)) 
  
  p[[i]] <- g
}



png('results_refined/Fig5_HNF4A_TPM_log_revision.png', width = 12, height = 8, units = 'in', res = 200)
print(p[[i]])
dev.off()



x = aggregate(df_tp, list(df_tp$tissue), mean)
as.character(x[x$TPM > 0.1, 1])

