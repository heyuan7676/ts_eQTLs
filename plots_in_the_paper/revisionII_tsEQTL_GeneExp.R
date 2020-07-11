setwd('/Users/Yuan/Desktop/plots/')

dat = read.table('data/revisionII/ts_eQTLs_ratio_geneExpLevel.txt', sep='\t', header=T,stringsAsFactors = F)
dat = dat[dat$has_ts_eQTL != 0, ]

g1 = ggplot(data = dat, aes(x = V5, y = has_ts_eQTL)) + 
  geom_boxplot() + 
  xlab('Ratio of gene expression in one tissue compared to other tissues')+
  ylab('Proportion of genes with ts-eQTLs') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90))


g = ggdraw() + 
  draw_plot(g1)
save_plot("results_refined/prop_genes_ts_eQTLs.png", g, base_width = 6, base_height = 4)



dat = read.table('data/revisionII/ts_eQTLs_TPM_thresholded.txt', sep='\t', header=T,stringsAsFactors = F)
dat$TPM_threshold = factor(dat$TPM_threshold, levels = unique(dat$TPM_threshold))

g1 = ggplot(data = dat, aes(x = TPM_threshold, y = ts_eQTL_proportion)) + 
  geom_boxplot() + 
  xlab('Lower threshold of genes \n with median TPM in the tissue') + 
  ylab('Proportion of ts-eQTLs\n among tested eQTLs')+
  theme_bw()

g2 = ggplot(data = dat, aes(x = TPM_threshold, y = ts_gene_proportion)) + 
  geom_boxplot() + 
  xlab('Lower threshold of genes \n with median TPM in the tissue') + 
  ylab('Proportion of genes with ts-eQTLs\n among tested genes')+
  theme_bw()


g = ggdraw() + 
  draw_plot(g1 , x = 0, y = 0, width = 0.45, height = 1) + 
  draw_plot(g2 , x = 0.5, y = 0, width = 0.45, height = 1) 
save_plot("results_refined/prop_ts_eQTLs_thre.png", g, base_width = 8, base_height = 3)


