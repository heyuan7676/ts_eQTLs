setwd('/Users/Yuan/Desktop/plots/')
source('scripts_refined/utils.R')

### Fig 1.

plot_bars <- function(fm, method){
  fm = fm / matrix(rep(apply(fm, 2, function(x) max(abs(x))), dim(fm)[1]), nrow = dim(fm)[1], byrow = T)
  nT_each_factor = data.frame(apply(fm, 2, function(x) sum((abs(x) > 0.01))))
  colnames(nT_each_factor) = 'V1'
  nT_each_factor$V2 = factor(as.character(seq(1, dim(nT_each_factor)[1])), levels = seq(1, dim(nT_each_factor)[1]))
  
  g1 = ggplot(data = nT_each_factor, aes(x = V2, y =V1)) + geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle(method) + 
    xlab('Factors') + 
    ylab('')
  
  
  nF_each_tissue = data.frame(apply(fm, 1, function(x) sum((abs(x) > 0.01))))
  colnames(nF_each_tissue) = 'V1'
  nF_each_tissue$V2 = factor(tissues, levels = tissues)
  
  g2 = ggplot(data = nF_each_tissue, aes(x = V2, y =V1)) + geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle(method) + 
    xlab('') + 
    ylab('')
  
  return(list(g1, g2))
}



factor_matrix = read.table('data/Fig1_factor_matrix.txt', sep='\t', stringsAsFactors = F, row.names = 1, header=T)
flashr_NN = read.table('~/Desktop/flashr/results/flashr_NN.txt', sep='\t', header=T, stringsAsFactors = F)
flashr = read.table('~/Desktop/flashr/results/flashr_default.txt', sep='\t', header=T, stringsAsFactors = F)
flashr_bf = read.table('~/Desktop/flashr/results/flashr_gd_bf.txt', sep='\t', header=T, stringsAsFactors = F)

sn_spMF = plot_bars(factor_matrix, 'sn_spMF')
flashr = plot_bars(flashr, 'flashr_default')
flashr_NN = plot_bars(flashr_NN, 'flashr_NN')
flashr_bf = plot_bars(flashr_bf, 'flashr_bf')


fig = ggdraw() + 
  draw_plot(sn_spMF[[1]], x = 0, y = 0.7, width = 0.5, height = 0.3) + 
  draw_plot(flashr[[1]], x = 0, y = 0.35, width = 0.5, height = 0.3) + 
  draw_plot(flashr_bf[[1]], x = 0.5, y = 0.7, width = 0.5, height = 0.3) + 
  draw_plot(flashr_NN[[1]], x = 0.5, y = 0.35, width = 0.5, height = 0.3)
save_plot("Supplmentary_plots/N_T_F_1.png", fig, base_width = 7, base_height = 10)



fig = ggdraw() + 
  draw_plot(sn_spMF[[2]] + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4)), x = 0, y = 0.5, width = 0.5, height = 0.5) + 
  draw_plot(flashr[[2]] + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4)), x = 0, y = 0, width = 0.5, height = 0.5) + 
  draw_plot(flashr_bf[[2]] + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4)), x = 0.5, y = 0.5, width = 0.5, height = 0.5) + 
  draw_plot(flashr_NN[[2]] + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4)), x = 0.5, y = 0, width = 0.5, height = 0.5)
save_plot("Supplmentary_plots/N_T_F_2.png", fig, base_width = 16, base_height = 16)





