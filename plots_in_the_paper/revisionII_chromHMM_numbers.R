
setwd('/Users/Yuan/Desktop/plots/')

plot_ratio <- function(status, suffix, title_text, xlab_text, color_pad = 'Blues'){
  data = read.table(paste0('data/revisionII//Number_of_tissues_',status,'_proportion_melt',suffix,'.txt'), 
                    sep = '\t', header=T, stringsAsFactors = F)
  data$variable = factor(data$variable, levels = paste0('Group', seq(0,22)))
  
  
  data_normalized_by_group0 = NULL
  for(k in seq(0,22)){
    x = data[data$Number_tissues == k, ]
    x$value = x$value / x$value[1]
    data_normalized_by_group0 = rbind(data_normalized_by_group0, x)
  }
  data_normalized_by_group0 = data_normalized_by_group0[data_normalized_by_group0$variable != 'Group0', ]
  
  mycolors <- colorRampPalette(brewer.pal(8, color_pad))(23)
  data_normalized_by_group0 = data_normalized_by_group0[data_normalized_by_group0$variable != 'Group0', ]
  data_normalized_by_group0 = data_normalized_by_group0[data_normalized_by_group0$Number_tissues > 0, ]
  data_normalized_by_group0$Number_tissues = factor(data_normalized_by_group0$Number_tissues)
  tp = sort(unique(data_normalized_by_group0$Number_tissues))
  g = ggplot(data = data_normalized_by_group0, aes(x = Number_tissues, y = value, fill = Number_tissues)) +
    scale_fill_manual(values = mycolors) + 
    geom_boxplot() + 
    theme_bw()  + 
    ggtitle(title_text) + 
    xlab(xlab_text) + 
    ylab('Fraction of ts-eQTLs variants/\n  fraction of u-eQTLs variants') + 
    geom_hline(yintercept=1, linetype="dashed", color = "grey") + 
    scale_x_discrete(breaks = tp[seq(1, length(tp), by = 2)])
  
  return(g)
}



g_TssA = plot_ratio(status = '1_TssA', suffix = '', title_text = '', xlab_text = 'Number of tissues with promoter overlapping eQTL variants')
g_TssA_rd = plot_ratio(status = '1_TssA', suffix = '_rd', title_text = '', xlab_text= 'Number of tissues with promoter overlapping random variants')

g_Enh = plot_ratio(status = '7_Enh', suffix = '', title_text = '', xlab_text = 'Number of tissues with enhancer overlapping eQTL variants', color_pad = 'Reds')
g_Enh_rd = plot_ratio(status = '7_Enh', suffix = '_rd', title_text = '', xlab_text= 'Number of tissues with enhancer overlapping random variants', color_pad = 'Reds')

g = ggdraw() + 
  draw_plot(g_TssA + theme(legend.position = 'none'), x = 0, y = 0.5, width = 0.5, height = 0.5)  + 
  draw_plot(g_TssA_rd + theme(legend.position = 'none', axis.title.y = element_blank()), x = 0.5, y = 0.5, width = 0.45, height = 0.5) +
  draw_plot(g_Enh + theme(legend.position = 'none'), x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(g_Enh_rd + theme(legend.position = 'none', axis.title.y = element_blank()), x = 0.5, y = 0, width = 0.45, height = 0.5)

save_plot("results_refined/Number_pro_enhs.png", g, base_width = 10, base_height = 7)

