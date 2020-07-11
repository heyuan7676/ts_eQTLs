setwd('/Users/Yuan/Desktop/plots/')
source('scripts_refined/utils.R')

flashr_signs = read.table('data/flashr_signs_df.txt', sep='\t', header = T)
flashr_signs$Factor = factor(flashr_signs$Factor, levels = sort(unique(flashr_signs$Factor)))
g1 = ggplot(data = flashr_signs, aes(x = Factor, y = value, fill = variable)) + 
    geom_bar(stat="identity", width=.5, position = "dodge") + 
    ylab('Number of eQTLs') + 
  theme_bw() + 
  ggtitle('flashr_default')




flashr_signs = read.table('data/flashr_signs_bf_df.txt', sep='\t', header = T)
flashr_signs$Factor = factor(flashr_signs$Factor, levels = sort(unique(flashr_signs$Factor)))
g2 = ggplot(data = flashr_signs, aes(x = Factor, y = value, fill = variable)) + 
  geom_bar(stat="identity", width=.5, position = "dodge") + 
  ylab('Number of eQTLs') + 
  theme_bw() + 
  ggtitle('flashr_bf')



fig_signs = ggdraw() + 
  draw_plot(g1, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(g2, x = 0, y = 0, width = 1, height = 0.5)
  
save_plot("results_refined/Fig_flashr_flashrbf_signs.png", fig_signs, base_width = 6, base_height = 5)
