library(plot.matrix)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(cowplot)


metrics = read.table('simulation/output/metrics/metrics.txt', sep=' ', header=T, stringsAsFactors = F)
metrics$method = factor(metrics$method, 
                           levels = c("sn_spMF", "flashr_nn", "flashr_backfit","flashr_default",
                                      "softImpute", "PCA", "SSVD", 
                                      "PMA_cv1", "PMA_cv2","NMF"))
metrics$tau = factor(metrics$tau)
metrics[is.na(metrics$precision), "precision"] = 0


evaluation_theme =  theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x = element_text(angle = 90)) +
  theme_bw() 


mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(12)

f_corr_plot = ggplot(data = metrics, aes(tau, f_corr, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Correlation with true factor values") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme

l_corr_plot = ggplot(data = metrics, aes(tau, l_corr, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Correlation with true loading values") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme

precision_plot = ggplot(data = metrics, aes(tau, precision, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Precision") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme

recall_plot = ggplot(data = metrics, aes(tau, recall, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Recall") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme



g = ggdraw() + 
  draw_plot(f_corr_plot + theme(legend.position = 'none'), x = 0, y = 0.5, width = 0.4, height = 0.5)  + 
  draw_plot(l_corr_plot, x = 0.45, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(recall_plot + theme(legend.position = 'none'), x = 0, y = 0, width = 0.4, height = 0.5) +
  draw_plot(precision_plot + theme(legend.position = 'none'), x = 0.45, y = 0, width = 0.4, height = 0.5)

save_plot("simulation/output/metrics/metrics.png", g, base_width = 12, base_height = 6)





