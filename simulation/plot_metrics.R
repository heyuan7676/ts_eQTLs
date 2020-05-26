suppressWarnings(library(plot.matrix))
suppressWarnings(library(ggplot2))
suppressWarnings(library(gplots))
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(cowplot))


metrics = read.table('simulation/output/metrics/metrics.txt', sep=' ', header=T, stringsAsFactors = F)
metrics = metrics[metrics$tau != 1, ]
metrics$method = factor(metrics$method, 
                           levels = c("sn_spMF", "flashr_nn", "flashr_backfit","flashr_default",
                                      "softImpute", "PCA", "SSVD", 
                                      "PMA_cv1", "PMA_cv2","NMF", "NBSFA"))
metrics$tau = factor(metrics$tau)
metrics[is.na(metrics$u_precision), "u_precision"] = 0
metrics[is.na(metrics$u_recall), "u_recall"] = 0
metrics[is.na(metrics$ts_precision), "ts_precision"] = 0
metrics[is.na(metrics$ts_recall), "ts_recall"] = 0


evaluation_theme =  theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x = element_text(angle = 90)) +
  theme_bw() 


mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(12)

f_corr_plot = ggplot(data = metrics, aes(tau, f_corr, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Correlation with true factor values") + 
  xlab("") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme

l_corr_plot = ggplot(data = metrics, aes(tau, l_corr, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Correlation with true loading values") + 
  xlab("") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme

u_precision_plot = ggplot(data = metrics, aes(tau, u_precision, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Precision of u-eQTLs") + 
  xlab("") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme

u_recall_plot = ggplot(data = metrics, aes(tau, u_recall, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Recall of u-eQTLs") + 
  xlab("") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme



ts_precision_plot = ggplot(data = metrics, aes(tau, ts_precision, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Precision of ts-eQTLs") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme

ts_recall_plot = ggplot(data = metrics, aes(tau, ts_recall, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Recall of ts-eQTLs") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme




g = ggdraw() + 
  draw_plot(f_corr_plot + theme(legend.position = 'none'), x = 0, y = 0.66, width = 0.4, height = 0.33)  + 
  draw_plot(l_corr_plot + theme(legend.position = 'none'), x = 0.45, y = 0.66, width = 0.4, height = 0.33) +
  draw_plot(u_recall_plot + theme(legend.position = 'none'), x = 0, y = 0.33, width = 0.4, height = 0.33) +
  draw_plot(u_precision_plot, x = 0.45, y = 0.33, width = 0.485, height = 0.33)+
  draw_plot(ts_recall_plot + theme(legend.position = 'none'), x = 0, y = 0, width = 0.4, height = 0.33)+
  draw_plot(ts_precision_plot + theme(legend.position = 'none'), x = 0.45, y = 0, width = 0.4, height = 0.33)

save_plot("simulation/output/metrics/metrics.png", g, base_width = 16, base_height = 10)





