suppressWarnings(library('plyr'))
suppressWarnings(library('dplyr'))
suppressWarnings(library('RColorBrewer'))
suppressWarnings(library('pheatmap'))
suppressWarnings(library('cowplot'))
suppressWarnings(library('reshape2'))
suppressWarnings(library('R.matlab'))
suppressWarnings(library('ggplot2'))
source('simulation/compare_methods.R')


wrap_multiple_runs <- function(tau){
	simulation_N = 1

	evaulation_i = data.frame()
	for(k in seq(1,simulation_N)){
  		print(paste0('Seed', k))
  		evaulation_i = rbind(evaulation_i, compare_methods(tau = tau,seed = k))
	}
	evaulation_i$tau = tau
	return(evaulation_i)
}


evaluation = list()
iii = 1
for(tau in c(1000, 500, 100, 50, 20, 10, 5)){
	evaluation[[iii]] = wrap_multiple_runs(tau = tau)
	iii = iii + 1
}


evaluation = do.call("rbind", evaluation)
evaluation$method = factor(evaluation$method,
                                   levels=c("True", "sn_spMF","NMF",
                                            "flashr_nn",'flashr_backfit','flashr_default',
                                            'NBSFA', "softImpute", "PCA",
                                            "SSVD", "PMA_cv1", "PMA_cv2"))
evaluation$tau = factor(evaluation$tau)
write.table(evaluation,'simulation/output/evaluation.txt', sep='\t')


evaluation_theme =  theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x = element_text(angle = 90)) +
  theme_bw() 



mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(12)

g1 = ggplot(data = evaluation, aes(tau, abs(V1), fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Correlation with true loading values") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme+ guides(fill=guide_legend(ncol=2))

g2 = ggplot(data = evaluation, aes(tau, abs(V2), fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Correlation with true factor values") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_manual(values = mycolors) + 
  evaluation_theme + guides(fill=guide_legend(ncol=2))

g3 = ggplot(data = evaluation, aes(tau, V3, fill = method)) + 
  geom_boxplot(lwd=0.2, outlier.size = 0.2) + 
  ggtitle("Relative root mean squared error") + 
  xlab("Precision of the error values") + 
  ylab("") + 
  scale_fill_brewer(palette="Set1") + 
  evaluation_theme


g = ggdraw() + 
  #draw_plot(g3, x = 0, y = 0.5, width = 0.665, height = 0.5)  + 
  draw_plot(g1, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(g2 + theme(legend.position = 'none'), x = 0, y = 0, width = 0.742, height = 0.5)
save_plot("Evaluation.png", g, base_width = 9, base_height = 6)

