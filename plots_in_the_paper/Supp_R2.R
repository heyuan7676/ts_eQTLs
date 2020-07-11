setwd('/Users/Yuan/Desktop/plots/')
source('scripts_refined/utils.R')

### gtex colors
tissues = read.table('data/tissues.txt', sep='\t', header=F, stringsAsFactors = F)
tissues = gsub("-","",tissues$V1)

gtex_col = read.table('data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
rownames(gtex_col) = gtex_col$tissue_id
gtex_col = gtex_col[tissues,]
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)



### Fig 1.

plot_factor_matrix <- function(fm, saveFile, factor_names_to_use = NA){
  fm = fm / matrix(rep(apply(fm, 2, function(x) max(abs(x))), dim(fm)[1]), nrow = dim(fm)[1], byrow = T)
  factor_matrix = data.frame(fm)
  
  rownames(factor_matrix) = tissues
  colnames(factor_matrix) = paste0("Factor", seq(1, dim(factor_matrix)[2]))
  
  factor_matrix$Tissue <- factor(rownames(factor_matrix), levels = gtex_col$tissue_id)
  factor_bar_df = melt(factor_matrix, id.vars = 'Tissue')
  
  if(sum(is.na(factor_names_to_use)) == 0){
    factor_bar_df$variable = as.character(factor_bar_df$variable)
    for(k in seq(1,23)){
      factor_bar_df[factor_bar_df$variable == paste0('Factor', k), "variable"] = paste0('Factor', k, '\n(', factor_names_to_use[k], ')')
    }
    factor_bar_df$variable = factor(factor_bar_df$variable, levels = unique(factor_bar_df$variable))
    
  }
  
  fig1_g = ggplot(data = factor_bar_df, aes(x = Tissue, y =value, fill = Tissue)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(values = gtex_col$tissue_color_hex)+
    facet_wrap(~variable, ncol = 4) + 
    xlab("") + 
    ylab("")+
    #theme(axis.text.x = element_text(angle = 90, hjust=0.97,vjust=0.2,size=7,family="Helvetica"))  +
    theme(axis.title = element_text(family="Helvetica",size=20)) +
    theme(plot.title = element_text(size = rel(1)), panel.background = element_blank()) +
    theme(axis.ticks.x =element_blank()) + 
    theme(axis.text.x = element_blank()) +
    theme(panel.grid = element_line(colour = 'grey', size = 0.2)) +
    #guides(col = FALSE) + 
    #theme(panel.background = element_rect(colour = "black")) +
    #theme(legend.justification = c(1, 1), legend.position = c(1, 0.9)) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2)) + 
    theme(legend.position = "none")
  
  
  pdf(paste0('results/',saveFile,'.pdf'),width=8,height=12)
  print(fig1_g) 
  dev.off()
  
}

factor_matrix = read.table('data/Fig1_factor_matrix.txt', sep='\t', stringsAsFactors = F, row.names = 1, header=T)
fm = factor_matrix
fm = fm / matrix(rep(apply(fm, 2, function(x) max(abs(x))), dim(fm)[1]), nrow = dim(fm)[1], byrow = T)
apply(fm, 2, function(x) tissues[which(abs(x) > 1e-5)])

g = ggplot(data = df, aes(x = ts_tissues, y =Freq)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(rotation = 90)) + 
  ggtitle('Number of ts-factors each tissue is in')

vif = 1 / (1 - summary(lm(as.matrix(factor_matrix[,1]) ~ 0 + as.matrix(factor_matrix[,seq(2,23)])))$r.squared)
plot_factor_matrix(factor_matrix, 'sn_spMF', factor_names)


load('data/Factor_matrices.RData')
pma_f = x[['pma_f']]
pma_f_cv2 = x[['pma_cv2_f']]
plot_factor_matrix(pma_f, 'pma_f')
plot_factor_matrix(pma_f_cv2, 'pma_f_cv2')

softimpute_f = read.table('data/softImpute_fm.txt', sep='\t', header=T, stringsAsFactors = F)
fm = softimpute_f
fm = fm / matrix(rep(apply(fm, 2, function(x) max(abs(x))), dim(fm)[1]), nrow = dim(fm)[1], byrow = T)
apply(fm, 2, function(x) "Adipose_Visceral_Omentum" %in% tissues[which(abs(x) > 0.01)])
plot_factor_matrix(softimpute_f, 'softimpute_f')


flashr_NN = read.table('~/Desktop/flashr/results/flashr_NN.txt', sep='\t', header=T, stringsAsFactors = F)
vif = 1 / (1 - summary(lm(as.matrix(flashr_NN[,1]) ~ 0 + as.matrix(flashr_NN[,seq(2,23)])))$r.squared)
summary(lm(as.matrix(flashr_NN[,1]) ~ 0 + as.matrix(flashr_NN[,c(3,5,6,8,10,12,14,17)])))$r.squared
plot_factor_matrix(flashr_NN, 'flashr_NN')

flashr = read.table('~/Desktop/flashr/results/flashr_default.txt', sep='\t', header=T, stringsAsFactors = F)
plot_factor_matrix(flashr, 'flashr')

flashr = read.table('~/Desktop/flashr/results/flashr_gd_bf.txt', sep='\t', header=T, stringsAsFactors = F)
fm = flashr
fm = fm / matrix(rep(apply(fm, 2, function(x) max(abs(x))), dim(fm)[1]), nrow = dim(fm)[1], byrow = T)
apply(fm, 2, function(x) tissues[which(abs(x) > 0.5)])
plot_factor_matrix(flashr, 'flashr_gd_bf')



### distribution of R2

methods = c()
r2 = c()
quantiles = c()
for(method in c("sn_spMF", "flashr_bf", "flashr_default","flashr_NN", 
                "softImpute", "PMD", "PMD_cv2")){
  m1_r2 = read.table(paste0('data/R2/',gsub('_default', '', method),'_Loadings_beta_BH_alpha0.05_corrected_R2.txt'), sep='\t')
  print(c(method, mean(m1_r2$V1 > 0), mean(m1_r2$V1 > 0.2),mean(m1_r2$V1 > 0.6)))
  
  m1_r2_quantiles = m1_r2 %>% mutate(quantile = ntile(V1, 10))
  methods = c(methods, rep(method,dim(m1_r2_quantiles)[1]))
  r2 = c(r2, m1_r2_quantiles$V1)
  quantiles = c(quantiles, m1_r2_quantiles$quantile)
}

m_r2 = data.frame(methods)
m_r2$r2 = r2
m_r2$quantile = quantiles

m_r2$methods = factor(m_r2$methods, 
                      levels = c("sn_spMF", 'flashr_bf',"flashr_default", "flashr_NN",
                                 "softImpute", "PMD", "PMD_cv2"))


m_r2$quantile = paste0("Quantile: ", (m_r2$quantile - 1)*10, '% - ', m_r2$quantile * 10, "%")

g  = ggplot(data = m_r2, aes(x=methods, y=r2, fill = methods)) + 
  geom_violin() + 
  facet_wrap(~quantile) + 
  scale_color_brewer(palette="Set2") + 
  xlab("")+ 
  ylab("R2")+ 
  ggtitle('Distribution of R2 across eQTLs') + 
  theme_bw() + 
  theme(axis.text.x = element_blank())

bn = bquote(''*R^2*'')
main = bquote('Distribution of '*R^2*' across eQTLs')

g  = ggplot(data = m_r2, aes(x=r2, color = methods)) + 
  geom_density() + 
  scale_color_brewer(palette="Set2") + 
  xlab("")+ 
  ylab('Density')+ 
  ggtitle(main) + 
  theme_bw()


save_plot("Supplmentary_plots//R2_distribution.png", g, base_width = 7.4, base_height=5)

