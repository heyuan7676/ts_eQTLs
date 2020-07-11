setwd('/Users/Yuan/Desktop/plots/')
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(lemon)


### gtex colors
gtex_col = read.table('data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id),]
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)
rownames(gtex_col) = gtex_col$tissue_id
head(gtex_col)

### factor matrix
fn = 'SparseMF_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a149_l1490'
load(paste0('data/R_code/', fn, '.RData'))
factor_matrix = as.data.frame(FactorM)
med_max = median(apply(factor_matrix, 2, max))
factor_matrix = factor_matrix / matrix(rep((apply(factor_matrix, 2, max) / med_max), 49), nrow = 49, byrow = T)


tissues_in_factors_0 = apply(factor_matrix, 2, function(x) tissues[x>0])
tissues_in_factors = rep('0', length(tissues_in_factors_0))
tissues_in_factors[1]  = 'Universal'
for(k in seq(2, length(tissues_in_factors_0))){
  tissues_in_factors[k] = paste(unlist(tissues_in_factors_0[k]), collapse =';')
}
tissues_in_factors[17] = 'Brain tissues'

factor_order = c(1, 17, 20, 3, 5, 10, 12, 4,
                 22, 13, 23, 14, 16, 6, 2, 19, 
                 18, 8, 21,11, 7, 9, 24, 15)
factor_matrix = factor_matrix[, factor_order]
tissues_in_factors = tissues_in_factors[factor_order]

rownames(factor_matrix) = gtex_col$tissue_id
colnames(factor_matrix) = paste0("Factor", seq(1, dim(factor_matrix)[2]))

factor_matrix$Tissue <- factor(rownames(factor_matrix), levels = gtex_col$tissue_id)
factor_bar_df = melt(factor_matrix, id.vars = 'Tissue')

factor_bar_df$variable_plot = as.character(factor_bar_df$variable)
for(k in seq(1,dim(factor_matrix)[2])){
  factor_bar_df[factor_bar_df$variable_plot == paste0('Factor', k), "variable_plot"] = paste0('Factor', k, '\n(', tissues_in_factors[k], ')')
}
factor_bar_df$variable_plot = factor(factor_bar_df$variable_plot, levels = unique(factor_bar_df$variable_plot))



fig1_g = ggplot(data = factor_bar_df, aes(x = Tissue, y =value, fill = Tissue)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = gtex_col$tissue_color_hex)+
  facet_rep_wrap(~variable_plot, ncol = 4) +  #repeat.tick.labels = 'left'
  xlab("Tissues") + 
  ylab("Values in the factors")+
  theme(strip.text.x = element_text(size = 9)) + 
  theme(axis.ticks.x =element_blank()) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y =element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15)) + 
  theme(legend.position = "none")  + 
  background_grid(major = 'xy', minor = 'none')  + 
  ylim(0,7.5)

png('Supplmentary_plots/Fig1_1.png', width = 8, height = 12, units = 'in', res = 200)
print(fig1_g) 
dev.off()



### R code
### factor matrix
fn = 'SparseMF_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a149_l1490'
load(paste0('data/R_code/', fn, '.RData'))
factor_matrix = as.data.frame(FactorM)
med_max = median(apply(factor_matrix, 2, max))
factor_matrix = factor_matrix / matrix(rep((apply(factor_matrix, 2, max) / med_max), 49), nrow = 49, byrow = T)


### matlab code
### factor matrix
factor_matrix = read.table('data/Fig1_factor_matrix.txt', sep='\t', stringsAsFactors = F, row.names = 1)
factor_matrix = factor_matrix[rownames(factor_matrix) != "", ]
rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
colnames(factor_matrix) = paste0("Factor", seq(1, dim(factor_matrix)[2]))





