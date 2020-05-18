suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(gridExtra))
suppressWarnings(library(cowplot))
suppressWarnings(library(lemon))


### gtex colors for the tissues
gtex_col = read_tsv('../data/gtex_colors.txt')

### 

plot_factor_matrix <- function(Factor_fn){
  factor_matrix = read.table(Factor_fn, sep='\t', stringsAsFactors = F, header = F, row.names = 1)
  rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
  nF = ncol(factor_matrix)
  
  med_max = median(apply(factor_matrix, 2, max))
  factor_matrix = factor_matrix / matrix(rep((apply(factor_matrix, 2, max) / med_max), nrow(factor_matrix)), 
                                         nrow = nrow(factor_matrix), byrow = T)

  ## assign names to each factor
  tissues_in_factors_0 = apply(factor_matrix, 2, function(x) tissues[x>1e-5])
  tissues_in_factors = rep('0', length(tissues_in_factors_0))
  tissues_in_factors[sapply(tissues_in_factors_0, function(x) length(x)) > 40]  = 'Universal'
  for(k in seq(2, length(tissues_in_factors_0))){
    tissues_in_factors[k] = paste(unlist(tissues_in_factors_0[k]), collapse =';')
  }

  rownames(factor_matrix) = gtex_col$tissue
  colnames(factor_matrix) = paste0("Factor", seq(1, dim(factor_matrix)[2]))

  factor_matrix$Tissue = factor(rownames(factor_matrix), levels = gtex_col$tissue)
  factor_bar_df = melt(factor_matrix, id.vars = 'Tissue')

  factor_bar_df$variable_plot = as.character(factor_bar_df$variable)
  for(k in seq(1,dim(factor_matrix)[2])){
    factor_bar_df[factor_bar_df$variable_plot == paste0('Factor', k), "variable_plot"] = paste0('Factor', k, '\n(', tissues_in_factors[k], ')')
  }
  factor_bar_df$variable_plot = factor(factor_bar_df$variable_plot, levels = unique(factor_bar_df$variable_plot))

  fig_factors = ggplot(data = factor_bar_df, aes(x = Tissue, y =value, fill = Tissue)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(values = gtex_col$color_hex)+
    facet_rep_wrap(~variable_plot, ncol = 4) +  #repeat.tick.labels = 'left'
    xlab("Tissues") + 
    ylab("Values in the factors")+
    theme(strip.text.x = element_text(size = 6)) + 
    theme(axis.ticks.x =element_blank()) + 
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y =element_text(size = 10)) +
    theme(axis.title.x = element_text(size = 15)) + 
    theme(axis.title.y = element_text(size = 15)) + 
    theme(legend.position = "none")  + 
    background_grid(major = 'xy', minor = 'none')  + 
    ylim(0,7.5)
  
  save_fn = gsub('.txt', '_factors.pdf', Factor_fn)
  pdf(save_fn, height = 1.8 * ceiling((ncol(factor_matrix)-1) / 4))
  print(fig_factors)
  dev.off()
}
