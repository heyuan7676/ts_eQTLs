suppressWarnings(library(cowplot))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(gridExtra))
suppressWarnings(library(colortools))
suppressWarnings(library(colorspace))
suppressWarnings(library(readr))
suppressWarnings(library(dplyr))
suppressWarnings(library(lemon))

### gtex colors for the tissues
gtex_col = read_tsv('data/gtex_colors.txt')
gtex_col = gtex_col[order(gtex_col$tissue), ]

### plot the factor matrix
plot_factor_matrix <- function(Factor_fn){
  factor_matrix = read.table(Factor_fn, sep='\t', stringsAsFactors = F, header = T, row.names = 1)
  rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
  gtex_col = gtex_col %>% filter(tissue %in% rownames(factor_matrix))

  nF = ncol(factor_matrix)
  
  # assign colors
  if(sum(!rownames(factor_matrix) %in% gtex_col$tissue) > 0){
    stop("Some features are not assigned colors to them")
  }
  factor_matrix_col = (abs(factor_matrix) > 1e-5 ) *1
  for (row in rownames(factor_matrix)){
    factor_matrix_col[row, factor_matrix_col[row, ]==1] = gtex_col %>% filter(tissue == row) %>% pull(color_hex)
    factor_matrix_col[row, factor_matrix_col[row, ]==0] = '#FFFFFF'
  }
  
  # tiles for plotting
  df = data.frame()
  for (i in seq(0, nrow(factor_matrix)-1)){
    dfi = data.frame(x1 = seq(1,nF) , x2 = seq(2, nF+1), 
                     y1 = rep(i,nF), y2 = rep(i+1, nF), 
                     col = factor_matrix_col[i+1, ])
    df = rbind(df, dfi)
  }
  
  
  main_theme = theme(axis.title.x = element_text(size = 5),
                     axis.title.y = element_text(size = 5),
                     axis.text.x = element_text(size = 5),
                     legend.key.size = unit(0.1, "inches"),
                     legend.position = 'top',
                     legend.margin=margin(0,0,0,1),
                     legend.box.margin=margin(-10,-10,-10,0),                        
                     legend.title = element_text(size = 5), 
                     legend.text  = element_blank()) +
    theme_bw()+
    background_grid(major = 'none', minor = 'none') 
  
  fig_factor_matrix = ggplot() +
    scale_x_continuous(name = 'Tissues', expand = c(0,0)) + 
    scale_y_continuous(name = 'Factors', expand = c(0,0)) + 
    coord_flip() + 
    geom_rect(data = df, mapping = aes(xmin =y1, xmax=y2, ymin=x1, ymax=x2), 
              fill = df$col,
              color = 'black', size = 0.01) + 
    xlab("Tissues") + 
    ylab("Factors")+
    main_theme + 
    theme_bw() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8)) + 
    theme(axis.ticks =element_blank()) +
    theme(axis.text = element_blank()) + 
    geom_hline(yintercept = seq(1,nF), color = "black", size=0.1) 
  
  
  save_fn = gsub('.txt', '.pdf', Factor_fn)
  pdf(save_fn, width = 0.1 * ncol(factor_matrix) + 1, 
      height = 0.1 * nrow(factor_matrix))
  print(fig_factor_matrix)
  dev.off()
}





### plot each factor in a panel
plot_factor <- function(Factor_fn, plot_panel_names = F){
  factor_matrix = read.table(Factor_fn, sep='\t', stringsAsFactors = F, header = T, row.names = 1)
  rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
  gtex_col = gtex_col %>% filter(tissue %in% rownames(factor_matrix))

  nF = ncol(factor_matrix)
  signs = sign(factor_matrix)
  
  factor_matrix = abs(factor_matrix)
  med_max = median(apply(factor_matrix, 2, max))
  factor_matrix = factor_matrix / matrix(rep((apply(factor_matrix, 2, max) / med_max), nrow(factor_matrix)), 
                                         nrow = nrow(factor_matrix), byrow = T)
  factor_matrix = factor_matrix * signs
  
  ## assign names to each factor
  tissues_in_factors_0 = apply(factor_matrix, 2, function(x) rownames(factor_matrix)[abs(x)>1e-5])
  tissues_in_factors = rep('0', length(tissues_in_factors_0))
  tissues_in_factors[sapply(tissues_in_factors_0, function(x) length(x)) > 40]  = 'Universal'
  for(k in seq(2, length(tissues_in_factors_0))){
    tissues_in_factors[k] = paste(unlist(tissues_in_factors_0[k]), collapse =';')
  }
  
  colnames(factor_matrix) = paste0("Factor", seq(1, dim(factor_matrix)[2]))
 
  factor_matrix$Tissue = factor(rownames(factor_matrix), levels = gtex_col$tissue)
  factor_bar_df = melt(factor_matrix, id.vars = 'Tissue')
  
  if(plot_panel_names){
    factor_bar_df$variable_plot = as.character(factor_bar_df$variable)
    for(k in seq(1,dim(factor_matrix)[2])){
      factor_bar_df[factor_bar_df$variable_plot == paste0('Factor', k), "variable_plot"] = paste0('Factor', k, '\n(', tissues_in_factors[k], ')')
    }
  }else{
    factor_bar_df$variable_plot = as.character(factor_bar_df$variable)
    for(k in seq(1,dim(factor_matrix)[2])){
      factor_bar_df[factor_bar_df$variable_plot == paste0('Factor', k), "variable_plot"] = paste0('Factor', k)
    }    
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
    ylim(min(0, min(factor_bar_df$value) - abs(min(factor_bar_df$value)) * 0.2),
         max(abs(factor_bar_df$value)) * 1.2)
  
  save_fn = gsub('.txt', '_factors.pdf', Factor_fn)
  pdf(save_fn, height = 1.8 * ceiling((ncol(factor_matrix)-1) / 4))
  print(fig_factors)
  dev.off()
}
