

### plot the input matrix
plot_X_matrix <- function(input_Xfn, N_eQTLs = 50,  
                                           high_color = 'red', low_color = 'blue', 
                                           nF = 23){
  factor_order_tp = sample(seq(1,nF), size= length(factor_order_tp), replace =F)
  row_to_use = as.vector(sapply(factor_order_tp, function(x) seq(100*(x-1)+1, 100 * (x-1)+N_eQTLs)))
  
  X_matrix = read.table(input_Xfn, sep='\t', stringsAsFactors = F, header = T)
  X_matrix = X_matrix[,seq(3, dim(X_matrix)[2])]
  X_matrix = X_matrix[rev(seq(1, dim(X_matrix)[1])), ]
  X_matrix = X_matrix[row_to_use, ]
  
  X_matrix = as.matrix(X_matrix)
  X_matrix[is.na(X_matrix)] = 0
  
  # do clustering for eQTLs
  dat.dendro <- as.dendrogram(hclust(d = dist(x = abs(X_matrix))))
  eQTL.order <- order.dendrogram(dat.dendro)
  
  X_matrix = as.data.frame(X_matrix)
  X_matrix$eQTL = seq(1, dim(X_matrix)[1])
  colnames(X_matrix) = c(gtex_col$tissue, "eQTL")
  
  X_plot  = melt(X_matrix, id.vars = 'eQTL')
  X_plot$variable = factor(X_plot$variable, levels = gtex_col$tissue)
  
  limit_value = 2
  X_plot[X_plot$value > limit_value, "value"] = limit_value
  X_plot[X_plot$value < -limit_value, "value"] = -limit_value
  
  
  
  main_theme = theme_classic()+
    theme(axis.text.x = element_blank()) + 
    theme(axis.text.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) + 
    theme(axis.ticks.y = element_blank()) +  
    theme(axis.line.x  = element_blank()) +
    theme(axis.line.y = element_blank()) +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_blank(), 
          legend.text  = element_text(size = 4),
          legend.key.size = unit(0.05, "inches"),
          legend.position = 'top')
  
  heatmap_Zscore = ggplot(data = X_plot, aes(x = variable, y = eQTL)) + 
    geom_tile(aes(fill = value), col = 'black', size=0.005) + 
    scale_fill_gradient2(high = high_color, mid = 'white', low = low_color, 
                         space ="Lab" ,
                         midpoint = 0, 
                         limits = c(-limit_value, limit_value)) + 
    ylab("") + 
    xlab("") +
    main_theme
  
  return(heatmap_Zscore)
}



