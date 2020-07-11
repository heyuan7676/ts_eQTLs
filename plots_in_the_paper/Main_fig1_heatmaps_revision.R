setwd('/Users/Yuan/Desktop/plots/')
library(ggplot2)
library(gplots)
library(reshape2)
library(gridExtra)
library(ComplexHeatmap)
library(cowplot)
library(lemon)
library(viridis)
library(dichromat)


### gtex colors
gtex_col = read.table('data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id),]
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)
rownames(gtex_col) = gtex_col$tissue_id

tissues = read.table('data/tissues.txt', sep='\t', header=F, stringsAsFactors = F)
tissues = tissues$V1

main_theme_Fig1 = theme_classic()+
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


plot_heatmap_Z_scores_B_values <- function(suffix, high_color, low_color, N_eQTLs = 50){
  factor_order_tp = seq(1,23)
  factor_order_tp = sample(factor_order_tp, size= length(factor_order_tp), replace =F)
  row_to_use = as.vector(sapply(factor_order_tp, function(x) seq(100*(x-1)+1, 100 * (x-1)+N_eQTLs)))
  
  nF_Y = read.table(paste0('data/Fig1_Y_',suffix,'.txt'), sep='\t', stringsAsFactors = F, header = T)
  nF_Y = nF_Y[,seq(3, dim(nF_Y)[2])]
  nF_Y = nF_Y[rev(seq(1, dim(nF_Y)[1])), ]
  nF_Y = nF_Y[row_to_use, ]
  
  nF_W = read.table(paste0('data/Fig1_W_',suffix,'.txt'), sep='\t', stringsAsFactors = F, header = T)
  nF_W = nF_W[,seq(3, dim(nF_W)[2])]
  nF_W = nF_W[rev(seq(1, dim(nF_W)[1])), ]
  nF_W = nF_W[row_to_use, ]
  
  nF_B = read.table(paste0('data/Fig1_B_',suffix,'.txt'), sep='\t', stringsAsFactors = F, header = T)
  nF_B = nF_B[,seq(3, dim(nF_B)[2])]
  nF_B = nF_B[rev(seq(1, dim(nF_B)[1])), ]
  nF_B = nF_B[row_to_use, ]
  
  #Z = as.matrix(nF_Y * nF_W)
  Z = as.matrix(nF_Y)
  Z[is.na(Z)] = 0
  
  # do clustering for eQTLs
  dat.dendro <- as.dendrogram(hclust(d = dist(x = abs(Z))))
  eQTL.order <- order.dendrogram(dat.dendro)

  Z = as.data.frame(Z)
  Z$eQTL = seq(1, dim(nF_Y)[1])
  colnames(Z) = c(tissues, "eQTL")
  
  Z_plot  = melt(Z, id.vars = 'eQTL')
  Z_plot$variable = factor(Z_plot$variable, levels = tissues)
  
  limit_value = 2
  Z_plot[Z_plot$value > limit_value, "value"] = limit_value
  Z_plot[Z_plot$value < -limit_value, "value"] = -limit_value
  
  heatmap_Zscore = ggplot(data = Z_plot, aes(x = variable, y = eQTL)) + 
    geom_tile(aes(fill = value), col = 'black', size=0.005) + 
    scale_fill_gradient2(high = high_color, mid = 'white', low = low_color, 
                         space ="Lab" ,
                         midpoint = 0, 
                         limits = c(-limit_value, limit_value)) + 
    ylab("") + 
    xlab("") + 
    main_theme_Fig1
  
  #tp = gsub("-", "", levels(Z_plot$variable))
  #color_bar = data.frame("x" = seq(1,49), "y" = rep(0,1), "color" = gtex_col[tp, "tissue_color_hex"])
  #heatmap_Zscore = heatmap_Zscore + geom_point(data = color_bar, aes(x = x, y = y), col = color_bar$color, size = 0.7)
  
  panel_width = unit(0.15, "npc")
  heatmap_Zscore = heatmap_Zscore + 
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,20,-20,-10)) +
    guides(fill= guide_colorbar(title = '', 
                                barwidth = panel_width,
                                title.position = 'top'))
  
  
  ### heatmap of B
  nF_B$eQTL = Z$eQTL
  B_plot = melt(nF_B, id.vars = 'eQTL')
  
  limit_value = max(c(abs(min(Z_plot$value)), max(Z_plot$value)))
  B_plot[B_plot$value < -limit_value, "value"] = -limit_value
  B_plot[B_plot$value > limit_value, "value"]  = limit_value
  
  heatmap_B = ggplot(data = B_plot, aes(x = variable, y = eQTL)) + 
    geom_tile(aes(fill = value), col = 'black', size = 0.005) + 
    scale_fill_gradient2(high = high_color, mid = 'white', low = low_color, 
                         midpoint = 0,
                         breaks = c(-0.2, 0, 0.2)) + 
    ylab("") + 
    xlab("") + 
    main_theme_Fig1
  
  panel_width = unit(0.075, "npc")
  heatmap_B = heatmap_B + 
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,20,-20,0)) + 
    guides(fill= guide_colorbar(title = '', 
                                barwidth = panel_width,
                                title.position = 'top'))
  
  
  ### heatmap of F
  factor_matrix = read.table('data/Fig1_factor_matrix.txt', sep='\t', stringsAsFactors = F, row.names = 1)
  factor_matrix = factor_matrix[rownames(factor_matrix) != "", ]
  rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
  factor_matrix$tissue = rownames(factor_matrix)
  
  F_plot = melt(factor_matrix, id.vars = 'tissue')
  F_plot$tissue = factor(as.character(F_plot$tissue), levels = gsub("-", "", levels(Z_plot$variable)))
  F_plot$variable = factor(F_plot$variable, levels = rev(as.character(unique(F_plot$variable))))
  #F_plot$variable = factor(F_plot$variable, levels = rev(levels(unique(F_plot$variable)))[factor_order_tp])
  
  heatmap_F = ggplot(data = F_plot, aes(x = tissue, y = variable)) + 
    geom_tile(aes(fill = value), col = 'black', size=0.005) + 
    scale_fill_gradient2(high = high_color, mid = 'white', low = low_color, midpoint = 0) + 
    ylab("") + 
    xlab("") + 
    scale_y_discrete(expand = c(0.07,0)) + 
    geom_hline(yintercept = seq(0,23) + 0.5, color = "black", size=0.08) + 
    main_theme_Fig1
  
  #tp = gsub("-", "", levels(F_plot$tissue))
  #color_bar = data.frame("x" = seq(1,49), "y" = rep(0,1), "color" = gtex_col[tp, "tissue_color_hex"])
  #heatmap_F = heatmap_F + geom_point(data = color_bar, aes(x = x, y = y), col = color_bar$color, size = 0.6)
  
  panel_width = unit(0.15, "npc")
  heatmap_F = heatmap_F + 
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,40,-5,0)) + 
    guides(fill= guide_colorbar(title = '', 
                                barwidth = panel_width,
                                title.position = 'top'))
  
  ## heatmap of weights
  #W = as.matrix(1/nF_W)
  W = nF_W
  W[is.na(W)] = 0
  W = as.data.frame(W)
  
  W$eQTL = seq(1, dim(nF_Y)[1])
  colnames(W) = c(tissues, "eQTL")
  
  W_plot  = melt(W, id.vars = 'eQTL')
  W_plot$variable = factor(W_plot$variable, levels = levels(Z_plot$variable))
  
  W_limit = 100
  W_plot[W_plot$value > W_limit, "value"] = W_limit
  heatmap_W = ggplot(data = W_plot, aes(x = variable, y = eQTL)) + 
    geom_tile(aes(fill = value), col = 'white', size=0.005) + 
    scale_fill_gradient2(high = high_color, 
                         mid = 'white', 
                         low = low_color, 
                         midpoint = 0) + 
    ylab("") + 
    xlab("") + 
    main_theme_Fig1
  
  #tp = gsub("-", "", levels(W_plot$variable))
  #color_bar = data.frame("x" = seq(1,49), "y" = rep(0,1), "color" = gtex_col[tp, "tissue_color_hex"])
  #heatmap_W = heatmap_W + geom_point(data = color_bar, aes(x = x, y = y), col = color_bar$color, size = 0.7)
  
  panel_width = unit(0.15, "npc")
  heatmap_W = heatmap_W + 
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,20,-20,-10)) + 
    guides(fill= guide_colorbar(title = '', 
                                barwidth = panel_width,
                                title.position = 'top'))
  
  
  return(list(heatmap_Zscore, heatmap_B, heatmap_F, heatmap_W))
}


gp = plot_heatmap_Z_scores_B_values('onefactor', 'red', 'blue', N_eQTLs = 50)

Fig1 = ggdraw() + 
  draw_plot(gp[[1]] + theme(axis.line = element_blank()), x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(gp[[2]] + theme(axis.line = element_blank()), x = 0.25, y = 0, width = 0.15, height = 1) + 
  draw_plot(gp[[3]] + theme(axis.line = element_blank()), x = 0.4, y = 0.4, width = 0.25, height = 0.35) + 
  draw_plot(gp[[4]] + theme(axis.line = element_blank()), x = 0.7, y = 0, width = 0.25, height = 1)

save_plot("results_refined/Fig1_heatmaps_revision_2.png", Fig1, base_width = 7.4, base_height = 4)
