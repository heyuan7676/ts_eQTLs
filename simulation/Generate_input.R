suppressWarnings(library('RColorBrewer'))
suppressWarnings(library('pheatmap'))
suppressWarnings(library('cowplot'))
suppressWarnings(library('reshape2'))
suppressWarnings(library('R.matlab'))


generate_input <- function(N=100, P = 10, K=5, saveplot = F, tau, savedir, seed){
  set.seed(seed)
  dir.create(savedir, showWarnings = F)
  
  ## generate loading values
  # probability that a data point has each of the factor
  loading_pattern = list()
  lp_idx = 1
  for(i in seq(1,K)){
    tp = utils::combn(K, i)
    for(col in seq(1, dim(tp)[2])){
      loading_pattern[[lp_idx]] = tp[,col]
      lp_idx = lp_idx + 1
    }
  }
  
  loading_pattern_prob = c(0.0, 
                           0.25, 0.2, 0.1, 0.15, 0.1, 
                           0.05, 0,0, 0,
                           0.05, 0.025, 0.025,0,
                           0, 0, 0, 0.05,
                           rep(0, length(loading_pattern)-(5+4+4+4)))
  
  l = matrix(rep(0, N*K), nrow = N)
  l_idx = sample(seq(1, length(loading_pattern)+1), N, replace = TRUE, prob = loading_pattern_prob)
  for(i in seq(1,N)){
    if(l_idx[i] > 1){
      for(ij in loading_pattern[[l_idx[i]-1]]){
        l[i, ij] = rnorm(n=1)
      }
    }
  }
  
  l = l[order(l_idx), ]
  colorramp = brewer.pal(11,"RdBu")
  applycolors = colorRampPalette(colorramp)
  
  
  l_melt = data.frame(l)
  colnames(l_melt) = seq(1,5)
  l_melt$tissue = rownames(l_melt)
  l_melt = melt(l_melt)
  l_melt$tissue = as.numeric(l_melt$tissue)
  
  if(saveplot){
    heatmap_g = ggplot(data = l_melt, aes(x = tissue, y = variable)) + 
      geom_tile(aes(fill = value)) +
      scale_fill_gradient2(high = 'orange', mid = 'white', low = 'darkblue', midpoint =  0, limits=c(min(l_melt$value),max(l_melt$value))) +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
      ylab("Loadings") + xlab("Data points")  + 
      background_grid(major = 'xy', minor = 'none') + 
      theme(legend.position="top") + 
      labs(fill = "") 
    
    panel_width = unit(0.8,"npc") - sum(ggplotGrob(heatmap_g)[["widths"]][-3]) - unit(1,"line")
    heatmap_g = heatmap_g + guides(fill= guide_colorbar(barwidth = panel_width))
    
    save_plot(paste0(savedir, "/Simulated_Loading_tau",tau,"_seed",seed, ".png"), 
              heatmap_g, base_width = 3.5, base_height = 8)
  }

  ## generate factor values
  f = matrix(rep(0, P * K), nrow = P)
  f[,1] = 1
  f[,2] = c(rep(0, 5), 1, rep(0,4))
  f[,3] = c(rep(0, 2), rep(1,2), rep(0,6))
  f[,4] = c(rep(0, 6), 1, rep(0,3))
  f[,5] = c(rep(0, 8), rep(1,2))
  
  f_melt = data.frame(f)
  colnames(f_melt) = seq(1,5)
  f_melt$tissue = rownames(f_melt)
  f_melt = melt(f_melt)
  f_melt$tissue = as.numeric(f_melt$tissue)
  
  if(saveplot){
  heatmap_g = ggplot(data = f_melt, aes(x = tissue, y = variable)) + 
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = '#91CF60', mid = 'white', low = 'white', midpoint =  0, limits=c(min(f_melt$value),max(f_melt$value))) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    ylab("Tissues") + xlab("Factors")  + 
    background_grid(major = 'xy', minor = 'none') + 
    theme(legend.position="top") + 
    labs(fill = "") 
  
  panel_width = unit(0.8,"npc") - sum(ggplotGrob(heatmap_g)[["widths"]][-3]) - unit(1,"line")
  heatmap_g = heatmap_g + guides(fill= guide_colorbar(barwidth = panel_width))
  
  save_plot(paste0(savedir, "/Simulated_Factors_tau",tau,"_seed",seed, ".png"), 
            heatmap_g, base_width = 5, base_height = 3)
  }
  
  
  ## generate errors
  E = matrix(nrow = N, ncol = P)
  for(i in seq(1,N)){
    for(j in seq(1,P)){
      E[i,j] = rnorm(1, sd = sqrt(1/tau))
    }
  }
  
  e_melt = data.frame(E)
  colnames(e_melt) = seq(1,10)
  e_melt$tissue = rownames(e_melt)
  e_melt = melt(e_melt)
  e_melt$tissue = as.numeric(e_melt$tissue)
  
  ## generate input matrix
  X = as.data.frame(l %*% t(f) + E)
  
  x_melt = data.frame(X)
  colnames(x_melt) = seq(1,10)
  x_melt$tissue = rownames(x_melt)
  x_melt = melt(x_melt)
  x_melt$tissue = as.numeric(x_melt$tissue)
  
  if(saveplot){
  heatmap_g = ggplot(data = x_melt, aes(x = tissue, y = variable)) + 
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = '#d13303', mid = 'white', low = '#4E84C4', midpoint =  0, limits=c(min(x_melt$value),max(x_melt$value))) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    ylab("Tissues") + xlab("Data points")  + 
    background_grid(major = 'xy', minor = 'none') + 
    theme(legend.position="top") + 
    labs(fill = "") 
  
  panel_width = unit(0.8,"npc") - sum(ggplotGrob(heatmap_g)[["widths"]][-3]) - unit(1,"line")
  heatmap_g = heatmap_g + guides(fill= guide_colorbar(barwidth = panel_width))
  save_plot(paste0(savedir, "/Simulated_Input_tau",tau,"_seed",seed, ".png"), 
            heatmap_g, base_width = 5, base_height = 8)
  
  
  heatmap_g = ggplot(data = e_melt, aes(x = tissue, y = variable)) + 
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = '#d13303', mid = 'white', low = '#4E84C4', midpoint =  0, limits=c(min(x_melt$value),max(x_melt$value))) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    ylab("Tissues") + xlab("Data points")  + 
    background_grid(major = 'xy', minor = 'none') + 
    theme(legend.position="top") + 
    labs(fill = "") 
  
  panel_width = unit(0.8,"npc") - sum(ggplotGrob(heatmap_g)[["widths"]][-3]) - unit(1,"line")
  heatmap_g = heatmap_g + guides(fill= guide_colorbar(barwidth = panel_width))
  
  save_plot(paste0(savedir, "/Simulated_Error_tau",tau,"_seed",seed, ".png"),  
            heatmap_g, base_width = 5, base_height = 8)
  }
  
  print(paste0('abs mean Error / abs mean X = ', mean(abs(E)) / mean(abs(l %*% t(f)))))
  print(paste0('var Error / var X = ', (1/tau) / var(c(as.matrix(X)))))
  
  W = matrix(1, ncol = ncol(X), nrow = nrow(X))
  

  X = cbind(paste0("Gene", seq(1, nrow(X))), paste0("SNP", seq(1, nrow(X))), X)
  W = cbind(paste0("Gene", seq(1, nrow(X))), paste0("SNP", seq(1, nrow(X))), W)
  colnames(X) = c("Gene", "SNP", paste0("Feature", seq(1, ncol(X)-2)))
  colnames(W) = c("Gene", "SNP", paste0("Feature", seq(1, ncol(X)-2)))


  tempfn=paste0(savedir, '/X_tau', tau, '_seed', seed, '.mat')
  print(tempfn)
  tp = X[, seq(3, ncol(X))]
  writeMat(paste0(savedir, '/X_tau', tau, '_seed', seed, '.mat'), X = tp)
  save(X, f, l, W, file = paste0(savedir, '/Input_tau', tau, '_seed', seed, '.RData'))
  write.table(X, paste0(savedir, '/Input_tau', tau, '_seed', seed, '_X.txt'), sep='\t', quote = F, row.names = F)
  write.table(W, paste0(savedir, '/Input_tau', tau, '_seed', seed, '_W.txt'), sep='\t', quote = F, row.names = F)
  write.table(f, paste0(savedir, '/Input_tau', tau, '_seed', seed, '_f.txt'), sep='\t', quote = F, row.names = F)
  write.table(l, paste0(savedir, '/Input_tau', tau, '_seed', seed, '_l.txt'), sep='\t', quote = F, row.names = F)
}





library(optparse)



option_list = list(
        make_option(c("-t", "--tau"), type = "double", default=10, help="Precision of error"),
        make_option(c("-s", "--seed"), type = "integer", default=1, help="Seed to generate random input"),
        make_option(c("-d", "--dir"), type = "character", default='simulation/input/', help="directory to save the generated input"))


opt = parse_args(OptionParser(option_list=option_list))
generate_input(tau = opt$tau, savedir = opt$dir, seed = opt$seed)
