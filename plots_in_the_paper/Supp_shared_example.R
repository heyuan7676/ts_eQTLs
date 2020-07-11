setwd('/Users/Yuan/Desktop/plots/')
library(ggplot2)
library(reshape2)
library(gridExtra)

### tissues
tissues = read.table('data/tissues.txt', sep='\t', header=F, stringsAsFactors = F)
tissues = tissues$V1
nT = length(tissues)


### gtex colors
gtex_col = read.table('data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id),]
rownames(gtex_col) = tissues
gtex_col$tissue_id = tissues
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)
head(gtex_col)


### 

readin <- function(group, suffix){
  Y_df = read.table(paste0('data/examples/Fig2_example_',group,'_Y_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  W_df = read.table(paste0('data/examples/Fig2_example_',group,'_W_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  W_df = 1 / W_df
  #P_df = Y_df
  #B_df = Y_df
  #P_df = read.table(paste0('data/examples/Fig2_example_',group,'_PV_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  #B_df = read.table(paste0('data/examples/Fig2_example_',group,'_B_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  return(list(Y_df, W_df)  )
}
plot_one_vec <- function(pr_data_point){
  
  fsize=8*3/3
  
  gtex.theme = theme_classic() + 
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_text(size = fsize),
          axis.title = element_text(size = fsize), 
          legend.text = element_text(size = fsize), 
          legend.title = element_text(size = fsize),
          strip.background = element_blank())
  
  pr_data_point = pr_data_point[seq(nrow(pr_data_point), 1), ]
  pd_g = ggplot() +
    geom_pointrange(data=pr_data_point, aes(x=tissues, y=effect_size, ymin=lower, ymax=higher,color=col), fatten=0.5) + 
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    #annotate("text", x = cis_eQTL_tis-0.5, y = pr_data_point[cis_eQTL_tis, 'higher']+0.1, label = "*") + 
    coord_flip() +  # flip coordinates (puts labels on y axis) 
    scale_color_identity()+
    xlab("") + ylab("Effect size") +
    gtex.theme  + 
    theme(axis.ticks.y = element_blank(), axis.line.y=element_blank(),
          axis.text.x = element_text(size = 10, colour = 'black'),
          axis.title.x = element_text(size = 12))
  
  return(pd_g)
  
}
plot_one_example <- function(group, suffix, pair_name){
  data = readin(group, suffix)
  Y_df = data[[1]][pair_name, ]
  W_df = data[[2]][pair_name, ]
  
  pr_data_point <- t(data.frame(rbind(Y_df, W_df)))
  colnames(pr_data_point) = c("effect_size", "standard_error")
  rownames(pr_data_point) = tissues
  pr_data_point = as.data.frame(pr_data_point)
    
  conf.level = 0.95
  ci.vals = -qnorm( ( 1 - conf.level ) / 2 )
  pr_data_point$lower = pr_data_point$effect_size - ci.vals * pr_data_point$standard_error
  pr_data_point$higher = pr_data_point$effect_size + ci.vals * pr_data_point$standard_error
  pr_data_point$tissue = rownames(pr_data_point)
  pr_data_point$col    = gtex_col[pr_data_point$tissue, "tissue_color_hex"]
    
  pd_g = plot_one_vec(pr_data_point)
  return(pd_g)
  #pdf(paste0('results/Fig2_example',gsub(":","_",pair_name),'.pdf'), width = 2, height = 4)
  #png(paste0('results/Fig2_example',gsub(":","_",pair_name),'.png'), width = 2, height = 4, units = 'in', res = 200)
  #print(pd_g)
  #dev.off()
    
}

plot_example <- function(group, suffix, numberPlots){
  data = readin(group, suffix)
  Y_df = data[[1]]
  W_df = data[[2]]
  #P_df = data[[3]]
  #B_df = data[[4]]
  
  for (i in seq(1, numberPlots)){
    #print(i)
    
    pr_data_point <- t(data.frame(rbind(Y_df[i, ], W_df[i, ])))
    colnames(pr_data_point) = c("effect_size", "standard_error")
    rownames(pr_data_point) = tissues
    pr_data_point = as.data.frame(pr_data_point)
    
    conf.level = 0.95
    ci.vals = -qnorm( ( 1 - conf.level ) / 2 )
    pr_data_point$lower = pr_data_point$effect_size - ci.vals * pr_data_point$standard_error
    pr_data_point$higher = pr_data_point$effect_size + ci.vals * pr_data_point$standard_error
    pr_data_point$tissue = rownames(pr_data_point)
    pr_data_point$col    = gtex_col[pr_data_point$tissue, "tissue_color_hex"]
    
    pd_g = plot_one_vec(pr_data_point) + theme(axis.text.x = element_text(size = 8))
    
    pdf(paste0('results/',suffix,'/', group, '_',suffix,'_', gsub(":","_", rownames(Y_df)[i]), '.pdf'), width = 2, height = 4)
    print(pd_g)
    dev.off()
    
  }
  
}



x = plot_one_example('group1','cisInOne', 'ENSG00000141564.13:chr17_80739799')
y = plot_one_example('group0','cisInMore', 'ENSG00000172062.16:chr5_70659368')
z = plot_one_example('group1','cisInMore', 'ENSG00000162971.10:chr2_199896221')


for(g in seq(1,23)){
  print(g)
  #plot_example(paste0("group", g), 'cisInMore', 100) 
  plot_example(paste0("group", g), 'InTwoThreeFactors', 20) 
}

plot_example(paste0("group", 0), 'sn_spMF_not_flashrNN', 100)

plot_example('groupBrain', 'heuristic_1_not2', 100)


x = plot_one_example('group0', 'sn_spMF_not_flashrNN', 'ENSG00000197134.11:chr19_21938491')


data = readin('Spleen', 'cisNotinFactor')

pr = 'ENSG00000213390.10:chr10_97289253'
Y_df = data[[1]][pr,]
W_df = data[[2]][pr,]
P_df = data[[3]][pr,]
B_df = data[[4]][pr,]



#### plot the factors

X = read.table('data/Fig1_factor_matrix.txt', sep='\t', stringsAsFactors = F, row.names = 1)
X = X[rownames(X) != "", ]
max_median = max(apply(X,2,median))
X = apply(X, 2, function(x) max_median/max(x)*x)

Y_fitted = X %*% t(B_df)

for(k in seq(1,1)){
  shared_factor = as.data.frame(X[,k])
  colnames(shared_factor) = c("effect_size")
  shared_factor$lower = shared_factor$effect_size
  shared_factor$higher = shared_factor$effect_size
  shared_factor$tissue = rownames(shared_factor)
  shared_factor$col = gtex_col$tissue_color_hex
  shared_factor_g = plot_one_vec(shared_factor)
  
  
  pdf(paste0('results/factors/Fig0_factor',k-1,'.pdf'), width = 1.5, height = 4)
  print(shared_factor_g)
  dev.off() 
}

