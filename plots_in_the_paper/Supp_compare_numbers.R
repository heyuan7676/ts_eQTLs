setwd('/Users/Yuan/Desktop/plots/')
library(cowplot)
library(ggplot2)
library(reshape2)
library(gridExtra)

### factor names

shared_factor_color = 'navy'
ts_factor_color = 'dodgerblue'
ts_factor_color_unmatched = 'lightskyblue'


### gtex colors
gtex_col = read.table('data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id),]
rownames(gtex_col) = gtex_col$tissue_id
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)


main_theme = theme(axis.title.y = element_text(size=6)) +
  theme(axis.title.x = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 5)) + 
  theme(axis.text.y = element_text(size = 5)) + 
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6))

### fig2 - A
sig_N_dat = read.table('data/Fig2_ts_shared_pairs_flashr.txt', sep='\t', header=T, stringsAsFactors = F)
sig_N_dat$Factor = factor_names
sig_N_dat$Factor_ordered = factor(factor_names, levels = factor_names)

fig2_A = ggplot(data =sig_N_dat, aes(x = Factor_ordered, y = X0, fill = factor(ifelse(Factor=="Shared", "shared", "tissue-specific")))) + 
  geom_bar(stat = 'identity') + 
  xlab("") + 
  ylab("Proportion of eQTLs")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = c(0.8,0.8)) +
  theme(legend.key.size = unit(0.03, "inches"))+
  scale_fill_manual(name = "eQTLs", values=c(shared_factor_color, ts_factor_color)) +
  background_grid(major = 'xy', minor = 'none') + 
  main_theme


### Fig2 - B: number of factors per eQTL
number_factors <- function(fn, title_text){
  N_Ftor = read.table(fn, sep='\t', header=T, stringsAsFactors = F)
  N_Ftor = N_Ftor[N_Ftor$Number_F > 0, ]
  N_Ftor = N_Ftor[N_Ftor$group != 'All', ]
  
  number_factors_ggplot = ggplot(data =N_Ftor, aes(x = Number_F, y = Freq, fill = group)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    xlab("Number of factors") + 
    ylab("Frequency") + 
    ggtitle(title_text) + 
    main_theme + 
    scale_fill_manual(name = "eQTLs", values=c(shared_factor_color, ts_factor_color))
 
  return(number_factors_ggplot) 
}

f1 = number_factors('data/Fig2_number_F_for_eQTL.txt', 'spMF_eQTL')
f2 = number_factors('data/Fig2_number_F_for_eQTL_flashr.txt', 'flashr_eQTL')
f3 = number_factors('data/Fig2_number_F_for_genes.txt', 'spMF_genes')
f4 = number_factors('data/Fig2_number_F_for_genes_flashr.txt', 'flashr_genes')

fig2_g = ggdraw() + 
  draw_plot(f1, x = 0, y = 0, width = 0.5, height = 0.5) + 
  draw_plot(f2, x = 0.5, y = 0, width = 0.5, height = 0.5) + 
  draw_plot(f3, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(f4, x = 0.5, y = 0.5, width = 0.5, height = 0.5)

save_plot("results_refined/Fig2_number_factors_temp.png", fig2_g, base_width = 7.4)



number_factors_compare_methods <- function(){
  
  N_Ftor_collect = data.frame()
  
  fn = 'data/Fig2_number_F_for_genes.txt'
  N_Ftor = read.table(fn, sep='\t', header=T, stringsAsFactors = F)
  N_Ftor = N_Ftor[N_Ftor$Number_F > 0, ]
  N_Ftor = N_Ftor[N_Ftor$group == 'ts_genes', ]
  N_Ftor$group = 'spMF'
  N_Ftor_collect = rbind(N_Ftor_collect, N_Ftor)
  
  fn = 'data/Fig2_number_F_for_genes_flashr.txt'
  N_Ftor = read.table(fn, sep='\t', header=T, stringsAsFactors = F)
  N_Ftor = N_Ftor[N_Ftor$Number_F > 0, ]
  N_Ftor = N_Ftor[N_Ftor$group == 'ts_genes', ]
  N_Ftor$group = 'flashr'
  N_Ftor_collect = rbind(N_Ftor_collect, N_Ftor)

  fn = 'data/Fig2_number_F_for_genes_thresholding.txt'
  N_Ftor = read.table(fn, sep='\t', header=T, stringsAsFactors = F)
  N_Ftor = N_Ftor[N_Ftor$Number_F > 0, ]
  N_Ftor = N_Ftor[N_Ftor$group == 'ts_genes', ]
  N_Ftor$group = 'thresholding'
  N_Ftor_collect = rbind(N_Ftor_collect, N_Ftor)
  
  fn = 'data/Fig2_number_F_for_genes_cbset.txt'
  N_Ftor = read.table(fn, sep='\t', header=T, stringsAsFactors = F)
  N_Ftor = N_Ftor[N_Ftor$Number_F > 0, ]
  N_Ftor = N_Ftor[N_Ftor$group == 'ts_genes', ]
  N_Ftor$group = 'cbset'
  N_Ftor_collect = rbind(N_Ftor_collect, N_Ftor)
  
  N_Ftor_collect = N_Ftor_collect[,c("group", "Number_F", "Freq")]
  N_Ftor_collect = dcast(N_Ftor_collect, Number_F~group)
  
  N_Ftor_collect[!complete.cases(N_Ftor_collect$spMF), "spMF"] = 0
  N_Ftor_collect[!complete.cases(N_Ftor_collect$flashr), "flashr"] = 0
  N_Ftor_collect[!complete.cases(N_Ftor_collect$thresholding), "thresholding"] = 0
  N_Ftor_collect[!complete.cases(N_Ftor_collect$cbset), "cbset"] = 0
  
  N_Ftor_collect = melt(N_Ftor_collect, id.vars = 'Number_F')
  
  compare_g = ggplot(data =N_Ftor_collect, aes(x = Number_F, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    xlab("Number of factors") + 
    ylab("Frequency of ts-genes") + 
    ggtitle("Compare three methods to genes in credible sets") + 
    main_theme + 
    scale_fill_manual(name = "eQTLs", values=tetradic("mediumblue"))
  
  save_plot("Supplmentary_plots//Compare_ts_genes.png", compare_g, base_width = 4, base_height = 2.4)
  
  return(N_Ftor_collect)
}
N_Ftor_collect = number_factors_compare_methods()


### Examples


readin <- function(group, suffix){
  Y_df = read.table(paste0('data/examples/Fig2_example_',group,'_Y_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  W_df = read.table(paste0('data/examples/Fig2_example_',group,'_W_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  W_df = 1 / W_df
  P_df = read.table(paste0('data/examples/Fig2_example_',group,'_PV_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  B_df = read.table(paste0('data/examples/Fig2_example_',group,'_B_',suffix,'.txt'), sep='\t', header=T, row.names = 1)
  return(list(Y_df, W_df, P_df, B_df)  )
}


#plot_one_example <- function(group, suffix, pair_name){
#data = readin(group, suffix)
#Y_df = data[[1]][pair_name, ]
#W_df = data[[2]][pair_name, ]
#P_df = data[[3]][pair_name, ]
#B_df = data[[4]][pair_name, ]
#pr_data_point <- t(data.frame(rbind(Y_df, W_df)))

plot_one_example <- function(genei, snpi){
  
  pr_data_point = read.table(paste0('data/examples/Fig3_example_',genei,'_', snpi, '.txt'), sep='\t', header=T, row.names = 1)
  colnames(pr_data_point) = c("effect_size", "standard_error")
  pr_data_point$standard_error = 1/pr_data_point$standard_error
  rownames(pr_data_point) = tissues
  pr_data_point = as.data.frame(pr_data_point)
  
  conf.level = 0.95
  ci.vals = -qnorm( ( 1 - conf.level ) / 2 )
  pr_data_point$lower = pr_data_point$effect_size - ci.vals * pr_data_point$standard_error
  pr_data_point$higher = pr_data_point$effect_size + ci.vals * pr_data_point$standard_error
  pr_data_point$tissue = rownames(pr_data_point)
  pr_data_point$col    = gtex_col[pr_data_point$tissue, "tissue_color_hex"]
  
  pd_g = ggplot() +
    geom_pointrange(data=pr_data_point, aes(x=tissues, y=effect_size, ymin=lower, ymax=higher,color=col), fatten=0.5) + 
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis) 
    scale_color_identity()+
    xlab("") + ylab("Effect size") +
    main_theme  + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank())
  
  return(pd_g)
}

x = plot_one_example('ENSG00000198626.15', 'chr1_237071945')
y = plot_one_example('ENSG00000057294.14', 'chr12_32835973')
z = plot_one_example('ENSG00000196091.13', 'chr12_101681700')


x = plot_one_example('group1','cisInOne', 'ENSG00000141564.13:chr17_80739799')
y = plot_one_example('group0','cisInMore', 'ENSG00000172062.16:chr5_70659368')
z = plot_one_example('group1','cisInMore', 'ENSG00000162971.10:chr2_199896221')

fig2_g = ggdraw() + 
  draw_plot(fig2_A, x = 0, y = .5, width = 0.4, height = 0.5) + 
  draw_plot(fig2_B, x = 0, y = 0, width = 0.4, height = 0.55) + 
  draw_plot(x, x = 0.4, y = 0, width = 0.2, height = 1) +
  draw_plot(y, x = 0.6, y = 0, width = 0.2, height = 1) +
  draw_plot(z, x = 0.8, y = 0, width = 0.2, height = 1) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), x = c(0,0,0.4,0.6,0.8), y =c(1,0.55,1,1,1), size = 8)

save_plot("results_refined/Fig2_temp.png", fig2_g, base_width = 7.4)

