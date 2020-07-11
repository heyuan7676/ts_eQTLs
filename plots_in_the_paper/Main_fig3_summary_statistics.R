setwd('/Users/Yuan/Desktop/plots/')
library(cowplot)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(colortools)
library(VennDiagram)


tissues = read.table('data/tissues.txt', sep='\t', stringsAsFactors = F)
tissues = tissues$V1

tissue_sample_size = read.table('data/summary_table.tsv', sep='\t', header=T, row.names = 1)
tissue_sample_size = tissue_sample_size[rownames(tissue_sample_size) %in% tissues, ]

shared_factor_color = 'navy'
ts_factor_color = 'dodgerblue'
ts_factor_color_onlyts = 'lightskyblue'

factor_tissues = list(tissues, tissues[grep("Brain", tissues)], "Pituitary", "Spleen", 
                      c("Colon_Sigmoid", "Small_Intestine_Terminal_Ileum"), "Thyroid", "Adipose_Visceral_Omentum", 
                      "Nerve_Tibial", c("Adipose_Subcutaneous", "Breast_Mammary_Tissue"),
                      "Whole_Blood", tissues[grep("Cere", tissues)], 
                      tissues[grep("Heart", tissues)], tissues[grep("Skin", tissues)], 
                      "Cells_EBV-transformed_lymphocytes", 
                      c("Esophagus_Gastroesophageal_Junction", "Esophagus_Muscularis", "Colon_Transverse"), 
                      "Liver", "Stomach", "Testis", tissues[grep("Artery", tissues)], 
                      "Muscle_Skeletal", "Pancreas", "Lung", "Esophagus_Mucosa")
factor_nSample = sapply(factor_tissues, function(x) sum(tissue_sample_size[x, "nSample"]))

main_theme_Fig2 = theme_bw() + 
                theme(axis.title.x = element_text(size = 8),
                        axis.title.y = element_text(size = 8),
                        legend.title = element_text(size = 6), 
                        legend.text  = element_text(size = 5),
                        axis.text = element_text(size = 5),
                        legend.key.size = unit(0.1, "inches"),
                        legend.position = 'top',
                       legend.justification="left",
                        legend.margin=margin(0,0,0,1),
                        legend.box.margin=margin(-10,-10,-10,0)) +
                  background_grid(major = 'none', minor = 'none') 



### fig2 - A
sig_N_dat = read.table('data/Fig2_sig_prop.txt', sep='\t', header=T, stringsAsFactors = F)
colnames(sig_N_dat) = c("Factor", "Proportion")
sig_N_dat$Factor = factor_names
#sig_N_dat$Factor_ordered = factor(factor_names, levels = factor_names[rev(order(factor_nSample))])
sig_N_dat$Factor_ordered = factor(factor_names, levels = factor_names)
sig_N_dat$color = ifelse(sig_N_dat$Factor=="Ubiquitous", "U-eQTLs", "TS-eQTLs")
sig_N_dat$color = factor(sig_N_dat$color, level = c("U-eQTLs", "TS-eQTLs"))

fig2_A = ggplot(data =sig_N_dat, aes(x = Factor_ordered, y = Proportion, fill = color))+
  geom_bar(stat = 'identity') + 
  xlab("") + 
  ylab("Fraction of eQTLs")+
  main_theme_Fig2 + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(color = 'black'), 
        panel.border = element_blank(),
        axis.line.x = element_line(), 
        axis.line.y = element_line())+
  scale_fill_manual(name = "", values=c(shared_factor_color, ts_factor_color))


### Fig2 - B, C: number of factors per eQTL
number_factors <- function(fn, title_text, N0){
  #N0 = 5301827
  # N0 = 17480
  N_Ftor = read.table(fn, sep='\t', header=T, stringsAsFactors = F)
  
  x = N_Ftor[N_Ftor$group == 'ts_pairs', ]
  y = N_Ftor[N_Ftor$group == 'shared_pairs', ]
  a = sum(x[x$Number_F >=2 ,"Freq"]) + sum(y[y$Number_F >=2, "Freq"])
  b = sum(x[x$Number_F >=1 ,"Freq"]) + sum(y[y$Number_F >=1, "Freq"])
  print(paste0('At least 1 ts-factor: ',a,' eQTLs; At least 2 ts-factors: ', b))
  
  # higher-level summary
  at_least_one_factor = N_Ftor[N_Ftor$group == 'All', ]
  no_factor_N = N0 - sum(at_least_one_factor$Freq)
  
  N_Ftor = N_Ftor[N_Ftor$group != 'All', ]
  
  shared_with_ts = N_Ftor[N_Ftor$group == 'shared_pairs', ]
  shared_with_ts = shared_with_ts[shared_with_ts$Number_F > 0, ]
  shared_with_ts = sum(shared_with_ts$Freq)
  
  N_Ftor_total = as.data.frame(rbind(c(no_factor_N, "other_eQTLs"),
                                     c(sum(N_Ftor[N_Ftor$group == 'ts_pairs', "Freq"])+shared_with_ts, "TS-eQTLs"), 
                                     c(sum(N_Ftor[N_Ftor$group == 'shared_pairs', "Freq"]), "U-eQTLs")))
  
  # venn plot
  area1 = sum(N_Ftor[N_Ftor$group == 'shared_pairs', "Freq"]) 
  area2 = sum(N_Ftor[N_Ftor$group == 'ts_pairs', "Freq"]) + shared_with_ts
  cross.area =shared_with_ts
  
  png('results_refined/Fig2_part2_venn.png', width = 2, height = 2, units = 'in', res = 400)
  draw.pairwise.venn(area1, area2, cross.area,
                     category = c('U-eQTLs', "TS-eQTLs"),
                     fill = c(shared_factor_color, ts_factor_color_onlyts),
                     cat.pos = c(190, 180),
                     cex = 0.65,
                     cat.cex = 0.8,
                     lwd = 1,
                     fontfamily = rep("sans",3),
                     cat.fontfamily = rep("sans", 2),
                     col = 'white',
                     alpha = 1,
                     label.col = 'white')
  dev.off()
  
  colnames(N_Ftor_total) = c("Freq", "eQTLs")
  N_Ftor_total$Freq = as.numeric(as.character(N_Ftor_total$Freq))
  N_Ftor_total$Freq = N_Ftor_total$Freq / N0
  
  N_Ftor_total$eQTLs = as.character(N_Ftor_total$eQTLs)
  N_Ftor_total = N_Ftor_total[N_Ftor_total$eQTLs != 'other_eQTLs', ]
  N_Ftor_total$eQTLs = factor(N_Ftor_total$eQTLs, levels = c("U-eQTLs", "TS-eQTLs", "other_eQTLs"))
  
  N_Ftor_total_plot = ggplot(data = N_Ftor_total, aes(x=eQTLs, y = Freq, fill = eQTLs)) + 
    geom_bar(stat = 'identity') + 
    xlab("") + 
    ylab("Fraction of eQTLs") + 
    scale_fill_manual(name = "", values=c(shared_factor_color, ts_factor_color, "grey")) + 
    theme(legend.position = 'none') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
          axis.text.y = element_text(color = 'black')) 
  
  ### for shared and ts each
  N_Ftor$Freq = as.numeric(N_Ftor$Freq) / N0
  N_Ftor$Number_F = as.numeric(N_Ftor$Number_F)
  N_Ftor$group = factor(N_Ftor$group)
  print(paste0("Fraction of eQTLs with >5 TS-factors:", sum(N_Ftor[N_Ftor$Number_F >=5, "Freq"]) / sum(N_Ftor$Freq)))
  
  N_Ftor[(N_Ftor$Number_F == 10) & (N_Ftor$group == 'ts_pairs'), "Freq"] = sum(N_Ftor[(N_Ftor$Number_F >=10) & (N_Ftor$group == 'ts_pairs'), "Freq"])
  N_Ftor[(N_Ftor$Number_F == 10) & (N_Ftor$group == 'shared_pairs'), "Freq"] = sum(N_Ftor[(N_Ftor$Number_F >=10) & (N_Ftor$group == 'shared_pairs'), "Freq"])
  
  N_Ftor = N_Ftor[N_Ftor$Number_F <=10, ]
  N_Ftor[N_Ftor$Number_F==10, "Number_F"] = ">=10"
  N_Ftor$Number_F = factor(N_Ftor$Number_F, levels = c(as.character(seq(0,9)), ">=10"))
  
  number_factors_ggplot = ggplot(data =N_Ftor, aes(x = Number_F, y = Freq, fill = group)) + 
    geom_bar(stat = 'identity', position = 'stack') +
    #facet_wrap(~group) + 
    xlab("Number of tissue-specific factors") + 
    ylab("Fraction of eQTLs") + 
    ggtitle(title_text) + 
    scale_fill_manual(name = "eQTLs", values=c(shared_factor_color, ts_factor_color_onlyts), 
                      labels = c("With U-factor", "With only TS-factors(s)")) + 
    ggtitle("")
  
  return(list(N_Ftor_total_plot, number_factors_ggplot) )
}

### Fig2 - number of tissues per eQTL
number_tissues <-function(fn, N0){
  nT_df = read.table(fn, sep='\t', stringsAsFactors = F, header=T)
  nT_df$group = "ts_eQTLs"
  nT_df[nT_df$N_tissues >=49, "group"] = 'shared_eQTLs'
  
  nT_df[nT_df$N_tissues == 49, "eQTLs"] = sum(nT_df[nT_df$N_tissues >= 49, "eQTLs"])
  nT_df = nT_df[nT_df$N_tissues <= 49, ]
  nT_df$eQTLs = nT_df$eQTLs / N0
  
  fig2_D = ggplot(data = nT_df, aes(x = N_tissues, y = eQTLs, fill = group)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(name = "eQTLs", values=c(shared_factor_color, ts_factor_color_onlyts), 
                      labels = c("With U-factor", "With only TS-factors(s)")) + 
    theme(legend.position = 'none') + 
    ylab("Fraction of eQTLs") + 
    xlab("Number of tissues") + 
    main_theme_Fig2
  
  return(fig2_D)
}



N0 = 5301827
fn = 'data/Fig2_number_F_for_eQTL.txt'
fig2_B = number_factors(fn, 'spMF_eQTLs', N0)[[1]] + 
  main_theme_Fig2
fig2_C = number_factors(fn, 'spMF_eQTLs',N0)[[2]] + 
  main_theme_Fig2

fn = 'data/Fig2_number_tissues_for_eQTL.txt'
fig2_D = number_tissues(fn, N0)


fig2_part2 = ggdraw() + 
  draw_plot(fig2_A, x = 0.04, y = .4, width = 0.45, height = 0.55) + 
  draw_plot(fig2_B + theme(legend.position = 'none'), x = 0.52, y = 0.4435, width = 0.25, height = 0.58) + 
  draw_plot(fig2_C, x = 0.045, y = 0, width = 0.45, height = 0.45) + 
  draw_plot(fig2_D, x = 0.515, y = 0.0019, width = 0.45, height = 0.55)
  #draw_plot_label(label = c("A", "B", "C", "D", "E"), x = c(0, 0.48, 0.72, 0, 0.48)+0.06, y =c(0.98,0.98, 0.98, 0.51, 0.51), size = 8)
save_plot("results_refined/Fig2_part2_revision.png", fig2_part2, base_width = 7.4)





fn = 'data/Fig2_number_F_for_eQTL_flashr.txt'
fig2_B_flashr = number_factors(fn, 'flashr_eQTLs', N0)[[1]] + 
  main_theme_Fig2
fig2_C_flashr = number_factors(fn, 'flashr_eQTLs',N0)[[2]] + 
  main_theme_Fig2
fn = 'data/Fig2_number_tissues_for_eQTL_flashr.txt'
fig2_D_flashr = number_tissues(fn, N0)



fn = 'data/Fig2_number_F_for_eQTL_thresholding.txt'
fig2_B_thresholding = number_factors(fn, 'thresholding', N0)[[1]] + 
  main_theme_Fig2
fig2_C_thresholding = number_factors(fn, 'thresholding',N0)[[2]] + 
  main_theme_Fig2
fn = 'data/Fig2_number_tissues_for_eQTL_thresholding.txt'
fig2_D_thresholding = number_tissues(fn, N0)




fig2_part2_flashr = ggdraw() + 
  draw_plot(g_flashr_compare + theme(plot.title=element_text(size=8)) + theme(axis.text.x=element_text(angle = 45)), x = 0, y = 0.45, width = 0.6, height = 0.6) + 
  draw_plot(fig2_B_flashr, x = 0.6, y = 0.4, width = 0.3, height = 0.6) + 
  draw_plot(fig2_C_flashr, x = 0.0, y = 0, width = 0.54, height = 0.65) + 
  draw_plot(fig2_D_flashr, x = 0.47, y = 0, width = 0.55, height = 0.58)
save_plot("results_refined//Fig2_par2_flashr.png", fig2_part2_flashr, base_width = 7.4, base_height = 4)


fig2_part2_thresholding = ggdraw() + 
  draw_plot(g_thresholding_compare + theme(plot.title=element_text(size=8)) + theme(axis.text.x=element_text(angle = 45)), x = 0, y = 0.45, width = 0.6, height = 0.6) 
  #draw_plot(fig2_B_thresholding, x = 0.6, y = 0.4, width = 0.3, height = 0.6) + 
  #draw_plot(fig2_C_thresholding, x = 0.0, y = 0, width = 0.54, height = 0.65) + 
  #draw_plot(fig2_D_thresholding, x = 0.47, y = 0, width = 0.55, height = 0.58)
save_plot("results_refined//Fig2_par2_thresholding.png", fig2_part2_thresholding, base_width = 7.4, base_height = 4)



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
    main_theme_Fig2 + 
    scale_fill_manual(name = "eQTLs", values=tetradic("mediumblue"))
  
  save_plot("Supplmentary_plots//Compare_ts_genes.png", compare_g, base_width = 4, base_height = 2.4)
  
  return(N_Ftor_collect)
}
#N_Ftor_collect = number_factors_compare_methods()

