setwd('/Users/Yuan/Desktop/plots/')
library(cowplot)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(colortools)
library(colorspace)

tissues = read.table('data/tissues.txt', sep='\t', stringsAsFactors = F)
tissues = tissues$V1
nT = length(tissues)

### gtex colors
gtex_col = read.table('data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id),]
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)
rownames(gtex_col) = gtex_col$tissue_id

main_theme_Fig2 = theme(axis.title.x = element_text(size = 5),
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



### Examples

### factor matrix

plot_one_example_yw <- function(genei, snpi){
  
  pr_data_point = read.table(paste0(local_dir, '/Fig3_example_yw_',genei,':', snpi, '.txt'), sep='\t', header=T, row.names = 1)
  colnames(pr_data_point) = c("effect_size", "standard_error")
  pr_data_point$standard_error = 1/pr_data_point$standard_error
  rownames(pr_data_point) = gsub("-", "", tissues)
  pr_data_point = as.data.frame(pr_data_point)
  
  pr_data_point[is.na(pr_data_point$effect_size), "effect_size"] = 0
  pr_data_point[is.na(pr_data_point$standard_error), "standard_error"] = 0
  
  conf.level = 0.95
  ci.vals = -qnorm( ( 1 - conf.level ) / 2 )
  pr_data_point$lower = pr_data_point$effect_size - ci.vals * pr_data_point$standard_error
  pr_data_point$higher = pr_data_point$effect_size + ci.vals * pr_data_point$standard_error
  pr_data_point$col    = gtex_col[rownames(pr_data_point), "tissue_color_hex"]
  
  pr_data_point = pr_data_point[seq(nrow(pr_data_point), 1), ]
  pr_data_point = pr_data_point[seq(nrow(pr_data_point), 1), ]
  
  pr_data_point = pr_data_point[seq(nrow(pr_data_point), 1), ]
  pr_data_point$tissue = factor(rownames(pr_data_point), levels = rownames(pr_data_point))
  
  pd_g = ggplot() +
    geom_pointrange(data=pr_data_point, aes(x=tissues, y=effect_size, ymin=lower, ymax=higher,color=col), fatten=0.5) + 
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis) 
    scale_color_identity()+
    xlab("") + ylab("Effect size") +
    main_theme_Fig2  + 
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.line.y=element_blank())
  
  return(pd_g)
}

plot_one_example_pb <-function(genei, snpi, breaks_vec = NULL, labels_vec = NULL){
  df_beta = read.table(paste0(local_dir, '/Fig3_example_pb_',genei,':', snpi, '.txt'), sep='\t', header=T)
  print(c(min(df_beta$beta), max(df_beta$beta)))
  print(factor_names[which(df_beta$beta_corrected !=0)])
  print(df_beta$beta_corrected)
  
  
  idx_x = abs(df_beta$beta) >= min(abs(df_beta[df_beta$beta_corrected!=0,"beta"]))
  idx_y = df_beta$beta_corrected == 0
  df_beta[which(idx_x * idx_y == 1), "beta"] = 0
  
  #df_beta[df_beta$pv == -1, "beta"] = 0
  df_beta$xmin = 0
  df_beta$xmax = 1
  df_beta$ymin = rev(seq(0,22))
  df_beta$ymax = rev(seq(1,23))
  
  df_beta$pv = as.numeric(df_beta$pv)
  df_beta$pv[(df_beta$pv < 0.05) & (df_beta$beta_corrected == 0)] = 1
  
  df_beta$beta = round(df_beta$beta,2)
  df_beta$beta_corrected = round(df_beta$beta_corrected,2)
  
  df_beta$X = rev(df_beta$X)
  fig2_E = ggplot(data = df_beta) + 
    #geom_rect(data = df_beta, 
    #          mapping = aes(xmin =xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
    #          fill = 'white') + 
    #geom_text(data = df_beta,
    #          aes(x = 0.3, y = ymin, label = beta),
    #          size = 1, vjust = -0.4, hjust = 0, check_overlap = TRUE) +
    geom_tile(aes(x = 1, y = X, fill = beta)) + 
    theme_bw() +
    scale_fill_gradient2(name = "loading", high = 'black', mid = 'white', low = 'black', 
                         midpoint = 0, breaks =breaks_vec, labels = labels_vec, ) +
    scale_x_continuous(name = "", expand = c(0,0)) + 
    scale_y_continuous(name = "", expand = c(0,0)) + 
    theme(axis.title = element_text(size=6))  + 
    main_theme_Fig2 + 
    theme(axis.ticks.x =element_blank()) +
    theme(axis.ticks.y =element_blank()) +
    theme(axis.text.x = element_blank()) + 
    theme(axis.text.y = element_blank()) + 
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())  + 
    geom_hline(yintercept = seq(0,22)+0.5, color = "black", size=0.05)
  
  
  df_beta$logP = rev(-log10(df_beta$pv))
  df_beta$factor = as.numeric(rownames(df_beta))
  
  fig2_F = ggplot(data =df_beta, aes(x = factor, y = logP)) + 
    geom_bar(stat = 'identity') + 
    xlab("Factors") + 
    ylab("-log10(BH-corrected p-value)") + 
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())  + 
    geom_hline(yintercept = -log10(0.05), color = "red", size=0.1) +  
    theme(axis.text.x = element_blank()) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.ticks.x = element_blank()) + 
    theme(axis.ticks.y = element_blank()) +
    main_theme_Fig2
  
  return(list(fig2_E, fig2_F))
}

plot_factor_matrix <- function(genei, snpi){
  factor_matrix = read.table('data/Fig1_factor_matrix.txt', sep='\t', stringsAsFactors = F, row.names = 1)
  factor_matrix = factor_matrix[rownames(factor_matrix) != "", ]
  rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
  factor_matrix = factor_matrix[seq(nrow(factor_matrix), 1), ]
  nF = dim(factor_matrix)[2]
  
  factor_matrix_col = (factor_matrix > 1e-5 ) *1
  for (row in rownames(factor_matrix)){
    factor_matrix_col[row, factor_matrix_col[row, ]==1] = gtex_col[row,"tissue_color_hex"]
    factor_matrix_col[row, factor_matrix_col[row, ]==0] = '#FFFFFF'
  }
  
  df_beta = read.table(paste0(local_dir, '/Fig3_example_pb_',genei,':', snpi, '.txt'), sep='\t', header=T, row.names = 1)
  factor_matrix_col[,df_beta$beta_corrected == 0] = lighten(factor_matrix_col[,df_beta$beta_corrected == 0], amount = 0.9)
  
  df = data.frame()
  for (i in seq(0, nT-1)){
    dfi = data.frame(x1 = seq(1,nF) , x2 = seq(2, nF+1), y1 = rep(i,nF), y2 = rep(i+1, nF), col = factor_matrix_col[i+1, ])
    df = rbind(df, dfi)
  }
  
  fig2_D = ggplot() +
    scale_x_continuous(name = 'Tissues', expand = c(0,0)) + 
    scale_y_continuous(name = 'Factors', expand = c(0,0)) + 
    coord_flip() + 
    geom_rect(data = df, mapping = aes(xmin =y1, xmax=y2, ymin=x1, ymax=x2), 
              fill = df$col,
              color = 'black', size = 0.01) + 
    xlab("Tissues") + 
    ylab("Factors")+
    main_theme_Fig2 + 
    theme_bw() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8)) + 
    theme(axis.ticks =element_blank()) +
    theme(axis.text = element_blank()) + 
    geom_hline(yintercept = seq(1,23), color = "black", size=0.1) 
  return(fig2_D)
}

one_genei_snpi <- function(genei, snpi, breaks_vec= NULL, labels_vec=NULL){
  Fig2C = plot_one_example_yw(genei, snpi)
  Fig2D = plot_factor_matrix(genei, snpi)
  Fig2E = plot_one_example_pb(genei, snpi, breaks_vec, labels_vec)[[1]] + 
    main_theme_Fig2 + 
    theme(legend.key.size = unit(0.05, "inches")) + 
    theme(legend.position="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-9,-9,-9,-9)) + 
    theme(axis.text.x = element_blank()) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.ticks.x = element_blank()) + 
    theme(axis.ticks.y = element_blank())
  
  return(list(Fig2C, Fig2D, Fig2E))  
}

local_dir = 'data/examples/'

# 1 TS
fn_name = 'Fig3_example_yw_ENSG00000151948.11:chr12_128918412.txt'
x = strsplit(gsub(".txt", "", gsub("Fig3_example_yw_", "", fn_name)), ":")[[1]]
pair1 = one_genei_snpi(x[1], x[2], c(0, 0.1), c("0", "0.1"))

# 1 Universal
fn_name = 'Fig3_example_yw_ENSG00000160201.11:chr21_43025075.txt'
## cis-eQTL in 6 tissues
x = strsplit(gsub(".txt", "", gsub("Fig3_example_yw_", "", fn_name)), ":")[[1]]
pair2 = one_genei_snpi(x[1], x[2], c(0, 0.03), c("0", "0.03"))

# 2 Ts
fn_name = 'Fig3_example_yw_ENSG00000068724.15:chr2_46965017.txt'
x = strsplit(gsub(".txt", "", gsub("Fig3_example_yw_", "", fn_name)), ":")[[1]]
pair3 = one_genei_snpi(x[1], x[2], c(-0.05, 0, 0.05), c("-0.05", " 0", " 0.05"))

# 1 Universal + 1 Ts
fn_name = 'Fig3_example_yw_ENSG00000170458.13:chr5_140721484.txt'
## cis-eQTL in 10 tissues including testis
x = strsplit(gsub(".txt", "", gsub("Fig3_example_yw_", "", fn_name)), ":")[[1]]
pair4 = one_genei_snpi(x[1], x[2], c(0, 0.04), c("0", "0.04"))



# brain
fn_name = 'Fig3_example_yw_ENSG00000125337.16:chr6_168782339.txt'
x = strsplit(gsub(".txt", "", gsub("Fig3_example_yw_", "", fn_name)), ":")[[1]]
pair5 = one_genei_snpi(x[1], x[2], c(-0.04, 0, 0.04 ), c("-0.04", " 0", "0.04"))


# brain + 1 Ts
fn_name = 'Fig3_example_yw_ENSG00000275700.4:chr17_37010841.txt'
## ciseqTL in: 4 brain tissues, nerve tibial
x = strsplit(gsub(".txt", "", gsub("Fig3_example_yw_", "", fn_name)), ":")[[1]]
pair6 = one_genei_snpi(x[1], x[2], c(0, 0.05 ), c("0", "0.05"))


fig2_part1 = ggdraw() + 
  draw_plot(pair1[[1]] + 
              theme(panel.border = element_blank(),
                    axis.line.x = element_line(size=0.5),
                    axis.text.x = element_text(size=5)) + 
              scale_y_continuous(breaks=c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) + 
              xlab("") + ylab(""),
            x = 0, y = 0.072, width = 0.15, height = 0.9) +
  draw_plot(pair1[[2]] + 
              theme(panel.border = element_rect(size=0.25, color = 'black'),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()), 
            x = 0.21, y = 0.15, width = 0.17, height = 0.82) + 
  draw_plot(pair1[[3]] + 
              theme(legend.text = element_text(size=4),
                    legend.title = element_text(size=5)), 
            x = 0.21 + 0.16, y = 0.33, width = 0.08, height = 0.38) + 
  
  draw_plot(pair6[[1]] + 
              theme(panel.border = element_blank(),
                    axis.line.x = element_line(size=0.5),
                    axis.text.x = element_text(size=5)) + 
              scale_y_continuous(breaks=c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) + 
              xlab("") + ylab(""),
            x = 0.5, y = 0.072, width = 0.15, height = 0.9) +
  draw_plot(pair6[[2]] + 
              theme(panel.border = element_rect(size=0.25, color = 'black'),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()), 
            x = 0.5 + 0.21, y = 0.15, width = 0.17, height = 0.82) + 
  draw_plot(pair6[[3]] + 
              theme(legend.text = element_text(size=4),
                    legend.title = element_text(size=5)), 
            x = 0.5 + 0.21 + 0.16, y = 0.33, width = 0.08, height = 0.38)
  #draw_plot_label(label = c("A", "B"), x = c(0.05,0.52), y =c(1,1), size = 8)
save_plot("results_refined/Fig2_part1_1_revision_II.png", fig2_part1, base_width = 7.4)



fig2_part2 = ggdraw() + 
  draw_plot(pair2[[1]] + 
              theme(panel.border = element_blank(),
                    axis.line.x = element_line(size=0.5),
                    axis.text.x = element_text(size=5)) + 
              scale_y_continuous(breaks=c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) + 
              xlab("") + ylab(""),
            x = 0, y = 0.072, width = 0.15, height = 0.9) +
  draw_plot(pair2[[2]] + 
              theme(panel.border = element_rect(size=0.25, color = 'black'),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()), 
            x = 0.21, y = 0.145, width = 0.17, height = 0.83) + 
  draw_plot(pair2[[3]] + 
              theme(legend.text = element_text(size=4),
                    legend.title = element_text(size=5)), 
            x = 0.21 + 0.16, y = 0.33, width = 0.08, height = 0.38) + 
  
  draw_plot(pair4[[1]] + 
              theme(panel.border = element_blank(),
                    axis.line.x = element_line(size=0.5),
                    axis.text.x = element_text(size=5)) + 
              scale_y_continuous(breaks=c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5")) + 
              xlab("") + ylab(""),
            x = 0.5, y = 0.072, width = 0.15, height = 0.9) +
  draw_plot(pair4[[2]] + 
              theme(panel.border = element_rect(size=0.25, color = 'black'),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank()), 
            x = 0.5 + 0.21, y = 0.145, width = 0.17, height = 0.83) + 
  draw_plot(pair4[[3]] + 
              theme(legend.text = element_text(size=4),
                    legend.title = element_text(size=5)), 
            x = 0.5 + 0.21 + 0.16, y = 0.33, width = 0.08, height = 0.38)
save_plot("results_refined/Fig2_part2_1_revision_II.png", fig2_part2, base_width = 7.4)



fig2_part1_3 = ggdraw() + 
  draw_plot(pair5[[1]], x = 0, y = 0.072, width = 0.23, height = 0.9) + 
  draw_plot(pair5[[2]], x = 0.20, y = 0.19, width = 0.18, height = 0.715) + 
  draw_plot(pair5[[3]], x = 0.38, y = 0.32, width = 0.08, height = 0.4) + 
  
  draw_plot(pair6[[1]], x = 0.45, y = 0.073, width = 0.23, height = 0.9) + 
  draw_plot(pair6[[2]], x = 0.20 + 0.45, y = 0.195, width = 0.2, height = 0.705) + 
  draw_plot(pair6[[3]], x = 0.41 + 0.45, y = 0.32, width = 0.08, height = 0.4) + 
  
  draw_plot_label(label = c("E", "F"), x = c(0.05,0.5), y =c(0.95,0.95), size = 8)

save_plot("results_refined/Fig2_part1_3_revision.png", fig2_part1_3, base_width = 7.4)





fn_name = 'Fig3_example_yw_ENSG00000203880.11:chr20_64239294.txt'
## ciseqTL in: 4 brain tissues, nerve tibial
x = strsplit(gsub(".txt", "", gsub("Fig3_example_yw_", "", fn_name)), ":")[[1]]
pair7 = one_genei_snpi(x[1], x[2], c(0, 0.05 ), c("0", "0.05"))
pair7_corrected = one_genei_snpi(x[1], x[2], c(0, 0.05 ), c("0", "0.05"))


fig2_flip = ggdraw() + 
  draw_plot(pair7[[1]], x = 0, y = 0.072, width = 0.23, height = 0.9) + 
  draw_plot(pair7[[2]], x = 0.20, y = 0.19, width = 0.18, height = 0.715) + 
  draw_plot(pair7[[3]], x = 0.38, y = 0.32, width = 0.08, height = 0.4) + 
  
  draw_plot(pair7_corrected[[1]], x = 0.45, y = 0.073, width = 0.23, height = 0.9) + 
  draw_plot(pair7_corrected[[2]], x = 0.20 + 0.45, y = 0.195, width = 0.2, height = 0.705) + 
  draw_plot(pair7_corrected[[3]], x = 0.41 + 0.45, y = 0.32, width = 0.08, height = 0.4)

save_plot("Supplmentary_plots/Flipped_sign_revision_II.png", fig2_flip, base_width = 7.4)

