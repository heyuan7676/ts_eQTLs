source('scripts_refined/utils.R')

pointSize = 1
textSize = 3
spaceLegend = 0.2

main_theme_Fig5 = theme_bw() + 
  background_grid(major = 'none', minor = 'none') + 
  theme(axis.title.x = element_text(size = 8,
                                    color = 'black'),
        axis.title.y = element_text(size = 8, color = 'black'),
        strip.text = element_text(size=4, color = 'black'),
        legend.title = element_text(size = 6, color = 'black'), 
        legend.text  = element_text(size = 5, color = 'black'),
        axis.text = element_text(size = 5, color = 'black'),
        legend.position = 'none',
        plot.title = element_text(face = "bold", size=8, color = 'black', hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(color = 'black'), 
        axis.line.y = element_line(color = 'black'))


## read in
readin_TF_result <- function(method, factor_names, shared_text = 0, add = 1, TPM = 0, HIT = 10, HIT_upper = 500){
  TssA_dat_sp = read.table(paste0('data/Fig5_TF_enrichment_',method, '_ROADMAP_1_TssA_iter2_noRes','.txt'), sep='\t', header=T, stringsAsFactors = F)
  Enh_dat_sp  = read.table(paste0('data/Fig5_TF_enrichment_',method, '_ROADMAP_7_Enh_iter2_noRes','.txt'), sep='\t', header=T, stringsAsFactors = F)
  TssA_dat_sp$region = 'TssA'
  Enh_dat_sp$region = 'Enh'
  thr_pv = rbind(TssA_dat_sp, Enh_dat_sp)

  thr_pv$background_nonTFBS = floor(thr_pv$background_nonTFBS/2) 
  thr_pv$background_TFBS = floor(thr_pv$background_TFBS /2)
  thr_pv = thr_pv[floor(thr_pv$background_TFBS) + thr_pv$ts_TFBS > HIT,]
  thr_pv = thr_pv[floor(thr_pv$background_TFBS) + thr_pv$ts_TFBS < HIT_upper,]
  
  remove_TFs = c('MEF2A', 'MEF2B', 'MEF2C','MEF2D',
                 'NFAT5', 'ZNF384', 'FOXK1', 'IRF1',
                 'FOXJ3', 'SRF', 'FOXD3',
                 'POU2F2', 'POU2F3', 'LIN54', 'E2F1')
  thr_pv = thr_pv[!thr_pv$TF %in% remove_TFs, ]
  
  thr_pv$TF_TPM = as.numeric(thr_pv$TF_TPM)
  thr_pv = thr_pv[thr_pv$TF_TPM > TPM, ]
  thr_pv = add_count(thr_pv, factor_names)
  
  if(!is.null(factor_names)){
    if(is.numeric(thr_pv$group[1])){
      thr_pv$Factor_name = factor_names[as.numeric(thr_pv$group)+add]  
    }else{
      thr_pv$Factor_name = thr_pv$group
      thr_pv[thr_pv$Factor_name == 'Shared', "Factor_name"] = "Ubiquitous"
    }
  }else{
    thr_pv$Factor_name = thr_pv$group
  }
  
  thr_pv$BH_pvalue = p.adjust(thr_pv$PV_background, method = 'BH')
  TssA_dat_sp = thr_pv[thr_pv$region == 'TssA', ]
  Enh_dat_sp  = thr_pv[thr_pv$region == 'Enh',]
  
  TssA_dat_sp_sig = TssA_dat_sp[TssA_dat_sp$BH_pvalue < alpha, ]
  Enh_dat_sp_sig = Enh_dat_sp[Enh_dat_sp$BH_pvalue < alpha, ]
  
  print(paste0("TssA: ", 'Shared-', length(unique(TssA_dat_sp_sig[TssA_dat_sp_sig$group == shared_text, "TF"])), 
          ' TS-', length(unique(TssA_dat_sp_sig[TssA_dat_sp_sig$group != shared_text, "TF"])),
          ' Total-', length(unique(TssA_dat_sp_sig$TF))))
  
  print(paste0("Enh: ", 'Shared-', length(unique(Enh_dat_sp_sig[Enh_dat_sp_sig$group == shared_text, "TF"])), 
          ' TS-', length(unique(Enh_dat_sp_sig[Enh_dat_sp_sig$group != shared_text, "TF"])),
          ' Total-', length(unique(Enh_dat_sp_sig$TF))))
  
  return(thr_pv)
}

add_count <- function(df, factor_names){
  TF_count = df %>% group_by(TF) %>% dplyr::summarise(count = n())
  TF_count = as.data.frame(TF_count)
  rownames(TF_count) = TF_count$TF
  
  sig_df = df[df$BH_pvalue < alpha, ]
  sig_TF_count = sig_df %>% group_by(TF) %>% dplyr::summarise(count = n())
  sig_TF_count = as.data.frame(sig_TF_count)
  rownames(sig_TF_count) = sig_TF_count$TF
  
  df$count = sig_TF_count[df$TF, "count"]
  df[is.na(df$count), "count"] = 0
  df$total_count = TF_count[df$TF, "count"]
  
  return(df)
}


TPM_threshold = 1
HIT_threshold = 10
alpha = 0.05


sn_spMF = readin_TF_result('sn_spMF', factor_names, TPM = TPM_threshold, HIT = HIT_threshold)

flashr = readin_TF_result('flashr', factor_names_flashr, TPM = TPM_threshold, HIT = HIT_threshold)
flashr_bf = readin_TF_result('flashr_bf', factor_names_flashr, TPM = TPM_threshold, HIT = HIT_threshold)
flashr_NN = readin_TF_result('flashr_NN', factor_names_flashr_NN, TPM = TPM_threshold, HIT = HIT_threshold)

si = readin_TF_result('softImpute', factor_seq_names, TPM = TPM_threshold, HIT = HIT_threshold)
pmd = readin_TF_result('PMD', factor_seq_names, TPM = TPM_threshold, HIT = HIT_threshold)
pmd2 = readin_TF_result('PMD_cv2', factor_seq_names, TPM = TPM_threshold, HIT = HIT_threshold)
heuristic_1 = readin_TF_result('heuristic_1', shared_text = 0, 
                               c("Ubiquitous", tissues), TPM = TPM_threshold, HIT = HIT_threshold)
heuristic_2 = readin_TF_result('heuristic_2', shared_text = 0, 
                               factor_names_thresholding, TPM = TPM_threshold, HIT = HIT_threshold)







## hits per factor
hit_per_factor_to_pv <- function(df, factor_names_use){
  TF_count = df %>% group_by(TF) %>% dplyr::summarise(count = n())
  TF_count = as.data.frame(TF_count)
  rownames(TF_count) = TF_count$TF
  
  # get the row per top pv
  tissue_top_pv = df %>% group_by(TF) %>% dplyr::summarise(group = group[which.min(pvalue)])
  tissue_top_pv = as.data.frame(tissue_top_pv)
  tissue_top_pv$count = TF_count[tissue_top_pv$TF, "count"]
  rownames(tissue_top_pv) = tissue_top_pv$TF
  top_count = as.data.frame(tissue_top_pv %>% group_by(group) %>% dplyr::summarise(count=n())) 
  
  # count numbers
  factor_hits = as.data.frame(df %>% group_by(group) %>% dplyr::summarise(count=n()))
  factor_hits = merge(factor_hits, top_count, by.x = 'group', by.y = 'group', all = T)
  colnames(factor_hits) = c("group", "count", "top_count")
  
  # format
  if(is.numeric(factor_hits$group[1])){
    factor_hits$group = factor_names_use[as.numeric(factor_hits$group)+1] 
  }
  factor_hits$group = factor(factor_hits$group , levels = factor_hits$group )
  factor_hits[is.na(factor_hits$top_count), "top_count"] = 0
  factor_hits$not_top_count = factor_hits$count - factor_hits$top_count
  factor_hits$count <- NULL
  
  factor_hits = melt(factor_hits)
  factor_hits[is.na(factor_hits$value), "value"] = 0
  
  return(factor_hits)
}


### compare TF for shared and ts in enhancer / promoter
numberTFs_shared_ts_eQTLs <- function(TssA_dat, Enh_dat, shared_text = 0, unique_TF = TRUE){
  TssA_dat_sig = TssA_dat[TssA_dat$BH_pvalue < 0.05, ]
  Enh_dat_sig  = Enh_dat[Enh_dat$BH_pvalue < 0.05, ]
  if(unique_TF){
    A = unique(TssA_dat_sig[TssA_dat_sig$group == shared_text, "TF"])
    B = unique(TssA_dat_sig[TssA_dat_sig$group != shared_text, "TF"])
    AB = intersect(A,B)
    C = length(unique(TssA_dat$TF)) - length(unique(TssA_dat_sig$TF))
    
    D = unique(Enh_dat_sig[Enh_dat_sig$group == shared_text, "TF"])
    E = unique(Enh_dat_sig[Enh_dat_sig$group != shared_text, "TF"])
    DE = intersect(D,E)
    G = length(unique(Enh_dat$TF)) - length(unique(Enh_dat_sig$TF))
  }else{
    A = (TssA_dat_sig[TssA_dat_sig$group == shared_text, "TF"])
    B = (TssA_dat_sig[TssA_dat_sig$group != shared_text, "TF"])
    AB = intersect(A,B)
    C = length((TssA_dat$TF)) - length((TssA_dat_sig$TF))
    
    D = (Enh_dat_sig[Enh_dat_sig$group == shared_text, "TF"])
    E = (Enh_dat_sig[Enh_dat_sig$group != shared_text, "TF"])
    DE = intersect(D,E)
    G = length((Enh_dat$TF)) - length((Enh_dat_sig$TF))
  }
  
  AB = length(AB)
  A = length(A) - AB
  B = length(B) - AB
  
  DE = length(DE)
  D = length(D) - DE
  E = length(E) - DE
  
  compare_TssA_Enh_dat = matrix(c(c(A,B,AB,C, D,E,DE,G), c("Shared", "ts", "both","notSig","Shared", "ts", "both","notSig"), c("TssA", "TssA", "TssA","TssA","Enh", "Enh", "Enh", "Enh")), ncol = 3)
  compare_TssA_Enh_dat = as.data.frame(compare_TssA_Enh_dat)
  colnames(compare_TssA_Enh_dat) = c("value", "eQTL", "chrom_status")
  compare_TssA_Enh_dat$value = as.numeric(as.character(compare_TssA_Enh_dat$value))
  compare_TssA_Enh_dat$eQTL = factor(as.character(compare_TssA_Enh_dat$eQTL), levels = c("both", "Shared", "ts", "notSig"))
  
  Fig_summary_numbers = ggplot(data = compare_TssA_Enh_dat, aes(x=chrom_status, y=value, fill=eQTL)) + 
    geom_bar(stat = 'identity') + 
    xlab("")+
    ylab("Number of TFs")+
    scale_fill_manual(values=c( "purple", shared_factor_color, ts_factor_color, 'grey'), 
                      name="Enriched in",
                      labels=c("Both","Only U-eQTLs", "Only TS-eQTLs", "Not enriched")) +
    guides(fill= guide_legend(title.position = 'top',
                              nrow=4,byrow=TRUE)) + 
    theme(legend.position = 'right',
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10)) + 
    main_theme_Fig5
  
  return(Fig_summary_numbers)
}

### TFs each factor
numberTFs_each_factor <- function(TssA_dat_sig, Enh_dat_sig, factor_names_use){
  
  TssA_top_count = hit_per_factor_to_pv(TssA_dat_sig, factor_names_use)
  TssA_top_count$annotation = 'TssA'
  TssA_top_count[TssA_top_count$group == 'Shared', "group"] = 'Ubiquitous'
  Enh_top_count = hit_per_factor_to_pv(Enh_dat_sig, factor_names_use)
  Enh_top_count$annotation = 'Enh'
  Enh_top_count[Enh_top_count$group == 'Shared', "group"] = 'Ubiquitous'
  
  bb_plot_df = rbind(TssA_top_count, Enh_top_count)
  bb_plot_df = ddply(bb_plot_df, c("group", "annotation"), function(x) colSums(x["value"]))
  bb_plot_df$annotation = factor(bb_plot_df$annotation, levels = c("TssA", "Enh"))
  
  fig_TFs = ggplot(bb_plot_df, aes(group)) + 
    geom_bar(data=bb_plot_df[bb_plot_df$annotation=='TssA',], aes(y=value, fill=annotation), stat = 'identity') + 
    geom_bar(data=bb_plot_df[bb_plot_df$annotation=='Enh',], aes(y=-value, fill=annotation), stat = 'identity') + 
    scale_fill_manual(name = "Regulatory region", 
                      values=c('TssA' ='red', 'Enh' = 'yellow'), 
                      labels = c("Enhancer", "Promoter")) + 
    xlab("") + 
    geom_hline(yintercept = 0,colour = "black") + 
    #theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
    scale_y_continuous("Number of enriched TFs",
                       breaks = c(-100, 0, 50), labels = c(100, 0, 50)) + 
    coord_flip() + 
    theme(legend.position = 'top') +
    main_theme_Fig5
  
  return(fig_TFs)
}










### Number of factors each TF is enriched in
number_factors_enriched_Tfs <- function(dat_sig, shared_text){
  dat_sig_shared = dat_sig[dat_sig$group == shared_text, ]
  dat_sig_shared = unique(dat_sig_shared[,c("TF", "count")])
  dat_sig_ts = dat_sig[dat_sig$group != shared_text, ]
  dat_sig_ts = dat_sig_ts[!dat_sig_ts$TF %in% dat_sig_shared$TF, ]
  dat_sig_ts = unique(dat_sig_ts[,c("TF", "count")])
  
  dat_sig_shared_count = as.data.frame(table(dat_sig_shared$count - 1))
  colnames(dat_sig_shared_count) = c("Number_of_factors", "Freq")
  dat_sig_shared_count$group = "Ubiquitous"
  
  dat_sig_ts_count = as.data.frame(table(dat_sig_ts$count))
  colnames(dat_sig_ts_count) = c("Number_of_factors", "Freq")
  dat_sig_ts_count$group = 'Tissue-specific'
  
  sig_count = rbind(dat_sig_shared_count,dat_sig_ts_count)
  sig_count$Number_of_factors = factor(as.character(sig_count$Number_of_factors), 
                                       levels = as.character(sort(unique(as.numeric(as.character(sig_count$Number_of_factors))))))
  sig_count$group = factor(sig_count$group, levels = c("Ubiquitous", "Tissue-specific"))
  
  fig_nF = ggplot(data = sig_count, aes(x = Number_of_factors, y = Freq, fill = group))+
    geom_bar(stat = 'identity') + 
    xlab("Number of TS-factors") + 
    ylab("Number of TFs")+
    theme(legend.position = c(0, 1), 
          legend.justification = c(0, 0),
          legend.direction = "horizontal") + 
    theme(legend.key.size = unit(0.03, "inches"))+
    theme(legend.position = 'top') + 
    scale_fill_manual(name = "Enriched in", 
                      values=c(shared_factor_color,ts_factor_color),
                      labels = c("U-eQTLs", "Only TS-eQTLs")) + 
    main_theme_Fig5
  return(fig_nF)
}


### plot individual TF
plot_bar_pvalue <- function(TF, df){
  tp = df[df$TF == TF, ]
  group_levels = tp$Factor_name[order(tp$pvalue, decreasing = F)]
  tp$Factor_name = factor(tp$Factor_name, levels = group_levels)
  tp$logP = -log10(tp$pvalue)
  tp$sig = factor(ifelse(tp$BH_pvalue < alpha, "sig", "unsig"))
  
  g = ggplot(data = tp, aes(x = Factor_name, y = logP, fill = sig))  + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    main_theme_Fig5 + 
    xlab("") + 
    ylab(paste0(expression('-log'[10])," (P value)")) + 
    theme(axis.text.x =element_text(angle = 45, hjust = 1, color = 'black'),
          axis.text.y = element_text(color = 'black'))  + 
    scale_fill_manual(values = c("sig" = "grey10", "unsig" = "lightgrey")) + 
    labs(title = TF)  +
    theme(legend.position = 'none') + 
    theme(axis.line.x = element_line(size = 0.4), 
          axis.line.y = element_line(size = 0.4))
  
  return(g)
}


TPM_threshold = 1
HIT_threshold = 10
alpha = 0.05

numberTFs_shared_ts_eQTLs <- function(TssA_dat, Enh_dat, shared_text = 0, unique_TF = TRUE){
  TssA_dat_sig = TssA_dat[TssA_dat$BH_pvalue < 0.05, ]
  Enh_dat_sig  = Enh_dat[Enh_dat$BH_pvalue < 0.05, ]
  if(unique_TF){
    A = unique(TssA_dat_sig[TssA_dat_sig$group == shared_text, "TF"])
    B = unique(TssA_dat_sig[TssA_dat_sig$group != shared_text, "TF"])
    AB = intersect(A,B)
    C = length(unique(TssA_dat$TF)) - length(unique(TssA_dat_sig$TF))
    
    D = unique(Enh_dat_sig[Enh_dat_sig$group == shared_text, "TF"])
    E = unique(Enh_dat_sig[Enh_dat_sig$group != shared_text, "TF"])
    DE = intersect(D,E)
    G = length(unique(Enh_dat$TF)) - length(unique(Enh_dat_sig$TF))
  }else{
    A = (TssA_dat_sig[TssA_dat_sig$group == shared_text, "TF"])
    B = (TssA_dat_sig[TssA_dat_sig$group != shared_text, "TF"])
    AB = intersect(A,B)
    C = length((TssA_dat$TF)) - length((TssA_dat_sig$TF))
    
    D = (Enh_dat_sig[Enh_dat_sig$group == shared_text, "TF"])
    E = (Enh_dat_sig[Enh_dat_sig$group != shared_text, "TF"])
    DE = intersect(D,E)
    G = length((Enh_dat$TF)) - length((Enh_dat_sig$TF))
  }
  
  AB = length(AB)
  A = length(A) - AB
  B = length(B) - AB
  
  DE = length(DE)
  D = length(D) - DE
  E = length(E) - DE
  
  compare_TssA_Enh_dat = matrix(c(c(A,B,AB,C, D,E,DE,G), c("Shared", "ts", "both","notSig","Shared", "ts", "both","notSig"), c("TssA", "TssA", "TssA","TssA","Enh", "Enh", "Enh", "Enh")), ncol = 3)
  compare_TssA_Enh_dat = as.data.frame(compare_TssA_Enh_dat)
  colnames(compare_TssA_Enh_dat) = c("value", "eQTL", "chrom_status")
  compare_TssA_Enh_dat$value = as.numeric(as.character(compare_TssA_Enh_dat$value))
  compare_TssA_Enh_dat$eQTL = factor(as.character(compare_TssA_Enh_dat$eQTL), levels = c("both", "Shared", "ts", "notSig"))
  
  Fig_summary_numbers = ggplot(data = compare_TssA_Enh_dat, aes(x=chrom_status, y=value, fill=eQTL)) + 
    geom_bar(stat = 'identity') + 
    xlab("")+
    ylab("Number of TFs")+
    scale_fill_manual(values=c( "purple", shared_factor_color, ts_factor_color, 'grey'), 
                      name="Enriched in",
                      labels=c("Both","Only U-eQTLs", "Only TS-eQTLs", "Not enriched")) +
    guides(fill= guide_legend(title.position = 'top',
                              nrow=4,byrow=TRUE)) + 
    theme(legend.position = 'right',
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10)) + 
    main_theme_Fig5
  
  return(Fig_summary_numbers)
}

use_unique_TFs = T
yh = 420
Fig5A_1 = numberTFs_shared_ts_eQTLs(sn_spMF[sn_spMF$region == 'TssA', ], sn_spMF[sn_spMF$region == 'Enh', ], unique_TF = use_unique_TFs) + ylim(0, yh)
Fig5A_2 = numberTFs_shared_ts_eQTLs(flashr_bf[flashr_bf$region == 'TssA', ], flashr_bf[flashr_bf$region == 'Enh', ], unique_TF = use_unique_TFs) + ylim(0, yh)
Fig5A_3 = numberTFs_shared_ts_eQTLs(heuristic_1[heuristic_1$region == 'TssA', ], heuristic_1[heuristic_1$region == 'Enh', ], unique_TF = use_unique_TFs) + 
  ylim(0, yh) + 
  theme(legend.position =  'right', 
        legend.key.size = unit(0.1, "inches")) + 
  theme(axis.line.x = element_line(), 
        axis.line.y = element_line())


### sn-spMF in more details
numberTFs_each_factor <- function(TssA_dat_sig, Enh_dat_sig, factor_names_use){
  TssA_top_count = hit_per_factor_to_pv(TssA_dat_sig, factor_names_use)
  TssA_top_count$annotation = 'TssA'
  TssA_top_count[TssA_top_count$group == 'Shared', "group"] = 'Ubiquitous'
  Enh_top_count = hit_per_factor_to_pv(Enh_dat_sig, factor_names_use)
  Enh_top_count$annotation = 'Enh'
  Enh_top_count[Enh_top_count$group == 'Shared', "group"] = 'Ubiquitous'
  
  bb_plot_df = rbind(TssA_top_count, Enh_top_count)
  bb_plot_df = ddply(bb_plot_df, c("group", "annotation"), function(x) colSums(x["value"]))
  bb_plot_df$annotation = factor(bb_plot_df$annotation, levels = c("TssA", "Enh"))
  
  fig_TFs = ggplot(bb_plot_df, aes(group)) + 
    geom_bar(data=bb_plot_df[bb_plot_df$annotation=='TssA',], aes(y=value, fill=annotation), stat = 'identity') + 
    geom_bar(data=bb_plot_df[bb_plot_df$annotation=='Enh',], aes(y=-value, fill=annotation), stat = 'identity') + 
    scale_fill_manual(name = "Regulatory region", 
                      values=c('TssA' ='red', 'Enh' = 'yellow'), 
                      labels = c("Enhancer", "Promoter")) + 
    xlab("") + 
    geom_hline(yintercept = 0,colour = "black") + 
    #theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
    scale_y_continuous("Number of enriched TFs",
                       breaks = c(-100, 0, 50), labels = c(100, 0, 50)) + 
    coord_flip() + 
    theme(legend.position = 'top') +
    main_theme_Fig5
  
  return(fig_TFs)
}
## hits per factor
hit_per_factor_to_pv <- function(df, factor_names_use){
  TF_count = df %>% group_by(TF) %>% dplyr::summarise(count = n())
  TF_count = as.data.frame(TF_count)
  rownames(TF_count) = TF_count$TF
  
  # get the row per top pv
  tissue_top_pv = df %>% group_by(TF) %>% dplyr::summarise(group = group[which.min(PV_background)])
  tissue_top_pv = as.data.frame(tissue_top_pv)
  tissue_top_pv$count = TF_count[tissue_top_pv$TF, "count"]
  rownames(tissue_top_pv) = tissue_top_pv$TF
  top_count = as.data.frame(tissue_top_pv %>% group_by(group) %>% dplyr::summarise(count=n())) 
  
  # count numbers
  factor_hits = as.data.frame(df %>% group_by(group) %>% dplyr::summarise(count=n()))
  factor_hits = merge(factor_hits, top_count, by.x = 'group', by.y = 'group', all = T)
  colnames(factor_hits) = c("group", "count", "top_count")
  
  # format
  if(is.numeric(factor_hits$group[1])){
    factor_hits$group = factor_names_use[as.numeric(factor_hits$group)+1] 
  }
  factor_hits$group = factor(factor_hits$group , levels = factor_hits$group )
  factor_hits[is.na(factor_hits$top_count), "top_count"] = 0
  factor_hits$not_top_count = factor_hits$count - factor_hits$top_count
  factor_hits$count <- NULL
  
  factor_hits = melt(factor_hits)
  factor_hits[is.na(factor_hits$value), "value"] = 0
  
  return(factor_hits)
}
sn_spMF_sig = sn_spMF[sn_spMF$BH_pvalue < 0.05, ]
Fig5B = numberTFs_each_factor(sn_spMF_sig[sn_spMF_sig$region == 'TssA', ], 
                              sn_spMF_sig[sn_spMF_sig$region == 'Enh', ], 
                              factor_names)


flashr_sig = flashr[flashr$BH_pvalue < 0.05, ]
TssA_dat_flashr_sig = flashr_sig[flashr_sig$region == 'TssA', ]
Enh_dat_flashr_sig = flashr_sig[flashr_sig$region == 'Enh', ]
Fig5B_flashr = numberTFs_each_factor(TssA_dat_flashr_sig, Enh_dat_flashr_sig, factor_seq_names)
Fig5B_thr = numberTFs_each_factor(TssA_dat_thr_sig, Enh_dat_thr_sig, tissues)
Fig5B_si = numberTFs_each_factor(TssA_dat_si_sig, Enh_dat_si_sig,  factor_seq_names)
save_plot("results_refined/Fig4_flashr_details_revision.png", Fig5B_flashr, base_width = 7.4)
save_plot("results_refined/Fig4_si_details_revision.png", Fig5B_si, base_width = 7.4)



### Number of factors each TF is enriched in
number_factors_enriched_Tfs <- function(dat_sig, shared_text){
  dat_sig_shared = dat_sig[dat_sig$group == shared_text, ]
  dat_sig_shared = unique(dat_sig_shared[,c("TF", "count")])
  dat_sig_ts = dat_sig[dat_sig$group != shared_text, ]
  dat_sig_ts = dat_sig_ts[!dat_sig_ts$TF %in% dat_sig_shared$TF, ]
  dat_sig_ts = unique(dat_sig_ts[,c("TF", "count")])
  
  dat_sig_shared_count = as.data.frame(table(dat_sig_shared$count - 1))
  colnames(dat_sig_shared_count) = c("Number_of_factors", "Freq")
  dat_sig_shared_count$group = "Ubiquitous"
  
  dat_sig_ts_count = as.data.frame(table(dat_sig_ts$count))
  colnames(dat_sig_ts_count) = c("Number_of_factors", "Freq")
  dat_sig_ts_count$group = 'Tissue-specific'
  
  sig_count = rbind(dat_sig_shared_count,dat_sig_ts_count)
  sig_count$Number_of_factors = factor(as.character(sig_count$Number_of_factors), 
                                       levels = as.character(sort(unique(as.numeric(as.character(sig_count$Number_of_factors))))))
  sig_count$group = factor(sig_count$group, levels = c("Ubiquitous", "Tissue-specific"))
  
  fig_nF = ggplot(data = sig_count, aes(x = Number_of_factors, y = Freq, fill = group))+
    geom_bar(stat = 'identity') + 
    xlab("Number of TS-factors") + 
    ylab("Number of TFs")+
    theme(legend.position = c(0, 1), 
          legend.justification = c(0, 0),
          legend.direction = "horizontal") + 
    theme(legend.key.size = unit(0.03, "inches"))+
    theme(legend.position = 'top') + 
    scale_fill_manual(name = "Enriched in", 
                      values=c(shared_factor_color,ts_factor_color),
                      labels = c("U-eQTLs", "Only TS-eQTLs")) + 
    main_theme_Fig5
  return(fig_nF)
}
Fig5C = number_factors_enriched_Tfs(sn_spMF_sig[sn_spMF_sig$region == 'Enh', ], shared_text = 0)

Enh_dat_sp = sn_spMF[sn_spMF$region == 'Enh', ]
Enh_dat_sp$Factor_name = factor_names[Enh_dat_sp$group + 1]
plot_bar_pvalue <- function(TF, df){
  tp = df[df$TF == TF, ]
  group_levels = tp$Factor_name[order(tp$PV_background, decreasing = F)]
  tp$Factor_name = factor(tp$Factor_name, levels = group_levels)
  tp$logP = -log10(tp$PV_background)
  tp$sig = factor(ifelse(tp$BH_pvalue < alpha, "sig", "unsig"))
  
  g = ggplot(data = tp, aes(x = Factor_name, y = logP, fill = sig))  + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    main_theme_Fig5 + 
    xlab("") + 
    ylab(bquote('-'*log[10]*' (P value)')) + 
    #ylab(bquote('-log ('[10]')')) + 
    theme(axis.text.x =element_text(angle = 45, hjust = 1, color = 'black'),
          axis.text.y = element_text(color = 'black'))  + 
    scale_fill_manual(values = c("sig" = "grey10", "unsig" = "lightgrey")) + 
    labs(title = TF)  +
    theme(legend.position = 'none') + 
    theme(axis.line.x = element_line(size = 0.4), 
          axis.line.y = element_line(size = 0.4))
  
  return(g)
}
TF_list = c("GATA4", "HNF4A", "FOSL2")
TFs_pvalue_plots = list()
for(tfi in TF_list){
  TFs_pvalue_plots[[tfi]] = plot_bar_pvalue(tfi, Enh_dat_sp)
}


vn1=expression(sn-spMF)
vn2=expression(flashr[bf])
vn3=expression(heuristic[1])


fig4_g = ggdraw() + 
  draw_plot(Fig5A_1 + ggtitle(vn1), x = 0, y = 0.75, width = 0.28, height = 0.22) +
  draw_plot(Fig5A_2 + ggtitle(vn2) + ylab(""), x = 0.28, y = 0.75, width = 0.28, height = 0.22) +
  draw_plot(Fig5A_3 + ggtitle(vn3) + ylab(""), x = 0.55, y = 0.75, width = 0.455, height = 0.22) +
  draw_plot(Fig5B + 
              theme(legend.position = 'top',
                    legend.justification="left", 
                    legend.key.size = unit(0.1, "inches"),
                    legend.margin=margin(0,0,0,1),
                    legend.box.margin=margin(-10,-10,-5,0)) + 
              theme(axis.line.x = element_line(), 
                    axis.line.y = element_line()), 
            x = 0, y = 0.53, width = 0.95, height = 0.2) +
  draw_plot(Fig5C + 
              theme(legend.position = 'top',
                    legend.justification="left", 
                    legend.key.size = unit(0.1, "inches"), 
                    legend.margin=margin(0,0,0,1),
                    legend.box.margin=margin(-10,-10,-10,0)) + 
              theme(axis.line.x = element_line(), 
                    axis.line.y = element_line()), x = 0.009, y = 0.32, width = 0.5, height = 0.18) +
  draw_plot(TFs_pvalue_plots[["FOSL2"]] + ggtitle(expression(FOSL2)), x = 0.5, y = 0.262, width = 0.45, height = 0.26)  + 
  draw_plot(TFs_pvalue_plots[["GATA4"]] + ggtitle(expression(GATA4)), x = 0, y = 0.057, width = 0.5, height = 0.22) +
  draw_plot(TFs_pvalue_plots[["HNF4A"]] + ggtitle(expression(HNF4A)), x = 0.5, y = 0.04, width = 0.45, height = 0.23)

#draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), x = c(0,0,0, 0.5, 0, 0.5) + 0.035, 
#y =c(1,0.75,0.55,0.55,0.33,0.33) - 0.02, size = 8)

save_plot("results_refined/Fig4_v6_revision_2.png", fig4_g, base_width = 5, base_height = 10)






fig4_g = ggdraw() + 
  draw_plot(Fig5A_1 + ggtitle('sn-spMF'), x = 0, y = 0.78, width = 0.23, height = 0.2) +
  draw_plot(Fig5A_2 + ggtitle('flashr_bf') + ylab(""), x = 0.22, y = 0.78, width = 0.23, height = 0.2) +
  draw_plot(Fig5A_3 + ggtitle('heuristic') + ylab(""), x = 0.44, y = 0.78, width = 0.23, height = 0.2) +
  draw_plot(Fig5B, x = 0, y = 0.56, width = 1, height = 0.21) +
  draw_plot(Fig5C, x = 0, y = 0.35, width = 0.5, height = 0.2) +
  draw_plot(TFs_pvalue_plots[["FOSL2"]], x = 0.5, y = 0.305, width = 0.46, height = 0.265)  + 
  draw_plot(TFs_pvalue_plots[["GATA4"]], x = 0, y = 0.07, width = 0.5, height = 0.25) +
  draw_plot(TFs_pvalue_plots[["HNF4A"]], x = 0.47, y = 0.052, width = 0.5, height = 0.27) +
  
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), x = c(0,0,0, 0.5, 0, 0.5) + 0.035, 
                  y =c(1,0.8,0.58,0.58,0.33,0.33) - 0.02, size = 8)

save_plot("results_refined/Fig4_v6.png", fig4_g, base_width = 5, base_height = 10)





### compare TFs from promoter and enhancer
compare_TF_pro_enh <- function(Enh_dat_sig, TssA_dat_sig, factor_names, add = 1){
  pro_enh_overlap = c()
  pro_enh_overlap_TFs = list()
  for(group in unique(Enh_dat_sig$group)){
    tssa_tfs = TssA_dat_sig[TssA_dat_sig$group == group, "TF"]
    enh_tfs  = Enh_dat_sig[Enh_dat_sig$group == group, "TF"]
    x = length(intersect(tssa_tfs, enh_tfs))
    pro_enh_overlap_TFs[[as.character(group)]] = intersect(tssa_tfs, enh_tfs)
    pro_enh_overlap = c(pro_enh_overlap, length(tssa_tfs)-x, length(enh_tfs)-x,x)
  }
  
  pro_enh_overlap = data.frame(matrix(pro_enh_overlap, ncol=3, byrow = T))
  colnames(pro_enh_overlap) = c("only_TssA", "only_Enh", "Both")
  pro_enh_overlap$group = unique(Enh_dat_sig$group)
  
  pro_enh_overlap = melt(pro_enh_overlap, id.vars = 'group')
  if(is.numeric(pro_enh_overlap$group[1])){
    pro_enh_overlap$Factor = factor_names[pro_enh_overlap$group + add]
    pro_enh_overlap$Factor = factor(pro_enh_overlap$Factor, levels = factor_names)
  }else{
    pro_enh_overlap$Factor = pro_enh_overlap$group
    pro_enh_overlap[pro_enh_overlap$Factor == "Shared", "Factor"] = "Ubiquitous"
  }
  return(pro_enh_overlap)
}



sn_spMF_sig = sn_spMF[sn_spMF$BH_pvalue < 0.05, ]
TssA_dat_sn_spMF_sig = sn_spMF_sig[sn_spMF_sig$region == 'TssA', ]
Enh_dat_sn_spMF_sig = sn_spMF_sig[sn_spMF_sig$region == 'Enh', ]
pro_enh_overlap = compare_TF_pro_enh(Enh_dat_sn_spMF_sig, TssA_dat_sn_spMF_sig, factor_names)

x = pro_enh_overlap[order(pro_enh_overlap$Factor), ]
y = x[seq(3, dim(x)[1],3), "value"] + x[seq(1, dim(x)[1],3), "value"] + x[seq(2, dim(x)[1],3), "value"]
min(x[seq(3, dim(x)[1],3), "value"] / y)
max(x[seq(3, dim(x)[1],3), "value"] / y)
median(x[seq(3, dim(x)[1],3), "value"] / y)

fig_overlap = ggplot(data = pro_enh_overlap,aes(x=Factor, y=value, fill = variable)) + 
  geom_bar(stat = 'identity') + 
  main_theme_Fig5 + 
  ylab("Number of enriched TFs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position="right",
        legend.margin=margin(0,0,0,-2)) +
  scale_fill_manual(values = c("red", "yellow", "orange")) +
  theme(axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size = 0.5))

save_plot("results_refined//O_overlap_TssA_Enh_TFs_sn_spMF.png", fig_overlap, base_width = 5, base_height = 3)




heuristic2_sig = heuristic_2[heuristic_2$BH_pvalue < 0.05, ]
TssA_dat_heuristic2_sig = heuristic2_sig[heuristic2_sig$region == 'TssA', ]
Enh_dat_heuristic2_sig = heuristic2_sig[heuristic2_sig$region == 'Enh', ]
pro_enh_overlap = compare_TF_pro_enh(Enh_dat_heuristic2_sig, TssA_dat_heuristic2_sig, 
                                     c("Ubiquitous", tissues), add = 1)

x = pro_enh_overlap[order(pro_enh_overlap$Factor), ]
y = x[seq(3, dim(x)[1],3), "value"] + x[seq(1, dim(x)[1],3), "value"] + x[seq(2, dim(x)[1],3), "value"]
min(x[seq(3, dim(x)[1],3), "value"] / y)
max(x[seq(3, dim(x)[1],3), "value"] / y)
median(x[seq(3, dim(x)[1],3), "value"] / y)

fig_overlap = ggplot(data = pro_enh_overlap,aes(x=Factor, y=value, fill = variable)) + 
  geom_bar(stat = 'identity') + 
  main_theme_Fig5 + 
  ylab("Number of enriched TFs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position="right",
        legend.margin=margin(0,0,0,-2)) +
  scale_fill_manual(values = c("red", "yellow", "orange")) + 
  theme(axis.line.x = element_line(),
        axis.line.y = element_line())

save_plot("results_refined//O_overlap_TssA_Enh_heuristic_2_TFs.png", fig_overlap, base_width = 5, base_height = 3)




### RUN this for one time!
#for(tfi in unique(Enh_dat$TF)){
tfi = 'GABPA'
png(paste0('results_refined/TF_pvalues/Enh_noRes//TF_',tfi,'.png' ) , width = 3, height = 3, units = 'in', res = 300)
fig_i = plot_bar_pvalue(tfi, Enh_dat_sp)
print(fig_i)
dev.off()

#}



TF_list = c("HIF1A", "GATA6", "CLOCK",
            "BARHL1", "ASCL1","LHX2",
            "MYOG","ARNTL",
            "NR5A2","HNF1A", "FOXA1","ID2", "CREB1")

TF_list = c("HNF1A", "CLOCK", "ID2",
            "NR5A2", "CREB1","GATA6",
            "BARHL1", "ASCL1",
            "MYOG","ARNTL", "FOXA1", "HIF1A")
TFs_pvalue_plots = list()
for(tfi in TF_list){
  TFs_pvalue_plots[[tfi]] = plot_bar_pvalue(tfi, Enh_dat_sp)
}

png('Supplmentary_plots/Q_pvalur_barplots_2_revision.png', width = 8, height = 8, units = 'in', res = 200)
grid.newpage()
grid.draw(
  cbind(
    rbind(ggplotGrob(TFs_pvalue_plots[[1]]), ggplotGrob(TFs_pvalue_plots[[5]]), ggplotGrob(TFs_pvalue_plots[[9]]), size = "first"),
    rbind(ggplotGrob(TFs_pvalue_plots[[2]]), ggplotGrob(TFs_pvalue_plots[[6]]), ggplotGrob(TFs_pvalue_plots[[10]]), size = "first"), 
    rbind(ggplotGrob(TFs_pvalue_plots[[3]]), ggplotGrob(TFs_pvalue_plots[[7]]), ggplotGrob(TFs_pvalue_plots[[11]]), size = "first"), 
    rbind(ggplotGrob(TFs_pvalue_plots[[4]]), ggplotGrob(TFs_pvalue_plots[[8]]), ggplotGrob(TFs_pvalue_plots[[12]]), size = "first"), 
    size = "first"))
dev.off()






## TFs that significance correlates with/without TPM 
plot_cor_tpm_pvalue <- function(tfi, plot = F){
  x = Enh_dat[Enh_dat$TF == tfi, ]
  a = log10(x$TF_TPM)
  b = -log10(x$pvalue)
  if(plot){
    plot(a,b, main = tfi, xlab = 'log10(TPM)', ylab = '-log10(p-value)')
  }
  return(summary(lm(a~b))$adj.r.squared)
}
r2 = sapply(TssA_dat_sig$TF, function(tfi) plot_cor_tpm_pvalue(tfi))
r2[is.na(r2)] = 0
r2[r2 < 0] = 0
TssA_dat_sig$r2 = r2

TssA_dat_sig_r2 = unique(TssA_dat_sig[,c("TF", "r2")])
print(paste0("Among ", dim(TssA_dat_sig_r2)[1], " TFs enriched in enhancer, there are ", 
             sum(TssA_dat_sig_r2$r2>0.2), " TFs with r2 > 0.2"))

for(tfi in TssA_dat_sig_r2[TssA_dat_sig_r2$r2 > 0.5, "TF"]){
  plot_cor_tpm_pvalue(tfi, plot = T)
}




## ASB of HNF4A 
compute_OR_SE_ASB <- function(vec){
  ### NOT same as the funct ion in fig3.R
  or = (vec['inTest'] / vec['notinTest']) / (vec['bg_inTest'] / vec['bg_notinTest'])
  logor = log(or)
  se_logor = sqrt(1/vec['inTest'] + 1/vec['notinTest'] + 1/vec['bg_inTest'] + 1/vec['bg_notinTest'])
  ci_logor = c(logor - 1.96 * se_logor, logor + 1.96 * se_logor)
  return(c(logor, se_logor, ci_logor))
}
compute_meta_OR_CI_ASB <- function(test_df){
  ### NOT same as the funct ion in fig3.R
  test_df = test_df[,c("inTest", "notinTest", "bg_inTest", "bg_notinTest")]
  test_df$inTest = as.numeric(test_df$inTest)
  test_df$notinTest = as.numeric(test_df$notinTest)
  test_df$bg_inTest = as.numeric(test_df$bg_inTest)
  test_df$bg_notinTest = as.numeric(test_df$bg_notinTest)
  testdf_or_se = t(apply(test_df, 1, function(x) compute_OR_SE_ASB(x)))
  colnames(testdf_or_se) = c("logOR", "logSE", "logCI_low", "logCI_high")
  #testdf_meta = rma(yi = logOR, sei = logSE, data = testdf_or_se)
  result = exp(testdf_or_se)
  colnames(result) = c("OR", "SE", "CI_low", "CI_high")
  return(result)
}
wrap_up_asb_multiple_TFs <- function(TFs, reads_filter){
  df_plot_tfs = data.frame()
  for(TF in TFs){
    dat_asb_tissues = read.table(paste0('data/ASB/ASB_',TF,'_Filtered',reads_filter,'_peaks1_Enrichment.txt'), sep='\t', header=T, stringsAsFactors = F, row.names = 1)
    dat_asb_tissues = dat_asb_tissues[dat_asb_tissues$group == "matchedRC_005", ]
    colnames(dat_asb_tissues) =  c("inTest", "notinTest", "bg_inTest", "bg_notinTest", "OR","pV","group", "factor")
    dat_asb_tissues = dat_asb_tissues[dat_asb_tissues$factor %in% Enh_dat[Enh_dat$TF == 'HNF4A', 'group'], ]
    
    dat_asb_tissues_plot = as.data.frame(compute_meta_OR_CI_ASB(dat_asb_tissues))
    dat_asb_tissues_plot$group = dat_asb_tissues$group
    dat_asb_tissues_plot$factor = dat_asb_tissues$factor
    dat_asb_tissues_plot$TF = TF
    dat_asb_tissues_plot$pvalue = dat_asb_tissues$pV
    dat_asb_tissues_plot$inTest_total = dat_asb_tissues$inTest + dat_asb_tissues$bg_inTest
    
    df_plot = dat_asb_tissues_plot
    df_plot$factor = factor_names[df_plot$factor+1]
    df_plot$factor = factor(df_plot$factor, levels = rev(unique(df_plot$factor)))
    
    df_plot$group = factor(df_plot$group) 
    
    df_plot_tfs = rbind(df_plot_tfs, df_plot)
  }
  
  #x = df_plot_tfs[df_plot_tfs$TF == 'HNF4A',]
  #df_plot_tfs$factor = factor(df_plot_tfs$factor, levels = rev(x$factor[order(x$OR)]))
  
  tp = Enh_dat[Enh_dat$TF == 'HNF4A', ]
  group_levels = tp$Factor_name[order(tp$pvalue, decreasing = F)]
  df_plot_tfs$factor = factor(df_plot_tfs$factor, group_levels)
  
  ggplot_df = ggplot() +
    geom_pointrange(data=df_plot_tfs, aes(x=factor(factor), y=OR, ymin=CI_low, ymax=CI_high, group = TF, col = factor(TF)), 
                    fatten = 0.03, position=position_dodge(width=0.5), show.legend = TRUE, lwd = 0.4) + 
    geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    #coord_flip() +  # flip coordinates (puts labels on y axis) 
    main_theme_Fig5 + 
    xlab("") + ylab("Odds ratio") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.title = element_blank()) + 
    theme(legend.position = c(0.9, 0.9)) + 
    guides(shape = guide_legend(override.aes = list(size = 0.1)),
           color = guide_legend(override.aes = list(size = 0.1))) +
    scale_color_manual(values=c("grey", "black"),  name="") + 
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) + 
    scale_y_continuous(limits = c(0.4, 1.75), breaks = seq(0.5, 2, by = 0.5))
  
  return(ggplot_df)
}

TFs = c("HNF4A", "CTCF")
Fig4G = wrap_up_asb_multiple_TFs(TFs, 10)



fig4_g = ggdraw() + 
  draw_plot(Fig4A, x = 0, y = 0, width = 0.6, height = 0.7) +
  draw_plot(Fig4B, x = 0.6, y = 0, width = 0.3, height = 0.7)
save_plot("results_refined/Fig4_v5_flashr.png", fig4_g, base_width = 7.4)






TssA_plot = plot_summary_TF(TssA_dat_sp,factor_names)
save_plot("Supplmentary_plots/O_number_factors_TssA.png",  TssA_plot[[2]], base_width = 7.4, base_height = 5)

Enh_plot = plot_summary_TF(Enh_dat_sp,factor_names)
save_plot("Supplmentary_plots/O_number_factors_Enh.png",  Enh_plot[[2]], base_width = 7.4, base_height = 5)


### 

interesting_hits = Enh_dat_sig[Enh_dat_sig$score > 0.1, c("TF", "score", "TF_TPM", "count", "Factor_name")]
interesting_hits = interesting_hits[interesting_hits$Factor_name!='Shared', ]
interesting_hits[order(interesting_hits$TF), ]

