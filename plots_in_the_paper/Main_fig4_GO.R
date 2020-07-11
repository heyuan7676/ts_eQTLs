source('scripts_refined/utils.R')

library(distances)
library(scales)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(40)
myCol <- colorRampPalette(c("dodgerblue", "aliceblue", "white"))(100)
myCol <- colorRampPalette(c("red", "indianred1","lightpink","mistyrose","mintcream","lightcyan","aliceblue", "azure","white"))(100)

file_used = 'c5.bp.v6.2.symbols.gmt.txt'


main_theme_Fig3 = theme(axis.title.x = element_text(size = 8),
                        axis.title.y = element_text(size = 8),
                        legend.title = element_text(size = 6), 
                        legend.text  = element_text(size = 5),
                        axis.text.x = element_text(size = 5),
                        axis.text.y = element_text(size = 5),
                        legend.key.size = unit(0.1, "inches"),
                        legend.position = 'top',
                        legend.margin=margin(0,0,0,1),
                        legend.box.margin=margin(-10,-10,-10,0)) +
  background_grid(major = 'none', minor = 'none')



##### SAVE -- used for plots !!!
higher = 500
lower = 20
hit = 10



plot_go_per_factor <- function(dat_GO, group, prefix){
  dat_GO_group = dat_GO[dat_GO$group == group, ]
  dat_GO_group = dat_GO_group[order(dat_GO_group$BH_p), ]
  dat_GO_group = head(dat_GO_group, n=30)
  
  dat_GO_group$geneSet = factor(dat_GO_group$geneSet, levels = rev(dat_GO_group$geneSet))
  go_group = ggplot(data = dat_GO_group) + 
    geom_point(aes(x=geneSet, y=logP, color = OR)) + 
    coord_flip() + 
    theme_bw() + 
    scale_colour_gradient2() + 
    main_theme_Fig3 + 
    #theme(axis.text.y = element_text(angle = 45)) + 
    xlab("") + 
    ylab("-log10(p-value)") + 
    theme(legend.position = c(0.8,0.7)) + 
    background_grid(major = 'xy', minor = 'xy') + 
    theme(panel.border = element_blank(),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size = 0.5))
  save_plot(paste0("Supplmentary_plots/GO_enrichment_per_factor/",prefix,"_group", group, ".png"), go_group, base_width = 6)
  return(go_group)
}

readin_GO_result <- function(fn, HIT = 50, HIGHER = 500, LOWER = 50, alpha = 0.05, shared_text = NULL){
  dat_GO = read.table(fn, sep='\t', stringsAsFactors = F, header=T)
  dat_GO = dat_GO[dat_GO$test_inset + dat_GO$background_inset > HIT, ]
  dat_GO['OR_adjusted'] = ((dat_GO[,2]+1) / (dat_GO[,3] + 1)) /  ((dat_GO[,4] + 1) / (dat_GO[,5] + 1))
  dat_GO = dat_GO[,c(1,2,6,14, 7,10,11,12)]
  colnames(dat_GO) = c("geneSet", "numberGenes","OR", "OR_adjusted","pvalue", "group", "geneInSet", "BH_p")
  dat_GO = dat_GO[dat_GO$geneInSet < HIGHER, ]
  dat_GO = dat_GO[dat_GO$geneInSet > LOWER, ]
  
  dat_GO$BH_p = p.adjust(dat_GO$pvalue, method = 'BH')
  dat_GO$logP = -log10(dat_GO$pvalue)
  dat_GO$geneSet = gsub("GO_", "", dat_GO$geneSet)
  
  sig_dat = dat_GO[dat_GO$BH_p < alpha, ]
  print(paste0("all pathway:", length(unique(sig_dat$geneSet))))
  
  if(!is.null(shared_text)){
    sig_dat = dat_GO[dat_GO$group != shared_text, ] 
    sig_dat = sig_dat[sig_dat$BH_p < alpha, ]
  }else{
    sig_dat = dat_GO[dat_GO$BH_p < alpha, ]
  }
  enriched_dat = sig_dat[sig_dat$OR > 1, ]
  depleted_dat = sig_dat[sig_dat$OR < 1, ]
  
  enriched_N = length(unique(enriched_dat$geneSet))
  depleted_N = length(unique(depleted_dat$geneSet))
  
  print(paste0("ts_pathway:", enriched_N, " enriched; ", depleted_N, " depleted."))
  
  return(dat_GO)
}

alpha = 0.05
spMF_go = readin_GO_result('data/Fig3_GSEA_spMF.txt', HIT = hit, HIGHER = higher, LOWER = lower, alpha = alpha, shared_text = 0)
for(g in seq(0,22)){
  temp = plot_go_per_factor(spMF_go, g, 'spMF')
}




############################################ compare to thresholding and flashr

alpha = 0.05
spMF_go = readin_GO_result('data/Fig3_GSEA_spMF.txt', 
                           HIT = hit, HIGHER = higher, LOWER = lower, 
                           alpha = alpha, shared_text = 0)
flashr_go_df = readin_GO_result('data/Fig3_GSEA_flashr_df.txt', 
                             HIT = hit, HIGHER = higher, LOWER = lower, 
                             alpha = alpha, shared_text = 0)
flashr_go_bf = readin_GO_result('data/Fig3_GSEA_flashr_bf.txt', 
                             HIT = hit, HIGHER = higher, LOWER = lower, 
                             alpha = alpha, shared_text = 0)
flashr_NN_go = readin_GO_result('data/Fig3_GSEA_flashr_NN.txt', 
                             HIT = hit, HIGHER = higher, LOWER = lower, 
                             alpha = alpha, shared_text = 0)
thresholding_go = readin_GO_result('data/Fig3_GSEA_thresholding.txt', 
                                   HIT = hit, HIGHER = higher, LOWER = lower, 
                                   alpha = alpha, shared_text = 'Shared')
thresholding_2_go = readin_GO_result('data/Fig3_GSEA_Thresholding_Refined.txt', 
                                     HIT = hit, HIGHER = higher, LOWER = lower, 
                                     alpha = alpha, shared_text = -1)
softimpute_go = readin_GO_result('data/Fig3_GSEA_softImpute.txt', 
                                 HIT = hit, HIGHER = higher, LOWER = lower,
                                 alpha = alpha, shared_text = 0)
pmd_go = readin_GO_result('data/Fig3_GSEA_PMD.txt', 
                                 HIT = hit, HIGHER = higher, LOWER = lower,
                                 alpha = alpha, shared_text = 0)
pmd2_go = readin_GO_result('data/Fig3_GSEA_PMD_cv2.txt', 
                                 HIT = hit, HIGHER = higher, LOWER = lower,
                                 alpha = alpha, shared_text = 0)

spMF_go = spMF_go[spMF_go$group !=0, ]
flashr_go = flashr_go[flashr_go$group!=0, ]
thresholding_go = thresholding_go[thresholding_go$group != 'Shared' ,]
thresholding_2_go = thresholding_2_go[thresholding_2_go$group != -1,]

pvalue_list = list("sn_spMF"   = spMF_go[spMF_go$OR > 1, "pvalue"], 
                   "heuristic_1" = thresholding_go[thresholding_go$OR > 1, "pvalue"], 
                   "heuristic_2" = thresholding_2_go[thresholding_2_go$OR > 1, "pvalue"])
go_full_set = qqunif.plot(pvalue_list,auto.key=list(corner=c(.95,.05)))


############################################ compare to thresholding and flashr
alpha = 0.05
spMF_go_stringent = readin_GO_result('data/Fig3_GSEA_spMF_stringent.txt', HIT = hit, HIGHER = higher, LOWER=lower, alpha = alpha, shared_text = 0)
flashr_go_stringent = readin_GO_result('data/Fig3_GSEA_flashr_NN_stringent.txt', HIT = hit, HIGHER = higher, LOWER=lower, alpha = alpha, shared_text = 0)
thresholding_go_stringent = readin_GO_result('data/Fig3_GSEA_Thresholding_stringent.txt', HIT = hit, HIGHER = higher, LOWER=lower, shared_text = 'Shared')
thresholding_2_go_stringent = readin_GO_result('data/Fig3_GSEA_Thresholding_Refined_stringent.txt', HIT = hit, HIGHER = higher, LOWER=lower, shared_text = -1)
softimpute_go_stringent = readin_GO_result('data/Fig3_GSEA_softImpute_stringent.txt', HIT = hit, HIGHER = higher, LOWER = lower,alpha = alpha, shared_text = 0)
pmd_go_stringent = readin_GO_result('data/Fig3_GSEA_PMD_stringent.txt', HIT = hit, HIGHER = higher, LOWER = lower,alpha = alpha, shared_text = 0)
pmd2_go_stringent = readin_GO_result('data/Fig3_GSEA_PMD_cv2_stringent.txt', HIT = hit, HIGHER = higher, LOWER = lower,alpha = alpha, shared_text = 0)


spMF_go_stringent = spMF_go_stringent[spMF_go_stringent$group !=0, ]
thresholding_go_stringent = thresholding_go_stringent[thresholding_go_stringent$group != 'Shared', ]
thresholding_2_go_stringent = thresholding_2_go_stringent[thresholding_2_go_stringent$group != -1, ]

pvalue_list = list("sn_spMF"   = spMF_go_stringent[spMF_go_stringent$OR > 1, "pvalue"], 
                   "heuristic_1" = thresholding_go_stringent[thresholding_go_stringent$OR > 1, "pvalue"],
                   "heuristic_2" = thresholding_2_go_stringent[thresholding_2_go_stringent$OR > 1, "pvalue"])

go_stringent_set = qqunif.plot(pvalue_list,auto.key=list(corner=c(.95,.05)))


png('results_refined/Fig4_GO_qqplots_ORlarger1_four_methods_revision.png', width = 6, height = 5, units = 'in', res = 300)
print(c(go_full_set, go_stringent_set))
dev.off()




###################################### plot heatmap

alpha = 0.1
spMF_go_stringent = readin_GO_result('data/Fig3_GSEA_SparseMF_coph_v8_cbset_95_allPairs_filteredGenes.ciseQTL_results.complete_filteredSNPs.LDblocks_0.2_topPair_K30_a11_l110_stringent.txt', 
                                     HIT = hit, HIGHER = higher, LOWER=lower, 
                                     alpha = alpha, shared_text = 0)
spMF_go_stringent$geneSet = sapply(spMF_go_stringent$geneSet, 
                                   function(x) tolower(gsub("_", " ", x)))
enriched_gos = sort(unique(spMF_go_stringent[(spMF_go_stringent$pvalue < 0.001) & 
                                               (spMF_go_stringent$OR > 1), "geneSet"]))
factor_names_used = factor_names

dat_heatmap = spMF_go_stringent
dat_heatmap$logP = -log10(dat_heatmap$pvalue)
dat_heatmap = dat_heatmap[which(dat_heatmap$geneSet %in% enriched_gos), ]

dat_bh_p = dat_heatmap[,c("geneSet", "group", "logP")]
dat_bh_p = dcast(dat_bh_p, geneSet ~ group)
rownames(dat_bh_p) = dat_bh_p$geneSet
dat_bh_p = dat_bh_p[,seq(2,dim(dat_bh_p)[2])]
dat_bh_p[is.na(dat_bh_p)] = 0

# Run clustering for pathways
dat.matrix <- as.matrix(dat_bh_p)
dat.dendro <- as.dendrogram(hclust(d = dist(x = dat.matrix)), method = 'median')
geneset.order <- order.dendrogram(dat.dendro)

# Run clustering for factors
dat.matrix <- t(as.matrix(dat_bh_p))
dat.dendro <- as.dendrogram(hclust(d = dist(x = dat.matrix)), method = 'centroid')
factor.order <- order.dendrogram(dat.dendro)

dat_heatmap$geneSet = factor(x = dat_heatmap$geneSet, levels = rownames(dat_bh_p)[geneset.order], ordered = TRUE)
dat_heatmap$group = sapply(dat_heatmap$group, function(x) factor_names_used[x+1])
dat_heatmap$group = factor(x = dat_heatmap$group, levels = factor_names_used[factor.order], ordered = TRUE)


###  plot
dt_plot = dat_heatmap
alpha = 0.1
enriched_gos = sort(unique(spMF_go_stringent[(spMF_go_stringent$BH_p< alpha) & (spMF_go_stringent$OR > 1), "geneSet"]))
dt_plot = dt_plot[which(dt_plot$geneSet %in% enriched_gos), ]
include_factors = unique(dt_plot[dt_plot$BH_p < alpha, "group"])
dt_plot = dt_plot[sapply(dt_plot$group, function(x) x%in% include_factors),]
dt_plot[dt_plot$logP > 8, "logP"] = 8


dt_plot_point_dt = dt_plot[dt_plot$BH_p < alpha, ]


panel_width = unit(0.2, "npc")
go_plot = ggplot(data = dt_plot, aes(x = geneSet, y = group)) + 
  geom_tile(aes(fill = logP, width = 1, height = 1)) +
  scale_fill_gradient2(high = 'deepskyblue1', low = 'white' , na.value = "black")+
  geom_text(data = dt_plot_point_dt,  aes(x=geneSet, y=group, label = round(OR,1)), size = 1.5)+ 
  coord_flip() +
  xlab("") + ylab("") +  
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
  main_theme_Fig3 + 
  guides(fill= guide_colorbar(title = paste0(expression('-log'[10]), '(P value)'), 
                              barwidth = panel_width,
                              title.position = 'right')) + 
  theme(panel.background = element_rect(fill = "grey90"), 
        panel.border = element_blank()) + 
  theme(axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black'))


png('results_refined//Fig4_GO_enrichment_revision.png', width = 5.6, height = 8, units = 'in', res = 300)
print(go_plot)
dev.off()

