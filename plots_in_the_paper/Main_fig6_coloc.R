setwd('/Users/Yuan/Desktop/plots/')
library(ggplot2)
library(reshape2)
library(gridExtra)

### Mahattan plot in ggplot2: https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html

Mahattan_plot <- function(location, pv, SNP_center, title, background_color = 'grey', SNP_color = 'orange'){
  df = data.frame(matrix(c(location, pv), ncol=2))
  colnames(df) = c("chromosome_location", "P")
  df$CHR = 8
  
  df$is_highlight = 'none'
  df[df$chromosome_location == SNP_center, "is_highlight"] = "yes"
  df[df$chromosome_location != SNP_center, "is_highlight"] = "no"
  
  axisdf = df %>% group_by(CHR) %>% summarize(center=( max(chromosome_location) + min(chromosome_location) ) / 2 )
  
  xlim_min = max(SNP_center - window_size, 0)
  xlim_max = min(SNP_center + window_size, max(location))
  g = ggplot(df, aes(x=chromosome_location, y=-log10(P))) +
    ggtitle(title) +
    xlim(xlim_min, xlim_max) + 
    # Show all points
    #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.7, col = background_color) +
    geom_point(alpha=0.8, size=1, col = background_color) +
    #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    # custom X axis:
    #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    geom_point(data=subset(df, is_highlight=="yes"), color=SNP_color, size=1.2) +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}



SNP_center = 9325848
window_size = 1e6
ciseQTL_fn = 'Liver_ENSG00000173273.15.txt'

## cis-eQTL result
ciseQTL_wrap <- function(fn, figTitle, savefn){
  beta_liver_dat = read.table(paste0('data/forplot/', fn), sep='\t', header=T, stringsAsFactors  = F)
  SNP_location = as.vector(sapply(beta_liver_dat$variant_id,function(x) as.numeric(strsplit(x, "_")[[1]][2])))
  pv = beta_liver_dat$pval_nominal
  
  idx = which((SNP_location > SNP_center - window_size) * (SNP_location < SNP_center + window_size) == 1)
  beta_liver_dat = beta_liver_dat[idx, ]
  
  SNP_location = as.vector(sapply(beta_liver_dat$variant_id,function(x) as.numeric(strsplit(x, "_")[[1]][2])))
  pv = beta_liver_dat$pval_nominal
  
  png(paste0('results/', savefn), width = 4, height = 2, units = 'in', res = 400)
  g = Mahattan_plot(SNP_location, pv, SNP_center, figTitle)
  print(g)
  dev.off()
  return(g)
  
}
g = ciseQTL_wrap(ciseQTL_fn, 'ciseQTL_Liver', 'Fig6_ciseQTL_TNKS.png')

####### GWAS result
gwas_wrap <- function(fn, figTitle, savefn){
  gwas_dat = read.table(paste0('data/forplot/', fn), sep='\t', header=F, stringsAsFactors = F)
  gwas_dat = gwas_dat[complete.cases(gwas_dat$V2), ]
  gwas_SNP_location = as.vector(sapply(gwas_dat$V2,function(x) as.numeric(strsplit(x, "_")[[1]][2])))
  gwas_pvalue = gwas_dat$V8
  
  png(paste0('results/',savefn), width = 4, height = 2.1, units = 'in', res = 400)
  g_gwas = Mahattan_plot(gwas_SNP_location, gwas_pvalue, SNP_center, figTitle, "skyblue")
  print(g_gwas)
  dev.off()
  
  return(g_gwas)
}

g1 = gwas_wrap("GLGC_Mc_HDL_chr8_9325848_1MB.txt", "GWAS results for Mc_HDL","Fig6_GWAS_HDL.png")
g2 = gwas_wrap("GLGC_Mc_LDL_chr8_9325848_1MB.txt", "GWAS results for Mc_LDL","Fig6_GWAS_LDL.png")
g3 = gwas_wrap("UKB_high_cholesterol_chr8_9325848_1MB.txt", "GWAS results for UKB_high_cholesterol","Fig6_GWAS_high_cholesterol.png")



### colocalization analysis
library(coloc)
coloc_analysis <- function(ciseQTL_fn, gwas_fn){
  dat1 = read.table(paste0('data/forplot/', ciseQTL_fn), sep='\t', stringsAsFactors = F, header = T)
  dat2 = read.table(paste0('data/forplot/', gwas_fn), sep='\t', stringsAsFactors = F, header = T)
  
  liver_N = 208
  dat1$MAF = sapply(dat1$maf, function(x) min(c(x,1-x)))
  dat2$MAF = sapply(dat2$effect_allele_freq, function(x) min(c(x,1-x)))
  
  my.res <- coloc.abf(dataset1=list(pvalues=dat1$pval_nominal, N=liver_N, MAF = dat1$MAF, 
                                    type="quant", snp=dat1$variant_id),
                      dataset2=list(pvalues=dat2$pvalue, N=dat2$sample_size, MAF=dat2$MAF, 
                                    type="quant",snp=dat2$gtex_variant_id))
  return(my.res)
  #dat = my.res$results
  #tail(dat[order(dat$SNP.PP.H4), c("snp", "pvalues.df1", "pvalues.df2", "SNP.PP.H4")])
}


hdl.coloc = coloc_analysis(ciseQTL_fn, 'HDL.test.ma')
whichP = names(which.max(hdl.coloc$summary[2:6]))
Pvalue = max(hdl.coloc$summary[2:6])
hit = hdl.coloc$results[hdl.coloc$results$snp == 'chr8_9325848_A_G_b38',]
Pvalue_snp = hit['SNP.PP.H4']
tail(hdl.coloc$results[order(hdl.coloc$results$SNP.PP.H4), c("snp", "pvalues.df1", "pvalues.df2", "SNP.PP.H4")])

g1 = g1 + annotate("text", x = SNP_center+1e6/2, y = 40, size = 3,
                   label = paste0(whichP, ": ", round(Pvalue,2), "\n SNP.PP.H4: ", round(Pvalue_snp,2)))

ldl.coloc = coloc_analysis(ciseQTL_fn, 'LDL.test.ma')
whichP = names(which.max(ldl.coloc$summary[2:6]))
Pvalue = max(ldl.coloc$summary[2:6])
hit = ldl.coloc$results[ldl.coloc$results$snp == 'chr8_9325848_A_G_b38',]
Pvalue_snp = hit['SNP.PP.H4']
tail(ldl.coloc$results[order(ldl.coloc$results$SNP.PP.H4), c("snp", "pvalues.df1", "pvalues.df2", "SNP.PP.H4")])

g2 = g2 + annotate("text", x = SNP_center+1e6/2, y = 20, size = 3,
                   label = paste0(whichP, ": ", round(Pvalue,2), "\n SNP.PP.H4: ", round(Pvalue_snp,2)))

colesterol.coloc = coloc_analysis(ciseQTL_fn, 'high_cholesterol_test.ma')
whichP = names(which.max(colesterol.coloc$summary[2:6]))
Pvalue = max(colesterol.coloc$summary[2:6])
hit = colesterol.coloc$results[colesterol.coloc$results$snp == 'chr8_9325848_A_G_b38',]
Pvalue_snp = hit['SNP.PP.H4']
tail(colesterol.coloc$results[order(colesterol.coloc$results$SNP.PP.H4), c("snp", "pvalues.df1", "pvalues.df2", "SNP.PP.H4")])


g3 = g3 + annotate("text", x = SNP_center+1e6/2, y = 7, size = 3,
                   label = paste0(whichP, ": ", round(Pvalue,2), "\n SNP.PP.H4: ", round(Pvalue_snp,2)))

width = 0.5
fig_coloc = ggdraw() + 
  draw_plot(g, x = 0, y = .5, width = width, height = 0.5) + 
  draw_plot(g1, x = 0, y = 0, width = width, height = 0.5)  +
  draw_plot(g, x = width, y = .5, width = width, height = 0.5) +
  draw_plot(g2, x = width, y = 0, width = width, height = 0.5)
#draw_plot_label(label = c("A", "B", "C", "D", "E"), x = c(0,0,0.4,0.6,0.8), y =c(1,0.55,1,1,1), size = 8)

save_plot("results/Fig6_Coloc_PPP1R3B.png", fig_coloc, base_width = 7.4)



#### GWAS after COJO

cojo_dat = read.table('data/forplot/cojo.out.cma.cojo', sep='\t', header=T, stringsAsFactors = F)

cojo_SNP_location = cojo_dat$bp
cojo_pvalue = cojo_dat$pC





####### aggregate information

tpm_hnf4a_tnks = read.table('data/forplot/HNF4A_TNKS_TPM.txt', sep='\t', header=T, stringsAsFactors = F, row.names = 1)
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '0|0', "Genotype"] = 'AA'
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '1|1', "Genotype"] = 'GG'
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '1|0', 'Genotype'] = 'GA'
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '0|1', 'Genotype'] = 'GA'


tpm_hnf4a_tnks = read.table('data/forplot/HNF4A_AC_TPM.txt', sep='\t', header=T, stringsAsFactors = F, row.names = 1)
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '0|0', "Genotype"] = 'AA'
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '1|1', "Genotype"] = 'GG'
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '1|0', 'Genotype'] = 'GA'
tpm_hnf4a_tnks[tpm_hnf4a_tnks$Genotype == '0|1', 'Genotype'] = 'GA'


png('results/Fig6_HNF4A_TNKS_TPM_Genotpyes.png', width = 5, height = 2, units = 'in', res = 200)
ggplot(tpm_hnf4a_tnks) + geom_point(aes(x = HNF4A, y = TNKS), data = tpm_hnf4a_tnks) + 
  geom_boxplot(aes(y=TNKS)) + 
  facet_wrap(~Genotype)
dev.off()
