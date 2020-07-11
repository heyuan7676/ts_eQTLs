setwd('/Users/Yuan/Desktop/plots/')
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gplots)
library(RColorBrewer) 
library(distances)
library(scales)
library(wesanderson)
library(coloc)

main_theme_Fig6 = theme(axis.title.x = element_text(size = 8),
                        axis.title.y = element_text(size = 8),
                        legend.title = element_text(size = 6), 
                        legend.text  = element_text(size = 5),
                        axis.text = element_text(size = 5),
                        legend.position = 'top',
                        legend.margin=margin(0,0,0,1),
                        legend.box.margin=margin(-10,-10,-10,0),
                        legend.key.size = unit(0.05, "cm")) +
  theme(strip.text = element_text(size=4)) + 
  theme(axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) + 
  theme_bw() + 
  background_grid(major = 'none', minor = 'none')   + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(colour = 'black', size = 0.5), 
        axis.line.y = element_line(colour = 'black', size = 0.5), 
        axis.ticks = element_line(size = 0.3))



### liver-specific eQTL
pr_data_point = read.table('data/forplot/Fig6_example_y_w_exp2.txt', sep='\t', header=T)
colnames(pr_data_point) = c("effect_size", "standard_error")
#pr_data_point$standard_error = 1/pr_data_point$standard_error
rownames(pr_data_point) = gsub("-", "", tissues)
pr_data_point = as.data.frame(pr_data_point)

conf.level = 0.95
ci.vals = -qnorm( ( 1 - conf.level ) / 2 )
pr_data_point$lower = pr_data_point$effect_size - ci.vals * pr_data_point$standard_error
pr_data_point$higher = pr_data_point$effect_size + ci.vals * pr_data_point$standard_error
pr_data_point$col    = gtex_col[rownames(pr_data_point), "tissue_color_hex"]
pr_data_point$tissue = factor(rownames(pr_data_point), levels = rownames(pr_data_point))

Fig6A = ggplot() +
  geom_pointrange(data=pr_data_point, aes(x=tissues, y=effect_size, ymin=lower, ymax=higher,color=col), fatten=0.5) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  #coord_flip() +  # flip coordinates (puts labels on y axis) 
  scale_color_identity()+
  xlab("") + ylab("Effect size") +
  main_theme_Fig6  + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(size = 0.5))



### ASB
asb_dat = read.table('data/forplot/ASB_liver_HNF4A_ChIP_seqAligned_chr8_9325848.txt', sep='\t', header=F, stringsAsFactors = F)
asb_reads = as.data.frame(table(asb_dat$V4))
colnames(asb_reads) = c("Allele", "Number of reads")
asb_reads$Allele = c("A", "G")
asb_reads = melt(asb_reads)


Fig6B = ggplot(data = asb_reads, aes(x = Allele, y = value)) + 
  geom_bar(stat = 'identity') + 
  ylab("Number of reads") + 
  main_theme_Fig6



### TNKS TPM
#tpm_dat = read.table('data/forplot/TPM_TNKS_Liver.txt', sep='\t', header=T, stringsAsFactors = F)
#colnames(tpm_dat) = sapply(colnames(tpm_dat), function(x) paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep='.'))

#tpm_dat = read.table('data/forplot/TPM_normalized_TNKS_Liver.txt', sep='\t', header=T, stringsAsFactors = F)

tpm_dat = read.table('data/forplot/TNKS_residuals.txt', sep='\t', header=T, stringsAsFactors = F)
vcf_dat = read.table("data/forplot/Genotype_chr8_9325848.vcf", sep='\t', header=T, stringsAsFactors = F)
vcf_dat = t(vcf_dat)

tpm_dat = as.data.frame(t(tpm_dat))
colnames(tpm_dat) = c("TPM")
sampleIDs = intersect(rownames(tpm_dat), rownames(vcf_dat))


tpm_dat = data.frame(tpm_dat[sampleIDs,])
tpm_dat$Genotype = vcf_dat[sampleIDs, ]
colnames(tpm_dat) = c("TPM", "Genotype")

tpm_dat[tpm_dat$Genotype == '1|1', "Genotype"] = "GG"
tpm_dat[tpm_dat$Genotype == '1|0', "Genotype"] = "AG"
tpm_dat[tpm_dat$Genotype == '0|1', "Genotype"] = "AG"
tpm_dat[tpm_dat$Genotype == '0|0', "Genotype"] = "AA"

tpm_dat$Dose = 0
tpm_dat[tpm_dat$Genotype == 'AG', "Dose"] = 1
tpm_dat[tpm_dat$Genotype == 'AA', "Dose"] = 2

summary(lm(TPM~Dose, data = tpm_dat))

tpm_dat = melt(tpm_dat[,c("TPM", "Genotype")])

Fig6C = ggplot(data = tpm_dat, aes(x = Genotype, y = value)) + 
  geom_boxplot() + 
  ylab("Normalized expression of TNKS") + 
  main_theme_Fig6



Fig6 = ggdraw() + 
  draw_plot(Fig6A, x = 0, y = 0.55, width = 0.65, height = 0.42) + 
  draw_plot(Fig6B, x = 0.005, y = 0, width = 0.16, height = 0.5) + 
  draw_plot(Fig6C, x = 0.2, y = 0, width = 0.2, height = 0.5)
  #draw_plot_label(label = c("A", "B","C", "D", "E", "F"), x = c(0,0,0.18,0.35, 0.58, 0.58) + 0.02,  y =c(1, 0.65, 0.65, 0.65, 1, 0.65) - 0.02, size = 8)

save_plot("results_refined/Fig6_part1_revision.png", Fig6, base_width = 7.4)




Fig6 = ggdraw() + 
  draw_plot(pair7[[1]], x = 0, y = 0.072, width = 0.23, height = 0.9) + 
  draw_plot(pair7[[2]], x = 0.20, y = 0.19, width = 0.18, height = 0.715) + 
  draw_plot(pair7[[3]], x = 0.38, y = 0.32, width = 0.08, height = 0.4) + 
  draw_plot(Fig6B, x = 0.45, y = 0.6, width = 0.22, height = 0.45) + 
  draw_plot(Fig6C, x = 0.45, y = 0.3, width = 0.22, height = 0.45) + 
  draw_plot_label(label = c("A", "B","C", "D", "E", "F"), x = c(0,0.45,0.45,0.45, 0.7, 0.7) + 0.02,  y =c(1, 1,0.7,0.4,0.7,0.35) - 0.02, size = 8)

save_plot("results_refined/Fig6_part1.png", Fig6, base_width = 7.4)


