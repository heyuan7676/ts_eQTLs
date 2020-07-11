setwd('/Users/Yuan/Desktop/plots/')

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(lemon)


factor_names = c("Ubiquitous", "Brain tissues", "Pituitary", "Spleen", "Colon; Small intestine", "Thyroid", 
                 "Adipose Visceral Omentum", "Nerve Tibial", "Adipose; Mammary", 
                 "Whole Blood", "Brain Cerebellum tissues", "Heart tissues", "Skin tissues", "LCL", 
                 "Colon; Esophagus", "Liver", "Stomach", "Testis", "Artery tissues", 
                 "Muscle Skeletal", "Pancreas", "Lung", "Esophagus Mucosa")


factor_names_flashr = c("Ubiquitous", "Brain tissues", "Skin;Esophagus", "Artery; Whole Blood",
                        "Adipose; Digestive; Whole Blood", "Whole Blood", "Heart; Muscle", 
                        "Brain Cerebellum", "Adipose; Artery; Esophagus", 
                        "Brain 2", "LCL", "Testis", "Muscle Skeletal", "Thyroid", 
                        "Nerve_Tibial", "Esophagus_Mucosa", "Pancreas", "Lung",
                        "Colon; Small intestine", "Spleen; Whole Blood", "Liver", 
                        "Adrenal Gland", "Brain 3")


factor_names_flashr_bf = c("Ubiquitous", "Brain tissues", "Skin;Whole Blood", "Skin; Whole Blood",
                        "LCL; Muscle Skeletal", "Whole Blood", "Heart; Muscle", 
                        "Brain Cerebellum", "Adipose; Artery; Esophagus", 
                        "Brain 2", "LCL", "Testis", "Muscle Skeletal", "Thyroid", 
                        "Nerve_Tibial", "Esophagus_Mucosa", "Pancreas", "Lung",
                        "Colon; Small intestine", "Spleen; Whole Blood", "Liver", 
                        "Adrenal Gland", "Brain 3")



factor_names_thresholding = c("Ubiquitous","Adipose", "Adrenal_Gland", "Artery", "Brain",
                              "Cells_EBV-transformed_lymphocytes", "Cells_Cultured_fibroblasts",
                              "Colon", "Esophagus", "Heart","Kidney_Cortex", "Liver",
                              "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", 
                              "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary",
                              "Prostate", "Skin", "Small_Intestine_Terminal_Ileum", "Spleen",
                              "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

### factor matrix
factor_matrix = read.table('data/Fig1_factor_matrix.txt', sep='\t', stringsAsFactors = F, row.names = 1)
factor_matrix = factor_matrix[rownames(factor_matrix) != "", ]
rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
colnames(factor_matrix) = paste0("Factor", seq(1, dim(factor_matrix)[2]))

factor_matrix = factor_matrix / matrix(rep(apply(factor_matrix, 2, function(x) max(abs(x))), dim(factor_matrix)[1]), nrow = dim(factor_matrix)[1], byrow = T)


factor_matrix$Tissue <- factor(rownames(factor_matrix), levels = gtex_col$tissue_id)
factor_bar_df = melt(factor_matrix, id.vars = 'Tissue')

factor_bar_df$variable_plot = as.character(factor_bar_df$variable)
for(k in seq(1,23)){
  factor_bar_df[factor_bar_df$variable_plot == paste0('Factor', k), "variable_plot"] = paste0('Factor', k, '\n(', factor_names[k], ')')
}
factor_bar_df$variable_plot = factor(factor_bar_df$variable_plot, levels = unique(factor_bar_df$variable_plot))



fig1_g = ggplot(data = factor_bar_df, aes(x = Tissue, y =value, fill = Tissue)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = gtex_col$tissue_color_hex)+
  facet_rep_wrap(~variable_plot, ncol = 4) +  #repeat.tick.labels = 'left'
  xlab("Tissues") + 
  ylab("Values in the factors")+
  theme(strip.text.x = element_text(size = 8)) + 
  theme(axis.ticks.x =element_blank()) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y =element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15))   + 
  theme_bw() + 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line())

png('Supplmentary_plots/Fig1_revision.png', width = 8, height = 12, units = 'in', res = 200)
print(fig1_g) 
dev.off()
