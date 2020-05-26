suppressWarnings(library(ggplot2))
suppressWarnings(library(roloc))
suppressWarnings(library(readr))
suppressWarnings(library(dplyr))

### gtex colors for the tissues
gtex_col = read_tsv('data/gtex_colors.txt')
gtex_col = gtex_col[order(gtex_col$tissue), ]

factor_matrix = read.table('output/sn_spMF_K17_a1100_l150/sn_spMF_FactorMatrix_K17_a1100_l150_Run1.txt', sep='\t', stringsAsFactors = F, header = T, row.names = 1)
rownames(factor_matrix) = gsub('-', '', rownames(factor_matrix))
gtex_col = gtex_col %>% filter(tissue %in% rownames(factor_matrix))


color_bar = data.frame("x" = gtex_col$tissue, "y" = rep(0.2,1), "color" = gtex_col$color_hex)
color_fig = ggplot() + 
  geom_point(data = color_bar, aes(x = x, y = y), col = color_bar$color, size = 2) + 
  ylim(0,1) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, size=5, hjust = 1, vjust = 0.5)) + 
  theme(axis.line=element_blank(),axis.text.y = element_blank()) +
  xlab("")  + ylab("") + 
  scale_x_discrete(expand = c(0, 5)) + 
  theme(panel.background = element_rect(fill = "transparent"))


png("output/color.png", width = 800, height = 400, res = 200, bg="transparent")
print(color_fig)
dev.off()
