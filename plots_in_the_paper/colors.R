library(roloc)

color_bar = data.frame("x" = gtex_col$tissue_name, "y" = rep(0,1), "color" = gtex_col[, "tissue_color_hex"])
pdf("results/color.pdf")
ggplot() + 
  geom_point(data = color_bar, aes(x = x, y = y), col = color_bar$color, size = 2) + 
  ylim(0,1) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, size=5, hjust = 0.95)) + 
  theme(axis.line=element_blank(),axis.text.y = element_blank()) +
  xlab("")  + ylab("") + 
  scale_x_discrete(expand = c(0, 5))
dev.off()

color_bar$x_rev = factor(color_bar$x, levels = rev(color_bar$x))
png("results/color.png", res = 1000, width = 3, height = 5, units = 'in')
ggplot() + 
  geom_point(data = color_bar, aes(x = x_rev, y = y), col = color_bar$color, size = 2) + 
  ylim(0,1) + 
  theme(axis.ticks = element_blank(), axis.text.y = element_text(size=5, hjust = 1)) + 
  theme(axis.line=element_blank(),axis.text.x = element_blank()) +
  xlab("")  + ylab("") + 
  coord_flip()
dev.off()
