### Legend

panel_width = unit(0.15, "npc")
plot_object = plot_object + 
  theme(legend.key.size = unit(0.1, "inches"),
        legend.position = 'top',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-20,0)) +
  guides(fill= guide_colorbar(title = 'Effect size', 
                              barwidth = panel_width,
                              title.position = 'top'))


## discrete:  + guides(fill= guide_legend(title.position = 'top'))
