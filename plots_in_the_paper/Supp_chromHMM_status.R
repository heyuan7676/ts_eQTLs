setwd('~/Desktop/plots/')
source('scripts_refined/utils.R')

main_theme_Fig4 = theme(axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        legend.title = element_text(size = 7), 
                        legend.text  = element_text(size = 7),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8),
                        legend.key.size = unit(0.2, "inches"),
                        legend.position = 'top',
                        legend.margin=margin(0,0,0,1),
                        legend.box.margin=margin(-5,-5,-5,0)) +
  background_grid(major = 'none', minor = 'none')  


LMfn = 'spMF'
background='random'

status_list = c("1_TssA", "7_Enh")
  
for(status in status_list){
    fn = paste0('data/chromHMM_DNase_enrich_consistent/', LMfn, '_Active_enrichment_ROADMAP_',status,'_to', background, '_WithInOpen.txt')
    testN = read.table(fn, sep=',', header=T, stringsAsFactors = F)
    tissues_with_ENCODE = testN$inTissue
    
    fn = paste0('data/chromHMM_DNase_enrich_consistent/', LMfn, '_Active_enrichment_N_ROADMAP_',status,'_to', background, '_WithInOpen.txt')
    testN = read.table(fn, sep=',', header=T, stringsAsFactors = F)
    testN$tissue = tissues_with_ENCODE
    testN$matched = testN$matched == 1
    testN = testN[testN$group!=0, ]
    
    g = ggplot(aes(x = tissue, y = inTest, fill = matched), data = testN) + 
      geom_bar(stat="identity", width=.5, position = "dodge")
}



