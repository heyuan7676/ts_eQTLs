setwd('/Users/Yuan/Desktop/plots/')
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gplots)
library(RColorBrewer) 
library(distances)
library(scales)
library(wesanderson)
library(VennDiagram)
library(RVenn)
library(RColorBrewer) 
library(preprocessCore)
library(latticeExtra)
library(gplots)
library(metafor)
library(colortools)
library(dplyr)
library(cowplot)

library(plyr)
library(dplyr)
tissues = read.table('data/tissues.txt', stringsAsFactors = F)$V1

shared_factor_color = 'navy'
ts_factor_color = 'dodgerblue'
ts_factor_color_unmatched = 'grey'


factor_seq_names = c("Ubiquitous", seq(2,23))


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

factor_names_flashr_NN = c('Ubiquitous','Brain', 'Spleen;Whole Blood', 
                           'Muscle_Skeletal', 'Esophagus;Skin', 'Cells_EBV-transformed_lymphocytes', 
                           'Testis', 'Artery', 'Thyroid', 'Pancreas;Stomach', 'Nerve_Tibial', 'Adipose;Breast_Mammary_Tissue', 
                           'Brain_Cerebellum', 'Colon;Esophagus', 'Esophagus_Mucosa', 'Heart', 'Colon;Spleen',
                           'Adrenal_Gland', 'Brain2', 'Cells_Cultured_fibroblasts', 'Lung', 'Pituitary', 'Liver')



factor_names_thresholding = c("Ubiquitous","Adipose", "Adrenal_Gland", "Artery", "Brain",
                              "Cells_EBV-transformed_lymphocytes", "Cells_Cultured_fibroblasts",
                              "Colon", "Esophagus", "Heart","Kidney_Cortex", "Liver",
                              "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", 
                              "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary",
                              "Prostate", "Skin", "Small_Intestine_Terminal_Ileum", "Spleen",
                              "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

#Generate custom QQ plot
library(lattice)
#https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R#Under_the_Hood:_Multiple_P-value_Lists
qqunif.plot <- function(pvalues, 
                        should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                        xlab=expression(paste("Expected (",-log[10], " p-value)")),
                        ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                        draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                        already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                        par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}


tissues = read.table('data/tissues.txt', sep='\t', header=F, stringsAsFactors = F)
tissues = tissues$V1

### gtex colors
gtex_col = read.table('data/gtex_colors.txt', sep='\t', header = T, stringsAsFactors = F)
gtex_col = gtex_col[order(gtex_col$tissue_id),]
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)
rownames(gtex_col) = gtex_col$tissue_id
head(gtex_col)
