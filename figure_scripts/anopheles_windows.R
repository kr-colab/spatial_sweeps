setwd("~/phd/research/spatial_sweeps")
library(tidyverse)
library(data.table)
library(MASS)
library(patchwork)
library(ggiraph)
library(scales)
library(quantreg)
library(viridis)
library(splines)
library(npreg)
load("../../locator/data/cntrymap.Rdata")

plasma_gen = viridis_pal(option='C')
plasma = plasma_gen(15)
mako_gen = viridis_pal(option='G')
mako = mako_gen(15)

library(ggnewscale)

# Fig. 4: ANOPHELES GENOMIC WINDOWS

# dataframe for genome-wide
allsnphist <- fread('anopheles/out/all_snps_windowed_genome_scan_annotated.txt')
allsnphist$INV_group <- allsnphist$INV

# dataframe for quantiles
recsnphist <- fread('anopheles/out/high_recombination_windowed_genome_scan_annotated.txt')
recsnphist$INV_group <- recsnphist$INV

Rbstart <- 18575300
Rbstop <- 26767588
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$mids <= Rbstart),]$INV_group <- '2Rb1'
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$breaks <= Rbstart),]$INV_group <- '2Rb1'
recsnphist[(recsnphist$chromosome=='2R' &
              recsnphist$mids <= Rbstart),]$INV_group <- '2Rb1'
recsnphist[(recsnphist$chromosome=='2R' &
              recsnphist$breaks <= Rbstart),]$INV_group <- '2Rb1'
Rcstart <- 26750000
Rcstop <- 31473000
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$mids >= Rcstop),]$INV_group <- '2Rc2'
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$breaks >= Rcstop),]$INV_group <- '2Rc2'
recsnphist[(recsnphist$chromosome=='2R' &
              recsnphist$mids >= Rcstop),]$INV_group <- '2Rc2'
recsnphist[(recsnphist$chromosome=='2R' &
              recsnphist$breaks >= Rcstop),]$INV_group <- '2Rc2'
Lastart <- 20524058
Lastop <- 42165532
allsnphist[(allsnphist$chromosome=='2L' &
              allsnphist$mids <= Lastart),]$INV_group <- '2La1'
allsnphist[(allsnphist$chromosome=='2L' &
              allsnphist$breaks <= Lastart),]$INV_group <- '2La1'
allsnphist[(allsnphist$chromosome=='2L' &
              allsnphist$mids >= Lastop),]$INV_group <- '2La2'
allsnphist[(allsnphist$chromosome=='2L' &
              allsnphist$breaks >= Lastop),]$INV_group <- '2La2'
recsnphist[(recsnphist$chromosome=='2L' &
              recsnphist$mids <= Lastart),]$INV_group <- '2La1'
recsnphist[(recsnphist$chromosome=='2L' &
              recsnphist$breaks <= Lastart),]$INV_group <- '2La1'
recsnphist[(recsnphist$chromosome=='2L' &
              recsnphist$mids >= Lastop),]$INV_group <- '2La2'
recsnphist[(recsnphist$chromosome=='2L' &
              recsnphist$breaks >= Lastop),]$INV_group <- '2La2'



allsnphist<- rbind(allsnphist[allsnphist$chromosome != '2L',], 
                   allsnphist[allsnphist$chromosome=='2L' & allsnphist$INV_group != '',])
allsnphist<- rbind(allsnphist[allsnphist$chromosome != '2R',], 
                   allsnphist[allsnphist$chromosome=='2R' & allsnphist$INV_group != '',])
allsnphist <- allsnphist[!is.na(allsnphist$chromosome),]
allsnphist[allsnphist$inv=='','INV'] <- 'none'
allsnphist <- allsnphist[allsnphist$chromosome != 'X',]

known_loci <- data.frame(list('locus'=c('VGSC','GSTE','Cyp6p','Ace1','Rdl','Tep1'),
                              'start'=c(2358158, 28591663, 28420677, 3484107, 25363652, 11202091),
                              'stop'=c(2431617, 28602280, 28505816, 3495790, 25434556, 11206882),
                              'chromosome'=c('2L','3R','2R','2R', '2L','3L')))
known_loci$mids <- rowMeans(known_loci[,c('start','stop')])

# dataframe for masked data (don't plot low-recombination areas)
msnphist <- allsnphist
msnphist$statistic[msnphist$mean_cMMb < 1.5] <- NA

# vignette locations
cyp <- fread('anopheles/out/AGAP008552_outliers.csv')
cyp$Vignette <- 'Cytochrome P450\nCYP4H2'
cub <- fread('anopheles/out/AGAP029542_outliers.csv')
cub$Vignette <- 'CUB domain-\ncontaining protein'
grs <- fread('anopheles/out/NS-GR_outliers.csv')
grs$position <- grs$position.x
grs$chromosome <- grs$chromosome.x
grs$nstat <- grs$nstat.x
grs$Vignette <- 'Gustatory receptor'
vgs <- rbind(cyp[,c('position', 'chromosome', 'nstat', 'Vignette')],
             cub[,c('position', 'chromosome', 'nstat', 'Vignette')],
             grs[,c('position', 'chromosome', 'nstat', 'Vignette')])

ggplot() + theme_minimal() +
  
  # background (including masked regions)
  geom_line(data=allsnphist,
            aes(x=mids, y=statistic, 
                group=INV_group), 
            color='lightgray',
            lwd=0.4) +
  
  # quantile
  geom_line(data=recsnphist,
            aes(x=mids, y=Q, 
                group=INV_group,
                color=INV),
            alpha=0.75) +
  scale_color_manual(name='Inversion', values=rev(c(mako[2], mako[8], mako[11], mako[13])),
                     guide='none') +
  new_scale_color() +
  
  # low-recombination regions masked
  geom_line(data=msnphist,
            aes(x=mids, y=statistic, 
                color=INV,
                group=INV_group)) +
  scale_color_manual(name='Inversion', values=rev(c(mako[1], mako[7], mako[10], mako[14])), 
                     na.value = 'transparent') +
  
  # labels
  scale_x_continuous(name='Position (Mbp)',
                     breaks=c(0, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7),
                     labels=c('0', '10','20','30','40','50','60')) +
  ylab('Z score') +
  facet_wrap(vars(chromosome), scales = "free_x") +
  
  scale_y_continuous(limits=c(-4, NA))+#, trans='pseudo_log') +
  
  # annotate
  new_scale_color() +
  geom_point(data=known_loci, aes(x=mids, y=-3), shape=17) +
  geom_text(data=known_loci, aes(x=mids, y=-3.75, label=locus), size=8/.pt) +
  geom_point(data=vgs, aes(x=position, y=nstat+2, fill=Vignette), pch=21) +
  scale_fill_viridis(option='C', discrete = TRUE)

ggsave('figures/aggenomicwindows_annotated.pdf', width=8, height=4, units='in')  

  
  
  
  