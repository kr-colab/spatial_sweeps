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
library(ggrepel)
load("../../locator/data/cntrymap.Rdata")

plasma_gen = viridis_pal(option='C')
plasma = plasma_gen(15)
mako_gen = viridis_pal(option='G')
mako = mako_gen(15)

library(ggnewscale)

## CUB DOMAIN MAP

# outlier snps...
cub <- fread('anopheles/out/AGAP029542_outliers.csv')
annotation <- fread('anopheles/data/cub-domain-annotation.gff', sep='\t', fill=TRUE, skip=2)

type=c()
start=c()
stop=c()
strand=c()
name=c()

for (i in seq(nrow(annotation))){
  if (annotation[i,1] == '###'){ next }
  type <- c(type, annotation[i, 3][[1]])
  start <- c(start, annotation[i, 4][[1]])
  stop <- c(stop, annotation[i, 5][[1]])
  strand <- c(strand, annotation[i, 7][[1]])
  res <- str_match(annotation[i, 9], "Name=\\s*(.*?)\\s*;")
  name <- c(name, res[,2])
}

annotation <- data.frame(list('type'=type, 
                              'start'=start, 
                              'stop'=stop, 
                              'strand'=strand, 
                              'name'=name))
annotation$gene <- str_remove(annotation$name, '-RA')
annotation$ypos <- -0.05
annotation[annotation$gene == 'AGAP013740',]$ypos <- -0.15
yp <- -0.15
for (a in unique(annotation[annotation$start > 6350000 & annotation$start < 6375000,] %>% 
                    arrange(start) %>% pull(gene))){
  if (a == 'AGAP029542') next
  annotation[annotation$gene == a,]$ypos <- yp
  yp <- yp-0.05
}
yp <- -0.15
for (a in unique(annotation[annotation$start > 6375000,] %>% 
                 arrange(start) %>% pull(gene))){
  if (a == 'AGAP029542') next
  annotation[annotation$gene == a,]$ypos <- yp
  yp <- yp-0.1
}


pt <- ggplot() + theme_minimal() + xlab('Position')

for (g in unique(annotation$gene)){
  if (g == 'AGAP029542') next # plot focal gene last...
  
  exn <- annotation[(annotation$gene == g & annotation$type == 'CDS'),] %>% arrange(start)
  exn$end_intron <- c(exn$start[-1]-1, NA)
  exn$mid_pt <- (exn$end_intron - exn$stop)/2 + exn$stop
  pt <- pt + # exon blocks
             geom_rect(data=exn, aes(xmin=start, xmax=stop,
                          ymin=ypos-0.025, ymax=ypos+0.025,
                          fill=gene), linewidth=0, show.legend = F) +
             # left half of intron
             geom_segment(data=exn, aes(x=stop, xend=mid_pt,
                                        y=ypos, yend=ypos+0.01,
                                        color=gene), show.legend = F) +
             # right half of intron
              geom_segment(data=exn, aes(x=mid_pt, xend=end_intron,
                                       y=ypos+0.01, yend=ypos,
                                       color=gene), show.legend = F) + 
                geom_segment(data=annotation[(annotation$gene == g & annotation$type == 'five_prime_UTR'),],
                               aes(x=start, xend=stop,
                                   y=ypos, yend=ypos,
                                   color=gene), show.legend = F) +
              geom_segment(data=annotation[(annotation$gene == g & annotation$type == 'three_prime_UTR'),],
                                  aes(x=start, xend=stop,
                                      y=ypos, yend=ypos,
                                      color=gene), show.legend = F)
             # label gene
              if (exn$start < 6325000 ||
                  exn$start > 6375000 & exn$start < 6387500){
              # right aligned
                pt <- pt + geom_text(data=exn, aes(x=max(stop), y=ypos,
                                          label=gene), size=8/.pt, nudge_x=1000, hjust=0)
              }
              # left aligned
              else {
                pt <- pt + geom_text(data=exn, aes(x=min(start), y=ypos,
                                          label=gene), size=8/.pt, nudge_x=-1000, hjust=1)
              }
}

exn <- annotation[(annotation$gene == 'AGAP029542' & annotation$type == 'CDS'),] %>% arrange(start)
exn$ypos <- -0.05
exn$end_intron <- c(exn$start[-1]-1, NA)
exn$mid_pt <- (exn$end_intron - exn$stop)/2 + exn$stop
print(exn)

dts <- function(data, m=0.2, b=0.05){
  return(data*m + b)
}

dtr <- function(data, m=0.2, b=0.05){
  return((data - b)/m)
}

pt <- pt + 
  geom_rect(data=exn, aes(xmin=start, xmax=stop,
                             ymin=ypos-0.075, ymax=ypos+0.075), 
               fill=mako[12], color=mako[10],
               linewidth=0.1, 
               show.legend = F) +
  geom_segment(data=exn, aes(x=stop, xend=mid_pt,
                             y=ypos, yend=ypos+0.05), 
               color=mako[10], 
               show.legend = F) +
  geom_segment(data=exn, aes(x=mid_pt, xend=end_intron,
                             y=ypos+0.05, yend=ypos), 
               color=mako[10], 
               show.legend = F) + 
  # 5 prime utr
  geom_segment(data=annotation[(annotation$gene == g & annotation$type == 'five_prime_UTR'),],
               aes(x=start, xend=stop,
                   y=ypos, yend=ypos), color=mako[10], show.legend = F) +
  # 3 prime utr
  geom_segment(data=annotation[(annotation$gene == g & annotation$type == 'three_prime_UTR'),],
               aes(x=start, xend=stop,
                   y=ypos, yend=ypos), color=mako[10], show.legend = F) +
  geom_text(data=exn, aes(x=min(start), y=ypos,
                          label=gene), size=8/.pt, nudge_x=-2000, hjust=1) +
  scale_color_viridis(option='G', discrete=TRUE) +
  scale_fill_viridis(option='G', discrete=TRUE) +
  new_scale_fill() +
  # add snp data
  geom_point(data=cub,
             aes(x=position, y=dts(frequency), fill=area), pch=21) +
  scale_fill_viridis(option='C', name=expression(Area~(km^2)), labels=label_scientific()) +
  # with new axis!
  scale_y_continuous(
    # Add a second axis and specify its features
    sec.axis = sec_axis(~dtr(.)), 
    name="Frequency", 
    breaks=seq(0.05, 0.25, length=5),
    labels=c(0.00, 0.25, 0.50, 0.75, 1.00)
    ) 


pt

ggsave('figures/CUB-domaining-snps.pdf', width=11, height=3, units='in')
