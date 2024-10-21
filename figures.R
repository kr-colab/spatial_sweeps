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
pal = c('#182706',
        '#27301D',
        '#32413A',
        '#384A4C', 
        '#375B69',
        '#2E6D74',
        '#24AA87',
        '#25C589',
        '#2CD38B',
        '#36E48D',
        '#54EC87',
        '#66ED87',
        '#7CEA89',
        '#99ED95',
        '#ADEB9E')
palgen = colorRampPalette(pal)
palgen_rev = colorRampPalette(rev(pal))
colgen = colorRamp(pal)


# Fig. 1: SIMULATION

## cartoon of simulation

getsim <- function(tick){
  filepath=paste0('simulation/sweep_simulation/', tick, '.tick.txt')
  df <- data.frame(t(fread(filepath)))
  colnames(df) <- df[1, ]
  df <- df[-1, ]
  df[] <- lapply(df, function(x) type.convert(as.numeric(x)))
  df$mutated <- df$muts > 0
  df$frequency <- (sum(df$muts)) / (nrow(df)*2)
  df$tick <- tick
  row.names(df) <- NULL
  return(df)
}

df40 <- getsim(40)
df100 <- getsim(100)
df160 <- getsim(160)
simdf <- rbind(df40, df100, df160)

simpt <- ggplot() + theme_minimal(base_size = 10) + 
  geom_point(data=simdf %>% arrange(mutated), aes(x=x, y=y, color=mutated), size=1) +
  scale_color_manual(values=c(mako[1], mako[13]), name='Allele state', labels=c('Ancestral','Derived')) +
  coord_fixed(ratio=1) +
  ggtitle('Tick') +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=10),
        strip.text = element_text(size=8)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  facet_wrap(vars(tick))
simpt
padf <- fread('simulation/out/power_analysis_processed.txt')
summ_df <- padf %>% group_by(scoef, Nw, tick) %>% dplyr::summarize(area = mean(area),
                                                                   distance=mean(distance),
                                                                   frequency=mean(frequency),
                                                                   N_mutations=mean(N_mutations),
                                                                   outliers=mean(outliers),
                                                                   outlier_SNPdetected=mean(outlier_SNPdetected),
                                                                   windowoutliers=mean(windowoutliers),
                                                                   windowoutlier_SNPdetected=mean(windowoutlier_SNPdetected),
                                                                   windowoutlier_nonoutlier_SNPdetected=mean(windowoutlier_nonoutlier_SNPdetected))

# joint distribution
jdist <- ggplot(summ_df %>% arrange(scoef)) + theme_minimal(base_size = 10) +
  geom_smooth(aes(x=frequency, y=area, color=scoef, group=scoef),
              method='loess', span=0.05, linewidth=0.75) +
  facet_wrap(vars(Nw)) +
  scale_color_viridis(option='G', trans='log10', end = 0.9, name = 'Selection \ncoefficient') +
  xlab('Frequency') + ylab('Area') +
  ggtitle('Neighborhood size')+
  theme(plot.title = element_text(hjust = 0.5, vjust=0, size=10),
        legend.key.width = unit(0.75, 'line'),
        strip.text = element_text(size=8))

simpt / jdist + plot_annotation(tag_levels = 'A')
ggsave('figures/sweepsimjointdistribution.pdf', width=8, height=6, units='in')

# Fig. 2: POWER ANALYSIS

## example of how it works
snpdf <- fread('simulation/out/simulation_1260_frequency_area.txt')
snpdf <- snpdf[snpdf$area > 0,]

N = 1000
snpdf <- snpdf %>% arrange(by = freq)
bins <- rep(seq(from = 0, to=as.integer(nrow(snpdf)/ N)), each=N)[1:nrow(snpdf)]
snpdf$bin <- bins
snpdf <- snpdf %>% dplyr::group_by(bin) %>% dplyr::mutate(Q = quantile(area, probs=0.1),
                                            bin_frequency = mean(freq))
snpdf$outlier <- FALSE
snpdf[snpdf$area <= snpdf$Q,]$outlier <- TRUE
pajdist <- ggplot() + theme_minimal(base_size = 10) + 
  geom_point(data=snpdf, aes(x=freq, y=area), color=mako[4], alpha=0.05, pch=20, size=0.1) + 
  geom_point(data=snpdf[snpdf$scoef>0,], aes(x=freq, y=area), pch=21, fill=mako[14], size=2, stroke=0.75) +
  geom_line(data=snpdf, aes(x=bin_frequency, y=Q), linewidth=0.5, color=mako[1]) + 
  xlab('Frequency') + ylab('Area') +
  scale_y_continuous(limits=c(0,max(snpdf$area))) +
  theme(axis.line = element_line(colour = "black"))
pajdist


# windowed genome analysis
outlierhs <- hist(snpdf[snpdf$outlier == TRUE,]$position, 
                  breaks=1000)
breaks <- outlierhs$breaks[1:length(outlierhs$breaks)-1]
outlierhs <- data.frame(breaks=breaks, mids=outlierhs$mids, outlier=outlierhs$counts)
nonoutlierhs <- hist(snpdf[snpdf$outlier == FALSE,]$position, 
                     breaks=1000)
breaks <- nonoutlierhs$breaks[1:length(nonoutlierhs$breaks)-1]
nonoutlierhs <- data.frame(breaks=breaks, mids=nonoutlierhs$mids, nonoutlier=nonoutlierhs$counts)
snphist <- merge(outlierhs, nonoutlierhs)

snphist1 <- snphist[(snphist$nonoutlier == 0) & (snphist$outlier > 0),]
snphist1$ratio <- snphist1$outlier
snphist3 <- snphist[(snphist$nonoutlier > 0) & (snphist$outlier > 0),]
snphist3$ratio <- snphist3$outlier / snphist3$nonoutlier

snphist <- rbind(snphist1, snphist3)
snphist <- snphist %>% mutate(Q = quantile(ratio, probs=0.90, na.rm=TRUE))

pawindows <- ggplot(snphist) + theme_minimal(base_size = 8) +
  geom_line(aes(x=mids, y=ratio), color=mako[4], linewidth=0.1) + 
  geom_line(aes(x=mids, y=Q), color=mako[1], alpha=1, linewidth=0.5, linetype='dotted') +
  xlab('Position') + ylab('SF outlier ratio') +
  geom_point(aes(x=1e8/2, y=0.3), fill=mako[14], pch=25, size=2, stroke=0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6))
pawindows
# power analysis results
out <- ggplot(summ_df %>% arrange(-scoef)) + theme_minimal(base_size = 9) +
  #geom_hline(aes(yintercept=mean(padf$outliers / padf$N_mutations))) + 
  geom_smooth(aes(x=frequency, y=outlier_SNPdetected, color=scoef, group=scoef), 
              linewidth=0.5, method='loess', span=0.05) +
  facet_grid(cols=vars(Nw), scales='free') +
  scale_color_viridis(option='G', trans='log10', end = 0.9, name = 'Selection \ncoefficient') +
  xlab('Frequency') + ylab('P(SF outlier detected)') +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.key.width = unit(0.75, 'line'),
        axis.text.x = element_text(size=6))
wout <- ggplot(summ_df %>% arrange(-scoef)) + theme_minimal(base_size = 9) +
  #geom_hline(aes(yintercept=mean(padf$windowoutliers / padf$N_mutations))) + 
  geom_smooth(aes(x=frequency, y=windowoutlier_SNPdetected, color=scoef, group=scoef), 
              linewidth=0.5, method='loess', span=0.05) +
  facet_grid(cols=vars(Nw), scales='free') +
  scale_color_viridis(option='G', trans='log10', end = 0.9, name = 'Selection \ncoefficient') +
  xlab('Frequency') + ylab('P(WSF outlier detected)') +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.key.width = unit(0.75, 'line'),
        axis.text.x = element_text(size=6))

(pajdist + 
  inset_element(pawindows, 
                0.2, 0.01, 1, 0.45)) + 
(out / wout + plot_layout(guides = "collect") & theme(legend.text=element_text(size=6),
                                                      legend.title=element_text(size=8))) +
  plot_layout(widths = c(1, 1.5)) + plot_annotation(tag_levels = 'A')

ggsave('figures/poweranalysisres.pdf', width=8, height=4, units='in')

# Fig. 3: ANOPHELES JOINT DISTRIBUTION

# sampling locs
metadata <- fread('anopheles/data/admixture_k1_metadata.txt')
metadata_count <- metadata %>% dplyr::group_by(x, y) %>% dplyr::count()
ag_map <- ggplot()+coord_map(projection = "mollweide",
                             xlim=c(min(na.omit(c(metadata$x)))-2,
                                    max(na.omit(c(metadata$x)))+5),
                             ylim=c(min(na.omit(c(metadata$y)))-5,
                                    max(na.omit(c(metadata$y)))+5))+
  theme_classic()+theme(axis.title = element_blank(),
                        legend.title = element_text(size=8),
                        legend.text=element_text(size=6),
                        axis.text=element_text(size=6),
                        legend.box = "horizontal",
  )+
  geom_polygon(data=fortify(map),aes(x=long,y=lat,group=group),fill='white',color='black',lwd=0.2)+
  geom_count(data=metadata_count,aes(x=x,y=y,fill=n,size=n),color='black',pch=21) +
  scale_fill_viridis(option='C', guide = "legend")

ag_map


# joint distribution
allsnpsplot <- fread('anopheles/out/all_snps_INVincluded_annotated.txt')
allsnpsplot$plot_annotation <- 'none'
allsnpsplot[(allsnpsplot$annotation == '3_prime_UTR_variant' |
               allsnpsplot$annotation == '5_prime_UTR_variant' |
               allsnpsplot$annotation == '5_prime_UTR_premature_start_codon_gain_variant' |
               allsnpsplot$annotation == 'intron_variant' |
               allsnpsplot$annotation == 'splice_region_variant&intron_variant' |
               allsnpsplot$annotation == 'downstream_gene_variant' |
               allsnpsplot$annotation == 'upstream_gene_variant'),'plot_annotation'] <- 'Noncoding variant'
allsnpsplot[(allsnpsplot$annotation == 'missense_variant' |
               allsnpsplot$annotation == 'stop_gained'),'plot_annotation'] <- 'Nonsynonymous variant'
allsnpsplot[(allsnpsplot$annotation == 'synonymous_variant' |
               allsnpsplot$annotation == 'splice_region_variant&synonymous_variant'),'plot_annotation'] <- 'Synonymous variant'
allsnpsplot$plot_annotation <- factor(allsnpsplot$plot_annotation, levels=c('none','Noncoding variant','Synonymous variant','Nonsynonymous variant'))

allsnpsplot$INV <- factor(allsnpsplot$INV, levels=c('2La','2Rb','2Rc','none'))


# dataframe for outlier SNPs within outlier windows (ignoring intergenic regions)
outlier_notIG <- allsnpsplot[(allsnpsplot$outlier == TRUE & allsnpsplot$outlierwindow == TRUE & allsnpsplot$annotation != 'intergenic_region'),]
outlier_notIG$impact <- factor(outlier_notIG$impact, levels=c('MODIFIER', 'LOW', 'MODERATE', 'HIGH'))
outlier_notIG$plot_annotation <- factor(outlier_notIG$plot_annotation, levels=c('Noncoding variant','Synonymous variant','Nonsynonymous variant'))

missenseoutliersplot <- fread('anopheles/out/missense_outliers_thermompnn.txt')
# dataframe for missense outlier SNPs that change protein structure/stability
#missenseoutliersplot <- missenseoutliers[missenseoutliers$`ddG (kcal/mol)` != 0,] %>% dplyr::select(c('frequency','area','impact','annotation','INV','position','chromosome'))
missenseoutliersplot$bin <- 0
#missenseoutliersplot <- merge(missenseoutliersplot, allsnpsplot)#[allsnpsplot$outlier == TRUE & allsnpsplot$outlierwindow == TRUE & allsnpsplot$annotation == "stop_gained",] %>% dplyr::select(c('frequency','area','impact','annotation','INV','position','chromosome')))
missenseoutliersplot$plot_annotation <- 'Nonsynonymous variant'

ag_jdist <- ggplot() + theme_classic() +
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation=='none' & allsnpsplot$outlier==FALSE,], 
             mapping=aes(x=frequency, y=area), alpha=0.01, size=0.5, pch=20) + 
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation!='none' & allsnpsplot$outlier==FALSE,], 
             mapping=aes(x=frequency, y=area, color=plot_annotation), alpha=0.025, size=0.5, pch=20) +
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation=='none' & allsnpsplot$outlier==TRUE,], 
             mapping=aes(x=frequency, y=area), alpha=0.05, size=0.5, pch=20) + 
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation!='none' & allsnpsplot$outlier==TRUE,], 
             mapping=aes(x=frequency, y=area, color=plot_annotation), alpha=0.05, size=1, pch=20) + 
  geom_line(data=allsnpsplot, aes(x=frequency, y=cutoff, group=INV, linetype=INV), linewidth=0.5) + 
  geom_point(data=outlier_notIG %>% dplyr::arrange(impact),
             aes(x=frequency, y=area, color=plot_annotation, size=impact)) +
  scale_color_manual(labels=c('Noncoding variant','Synonymous variant','Nonsynonymous variant'),
                     values=c(mako[4],mako[10], mako[14])) +
  geom_point(data=missenseoutliersplot, aes(x=frequency, y=area, size=impact), shape=21, color='black', fill=mako[14]) +
  scale_size_manual(values=c(0.5, 1, 2, 3)) +
  scale_linetype_manual(values=c('dashed','dotdash','dotted','solid')) +
  scale_y_continuous(limits=c(0, NA))
ag_jdist

ag_jdist + 
  inset_element(ag_map +
                  theme(legend.box.spacing=margin(0),
                        legend.background = element_blank(), 
                        plot.background = element_blank()), 
                0.45, 0, 1, 0.45)
ggsave('figures/agjointdistribution.pdf', width=8, height=5, units='in')

# Fig. 4: ANOPHELES GENOMIC WINDOWS

allsnphist <- fread('anopheles/out/all_snp_histogram.txt')
allsnphist$INV_group <- allsnphist$INV

Rbstart <- 18575300
Rbstop <- 26767588
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$mids <= Rbstart),]$INV_group <- '2Rb1'
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$breaks <= Rbstart),]$INV_group <- '2Rb1'
Rcstart <- 26750000
Rcstop <- 31473000
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$mids >= Rcstop),]$INV_group <- '2Rc2'
allsnphist[(allsnphist$chromosome=='2R' &
              allsnphist$breaks >= Rcstop),]$INV_group <- '2Rc2'
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

allsnphist<- rbind(allsnphist[allsnphist$chromosome != '2L',], 
                   allsnphist[allsnphist$chromosome=='2L' & allsnphist$INV_group != '',])
allsnphist<- rbind(allsnphist[allsnphist$chromosome != '2R',], 
                   allsnphist[allsnphist$chromosome=='2R' & allsnphist$INV_group != '',])
allsnphist <- allsnphist[!is.na(allsnphist$chromosome),]
allsnphist[allsnphist$inv=='','INV'] <- 'none'
allsnphist <- allsnphist[allsnphist$chromosome != 'X',]

missenseoutliersplot$y <- 0

for (i in seq(nrow(missenseoutliersplot))){
  chr <- missenseoutliersplot[i,]$chromosome
  pos <- missenseoutliersplot[i,]$position
  offset <- max(allsnphist[allsnphist$chromosome=='2L','ratio'])/5
  missenseoutliersplot[i,'y'] <- (allsnphist[(allsnphist$breaks <= pos & (allsnphist$breaks + 100000) >= pos & allsnphist$chromosome==chr),]$ratio) + offset
}

window_plot <- ggplot() + theme_minimal() +
  geom_line(data = allsnphist, aes(x=mids, y=Q, group=INV_group, color=INV), 
            alpha=0.9, linewidth=0.3) + 
  scale_color_manual(values=rev(mako(8)[c(TRUE, FALSE)])) +
  new_scale_color() +
  geom_line(data = allsnphist, aes(x=mids, y=ratio, group=INV_group, color=INV), linewidth=0.3) + 
  scale_color_manual(values=rev(mako(8)[c(FALSE, TRUE)])) +
  
  
  #scale_linetype_manual(values=c('dashed','dotdash','dotted','solid')) +
  new_scale_color() +
  geom_point(data = missenseoutliersplot, aes(x=position, 
                                              y=y,
                                              size=frequency,
                                              fill=area), 
             pch=25,) +
  scale_size_continuous(range=c(0.3, 2), guide='legend', ) +
  scale_fill_viridis(option='C') +
  facet_wrap(vars(chromosome), scales = "free")#+
  #theme(legend.position = 'none')

window_plot

ggsave('figures/aggenomicwindows_annotated.pdf', width=8, height=4, units='in')


# Fig. 5: NONSYNONYMOUS WSF OUTLIERS

library(ggrepel)
library(ggnewscale)

mstex <- fread('anopheles/out/missense_outliers_manual_annotation.tsv')
mstex$`SNP frequency` <- as.numeric(mstex$`SNP frequency`)
mstex$`ddG (kcal/mol)` <- as.numeric(mstex$`ddG (kcal/mol)`)
mstex$`SNP area (km2)` <- as.numeric(mstex$`SNP area (km2)`)
mstex$Function <- factor(mstex$Function, levels=c('insecticide target','immune response','transcriptional regulation','DNA repair','sensory','cellular division','unannotated'))


mstex$xvar <- seq(nrow(mstex))
mstex$yvar <- mstex$`ddG (kcal/mol)`
mstex$yvar[mstex$yvar > 0] <- log10(mstex$yvar[mstex$yvar > 0])
mstex$yvar[mstex$yvar < 0] <- -1 * log10(abs(mstex$yvar[mstex$yvar < 0]))

ggplot() + theme_minimal() +
  geom_hline(aes(yintercept=0), alpha=0.5) +
  geom_text_repel(seed = 420,
                  data=mstex,
                  aes(x=`SNP frequency`, y=`ddG (kcal/mol)`, label=stringr::str_wrap(`Annotation`, 22)),
                  #xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
                  direction = "both",
                  lineheight = .75,
                  min.segment.length = 0,
                  size=2,
                  hjust='left',
                  segment.size=0.25,
                  segment.alpha=0.9,
                  force_pull=10,
                  force=0.5,
                  nudge_x=0.1,
                  nudge_y=0.025) +
  scale_size_continuous(trans='log10') +
  new_scale('size') +
  geom_point(data=mstex %>% arrange(-`ddG (kcal/mol)`),
             aes(x=`SNP frequency`, y=`ddG (kcal/mol)`, size=`SNP area (km2)`, color=`Function`)) +
  scale_color_manual(values=c(rev(mako_gen(6))[1:3], rev(plasma_gen(5))[1:3], 'lightgray')) +
    
    #'#24AA87','#54EC87', '#bcf0af', '#375B69','#808080','#A9A9A9','#DCDCDC')) +
  
  #scale_color_continuous(trans='log10') +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10, sigma=0.1)) +
  scale_x_continuous(trans='log10', breaks=c(0.1, 0.25, 0.5, 1))

ggsave('figures/agmissenseoutliers.pdf', width=8, height=4, units='in')

# SUPPLEMENT

# Fig. X: proportion of SNPs in power analysis

outlierprop <- ggplot(summ_df) + theme_minimal() +
  geom_smooth(aes(x=frequency, y=outliers/N_mutations, group=scoef, color=scoef), 
              method='loess', span=0.05) +
  facet_wrap(vars(Nw)) +
  scale_color_viridis(option='G', trans='log10', end = 0.9, name = 'Selection \ncoefficient')
woutlierprop <- ggplot(summ_df) + theme_minimal() +
  geom_smooth(aes(x=frequency, y=windowoutliers/N_mutations, group=scoef, color=scoef), 
              method='loess', span=0.05) +
  facet_wrap(vars(Nw)) +
  scale_color_viridis(option='G', trans='log10', end = 0.9, name = 'Selection \ncoefficient')
outlierprop / woutlierprop + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')


ggsave('figures/SUPP-poweranalysisoutlierprop.pdf', width=8, height=4, units='in')

# Fig. X: joint frequency-area distribution, split by inversion
allsnpsplot <- within(allsnpsplot, 
                      INV <- factor(INV, levels=c('none','2La','2Rb','2Rc')))
outlier_notIG <- within(outlier_notIG, 
                      INV <- factor(INV, levels=c('none','2La','2Rb','2Rc')))
missenseoutliersplot <- within(missenseoutliersplot, 
                        INV <- factor(INV, levels=c('none','2La','2Rb','2Rc')))


ag_jdist <- ggplot() + theme_minimal() +
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation=='none' & allsnpsplot$outlier==FALSE,], 
             mapping=aes(x=frequency, y=area), alpha=0.01, size=0.5, pch=20) + 
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation!='none' & allsnpsplot$outlier==FALSE,], 
             mapping=aes(x=frequency, y=area, color=plot_annotation), alpha=0.025, size=0.5, pch=20) +
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation=='none' & allsnpsplot$outlier==TRUE,], 
             mapping=aes(x=frequency, y=area), alpha=0.05, size=0.5, pch=20) + 
  geom_point(data=allsnpsplot[allsnpsplot$plot_annotation!='none' & allsnpsplot$outlier==TRUE,], 
             mapping=aes(x=frequency, y=area, color=plot_annotation), alpha=0.05, size=1, pch=20) + 
  geom_line(data=allsnpsplot, aes(x=frequency, y=cutoff), linewidth=0.5) + 
  geom_point(data=outlier_notIG %>% dplyr::arrange(impact),
             aes(x=frequency, y=area, color=plot_annotation, size=impact)) +
  scale_color_manual(labels=c('Noncoding variant','Synonymous variant','Nonsynonymous variant'),
                     values=c('#375B69','#2CD38B', '#ADEB9E')) +
  geom_point(data=missenseoutliersplot, aes(x=frequency, y=area, size=impact), shape=21, color='black', fill='#ADEB9E') +
  scale_size_manual(values=c(0.5, 1, 2, 3)) +
  scale_y_continuous(limits=c(0, NA)) +
  facet_wrap(vars(INV))
ag_jdist              
              
ggsave('figures/SUPP-agjdistinversion.pdf', width=8, height=4, units='in')
