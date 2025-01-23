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
  geom_point(data=snpdf, aes(x=freq, y=area), color=mako[4], alpha=0.01, pch=20, size=0.1) + 
  geom_point(data=snpdf[snpdf$scoef>0,], aes(x=freq, y=area), pch=21, fill=mako[13], size=2, stroke=0.75) +
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
snphist$snps <- snphist$outlier + snphist$nonoutlier
snphist$statistic <- NA
emp_q <- sum(snphist$outlier) / sum(snphist$snps)
st <- ((snphist$outlier / snphist$snps) - emp_q) / (sqrt((emp_q * (1-emp_q))/snphist$snps))
snphist$statistic <- st
snphist <- snphist %>% mutate(Q = quantile(statistic, probs=0.90, na.rm=TRUE))

pawindows <- ggplot(snphist) + theme_minimal(base_size = 8) +
  geom_line(aes(x=mids, y=statistic), color=mako[5], linewidth=0.1, alpha=0.9) + 
  geom_line(aes(x=mids, y=Q), color=mako[1], alpha=1, size=0.5, lty='11') +
  ylab('Z score') +
  geom_point(aes(x=1e8/2, y=5), fill=mako[13], pch=25, size=2, stroke=0.5) +
  scale_x_continuous(name='Position (Mbp)', breaks=c(0, 2.5e7, 5.0e7, 7.5e7, 10e8),
                     labels=c('0','25','50','75','100')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=6))
pawindows

summ_df <- fread('simulation/power_analysis_cutoff_summary.txt')

outliersplot <- ggplot(summ_df) + theme_minimal(base_size = 8) + 
  geom_line(aes(x = cutoff, y=outlier_SNPdetected, color=scoef, group=scoef)) +
  geom_line(aes(x = cutoff, y = outliers/N_mutations), color='gray', lty='21') +
  scale_color_viridis(name='Selection \ncoefficient', option='G', trans='log10', end = 0.9) +
  facet_grid(cols=vars(Nw)) +
  scale_x_reverse() +
  xlab('') +
  ylab('P(SF outlier detected)') +
  ggtitle("Neighborhood Size") +
  theme(plot.title = element_text(hjust = 0.5, size=8),
        axis.title.y = element_text(size=8))

windowsplot <- ggplot(summ_df) + theme_minimal(base_size = 8) +
  geom_line(aes(x = cutoff, y=windowoutlier_SNPdetected, color=scoef, group=scoef)) +
  geom_line(aes(x = cutoff, y = windowoutliers/N_mutations), color='gray', lty='21') +
  scale_color_viridis(name='Selection \ncoefficient', option='G', trans='log10', end=0.9) + 
  scale_y_continuous(limits=c(0,1)) +
  facet_grid(cols=vars(Nw)) +
  scale_x_reverse() +
  xlab('Cutoff') +
  ylab('P(WSF outlier detected)') +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size=8))

outliersplot / windowsplot

(pajdist + 
    inset_element(pawindows, 
                  0.2, 0.01, 1, 0.45))  + plot_annotation(tag_levels = 'A') | 
((outliersplot + theme(panel.spacing = unit(0.5, "cm"))) / 
      (windowsplot+ theme(panel.spacing = unit(0.5, "cm")))) + 
         plot_layout(guides='collect') + theme(panel.spacing = unit(0.5, "cm")) +
  plot_layout(widths = c(1, 1.5))




ggsave('figures/poweranalysisres.pdf', width=8, height=4, units='in')
