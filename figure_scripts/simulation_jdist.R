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
  df$hull <- FALSE
  row.names(df) <- NULL
  
  unmutated_df <- df[df$mutated==FALSE,]
  mutated_df <- df[df$mutated==TRUE,]
  row.names(mutated_df) <- NULL
  area_hull <- chull(mutated_df$x, mutated_df$y)
  area_hull <- c(area_hull, area_hull[1])
  area_hull_coords <- mutated_df[area_hull,c('x','y')]
  area_hull_poly <- Polygon(area_hull_coords, hole=F)
  print(c(tick, unique(df$frequency), area_hull_poly@area))
  mutated_df[area_hull,]$hull <- TRUE
  # some weirdness to get hulls graphed right...
  mutated_df_ <- mutated_df[mutated_df$hull==FALSE,]
  coords <- coordinates(area_hull_poly)
  for (i in (as.integer(rownames(coords)))){
    mutated_df_ <- rbind(mutated_df_, mutated_df[i,])
  }
  
  df <- rbind(unmutated_df, mutated_df_)
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
  geom_polygon(data=simdf[simdf$hull,], aes(x=x, y=y), color=mako[7], fill=mako[7], alpha=0.75, linewidth=0.5) +
  geom_point(data=simdf[simdf$mutated,], aes(x=x, y=y), color=mako[13], size=1, inherit.aes = FALSE) +
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

(simpt + theme(panel.spacing = unit(0.5, "cm"))) / (jdist  + theme(panel.spacing = unit(1, "cm"))) + plot_annotation(tag_levels = 'A')
ggsave('figures/sweepsimjointdistribution.pdf', width=8, height=6, units='in')
