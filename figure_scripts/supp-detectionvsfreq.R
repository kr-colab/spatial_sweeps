setwd("~/phd/research/spatial_sweeps")
library(tidyverse)
library(data.table)
library(MASS)
library(patchwork)
library(ggiraph)
library(scales)
library(quantreg)
library(viridis)

plasma_gen = viridis_pal(option='C')
plasma = plasma_gen(15)
mako_gen = viridis_pal(option='G')
mako = mako_gen(15)

library(ggnewscale)

summ_df <- fread('simulation/power_analysis_sliding_cutoff.txt')

ggplot(summ_df) + theme_bw() +
  geom_smooth(aes(x=frequency, y=outlier_SNPdetected, color=scoef)) +
  facet_grid(rows=vars(cutoff), cols=vars(Nw))
