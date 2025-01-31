setwd("~/phd/research/spatial_sweeps/")
library(tidyverse)
library(data.table)
library(patchwork)
library(scales)
library(quantreg)
library(viridis)

plasma_gen = viridis_pal(option='C')
plasma = plasma_gen(15)
mako_gen = viridis_pal(option='G')
mako = mako_gen(15)


Neutral <- c(rep(NA, 40), rescale(rgamma(20, shape=10))/7, rep(NA, 40))
Beneficial <- c(rep(NA, 45), seq(0, 0.5, length.out = 4), 0.6, 0.6, seq(0.5, 0, length.out = 4), rep(NA, 45))
xpos <- seq(100)
df <- data.frame(xpos, Neutral, Beneficial)
df <- df %>% pivot_longer(cols = c('Neutral','Beneficial'), names_to="Allele")
df$Allele <- factor(df$Allele, levels=c('Beneficial', 'Neutral'))
origin <- ggplot(df) + theme_bw() + geom_area(aes(x=xpos, y=value, color=Allele, fill=Allele, group=Allele), position='identity', alpha=0.75) +
  #facet_wrap(~Allele) + 
  scale_color_manual(values=c(mako[10], plasma[1])) +
  scale_fill_manual(values=c(mako[13], plasma[3])) +
  theme(
    strip.background = element_blank(),
    #legend.position='none',
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits=c(0,100))+
  labs(y = 'Frequency')

Neutral <- c(rep(NA, 20), rescale(rgamma(60, shape=10))/7, rep(NA, 20))
Beneficial <- c(rep(NA, 40), seq(0, 1, length.out = 10), seq(1, 0, length.out = 10), rep(NA, 40))
xpos <- seq(100)
df <- data.frame(xpos, Neutral, Beneficial)
df <- df %>% pivot_longer(cols = c('Neutral','Beneficial'), names_to="Allele")
df$Allele <- factor(df$Allele, levels=c('Beneficial', 'Neutral'))
nexts <- ggplot(df) + theme_bw() + geom_area(aes(x=xpos, y=value, color=Allele, fill=Allele, group=Allele), position='identity', alpha=0.75) +
  #facet_wrap(~Allele) + 
  scale_color_manual(values=c(mako[10], plasma[1])) +
  scale_fill_manual(values=c(mako[13], plasma[3])) +
  theme(
    strip.background = element_blank(),
    #legend.position='none',
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits=c(0,100))+
  labs(y = 'Frequency')



Neutral <- rescale(rgamma(100, shape=10))/7
Beneficial <- c(rep(NA, 25), seq(0, 1, length.out = 10), rep(1, 30), seq(1, 0, length.out = 10), rep(NA, 25))
xpos <- seq(100)
df <- data.frame(xpos, Neutral, Beneficial)
df <- df %>% pivot_longer(cols = c('Neutral','Beneficial'), names_to="Allele")
df$Allele <- factor(df$Allele, levels=c('Beneficial', 'Neutral'))
later <- ggplot(df) + theme_bw() + geom_area(aes(x=xpos, y=value, color=Allele, fill=Allele, group=Allele), position='identity', alpha=0.75) +
  #facet_wrap(~Allele) + 
  scale_color_manual(values=c(mako[10], plasma[1])) +
  scale_fill_manual(values=c(mako[13], plasma[3])) +
  theme(
    strip.background = element_blank(),
   # legend.position='none',
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = 'Frequency', x='Position')

origin/nexts/later + plot_annotation(tag_levels = 'A') + plot_layout(guides='collect')
ggsave('figures/allele_conceptviz.pdf', width=8, height=5, units='in')
