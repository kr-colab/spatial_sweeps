setwd("~/phd/research/spatial_sweeps")
library(ggplot2)
library(data.table)
library(scales)
library(patchwork)
library(viridis)
library(ggnewscale)

plasma_gen = viridis_pal(option='C')
plasma = plasma_gen(15)
mako_gen = viridis_pal(option='G')
mako = mako_gen(15)

# joint distribution
allsnps <- fread('anopheles/out/all_snps_genome_scan_annotated.txt')
# annotate plot based on 
allsnps$plot_annotation <- 'none'
allsnps[allsnps$annotation == 'intergenic_region' | 
          allsnps$annotation == 'non_coding_transcript_exon_variant' |
          allsnps$annotation == 'intragenic_variant','plot_annotation'] <- 'Intergenic variant'
allsnps[(allsnps$annotation == '3_prime_UTR_variant' |
           allsnps$annotation == '5_prime_UTR_variant' |
           allsnps$annotation == 'intron_variant' |
           allsnps$annotation == 'downstream_gene_variant' |
           allsnps$annotation == 'upstream_gene_variant'),'plot_annotation'] <- 'Noncoding variant'
allsnps[(allsnps$annotation == 'missense_variant' |
           allsnps$annotation == '5_prime_UTR_premature_start_codon_gain_variant' |
           allsnps$annotation == 'splice_region_variant&intron_variant' |
           allsnps$annotation == 'stop_gained' |
           allsnps$annotation == 'splice_region_variant&synonymous_variant' |
           allsnps$annotation == 'missense_variant&splice_region_variant' |
           allsnps$annotation == 'splice_region_variant' |
           allsnps$annotation == 'stop_lost&splice_region_variant' |
           allsnps$annotation == 'splice_acceptor_variant&intron_variant' |
           allsnps$annotation == 'splice_donor_variant&intron_variant' |
           allsnps$annotation == 'stop_lost' |
           allsnps$annotation == 'start_lost' |
           allsnps$annotation == 'splice_region_variant&stop_retained_variant' | 
           allsnps$annotation == 'stop_gained&splice_region_variant'),'plot_annotation'] <- 'Nonsynonymous variant'
allsnps[(allsnps$annotation == 'synonymous_variant'),'plot_annotation'] <- 'Synonymous variant'
allsnps$plot_annotation <- factor(allsnps$plot_annotation, levels=rev(c('Intergenic variant','Noncoding variant','Synonymous variant','Nonsynonymous variant')))

allsnps$INV <- factor(allsnps$INV, levels=c('2La','2Rb','2Rc','none'))

allsnps$alpha <- (1/sqrt(abs(0.75-allsnps$frequency)**2 + abs(0.25-rescale(allsnps$area))**2))
allsnps$size <- allsnps$alpha
allsnps[allsnps$plot_annotation == 'Nonsynonymous variant',]$alpha <- allsnps[allsnps$plot_annotation == 'Nonsynonymous variant',]$alpha * 2

cyp <- fread('anopheles/out/AGAP008552_outliers.csv')
cyp$Vignette <- 'Cytochrome P450 CYP4H2'
cub <- fread('anopheles/out/AGAP029542_outliers.csv')
cub$Vignette <- 'CUB domain-containing protein'
grs <- fread('anopheles/out/NS-GR_outliers.csv')
grs$frequency <- grs$frequency.x
grs$area <- grs$area.x
grs$Vignette <- 'Gustatory receptor'
vgs <- rbind(cyp[,c('frequency','area','Vignette')],
             cub[,c('frequency','area','Vignette')],
             grs[,c('frequency','area','Vignette')])

allsnps <- allsnps[allsnps$INV != 'none',]

for (i in unique(allsnps$INV)){

  assign(paste0('jdist_',i), ggplot() + theme_classic() +
    
    # plot all snps
    geom_point(data=allsnps[allsnps$INV == i,],
               aes(x=frequency, y=area,
                   alpha=alpha,
                   color=plot_annotation,
                   size=size)) +
    scale_alpha_continuous(range=c(0.01, 1)) +
    scale_size_continuous(range=c(0.01, 1)) +
    scale_color_manual(name='Annotation', values=rev(c(mako[1], mako[7], mako[10], mako[14]))) + 
    guides(alpha = "none", size="none") +
    
    # plot cutoff lines
    geom_line(data=allsnps[allsnps$INV==i,], aes(x=frequency, y=cutoff, group=INV), linewidth=0.5, color='white') + 
    geom_line(data=allsnps[allsnps$INV==i,], aes(x=frequency, y=cutoff, group=INV), linewidth=0.5) + 
    #scale_linetype_manual(name='Inversion', values=c('11', '31', '62', 'solid')) +
    new_scale('alpha') +
    new_scale('size') +
    
    # plot NS outliers
    geom_point(data=allsnps[allsnps$outlier & allsnps$plot_annotation=='Nonsynonymous variant' & allsnps$INV == i,],
               aes(x=frequency, y=area,
                   alpha=0.9,
                   size=size*2),
               color=mako[14]) +
    # scale_alpha_continuous(range=c(0.1, 1)) +
    scale_size_continuous(range=c(0.1, 2), trans='log') +
    #scale_fill_manual(values=rev(c(mako[4], mako[6], mako[10], mako[14]))) +
    guides(size="none", alpha="none") +
    scale_y_continuous(limits=c(0, NA)) +
    
    new_scale_color() +
    # plot vignettes!
    geom_point(data=vgs[vgs$INV == i,],
               aes(x=frequency, y=area, fill=Vignette),
               size=2, pch=21) +
    scale_fill_viridis(option='C', discrete = TRUE) +
    
    # axis labels
    xlab('Frequency') +
    ylab(expression(Area~(km^2),)))
}


jdist_2La + jdist_2Rb + jdist_2Rc + 
  plot_annotation(tag_levels = c('A')) + 
  plot_layout(guides = "collect", axis_titles = "collect") & theme(legend.position = 'bottom')
ggsave('figures/supplement-inversionjdist.pdf', width=8, height=4, units='in')

