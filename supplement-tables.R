setwd("~/phd/research/spatial_sweeps/")
library(tidyverse)
library(data.table)
library(stats)

# outlier snps table
allsnps <- fread('anopheles/out/all_snps_genome_scan_annotated.txt')
allsnps <- allsnps[allsnps$outlier,]
outliersnps <- allsnps %>% select(c('chromosome',
                     'position',
                     'alternate',
                     'frequency',
                     'area',
                     'annotation',
                     'impact',
                     'Gene ID',
                     'Product Description',
                     'INV',
                     'nstat_outlier',
                     'nstat')) %>%
                rename(!!c('alternate_allele'='alternate',
                     'SNPEff_annotation'='annotation',
                     'SNPEff_impact'='impact',
                     'Ensembl_gene_id'='Gene ID',
                     'product_description'='Product Description',
                     'inversion'='INV',
                     'window_outlier'='nstat_outlier',
                     'window_zscore'='nstat'))
fwrite(outliersnps, 'supp-data/outlier_snps.txt', sep='\t')

# missense snps table
nssnps <- fread('anopheles/out/missense_snps_stabilitychange.txt')
nssnps <- nssnps  %>% select(c('chromosome',
                               'position',
                               'alternate',
                               'frequency',
                               'area',
                               'annotation',
                               'impact',
                               'Gene ID',
                               'source_id',
                               'Product Description',
                               'ancestral_allele',
                               'derived_allele',
                               'Mutation',
                               'ddG (kcal/mol)',
                               'UniProt ID(s)',
                               'INV',
                               'nstat_outlier',
                               'nstat')) %>%
                      rename(!!c('alternate_allele'='alternate',
                                 'SNPEff_annotation'='annotation',
                                 'SNPEff_impact'='impact',
                                 'Ensembl_gene_id'='Gene ID',
                                 'transcript'='source_id',
                                 'product_description'='Product Description',
                                 'wt_NT' = 'ancestral_allele',
                                 'mut_NT' = 'derived_allele',
                                 'AA_change' = 'Mutation',
                                 'ddG_kcal_mol' = 'ddG (kcal/mol)',
                                 'uniprot_id' = 'UniProt ID(s)',
                                 'inversion'='INV',
                                 'window_outlier'='nstat_outlier',
                                 'window_zscore'='nstat'))
fwrite(nssnps, 'supp-data/nonsynonymous_snps.txt', sep='\t')
