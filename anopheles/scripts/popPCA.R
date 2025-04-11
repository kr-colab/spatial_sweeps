library(pcadapt)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1] # path to bed file

# run pcadapt on a bed file
########################

filename <- read.pcadapt(infile, type='bed')
x <- pcadapt(filename, K=4) 
positions <- read.table(str_replace(infile, '.bed', '.positions'))

pvalues <- x$pvalues
snpPass <- x$pass
pvalues <- pvalues[snpPass]
positions <- positions[snpPass,]
padj <- p.adjust(pvalues, method="bonferroni")

write.table(data.frame(position = positions, pval = pvalues, padj = padj), 
            file = paste0('out/', str_replace(basename(infile), '.bed', '.pcAdapt')), 
            sep='\t', row.names=FALSE)

filename <- read.pcadapt(infile, type='bed')
x <- pcadapt(filename, K=4, LD.clumping = list(size = 500, thr = 0.1)) 
positions <- read.table(str_replace(infile, '.bed', '.positions'))

pvalues <- x$pvalues
snpPass <- x$pass
pvalues <- pvalues[snpPass]
positions <- positions[snpPass,]
padj <- p.adjust(pvalues, method="bonferroni")

write.table(data.frame(position = positions, pval = pvalues, padj = padj), 
            file = paste0('out/', str_replace(basename(infile), '.bed', '.thinned.pcAdapt')), 
            sep='\t', row.names=FALSE)