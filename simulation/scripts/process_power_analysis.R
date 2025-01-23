library(scales)
#library(quantreg)
library(plyr)
library(tidyverse)
library(data.table)

padf <- fread('power_analysis_sliding_cutoff.txt', fill=TRUE)
padf <- data.table(padf)
# sort weird R printed output into separate columns
columns <- colnames(padf)[7:length(padf)]
padf <- padf %>% dplyr::select(colnames(padf)[1:9])
padf <- padf %>% tidyr::separate(frequency, into = columns, sep=" ")
padf <- padf %>% drop_na()
padf <- mutate_all(padf, function(x) as.numeric(as.character(x)))

##line up ticks across simulations (sweep started within ~3 ticks)
roundUp <- function(x,m) m*ceiling(x / m)
align_tick <- function(scoef, treetick, Nw, padf){
  sims <- unique(padf[(padf$scoef == scoef & padf$Nw == Nw),]$sim)
  for (s in sims){
      tix <- padf[(padf$scoef == scoef & padf$Nw == Nw & padf$sim == s),]$tick
      padf[(padf$scoef == scoef & padf$Nw == Nw & padf$sim == s),]$tick <- roundUp(tix, treetick)
    }
  return(padf)
}
for (Nw in c(10,100,1000)){
  padf <- align_tick(0.01, 2, Nw, padf)
  padf <- align_tick(0.001, 5, Nw, padf)
  padf <- align_tick(0.0001, 10, Nw, padf)
}

# summarize by cutoff
summ_df <- padf %>% group_by(scoef, Nw, cutoff) %>% dplyr::summarize(area = mean(area),
                                                                        distance=mean(distance),
                                                                   frequency=mean(frequency),
                                                                   N_mutations=mean(N_mutations),
                                                                   outliers=mean(outliers),
                                                      outlier_SNPdetected=mean(outlier_SNPdetected),
                                                      windowoutliers=mean(windowoutliers),
                                                      windowoutlier_SNPdetected=mean(windowoutlier_SNPdetected),
                                                      windowoutlier_nonoutlier_SNPdetected=mean(windowoutlier_nonoutlier_SNPdetected))

fwrite(summ_df, 'power_analysis_cutoff_summary.txt')

# summarize power across different frequencies
summ_df <- padf %>% group_by(scoef, Nw, tick, cutoff) %>% dplyr::summarize(area = mean(area),
                                                                        distance=mean(distance),
                                                                   frequency=mean(frequency),
                                                                   N_mutations=mean(N_mutations),
                                                                   outliers=mean(outliers),
                                                      outlier_SNPdetected=mean(outlier_SNPdetected),
                                                      windowoutliers=mean(windowoutliers),
                                                      windowoutlier_SNPdetected=mean(windowoutlier_SNPdetected),
                                                      windowoutlier_nonoutlier_SNPdetected=mean(windowoutlier_nonoutlier_SNPdetected))
for (ct in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)){
    fwrite(summ_df[summ_df$cutoff == ct,], paste0('power_analysis_frequency_cutoff',ct,'.txt'))
}