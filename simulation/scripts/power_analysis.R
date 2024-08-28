suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(tidyverse)))
suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(npreg)))
suppressMessages(suppressWarnings(require(splines)))
parser <- argparse::ArgumentParser()
parser$add_argument('--infile')
args <- parser$parse_args()
# what's the power to detect alleles under positive selection in our sims?


snpdf <- fread(args$infile)

#snpdf <- snpdf[snpdf$freq > min(snpdf$freq),]
snpdf <- snpdf[snpdf$area > 0,]

N = 1000
snpdf <- snpdf %>% arrange(by = freq)
bins <- rep(seq(from = 0, to=as.integer(nrow(snpdf)/ N)), each=N)[1:nrow(snpdf)]
snpdf$bin <- bins
snpdf <- snpdf %>% group_by(bin) %>% mutate(Q = quantile(area, probs=0.1),
                                            bin_frequency = mean(freq))
#snpdf <- snpdf %>% group_by(bin_frequency) %>% mutate(Q = mean(Q))
#qb <- snpdf %>% reframe(Q = Q, bin_frequency=bin_frequency) %>% distinct()
#mt <- bs(qb$bin_frequency, df=10)
#lt <- lm(qb$Q ~ mt)
#qb$cutoff <- predict(lt, data.frame(qb$bin_frequency))

#snpdf <- left_join(snpdf, qb)
snpdf$outlier <- FALSE
snpdf[snpdf$area <= snpdf$Q,]$outlier <- TRUE
#ggplot(snpdf) + geom_point(aes(x=freq, y=area, size=scoef, color=outlier)) + 
#  geom_line(aes(x=bin_frequency, y=cutoff)) + scale_y_continuous(limits=c(0,max(snpdf$area)))

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

#ggplot(snphist) + geom_line(aes(x=breaks, y=ratio)) + geom_line(aes(x=breaks, y=Q)) +
#  geom_point(x=1e8/2, y=0.4)

snphist$outlierwindow <- FALSE
snphist[snphist$ratio >= snphist$Q,'outlierwindow'] <- TRUE

snpdf$outlierwindow <- FALSE
for (i in seq(nrow(snphist))){
  startpos <- snphist[i, 'breaks']
  endpos <- startpos + (snphist[1, 'mids'] * 2)
  if (snphist[i, 'outlierwindow']){
    snpdf[snpdf$position >= startpos & snpdf$position <= endpos, 'outlierwindow'] <- TRUE
  }
}

print(paste(nrow(snpdf),
            snpdf[snpdf$scoef>0,]$freq,
            snpdf[snpdf$scoef>0,]$area,
            snpdf[snpdf$scoef>0,]$distance,
            sum(snpdf$outlier == TRUE),
            sum(snpdf$outlier == TRUE & snpdf$scoef > 0),
            sum(snpdf$outlierwindow == TRUE & snpdf$outlier == TRUE),
            sum(snpdf$outlierwindow == TRUE & snpdf$outlier==TRUE & snpdf$scoef > 0),
            sum(snpdf$outlierwindow == TRUE & snpdf$scoef > 0)))




