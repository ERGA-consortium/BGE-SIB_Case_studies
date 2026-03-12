library(tidyverse)
library(gridExtra)


fin_eu_fst <- read.table('../FST/LE_fin.LE_par.windowed.weir.fst',
                         sep = '\t', header = T)
fin_eu_fst$WEIGHTED_FST[fin_eu_fst$WEIGHTED_FST < 0] <- 0
chroms <- unique(fin_eu_fst$CHROM)[c(1:21,26,27)]

lepus_cordsum <- function(df, binsize){
  endcoord = 0
  for(i in chroms){
    df$sumcoord[df$CHROM == i] <- 
      df$BIN_START[df$CHROM == i] + endcoord
    endcoord <- endcoord + max(df$BIN_START[df$CHROM == i]) + binsize
  }
  autosomes <- df[df$CHROM %in% chroms, ]
  return(autosomes)
}

autosomes <- lepus_cordsum(fin_eu_fst, 100000)
color_mapping <- rep(c("#762a83", "#1b7837"),23)
aut_p <- ggplot(autosomes, aes(x=sumcoord, y=WEIGHTED_FST, color=CHROM, alpha = 0.6)) +
  geom_point() +
  scale_color_manual(values = color_mapping,
                     guide = "none") +
  theme_minimal()
aut_p

introg <- read.table('introgression_perc95.LE.LT.bed',
                     sep = '\t', header = F)
colnames(introg) <- c('chrom', 'startpos', 'endpos')
fai <- read.table('mLepEur1.pri_genomic.fa.fai',
                  sep = '\t', header = F)
colnames(fai) <- c('CHROM', 'LEN', 'V3', 'V4', 'V5')

nchrom <- 1:23
endcoord = 0
for(i in 1:23){
  introg$startsumcoord[introg$chrom == i] <- 
    introg$startpos[introg$chrom == i] + endcoord
  endcoord <- endcoord + fai$LEN[fai$CHROM == fai$CHROM[i]]
}



chrbreaks <- c(0)
chrstart=0
for(i in chroms){
  chrstart <- (chrstart + fai$LEN[fai$CHROM == i])
  chrbreaks[i] <- chrstart
}


introg$startsumcoord <- introg$startpos + chrbreaks[introg$chrom]
n = 0
for(i in chroms) {
  autosomes$chromnum[autosomes$CHROM == i] <- 1 + n
  n = n+1
}
autosomes$sumcoord <- autosomes$BIN_START + 25000 + chrbreaks[autosomes$chromnum]

ggplot() +
  geom_vline(data=introg, aes(xintercept=startsumcoord, alpha = 0.4), color = '#FF7010') +
  geom_point(data=autosomes, aes(x=sumcoord, y=WEIGHTED_FST, color=CHROM, alpha = 0.6)) +
  scale_color_manual(values = color_mapping,
                     guide = "none") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab('Genomic position') +
  scale_x_continuous(breaks=chrbreaks[1:23],
                     labels=c(1:23))
