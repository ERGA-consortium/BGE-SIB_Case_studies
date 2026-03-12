library(tidyverse)

setwd('D:/Work/lepushybrids/hybridhR/')
yreads <- read.table('../chrY/Lepus_all.Y.txt', sep ='\t', header = F)
head(yreads)
ylen <- 33560193
xlen <- 135183192
readlen <- 150

colnames(yreads) <- c('Sample', 'Readcount', 'Coverage', 'Spbypca')
yreads$Coverage <- (yreads$Readcount * readlen)/ ylen

ggplot(data = yreads) +
  geom_point(aes(x = Spbypca, y = Coverage))  

ggplot(data = yreads) +
  geom_point(aes(x = Spbypca, y = Readcount))  


xyreads <- read.table('../chrY/X_Y_reads_complete.txt', sep ='\t', header = F)
colnames(xyreads) <- c('Sample', 'X_reads', 'Y_reads', 'Species')

xyreads$Read_ratio <- xyreads$Y_reads / xyreads$X_reads
xyreads$Sex <- "F"
xyreads$Sex[xyreads$Read_ratio > 0.25] <- "M"
ggplot(data = xyreads) +
  geom_point(aes(x = Species, y = Read_ratio, color = Sex)) +
  theme_minimal()



xyreads$Y_cov <- xyreads$Y_reads * 150 / ylen
xyreads$X_cov <- xyreads$X_reads * 150 / xlen
xyreads$Cov_ratio <- xyreads$X_cov / xyreads$Y_cov
ggplot(data = xyreads) +
  geom_point(aes(x = Species, y = Cov_ratio, color = Sex)) +
  theme_minimal()


nrow(xyreads)
malehares <- nrow(xyreads[xyreads$Sex == 'M',])
femalehares <- nrow(xyreads[xyreads$Sex == 'F', ])
malehares/nrow(xyreads)
femalehares/nrow(xyreads)
