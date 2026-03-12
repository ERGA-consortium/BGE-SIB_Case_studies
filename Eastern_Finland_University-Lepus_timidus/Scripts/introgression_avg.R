setwd('D:/Work/lepushybrids/hybridhR/')
library(tidyverse)
library(RColorBrewer)

#FUNCTIONS
#This converts the raw input to something plottable
converter <- function(indata, coords){
  df <- as.data.frame(t(as.matrix(indata[ ,2:ncol(indata)])))
  samples <- indata[,1]
  colnames(df) <- c(samples)
  df <- df[1:(nrow(df)-1), ]
  df$pos <- coords
  return(df)
}

#Another converter step, to a long format
#Reason and goal: to easily plot all samples at once
create_dosage_data <- function(tbl, samples) {
  admix.data <- data.frame(ID = character(0), dosage = numeric(0), pos = numeric(0))
  datalist = list()
  for(i in 1:(ncol(tbl)-1)) {
    ID <- rep(samples[i], each=nrow(tbl))
    pos <- tbl[ , ncol(tbl)]
    dosage <- tbl[ , i]
    datalist[[i]] <- data.frame(ID, dosage, pos)
  }
  admix.data <- do.call(rbind, datalist)
  return(admix.data)
}

#Mean dosage calculation
calc_chr_mean <- function(chr.data){
  chrm <- c()
  for(i in 2:ncol(chr.data)){
    m <- mean(chr.data[ ,i])
    chrm <- c(chrm, m)
  }
  return(chrm)
}
calc_chr_mean(chr.data)


#CALC
#Do it for all chr and gather into one dataframe
gen <- 100
datalist = list()
for(i in 1:23){
  print(i)
  #Get data
  chr.data <-  read.table(paste0("../ELAI/",i,".LE.LT.gen",gen,".lc10.uC2.LE_dosages.txt"),
                                 sep = ' ', header = F)
  pos.data <- read.table(paste0('../ELAI/snpinfo/',i,'.LE.LT.gen',gen,'.lc10.uC2.snpinfo.txt'),
                                    sep = '\t', header = T)
  coords <- as.vector(pos.data$pos)
  
  #Columns for a table
  dosage <- calc_chr_mean(chr.data)
  dosage <- dosage[1:(length(dosage)-1)]
  chrom <- rep(i, length(dosage))
  pos <- coords
  #Gather all data
  datalist[[i]] <- data.frame(chrom, pos, dosage)
}
#Combine all data into single dataframe
avg.data <- do.call(rbind, datalist)
head(avg.data)
#write_tsv(avg.data, "mean_dosages.LT.LE.tsv")


#Plot it
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1.5, 2))
avg.data$ltdosage <- 2-(avg.data$dosage)
ggplot(avg.data) +
  geom_point(aes(pos, chrom, col=dosage), alpha=0.8) +
  theme_dark() +
  ggtitle(paste0("LT to LE introgression summary, gen",gen)) +
  sc +
  theme(panel.background = element_rect(fill = "lightgrey"))
#ggsave(paste0('introgression_avgs.LE.LT.chrall.gen',gen,'.pdf'), width = 400, height = 100, units = "mm", dpi = 600)

#Color by top 5%
#LE
perc <- sort(avg.data$dosage)[1:(0.05*nrow(avg.data))]
#LT
#perc <- sort(avg.data$dosage, decreasing = T)[1:(0.01*nrow(avg.data))]
percend <- perc[length(perc)-1]
avg.data$percentile <- 'Y'
avg.data$percentile[avg.data$dosage >= percend] <- 'N'

ggplot(avg.data) +
  geom_point(aes(pos, chrom, col=percentile), alpha=0.8) +
  theme_minimal() +
  ggtitle(paste0("LT to LE introgression summary, gen",gen)) +
  scale_color_manual(values = c('blue', 'darkorange')) +
  theme(legend.position = "none")
ggsave(paste0('introgression_perc95.LE.LT.chrall.gen',gen,'.pdf'), width = 400, height = 100, units = "mm", dpi = 600)

top <- avg.data[avg.data$percentile == 'Y',]
top$pos2 <- top$pos+1

write_tsv(top, paste0("introgression_perc95.gen",gen,".LE.LT.tsv"))




#Complete data re-read-in, if needed, just leaving this here
c.le.data <- read.table("sample.dataparental_inc.gen1.csv", sep = ',', header = T)
