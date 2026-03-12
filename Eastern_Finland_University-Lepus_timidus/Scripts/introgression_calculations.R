setwd('D:/Work/lepushybrids/hybridhR/')
library(tidyverse)
library(RColorBrewer)
library(ggforce)
library(ggrepel)

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

#1 sample all chr
sample.data <- data.frame(ID = character(0), dosage = numeric(0), pos = numeric(0), chrom = numeric(0))
#sample <- 'LE_Joe_584'

#Combine for LE
gen <- '1'
ver <- 'parental_inc'
for(chr in 1:23){
  le.dos <- read.table(paste0('../ELAI/',chr,'.LE.',ver,'.gen',gen,'.lc10.uC2.LE_dosages.txt'),
                       sep = ' ', header = F)
  snpinfo <- read.table(paste0('../ELAI/snpinfo/',chr,'.LE.LT.gen',gen,'.lc10.uC2.snpinfo.txt'),
                        sep = '\t', header = T)
  coords <- as.vector(snpinfo$pos)
  samples <- le.dos$V1

  #Format
  le_fin <- converter(le.dos, coords)
  #Reformat to long
  tbl <- le_fin[complete.cases(le_fin), ]
  le.long <- create_dosage_data(tbl, samples)
  
  for(sample in samples){
    le.tobind <- le.long[le.long$ID == sample, ]
    le.tobind$chrom <- chr
    sample.data <- rbind(sample.data, le.tobind)
  }
  print(chr)
}

write_csv(sample.data, paste0('sample.data',ver,'.gen',gen,'.csv'))

smean.data <- data.frame(sample = character(0), dosage = numeric(0))
for(s in samples){
    sample <- s
    dosage <- mean(sample.data[sample.data$ID == s, ]$dosage)
    g.data <- data.frame(sample, dosage)
    smean.data <- rbind(smean.data, g.data)
}
write_csv(smean.data, paste0('smean.data',ver,'.gen',gen,'.csv'))

#COmbine for LT
gen <- 1
ver <- 'parental_inc'
sample.data <- data.frame(ID = character(0), dosage = numeric(0), pos = numeric(0), chrom = numeric(0))
for(chr in 1:23){
  lt.dos <- read.table(paste0('../ELAI/',chr,'.LT.',ver,'.gen',gen,'.lc10.uC2.LE_dosages.txt'),
                       sep = ' ', header = F)
  snpinfo <- read.table(paste0('../ELAI/snpinfo/',chr,'.LT.LE.gen',gen,'.lc10.uC2.snpinfo.txt'),
                        sep = '\t', header = T)
  coords <- as.vector(snpinfo$pos)
  samples <- lt.dos$V1
  
  #Format
  le_fin <- converter(lt.dos, coords)
  #Reformat to long
  tbl <- le_fin[complete.cases(le_fin), ]
  le.long <- create_dosage_data(tbl, samples)
  
  for(sample in samples){
    le.tobind <- le.long[le.long$ID == sample, ]
    le.tobind$chrom <- chr
    sample.data <- rbind(sample.data, le.tobind)
  }
  print(chr)
}
write_csv(sample.data, paste0('sample.data.LT.',ver,'.gen',gen,'.csv'))

ltsmean.data <- data.frame(sample = character(0), dosage = numeric(0))
for(s in samples){
  sample <- s
  dosage <- mean(sample.data[sample.data$ID == s, ]$dosage)
  g.data <- data.frame(sample, dosage)
  ltsmean.data <- rbind(ltsmean.data, g.data)
}

write_csv(ltsmean.data, paste0('LT.smean.data',ver,'.gen',gen,'.csv'))


#Combine with PCA data
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y.v0702.tsv", header = T, sep = '\t')
xaut05.data <- read.table("../autosomes_prefilt/ALL.strict.pca.eigenvec")
xaut05.bound <- left_join(xaut05.data, hy.meta, by = c('V1' = 'Sample'))

ltle.smean.data <- rbind(ltsmean.data, smean.data)
ltle.b <- left_join(xaut05.bound, ltle.smean.data, by = c('V1' = 'sample'))

ltle.b[is.na(ltle.b$dosage), "dosage"] <- 0
ltle.b$dosage[ltle.b$dosage == 0 & ltle.b$Species == "LE"] <- 2
ltle.b$dosage[ltle.b$dosage == 0 & ltle.b$Species == "LE_parental"] <- 2

le.b <- ltle.b[ltle.b$Species == 'LE' | ltle.b$Species == 'LE_F1', ] # | ltle.b$Species == 'LE_parental', ] # | ltle.b$Species == 'LE_F1', ]
lt.b <- ltle.b[ltle.b$Species == 'LT', ]
myPalette <- colorRampPalette(brewer.pal(11, "YlOrRd"))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1.6, 2))
ggplot(data = lt.b) +
  geom_point(mapping = aes(x = V2, y = V3, col = dosage), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  #sc +
  ylab("PC2") + #(2,24%)") +
  xlab("PC1") + #(29,25 %)") +
  geom_text_repel(mapping = aes(x = V2, y = V3, col = dosage, label = V1)) +
  #coord_cartesian(xlim = c(0.065,0.17), ylim = c(-0.3,0.15)) +
  # geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
  #                        label = Group),
  #                    expand = unit(0.5, "mm")) +
  theme(legend.position="bottom")


ggplot(data = le.b) +
  geom_point(mapping = aes(x = V2, y = V3, col=dosage), size = 4, alpha = 0.5)  +
  theme_minimal() +
  scale_colour_stepsn(colours=c("red", "red", "darkorange", "orange", "orange", "black"),
                      breaks=c(1.5,1.7,1.9,1.999),
                      limits=c(1.45,2)) +
  xlab('PC1') +
  ylab('PC2')



ggplot(data = lt.b) +
  geom_point(mapping = aes(x = V2, y = V3, col=dosage), size = 4, alpha = 0.5)  +
  theme_minimal() +
  scale_colour_stepsn(colours=c("black", "orange", "red", "red"), #, "darkorange", "orange", "orange", "black"),
                      breaks=c(0,0,0.01, 0.1,0.3,0.5),
                      limits=c(0,0.6)) +
  xlab('PC1') +
  ylab('PC2') +
  theme(legend.position="bottom")



##############


#Plot
ggplot(allchr.data) +
  geom_point(aes(pos, chrom, col=dosage), alpha=0.8) +
  theme_dark() +
  ggtitle(paste0("LT to LE introgression: ",sample,"")) +
  sc +
  theme(panel.background = element_rect(fill = "lightgrey"))
ggsave(paste0('introgression_',sample,'_gen',gen,'.pdf'), width = 400, height = 100, units = "mm", dpi = 600)

##X-related stuff
x.dos <- read.table(paste0('../ELAI/X.LE.LT.gen',gen,'.lc10.uC2.LE_dosages.txt'),
sep = ' ', header = F)
snpinfo <- read.table(paste0('../ELAI/snpinfo/X.LE.LT.gen',gen,'.lc10.uC2.snpinfo.txt'),
                      sep = '\t', header = T)
coords <- as.vector(snpinfo$pos)
x_fin <- converter(x.dos, coords)
#Reformat to long
tbl <- x_fin[complete.cases(x_fin), ]
x.long <- create_dosage_data(tbl, samples)
x.tobind <- x.long[x.long$ID == sample, ]
x.tobind$chrom <- 'X'
samples <- x.dos$V1
allchr.data <- rbind(sample.data, x.tobind)
