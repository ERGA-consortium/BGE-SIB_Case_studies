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

#DATA
#Read in data
for(i in 1:23) {
  chr <- i
  gen <- '1'
  le.dos <- read.table(paste0('../ELAI/',chr,'.LE.LT.gen',gen,'.lc10.uC2.LE_dosages.txt'),
                       sep = ' ', header = F)
  lt.dos <- read.table(paste0('../ELAI/',chr,'.LE.LT.gen',gen,'.lc10.uC2.LT_dosages.txt'),
                       sep = ' ', header = F)
  snpinfo <- read.table(paste0('../ELAI/snpinfo/',chr,'.LE.LT.gen',gen,'.lc10.uC2.snpinfo.txt'),
                        sep = '\t', header = T)
  coords <- as.vector(snpinfo$pos)
  samples <- le.dos$V1
  
  #Format
  le_fin <- converter(le.dos, coords)
  lt_fin <- converter(lt.dos, coords)
  #Reformat to long
  tbl <- le_fin[complete.cases(le_fin), ]
  le.long <- create_dosage_data(tbl, samples)
  tbl <- lt_fin[complete.cases(lt_fin), ]
  lt.long <- create_dosage_data(tbl, samples)
  #Plots
  #Plot all samples at once
  
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 2))
  ggplot(le.long) +
    geom_point(aes(pos, ID, col=dosage), alpha=0.8) +
    theme_dark() +
    ggtitle(paste0("LT to LE introgression: chr ",chr," gen ",gen,"")) +
    sc +
    theme(panel.background = element_rect(fill = "lightgrey"))
  ggsave(paste0('introgression_chr',chr,'_gen',gen,'.pdf'), width = 400, height = 100, units = "mm", dpi = 600)
}


#1 sample all chr
sample.data <- data.frame(ID = character(0), dosage = numeric(0), pos = numeric(0), chrom = numeric(0))
sample <- 'LE_Joe_1027'
gen <- '1'
for(chr in 1:23){
  le.dos <- read.table(paste0('../ELAI/',chr,'.LE.LT.gen',gen,'.lc10.uC2.LE_dosages.txt'),
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
  le.tobind <- le.long[le.long$ID == sample, ]
  le.tobind$chrom <- chr
  sample.data <- rbind(sample.data, le.tobind)
}
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
write_tsv(allchr.data, 'introgression.dosage_data.combined.tsv')

ggplot(allchr.data) +
  geom_point(aes(pos, chrom, col=dosage), alpha=0.8) +
  theme_dark() +
  ggtitle(paste0("LT to LE introgression: ",sample,"")) +
  sc +
  theme(panel.background = element_rect(fill = "lightgrey"))
ggsave(paste0('introgression_',sample,'_gen',gen,'.pdf'), width = 400, height = 100, units = "mm", dpi = 600)


ggplot(x.long) +
  geom_point(aes(pos, ID, col=dosage), alpha=0.8) +
  theme_dark() +
  ggtitle(paste0("LE to LT introgression: Chr X")) +
  sc +
  theme(panel.background = element_rect(fill = "lightgrey"))

#############################
#Adding fst track
diff <- read.table("../ELAI/weir_fst/2.nonan.fst", sep = '\t', header = T)
w.dif <- read.table("../GENERAL/1.LELT.sites100.o25.csv", sep = ',', header = T)
colnames(diff) <- c ('chr', 'pos', 'fst')

ggplot(le.long) +
  geom_point(aes(pos, ID, col=dosage), alpha=0.8) +
  scale_color_distiller(palette = "YlOrRd") +
  theme_dark() +
  theme(panel.background = element_rect(fill = "lightgrey"))
ggplot(diff) +
  geom_col(aes(pos, chr, col=fst)) +
  scale_color_distiller(palette = "RdYlBu", direction = -1) +
  theme_minimal()
ggplot(w.dif) +
  geom_line(aes(x = mid, y = Fst_LE_LT, col = Fst_LE_LT), linewidth = 2) +
  theme_minimal() +
  scale_color_distiller(palette = "RdYlBu")


########################
#Plot only some samples
ggplot(bound)+
  geom_point(aes(pos, LE_Hel_477), col = 'green') +
  geom_point(aes(pos, LE_Joe_1027), col = 'blue') +
  geom_point(aes(pos, LE_Joe_228), col = 'yellow') +
  geom_point(aes(pos, LE_Val_34), col = 'orange') +
  geom_point(aes(pos, LE_Joe_377), col = 'cyan') +
  geom_point(aes(pos, LE_Joe_38), col = 'tan') +
  theme_minimal()


#Loop plotting all samples one by one
colNames <- samples
for(i in colNames){
  plt <- ggplot(le_fin, aes(x=pos, y=.data[[i]])) +
    geom_point(color="#B20000", alpha=0.5)
  #geom_hline(yintercept=0, size=0.06, color="black") + 
  #geom_smooth(method=lm, alpha=0.25, color="black", fill="black")
  print(plt)
}



#######
#Multigen comparison 
#Data
le.dos1 <- read.table('../ELAI/2.LE.LT.gen1.uC2.LE_dosages.txt', sep = ' ', header = F)
le.dos2 <- read.table('../ELAI/2.LE.LT.gen2.uC2.LE_dosages.txt', sep = ' ', header = F)
le.dos5 <- read.table('../ELAI/2.LE.LT.gen5.uC2.LE_dosages.txt', sep = ' ', header = F)
le.dos10 <- read.table('../ELAI/2.LE.LT.gen10.uC2.LE_dosages.txt', sep = ' ', header = F)
#le.dos20 <- read.table('../ELAI/2.LE.LT.gen20.uC2.LE_dosages.txt', sep = ' ', header = F)
snpinfo <- read.table('../ELAI/snpinfo/2.LE.LT.gen1.uC2.snpinfo.txt', sep = '\t', header = T)
coords <- as.vector(snpinfo$pos)

#Conversion
cle.dos1 <- converter(le.dos1, coords)
cle.dos2 <- converter(le.dos2, coords)
cle.dos5 <- converter(le.dos5, coords)
cle.dos10 <- converter(le.dos10, coords)

cled.dos1 <- create_dosage_data(cle.dos1, samples)
cled.dos2 <- create_dosage_data(cle.dos2, samples)
cled.dos5 <- create_dosage_data(cle.dos5, samples)
cled.dos10 <- create_dosage_data(cle.dos10, samples)

cled.dos1$gencount <- "1"
cled.dos2$gencount <- "2"
cled.dos5$gencount <- "5"
cled.dos10$gencount <- "10"

cled.bound <- rbind(cled.dos1, cled.dos2, cled.dos5, cled.dos10)


#Plot
ggplot()+
  geom_point(data=cle.dos1, aes(pos, LE_Oulu_457), col = 'cyan', alpha = 0.5) +
  geom_point(data=cle.dos2, aes(pos, LE_Oulu_457), col = 'lightblue', alpha = 0.5) +
  geom_point(data=cle.dos5, aes(pos, LE_Oulu_457), col = 'blue', alpha = 0.5) +
  geom_point(data=cle.dos10, aes(pos, LE_Oulu_457), col = 'darkblue', alpha = 0.5) +
  theme_minimal()


cled.bound <- cled.bound[cled.bound$ID == 'LE_Oulu_457',]
ggplot(cled.bound) +
  geom_point(aes(pos, gencount, col=dosage), alpha=0.3, size = 4) +
  #scale_color_distiller(palette = "") +
  theme_minimal() #+
  theme(panel.background = element_rect(fill = "lightgrey"))

#Lower cluster count comparison
le.dos1 <- read.table('../ELAI/1.LE.LT.gen1.uC2.LE_dosages.txt', sep = ' ', header = F)
le.dos2 <- read.table('../ELAI/1.LE.LT.gen1.lc10.uC2.LE_dosages.txt', sep = ' ', header = F)
snpinfo <- read.table('../ELAI/snpinfo/1.LE.LT.gen1.uC2.snpinfo.txt', sep = '\t', header = T)
coords <- as.vector(snpinfo$pos)

#Conversion
cle.dos1 <- converter(le.dos1, coords)
cle.dos2 <- converter(le.dos2, coords)

#Plot
ggplot()+
  geom_point(data=cle.dos1, aes(pos, LE_Oulu_457), col = 'blue', alpha=0.2) +
  geom_point(data=cle.dos2, aes(pos, LE_Oulu_457), col = 'red') +
  theme_minimal()

for(i in samples){
  plt <- ggplot() +
    geom_point(data=cle.dos1, aes(x=pos, y=.data[[i]]),
               color="#B20000", alpha=0.2) +
    geom_point(data=cle.dos2, aes(x=pos, y=.data[[i]]),
               color="blue", alpha=0.2) +
    theme_minimal()
  print(plt)
}
