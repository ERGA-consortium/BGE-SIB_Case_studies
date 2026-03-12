library(tidyverse)
library(dplyr)
library(tibble)
library(purrr)
library(ggthemes)
library(patchwork)
library(forcats)

setwd('D:/Work/lepushybrids/hybridhR/')
#The data
tbl2 <- read.table("../autosomes_prefilt/ALL.b.strict.4.Q") #/All.b.strict.4.Q")
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y.v0702.tsv", header = T, sep = '\t')


pop.names = data.frame(c("LE", "LT"))
colnames(pop.names) <- "V7"

famfile <- read.table("../autosomes_prefilt/ALL.b.strict.pruned.fam")
#barplot(t(as.matrix(tbl2)), col=rainbow(4),
#        xlab="Individual #", ylab="Ancestry", border=NA)

#Format it
fameta <- left_join(famfile, hy.meta, by = c('V2' = 'Sample'))
create_admix_data <- function(tbl, famfile) {
  admix.data <- data.frame(ID = character(0), value = numeric(0))
  datalist = list()
  for(i in 1:ncol(tbl)) {
    ID <- famfile$V2
    value <- tbl[ , i]
    pop.group <- i-1
    genot <- famfile$V1 #this line is modified from original
    datalist[[i]] <- data.frame(ID, value, pop.group, genot)
  }
  admix.data <- do.call(rbind, datalist)
  return(admix.data)
}

K2 <- create_admix_data(tbl2, fameta)
K2 <- left_join(K2, hy.meta, by = c('ID' = 'Sample'))
head(K2)

#K2$pop.group[K2$Species == 'LT' | K2$Species == 'LT_parental'] <- 2

#Plotting
#Palette
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")

#k2
k2.le <- K2[K2$Species == 'LE' | K2$Species == 'LE_parental' | K2$Species == 'LE_F1', ]
k2.lt <- K2[K2$Species == 'LT' | K2$Species == 'LT_parental', ]

ggplot(k2.le, aes(factor(ID), value, fill = factor(pop.group))) +
  geom_col(color = "gray", linewidth = 0.1) +
  theme_minimal() +
  labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_text(angle = 90),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank(),
    strip.text.x = element_text(angle = 90, size = 12)
  ) +
  scale_fill_manual(values = tol7qualitative,
                    #labels = c('LT_parental', 'LE_FIN', 'LE_parental', 'LT_FIN'))
                    labels = c('LE', 'LT'))


#_________________________________
#Individual ancestry checks
#_________________________________

#For Q4
tbl2 <- read.table("../autosomes_prefilt/ALL.b.strict.4.Q") #/All.b.strict.4.Q")
admixdata <- cbind(famfile, tbl2)[ ,c(2,7,8,9,10)]
colnames(admixdata) <- c('Sample', 'LT_parental', 'LE_FIN', 'LE_parental', 'LT_FIN')
admixdata[admixdata$LE_FIN > 0.999, ]

#For Q2
tbl2 <- read.table("../autosomes_prefilt/ALL.b.strict.2.Q") #/All.b.strict.4.Q")
admixdata <- cbind(famfile, tbl2)[ ,c(2,7,8)]
colnames(admixdata) <- c('Sample', 'LE', 'LT')
admixdata[admixdata$LE > 0.999, ][ ,1]
admixdata[admixdata$LT > 0.999, ][ ,1]
admixdata[admixdata$LT < 0.999 & admixdata$LE < 0.999, ][ ,1]



admixdata[admixdata$Sample == 'LE_Joe_584', ]


