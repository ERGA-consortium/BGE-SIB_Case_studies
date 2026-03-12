library(tidyverse)
library(dplyr)
library(tibble)
library(purrr)
library(ggthemes)
library(patchwork)
library(forcats)
library(ggrepel)
library(ggforce)

setwd('D:/Work/lepushybrids/hybridhR/')
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y.v0702.tsv", header = T, sep = '\t')
#This is the original test:
xaut05.data <- read.table("../pca/hybrid.fam/All.merge.r16.g025.pca.eigenvec")

#This is from the prefiltered vcf:
xaut05.data <- read.table("../autosomes_prefilt/ALL.strict.pca.eigenvec")
xaut05.bound <- left_join(xaut05.data, hy.meta, by = c('V1' = 'Sample'))
xaut05.bound$Species[xaut05.bound$V1 == 'LE_Joe_584'] <- "LE_F1"


xaut05.lt <- xaut05.bound[xaut05.bound$Species == 'LT' | xaut05.bound$Species == 'LT_parental', ]
xaut05.le <- xaut05.bound[xaut05.bound$Species == 'LE' , ] # | xaut05.bound$Species == 'LE_parental', ]

ggplot(data = xaut05.lt) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") + #(2,24%)") +
  xlab("PC1") + #(29,25 %)") +
  #geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  #coord_cartesian(xlim = c(0.065,0.17), ylim = c(-0.3,0.15)) +
  theme(legend.position="bottom")
#ggsave(paste0('PCA_',sample,'.png'), width = 600, height = 800, units = "mm", dpi = 600)

ggplot(data = xaut05.bound) +
  geom_point(mapping = aes(x = V3, y = V4, col = Species, shape = Species), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC3 (0,8 %)") +
  xlab("PC2 (2,24%)") +
  geom_text_repel(mapping = aes(x = V3, y = V4, col = Group, label = V1)) +
  #coord_cartesian(xlim = c(0.065,0.17), ylim = c(-0.3,0.15)) +
  #geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
  #                        label = Group),
  #                    expand = unit(0.5, "mm")) +
  theme(legend.position="bottom")


ltaut05.data <- read.table("../autosomes/LT.r9.g025.pca.eigenvec")
ltaut05.bound <- left_join(ltaut05.data, hy.meta, by = c('V1' = 'Sample'))
leaut05.data <- read.table("../autosomes/LE.r8.g025.pca.eigenvec")
leaut05.bound <- left_join(leaut05.data, hy.meta, by = c('V1' = 'Sample'))

lefin <- leaut05.bound[leaut05.bound$Species == "LE", ]
ltfin <- ltaut05.bound[ltaut05.bound$Species == 'LT', ]

ggplot(data = lefin) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  #coord_cartesian(xlim = c(0.065,0.17), ylim = c(-0.3,0.15)) +
  #geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
  #                        label = Group),
  #                    expand = unit(0.5, "mm")) +
  theme(legend.position="bottom")


r16g025.eigval <- read.table("All.merge.r16.g025.repca.eigenval", sep = '\t', header = F)
r16g025.eigval[1, ] / sum(r16g025.eigval)
r16g025.eigval[2, ] / sum(r16g025.eigval)
r16g025.eigval[3, ] / sum(r16g025.eigval)
lt.eigval <- read.table("LT.repca.eigenval", sep = '\t', header = F)
lt.eigval[1, ] / sum(lt.eigval)
lt.eigval[2, ] / sum(lt.eigval)
le.eigval <- read.table("LE.repca.eigenval", sep = '\t', header = F)
le.eigval[1, ] / sum(le.eigval)
le.eigval[2, ] / sum(le.eigval)
