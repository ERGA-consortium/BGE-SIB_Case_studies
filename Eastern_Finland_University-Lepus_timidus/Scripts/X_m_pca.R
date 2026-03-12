library(dplyr)
library(tibble)
library(purrr)
library(ggthemes)
library(ggrepel)
library(ggforce)

setwd('D:/Work/lepushybrids/hybridhR/')
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y.v0702.tsv", header = T, sep = '\t')
xm05.data <- read.table("../chrx/X.m.r9.g025.pca.eigenvec")

xm05.bound <- left_join(xm05.data, hy.meta, by = c('V1' = 'Sample'))

ggplot(data = xm05.bound) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  #  geom_mark_ellipse(aes(x = V2, y = V3, col = V1,
  #                        label = V1),
  #                    expand = unit(0.5, "mm")) +
  theme(legend.position="bottom")


#Plot per species
xm05.lt <- xm05.bound[xm05.bound$Species == 'LT_parental' | xm05.bound$Species == 'LT', ]
xm05.le <- xm05.bound[xm05.bound$Species == 'LE_parental' | xm05.bound$Species == 'LE' |
                        xm05.bound$Species == 'LE_F1', ]

xm05.oulu <- xm05.bound[xm05.bound$Group == 'LE_Oulu', ]
ggplot(data = xm05.le) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  theme(legend.position="bottom")
