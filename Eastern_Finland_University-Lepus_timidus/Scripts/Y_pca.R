library(tidyverse)
library(ggrepel)
library(ggforce)

setwd('D:/Work/lepushybrids/hybridhR/')
#PCA
hy.data <- read.table("../chrY/Y2.unfiltered.pca.eigenvec", header=F, sep='\t')
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y_cov.v2.tsv", header = T, sep = '\t')

hy.bound <- left_join(hy.data, hy.meta, c('V1' = 'Sample'))


#fullplot <- 
ggplot(data = hy.bound) +
  geom_point(mapping = aes(x = V2, y = V3, col = Species), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  #geom_text_repel(mapping = aes(x = V2, y = V3, col = Species, label = V1)) +
  #  geom_mark_ellipse(aes(x = V2, y = V3, col = V1,
  #                        label = V1),
  #                    expand = unit(0.5, "mm")) +
  theme(legend.position="bottom")
ggsave('pca_Y.pdf', width =140, height = 100, units = "mm", dpi = 300)

hy.le <- hy.bound[hy.bound$Species == 'LE' | hy.bound$Species == 'LE_parental',]
ggplot(data = hy.le) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  #geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  theme(legend.position="bottom") #+
  geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
                          label = Group),
                      expand = unit(0.5, "mm"))

hy.lt <- hy.bound[hy.bound$Species == 'LT' | hy.bound$Species == 'LT_parental',]
ggplot(data = hy.lt) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  theme(legend.position="bottom") +
  ylim(c(-0.1,0.1)) #+
  #coord_cartesian(ylim = c(-0.1,0.1))#+
  #geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
  #                      label = Group),
  #                  expand = unit(0.5, "mm"))


#___________________________________________
#PCA LE separate
#___________________________________________

hy.data <- read.table("../chrY/Y_LE.geno01.pca.eigenvec", header=F, sep='\t')
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y_cov.v2.tsv", header = T, sep = '\t')

hy.bound <- left_join(hy.data, hy.meta, c('V1' = 'Sample'))


ggplot(data = hy.bound) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = Group), size = 3) +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  theme(legend.position="bottom") #+
  geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
                        label = Group),
                    expand = unit(0.5, "mm"))  

  
#___________________________________________
#PCA LT separate
#___________________________________________
  
hy.data <- read.table("../chrY/Y_LT.geno05.pca.eigenvec", header=F, sep='\t')
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y_cov.v2.tsv", header = T, sep = '\t')
  
hy.bound <- left_join(hy.data, hy.meta, c('V1' = 'Sample'))
  
  
ggplot(data = hy.bound) +
    geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
    geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = Group), size = 3) +
    theme_minimal() +
    ylab("PC2") +
    xlab("PC1") +
    theme(legend.position="bottom") #+
    geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
                        label = Group),
                    expand = unit(0.5, "mm"))  
  