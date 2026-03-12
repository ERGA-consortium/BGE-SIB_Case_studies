library(tidyverse)
library(dplyr)
library(tibble)
library(purrr)
library(ggthemes)
library(ggrepel)
library(ggforce)

setwd('D:/Work/lepushybrids/hybridhR/')
hy.meta <- read.table("../hybridhR/Hybrid.sex_and_Y.v0702.tsv", header = T, sep = '\t')
xf05.data <- read.table("../chrx/X.fem.r7.g025.pca.eigenvec")

xf05.bound <- left_join(xf05.data, hy.meta, by = c('V1' = 'Sample'))

ggplot(data = xf05.bound) +
  geom_point(mapping = aes(x = V2, y = V3, col = Species), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  #geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  #  geom_mark_ellipse(aes(x = V2, y = V3, col = V1,
  #                        label = V1),
  #                    expand = unit(0.5, "mm")) +
  theme(legend.position="bottom")



xf75.data <- read.table("../chrx/X.fem.075.pca.eigenvec")
xf75.bound <- left_join(xf05.data, hy.meta, by = c('V1' = 'Sample'))
xf05.bound <- xf05.bound[xf05.bound$Species == 'LT_parental' | xf05.bound$Species == 'LT' |
                         xf05.bound$Species == 'LE_parental' | xf05.bound$Species == 'LE', ]

ggplot(data = xf05.bound) +
  geom_point(mapping = aes(x = V2, y = V3, col = Species), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  theme(legend.position="bottom")
ggsave('pca_X.pdf', width = 140, height = 100, units = "mm", dpi = 300)

xf05.lt <- xf05.bound[xf05.bound$Species == 'LT_parental' | xf05.bound$Species == 'LT', ]
xf05.le <- xf05.bound[xf05.bound$Species == 'LE_parental' | xf05.bound$Species == 'LE' ]

ggplot(data = xf05.bound) +
  geom_point(mapping = aes(x = V2, y = V3, col = Group), size = 2.5, alpha = 0.5)  +
  theme_minimal() +
  ylab("PC2") +
  xlab("PC1") +
  #geom_mark_ellipse(aes(x = V2, y = V3, col = Group,
  #                  label = Group),
  #                  expand = unit(0.5, "mm")) +
  #geom_text_repel(mapping = aes(x = V2, y = V3, col = Group, label = V1)) +
  theme(legend.position="bottom")
  #theme(legend.position="none")




#____________________________________

xf.tbl <- read.table("../chrx/X.fem.biall.05.rem.2.Q")

pop.names = data.frame(c("LE", "LT"))
colnames(pop.names) <- "V7"

famfile <- read.table("../chrx/X.fem.biall.05.rem.fam")

#Format it
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

K2 <- create_admix_data(xf.tbl, famfile)
K2 <- left_join(K2, hy.meta, by = c('ID' = 'Sample'))
head(K2)

tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")


ggplot(K2, aes(factor(ID), value, fill = factor(pop.group))) +
  geom_col(color = "gray", linewidth = 0.1) +
#  facet_grid(~fct_relevel(genot,"LER", "LE", "LT", "LTM"),
#             switch = "x", scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_text(angle = 90),
    panel.grid = element_blank(),
    legend.position = 'none',
    strip.text.x = element_text(angle = 90, size = 12)
  ) +
  scale_fill_manual(values = tol7qualitative)
