library(tidyverse)
library(pheatmap)

df <- read.table("../FST/ALL_subpop_fst_summary.tsv",
                 sep = '\t', header = F)
df$V4[df$V4 < 0] <- 0

#Exclude LE_Pun (only 1 indiv in subpop)
df <- df[!(df$V1 %in% "LE_Pun"),]
df <- df[!(df$V2 %in% "LE_Pun"),]

allp_m <- df %>%
  group_by(V1, V2) %>%
  summarise(V4 = mean(V4)) %>%
  ungroup() %>%

  spread(V2, V4, fill = 0) %>%
  
  # convert to matrix, X in row names, Y in column names
  remove_rownames() %>%
  column_to_rownames("V1") %>%
  as.matrix()

pdf('subpop_heatmap.pdf', width = 10, height = 10)
pheatmap(allp_m)
dev.off()
