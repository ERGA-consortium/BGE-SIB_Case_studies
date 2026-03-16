setwd("~/Public/ownCloud/Arnica_EUMS/09_Analysis/R_Katja")

library("poppr")
library("ggplot2")
library("reshape2")

### read in the data & view its content
#ArnicaPop <- read.genalex("20250320GenotypicTableWithoutLetters.csv", sep="\t") # KATJA: default from my Excel
#ArnicaPop <- read.genalex("20250221ArnicaSSRlenGenAlEx_NewlyTranslated.csv", sep=";")

Arnica_S <- read.genalex("20250321_GenotypicTable_With.csv", sep=";")
Arnica_L <- read.genalex("20250321_GenotypicTable_Without.csv", sep=";")
Arnica_F <- read.genalex("20250321_GenotypicTable_Fragments.csv", sep=";")

### collect all results in a single dataframe for comparative plotting
ArnicaPlotdata <- data.frame()
counter <- 0
for (k in c(Arnica_S, Arnica_L, Arnica_F)) {
  ArnicaPop <- k
  counter <- counter +1
  #View(ArnicaPop$tab) # crosstable of samples x alleles
  #View(ArnicaPop) # overview of the sub-variables of the ArnicaPop dataframe
  
  ### reformat the data to show presence/absence of each allele
  ArnicaAlleles <- data.frame(ArnicaPop$tab)
  ArnicaAlleles[is.na(ArnicaAlleles)] <- 0
  ArnicaAlleles[ArnicaAlleles == 2] <- 1 # comment out for allele-based stats
  #View(ArnicaAlleles)
  
  ### aggregate data by population
  ArnicaPopAlleles <- aggregate(ArnicaAlleles, list(ArnicaPop$pop), FUN=sum) 
  #View(ArnicaPopAlleles)
  
  ### output allele presence/absence table for the "with letters"/S dataset only
  if (counter == 1) {
    write.table(ArnicaPopAlleles, "ArnicaSSR_4GIS.csv", sep=",", quote=FALSE, row.names=FALSE)
  }
  
  ### reformat the data from crosstable to observation list
  ArnicaList2 <- melt(cbind(ArnicaPopAlleles))
  colnames(ArnicaList2) <- c("Pop", "Allele", "Count")
  #View(ArnicaList2)
  
  ### add a column which specifies the Locus
  ArnicaLoci2 <- list()
  for (i in ArnicaList2$Allele) {
    if (length(unlist(strsplit(i, ".", fixed = TRUE))) ==2) {
      ArnicaLoci2 <- c(ArnicaLoci2,unlist(strsplit(i, ".", fixed = TRUE))[1])
    }
    else {
      parts <- length(unlist(strsplit(i, ".", fixed = TRUE)))
      locname <- paste(unlist(strsplit(i, ".", fixed = TRUE))[1:parts-1], collapse="-")
      #print(locname)
      ArnicaLoci2 <- c(ArnicaLoci2,locname)
    }
  }
  #ArnicaLoci2[1:15]
  ArnicaList2$Locus <- factor(ArnicaLoci2, levels=unique(ArnicaLoci2))
  
  ### cut off the locus names from the allele names
  ArnicaAlls2 <- list()
  for (i in ArnicaList2$Allele) {
    parts <- length(unlist(strsplit(i, ".", fixed = TRUE)))
    ArnicaAlls2 <- c(ArnicaAlls2,unlist(strsplit(i, ".", fixed = TRUE))[parts])
  }
  #ArnicaAlls2[1:15]
  ArnicaList2$Allele <- factor(ArnicaAlls2, levels=unique(ArnicaAlls2))
  
  ### column to specify which study (and dataset: with/without letters, fragments) the data is from
  ArnicaStudies2 <- list()
  for (i in ArnicaList2$Pop) {
    ArnicaStudies2 <- c(ArnicaStudies2,unlist(strsplit(i, "-", fixed = TRUE))[1])
  }
  #View(ArnicaStudies2)
  ArnicaList2$Study <- factor(ArnicaStudies2, levels=unique(ArnicaStudies2))
  
  ### reformat the Pop column
  ArnicaPops2 <- list()
  for (i in ArnicaList2$Pop) {
    ArnicaPops2 <- c(ArnicaPops2,unlist(strsplit(i, "-", fixed = TRUE))[2])
  }
  #View(ArnicaPops2)
  ArnicaList2$Pop <- factor(ArnicaPops2, levels=unique(ArnicaPops2))
  
  ### add the data from the current dataset to the common dataframe
  ArnicaPlotdata <- rbind(ArnicaPlotdata, ArnicaList2)
}
#View(ArnicaPlotdata)

### sort the alleles according to size (! - no leading zeros, so normal sort does not fully work)
AllelesSorted <- c(sort(levels(ArnicaPlotdata$Allele))[236:277], 
                   sort(levels(ArnicaPlotdata$Allele))[1:235])
ArnicaPlotdata$Allele <- factor(ArnicaPlotdata$Allele, levels=AllelesSorted)

### sort the populations by geography & bind AM_09/10 and AM_31/32 together again
PopsSorted <- c("NOR03","NOR01", "DEN03", "SWE01", "AM09", "AM10", "LIT02", "GER07",
                "AM46", "AM30", "NET01", "AM15", "AM304", "GER06", "CZE01", "AM32", "AM31",
                "AM37", "LAG", "AM42", "CAR", "AM40", "AUS01", "ROM01", "ITA05", "AM43",
                "FRA02", "SALD", "PIND", "POZO", "FRA01", "COUT", "RIPO") # by Latitude
PopsLabeled <- c("NOR03","NOR01", "DEN03", "SWE01", "AM10", "AM10", "LIT02", "GER07",
                 "AM46", "AM30", "NET01", "AM15", "AM304", "GER06", "CZE01", "AM32", "AM32", 
                 "AM37", "LAG", "AM42", "CAR", "AM40", "AUS01", "ROM01", "ITA05", "AM43",
                 "FRA02", "SALD", "PIND", "POZO", "FRA01", "COUT", "RIPO")
ArnicaPlotdata$Pop <- factor(ArnicaPlotdata$Pop, levels=PopsSorted, labels=PopsLabeled)

### sort the loci by source (Berlin > Belgium)
LociSorted <- c(sort(levels(ArnicaPlotdata$Locus))[6:19], sort(levels(ArnicaPlotdata$Locus))[1:5], sort(levels(ArnicaPlotdata$Locus))[20:23])
ArnicaPlotdata$Locus <- factor(ArnicaPlotdata$Locus, levels=LociSorted)


### make a barplot of allele presence/absence in individuals of each population
colorlist <- c("#057E53", "#05AE53", "#05BE53", "#05EE53", 
               "#0500EE", "#0555EE", "#05AAEE", "#05EEEE",
               "#770053", "#AA0053", "#CC0053", "#FF0053")
ArnicaPlot <- ggplot(ArnicaPlotdata, aes(fill=Study, x=Allele, y=Count)) +
  geom_bar(position = "dodge", stat = "identity", width=0.3) +
  labs(x="alleles", y="individuals with allele") + theme_light() +
  scale_fill_manual(values=colorlist) +
  facet_grid(vars(Pop), vars(Locus), scales="free_x", axes="all_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4),
        axis.text.y = element_text(size=6))

pdf("ArnicaSSR_Barplots.pdf", width = 34*1.5, height = 24*1.5)
print(ArnicaPlot) # this plot is huge - printing takes ca. 1-2 minutes
dev.off()

### alternative plot saving method if ggplot is not stored in "ArnicaPlot"
#ggsave("SSRseq_Barplots.png", dpi=600, height=14, width=18, units="cm")

### plot each facet on an individual page in a pdf - does not work for now
#install.packages("ggforce")
#library(ggforce)

#ArnicaPlot + facet_grid_paginate(Pop ~ Locus, scales="free_x", ncol=1, nrow=1, page=1)

#pdf("SSRseq_Barplots_single.pdf", width = 5.5, height = 8)
#for (i in seq(length(unique(labels(ArnicaPlotdata$Pop)))*
#              #length(labels(ArnicaPlotdata$Locus)))){
#  ArnicaPage <- ArnicaPlot + facet_grid_paginate(Pop ~ Locus, 
#                                   scales="free_x", ncol=1, nrow=1, page=2)
#  print(ArnicaPage)
#}
#dev.off()
