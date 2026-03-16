


library("adegenet") 
  library("ape")
library("devtools")
  library("dplyr")
library("ggplot2")
library("gplots")
library("hierfstat")
library("lattice")
library("magrittr")
  library("pegas")
library("pheatmap")
library("poppr")
library("reshape2")
  library("tidyverse")

package?poppr
citation("poppr")

citation("tidyverse")

setwd("C:/Users/Elisabeet/Desktop/BA/R")

Arnica_D <- read.genalex("DATEN.csv", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = ";", recode = FALSE)
Arnica_L <- read.genalex("OhneFULL.csv", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = ";", recode = FALSE)
Arnica_D <- read.genalex("DuweEtAl2017_ArnicaMontana.csv", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = ";", recode = FALSE)
Arnica_S18 <- read.genalex("Mit18.csv", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = ";", recode = FALSE)
Arnica_L18 <- read.genalex("Ohne18.csv", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = ";", recode = FALSE)

# Private Alleles ---------------------------------------------------------

private_alleles(
  Arnica_D,
  form = alleles ~ .,
  report = "table",
  level = "population",
  count.alleles = TRUE,
  drop = FALSE
)

arnica_private2 <- private_alleles(Arnica_D, report = "data.frame")

arnica_private2 <- arnica_private2 %>% arrange(desc(allele))

data.p <- data.frame(
               S = c(30, 24, 21, 27, 10, 12, 13, 9, 12, 24, 29, 4, 11, 2, 2, 4, 8, 2, 1, 1, 5, 9, 12, 8, 0, 0, 0, 1, 2, 0, 0, 0, 3, 0, 1),
               L = c(15, 15, 14, 11, 8, 8, 6, 6, 6, 6, 6, 4, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

data_plot.p <- data.p %>%
  pivot_longer(cols = c(S, L), names_to = "Polymorphismus", values_to = "Allelanzahl")

ggplot(data_plot.p, aes(x = Allelanzahl, y = Population, fill = Polymorphismus)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Vergleich der Anzahl privater Allele (S vs. L)", x = "Anzahl private Allele", y = "Population", fill = "Polymorphismus", size = 0.5) +
  scale_fill_manual(values = c("#ED5564", "#edae49"),
                    labels = c("Länge", "Sequenz")) +
  scale_y_discrete(labels = c("1" = "SB-SALD",
                              "2" = "SB-PIND",
                              "3" = "SC-NOR03",
                              "4" = "SC-ITA05",
                              "5" = "SD-AM46",
                              "6" = "SB-RIPO",
                              "7" = "SD-AM43",
                              "8" = "SD-AM42",
                              "9" = "SD-AM40",
                              "10" = "SC-ROM01",
                              "11" = "SC-AUS01",
                              "12" = "SD-AM37",
                              "13" = "SC-LIT02",
                              "14" = "SD-AM30",
                              "15" = "SD-AM10",
                              "16" = "SC-SWE01",
                              "17" = "SC-ITA08",
                              "18" = "SC-GER07", "19" = "SR-LAG", "20" = "SR-ELS", "21" = "SC-ITA01", "22" = "SB-POZO",
                              "23" = "SB-COUT", "24" = "SR-CAR", "25" = "SD-AM32", "26" = "SD-AM21", "27" = "SD-AM15", "28" = "SC-NOR01", "29" = "SC-NET01",
                              "30" = "SC-GER06", "31" = "SC-GER03", "32" = "SC-FRA02", "33" = "SC-FRA01", "34" = "SC-DEN03", "35" = "SC-CZE01"))+
  theme_bw() +
  theme(legend.position="bottom")




ggplot(arnica_private2, aes(x = reorder(population, count), y = count)) +
  geom_col(fill = "steelblue") +
  coord_flip() + 
  theme_minimal(base_size = 10) +
  labs(x = "Population", y = "Anzahl private Allele",
       title = "Private Allele pro Population (Alte Daten)") 


data <- data.frame(
        S = c(25, 17, 12, 16, 19, 26, 22, 12, 17, 9, 6, 14, 27, 16, 28, 20, 11, 17),
        L = c(21, 15, 8, 16, 7, 23, 16, 8, 14, 5, 6, 9, 16, 11, 10, 13, 8, 7)
)

data.p$Population <- rownames(data.p)

data_plot <- data %>%
  pivot_longer(cols = c(S, L), names_to = "Polymorphismus", values_to = "Allelanzahl")
legend_labels <- c(
  "1" = "Arm01",
  "2" = "Arm02",
  "3" = "Arm03",
  "4" = "Arm04",
  "5" = "Arm05",
  "6" = "Arm06",
  "7" = "Arm08",
  "8" = "Arm09",
  "9" = "Arm10",
  "10" = "Arm11",
  "11" = "Armo03",
  "12" = "Am-AG-1",
  "13" = "Am-AG-10",
  "14" = "Am-AG-4B",
  "15" = "Am-AG-11",
  "16" = "Am-ATC-3",
  "17" = "Am-CT-2",
  "18" = "Am-AG-2B",
  "S" = "Sequenzpolymorphismus",
  "L" = "Längenpolymorphismus"
  )

ggplot(data_plot, aes(x = loci, y = Allelanzahl, fill = Polymorphismus)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Vergleich der Allelanzahl (S vs. L)", x = "Locus", y = "Allelanzahl", fill = "Polymorphismus", size = 0.5) +
  scale_fill_manual(values = c("#ED5564", "#edae49"),
                    labels = c("Länge", "Sequenz")) +
  scale_x_discrete(labels = c("1" = "Arm01",
                              "2" = "Arm02",
                              "3" = "Arm03",
                              "4" = "Arm04",
                              "5" = "Arm05",
                              "6" = "Arm06",
                              "7" = "Arm08",
                              "8" = "Arm09",
                              "9" = "Arm10",
                              "10" = "Arm11",
                              "11" = "Armo03",
                              "12" = "Am-AG-1",
                              "13" = "Am-AG-10",
                              "14" = "Am-AG-4B",
                              "15" = "Am-AG-11",
                              "16" = "Am-ATC-3",
                              "17" = "Am-CT-2",
                              "18" = "Am-AG-2B"))+
  theme_bw() +
  theme(legend.position="bottom")


library(tidyr)
library(dplyr)

# Beispiel: dein Datensatz heißt df
# df <- read.csv("deine_datei.csv")

df_wide <- arnica_private %>%
  pivot_wider(names_from = allele, values_from = count, values_fill = 0)

head(df_wide)

#Ende ----------------------------------------------------------------------------------------------------------------------------------------




#Subsets für die Pop von Sequenzpolymorphismus
AUS01_subset_s <- popsub(Arnica_S, sublist = "SC-AUS01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CZE01_subset_s <- popsub(Arnica_S, sublist = "SC-CZE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
DEN03_subset_s <- popsub(Arnica_S, sublist = "SC-DEN03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA01_subset_s <- popsub(Arnica_S, sublist = "SC-FRA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA02_subset_s <- popsub(Arnica_S, sublist = "SC-FRA02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER03_subset_s <- popsub(Arnica_S, sublist = "SC-GER03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER06_subset_s <- popsub(Arnica_S, sublist = "SC-GER06", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER07_subset_s <- popsub(Arnica_S, sublist = "SC-GER07", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA01_subset_s <- popsub(Arnica_S, sublist = "SC-ITA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA05_subset_s <- popsub(Arnica_S, sublist = "SC-ITA05", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA08_subset_s <- popsub(Arnica_S, sublist = "SC-ITA08", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LIT02_subset_s <- popsub(Arnica_S, sublist = "SC-LIT02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ROM01_subset_s <- popsub(Arnica_S, sublist = "SC-ROM01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SWE01_subset_s <- popsub(Arnica_S, sublist = "SC-SWE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NET01_subset_s <- popsub(Arnica_S, sublist = "SC-NET01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR01_subset_s <- popsub(Arnica_S, sublist = "SC-NOR01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR03_subset_s <- popsub(Arnica_S, sublist = "SC-NOR03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
COUT_subset_s <- popsub(Arnica_S, sublist = "SB-COUT", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
PIND_subset_s <- popsub(Arnica_S, sublist = "SB-PIND", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
POZO_subset_s <- popsub(Arnica_S, sublist = "SB-POZO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
RIPO_subset_s <- popsub(Arnica_S, sublist = "SB-RIPO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SALD_subset_s <- popsub(Arnica_S, sublist = "SB-SALD", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM10_subset_s <- popsub(Arnica_S, sublist = "SD-AM10", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM15_subset_s <- popsub(Arnica_S, sublist = "SD-AM15", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM21_subset_s <- popsub(Arnica_S, sublist = "SD-AM21", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM30_subset_s <- popsub(Arnica_S, sublist = "SD-AM30", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM32_subset_s <- popsub(Arnica_S, sublist = "SD-AM32", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM37_subset_s <- popsub(Arnica_S, sublist = "SD-AM37", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM40_subset_s <- popsub(Arnica_S, sublist = "SD-AM40", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM42_subset_s <- popsub(Arnica_S, sublist = "SD-AM42", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM43_subset_s <- popsub(Arnica_S, sublist = "SD-AM43", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM46_subset_s <- popsub(Arnica_S, sublist = "SD-AM46", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CAR_subset_s <- popsub(Arnica_S, sublist = "SR-CAR", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ELS_subset_s <- popsub(Arnica_S, sublist = "SR-ELS", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LAG_subset_s <- popsub(Arnica_S, sublist = "SR-LAG", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)

#für den Vergleich mit <50% Missing Rate
AUS01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-AUS01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CZE01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-CZE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
DEN03_subset_s18 <- popsub(Arnica_S18, sublist = "SC-DEN03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-FRA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA02_subset_s18 <- popsub(Arnica_S18, sublist = "SC-FRA02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER03_subset_s18 <- popsub(Arnica_S18, sublist = "SC-GER03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER06_subset_s18 <- popsub(Arnica_S18, sublist = "SC-GER06", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER07_subset_s18 <- popsub(Arnica_S18, sublist = "SC-GER07", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-ITA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA05_subset_s18 <- popsub(Arnica_S18, sublist = "SC-ITA05", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA08_subset_s18 <- popsub(Arnica_S18, sublist = "SC-ITA08", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LIT02_subset_s18 <- popsub(Arnica_S18, sublist = "SC-LIT02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ROM01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-ROM01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SWE01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-SWE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NET01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-NET01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR01_subset_s18 <- popsub(Arnica_S18, sublist = "SC-NOR01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR03_subset_s18 <- popsub(Arnica_S18, sublist = "SC-NOR03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
COUT_subset_s18 <- popsub(Arnica_S18, sublist = "SB-COUT", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
PIND_subset_s18 <- popsub(Arnica_S18, sublist = "SB-PIND", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
POZO_subset_s18 <- popsub(Arnica_S18, sublist = "SB-POZO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
RIPO_subset_s18 <- popsub(Arnica_S18, sublist = "SB-RIPO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SALD_subset_s18 <- popsub(Arnica_S18, sublist = "SB-SALD", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM10_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM10", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM15_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM15", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM21_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM21", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM30_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM30", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM32_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM32", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM37_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM37", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM40_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM40", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM42_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM42", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM43_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM43", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM46_subset_s18 <- popsub(Arnica_S18, sublist = "SD-AM46", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CAR_subset_s18 <- popsub(Arnica_S18, sublist = "SR-CAR", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ELS_subset_s18 <- popsub(Arnica_S18, sublist = "SR-ELS", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LAG_subset_s18 <- popsub(Arnica_S18, sublist = "SR-LAG", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)

#Für Längenpolymorphismus
AUS01_subset_l <- popsub(Arnica_D, sublist = "SC-AUS01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CZE01_subset_l <- popsub(Arnica_D, sublist = "SC-CZE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
DEN03_subset_l <- popsub(Arnica_D, sublist = "SC-DEN03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA01_subset_l <- popsub(Arnica_D, sublist = "SC-FRA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA02_subset_l <- popsub(Arnica_D, sublist = "SC-FRA02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER03_subset_l <- popsub(Arnica_D, sublist = "SC-GER03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER06_subset_l <- popsub(Arnica_D, sublist = "SC-GER06", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER07_subset_l <- popsub(Arnica_D, sublist = "SC-GER07", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA01_subset_l <- popsub(Arnica_D, sublist = "SC-ITA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA05_subset_l <- popsub(Arnica_D, sublist = "SC-ITA05", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA08_subset_l <- popsub(Arnica_D, sublist = "SC-ITA08", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LIT02_subset_l <- popsub(Arnica_D, sublist = "SC-LIT02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ROM01_subset_l <- popsub(Arnica_D, sublist = "SC-ROM01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SWE01_subset_l <- popsub(Arnica_D, sublist = "SC-SWE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NET01_subset_l <- popsub(Arnica_D, sublist = "SC-NET01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR01_subset_l <- popsub(Arnica_D, sublist = "SC-NOR01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR03_subset_l <- popsub(Arnica_D, sublist = "SC-NOR03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
COUT_subset_l <- popsub(Arnica_D, sublist = "SB-COUT", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
PIND_subset_l <- popsub(Arnica_D, sublist = "SB-PIND", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
POZO_subset_l <- popsub(Arnica_D, sublist = "SB-POZO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
RIPO_subset_l <- popsub(Arnica_D, sublist = "SB-RIPO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SALD_subset_l <- popsub(Arnica_D, sublist = "SB-SALD", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM10_subset_l <- popsub(Arnica_D, sublist = "SD-AM10", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM15_subset_l <- popsub(Arnica_D, sublist = "SD-AM15", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM21_subset_l <- popsub(Arnica_D, sublist = "SD-AM21", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM30_subset_l <- popsub(Arnica_D, sublist = "SD-AM30", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM32_subset_l <- popsub(Arnica_D, sublist = "SD-AM32", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM37_subset_l <- popsub(Arnica_D, sublist = "SD-AM37", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM40_subset_l <- popsub(Arnica_D, sublist = "SD-AM40", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM42_subset_l <- popsub(Arnica_D, sublist = "SD-AM42", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM43_subset_l <- popsub(Arnica_D, sublist = "SD-AM43", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM46_subset_l <- popsub(Arnica_D, sublist = "SD-AM46", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CAR_subset_l <- popsub(Arnica_D, sublist = "SR-CAR", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ELS_subset_l <- popsub(Arnica_D, sublist = "SR-ELS", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LAG_subset_l <- popsub(Arnica_D, sublist = "SR-LAG", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)

# Vergleich mit <50% Missing Rate
AUS01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-AUS01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CZE01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-CZE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
DEN03_subset_l18 <- popsub(Arnica_L18, sublist = "SC-DEN03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-FRA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA02_subset_l18 <- popsub(Arnica_L18, sublist = "SC-FRA02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER03_subset_l18 <- popsub(Arnica_L18, sublist = "SC-GER03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER06_subset_l18 <- popsub(Arnica_L18, sublist = "SC-GER06", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER07_subset_l18 <- popsub(Arnica_L18, sublist = "SC-GER07", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-ITA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA05_subset_l18 <- popsub(Arnica_L18, sublist = "SC-ITA05", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA08_subset_l18 <- popsub(Arnica_L18, sublist = "SC-ITA08", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LIT02_subset_l18 <- popsub(Arnica_L18, sublist = "SC-LIT02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ROM01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-ROM01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SWE01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-SWE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NET01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-NET01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR01_subset_l18 <- popsub(Arnica_L18, sublist = "SC-NOR01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR03_subset_l18 <- popsub(Arnica_L18, sublist = "SC-NOR03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
COUT_subset_l18 <- popsub(Arnica_L18, sublist = "SB-COUT", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
PIND_subset_l18 <- popsub(Arnica_L18, sublist = "SB-PIND", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
POZO_subset_l18 <- popsub(Arnica_L18, sublist = "SB-POZO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
RIPO_subset_l18 <- popsub(Arnica_L18, sublist = "SB-RIPO", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SALD_subset_l18 <- popsub(Arnica_L18, sublist = "SB-SALD", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM10_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM10", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM15_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM15", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM21_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM21", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM30_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM30", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM32_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM32", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM37_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM37", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM40_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM40", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM42_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM42", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM43_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM43", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
AM46_subset_l18 <- popsub(Arnica_L18, sublist = "SD-AM46", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CAR_subset_l18 <- popsub(Arnica_L18, sublist = "SR-CAR", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ELS_subset_l18 <- popsub(Arnica_L18, sublist = "SR-ELS", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LAG_subset_l18 <- popsub(Arnica_L18, sublist = "SR-LAG", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)



nLoc(Arnica_D)
nPop(Arnica_D)
nLoc(Arnica_L)
nPop(Arnica_L)

poppr(Arnica_D)
poppr(Arnica_L18)

basic.stats(Arnica_D)
basic.stats(Arnica_L18)
sd(basic.stats(Arnica_S$Fis), na.rm = TRUE)



basic.stats(AUS01_subset_s18)
basic.stats(CZE01_subset_s18)
basic.stats(AM10_subset_s18)
basic.stats(AM15_subset_s18)
basic.stats(AM21_subset_s18)
basic.stats(AM30_subset_s18)
basic.stats(AM32_subset_s18)
basic.stats(AM37_subset_s18)
basic.stats(AM40_subset_s18)
basic.stats(AM42_subset_s18)
basic.stats(AM43_subset_s18)
basic.stats(AM46_subset_s18)
basic.stats(CAR_subset_s18)
basic.stats(COUT_subset_s18)
basic.stats(DEN03_subset_s18)
basic.stats(ELS_subset_s18)
basic.stats(FRA01_subset_s18)
basic.stats(FRA02_subset_s18)
basic.stats(GER03_subset_s18)
basic.stats(GER06_subset_s18)
basic.stats(GER07_subset_s18)
basic.stats(ITA01_subset_s18)
basic.stats(ITA05_subset_s18)
basic.stats(ITA08_subset_s18)
basic.stats(LAG_subset_s18)
basic.stats(LIT02_subset_s18)
basic.stats(NET01_subset_s18)
basic.stats(NOR01_subset_s18)
basic.stats(NOR03_subset_s18)
basic.stats(PIND_subset_s18)
basic.stats(POZO_subset_s18)
basic.stats(RIPO_subset_s18)
basic.stats(ROM01_subset_s18)
basic.stats(SALD_subset_s18)
basic.stats(SWE01_subset_s18)

basic.stats(AUS01_subset_l)
summary(AUS01_subset_l)
basic.stats(CZE01_subset_l)
summary(CZE01_subset_l)
basic.stats(AM10_subset_l)
summary(AM10_subset_l)
basic.stats(AM15_subset_l)
summary(AM15_subset_l)
basic.stats(AM21_subset_l)
summary(AM21_subset_l)
basic.stats(AM30_subset_l)
summary(AM30_subset_l)
basic.stats(AM32_subset_l)
summary(AM32_subset_l)
basic.stats(AM37_subset_l)
summary(AM37_subset_l)
basic.stats(AM40_subset_l)
summary(AM40_subset_l)
basic.stats(AM42_subset_l)
summary(AM42_subset_l)
basic.stats(AM43_subset_l)
summary(AM43_subset_l)
basic.stats(AM46_subset_l)
summary(AM46_subset_l)
basic.stats(CAR_subset_l)
summary(CAR_subset_l)
basic.stats(COUT_subset_l)
summary(COUT_subset_l)
basic.stats(DEN03_subset_l)
summary(DEN03_subset_l)
basic.stats(ELS_subset_l)
summary(ELS_subset_l)
basic.stats(FRA01_subset_l)
summary(FRA01_subset_l)
basic.stats(FRA02_subset_l)
summary(FRA02_subset_l)
basic.stats(GER03_subset_l)
summary(GER03_subset_l)
basic.stats(GER06_subset_l)
summary(GER06_subset_l)
basic.stats(GER07_subset_l)
summary(GER07_subset_l)
basic.stats(ITA01_subset_l)
summary(ITA01_subset_l)
basic.stats(ITA05_subset_l)
summary(ITA05_subset_l)
basic.stats(ITA08_subset_l)
summary(ITA08_subset_l)
basic.stats(LAG_subset_l)
summary(LAG_subset_l)
basic.stats(LIT02_subset_l)
summary(LIT02_subset_l)
basic.stats(NET01_subset_l)
summary(NET01_subset_l)
basic.stats(NOR01_subset_l)
summary(NOR01_subset_l)
basic.stats(NOR03_subset_l)
summary(NOR03_subset_l)
basic.stats(PIND_subset_l)
summary(PIND_subset_l)
basic.stats(POZO_subset_l)
summary(POZO_subset_l)
basic.stats(RIPO_subset_l)
summary(RIPO_subset_l)
basic.stats(ROM01_subset_l)
summary(ROM01_subset_l)
basic.stats(SALD_subset_l)
summary(SALD_subset_l)
basic.stats(SWE01_subset_l)
summary(SWE01_subset_l)



fis1 <- basic.stats(ELS_subset_l)$Fis
fis1 <- fis1[,-c(2)] # don't do this if working with all pops, only within a subset
sd(fis1, na.rm = T)
view(fis1)
summary(Arnica_L18)

AM01_subset <- popsub(Arnica_D, sublist = "AM01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
basic.stats(AM01_subset)

Arnica_D@loc.n.all
Arnica_D@all.names
Arnica_S18@loc.n.all

AUS01_subset_l18@loc.n.all
AUS01_subset_l18@loc.n.all 
CZE01_subset_l18@loc.n.all 
DEN03_subset_l18@loc.n.all 
FRA01_subset_l18@loc.n.all 
FRA02_subset_l18@loc.n.all 
GER03_subset_l18@loc.n.all 
GER06_subset_l18@loc.n.all 
GER07_subset_l18@loc.n.all 
ITA01_subset_l18@loc.n.all 
ITA05_subset_l18@loc.n.all 
ITA08_subset_l18@loc.n.all 
LIT02_subset_l18@loc.n.all 
ROM01_subset_l18@loc.n.all 
SWE01_subset_l18@loc.n.all 
NET01_subset_l18@loc.n.all 
NOR01_subset_l18@loc.n.all 
NOR03_subset_l18@loc.n.all 
COUT_subset_l18@loc.n.all 
PIND_subset_l18@loc.n.all 
POZO_subset_l18@loc.n.all 
RIPO_subset_s18@loc.n.all 
SALD_subset_l18@loc.n.all 
AM10_subset_l18@loc.n.all 
AM15_subset_l18@loc.n.all 
AM21_subset_l18@loc.n.all 
AM30_subset_l18@loc.n.all 
AM32_subset_l18@loc.n.all 
AM37_subset_l18@loc.n.all 
AM40_subset_l18@loc.n.all 
AM42_subset_l18@loc.n.all 
AM43_subset_l18@loc.n.all 
AM46_subset_l18@loc.n.all 
CAR_subset_l18@loc.n.all 
ELS_subset_l18@loc.n.all 
LAG_subset_l18@loc.n.all 











HWE.test(Arnica_S,pop=NULL,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix"))
HWE.test.genind(Arnica_S,pop=NULL,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix"))
hw.test(Arnica_S)
Arnica_S18@all.names
# tells you all allele names per locus
# useful to see bp size range

summary(Arnica_L18)
summary(Arnica_D)
summary(Arnica_L)
summary(CZE01_subset_s18)
summary(AM10_subset_s)
summary(COUT_subset_s)
Na <- summary(CZE01_subset_l)$Numberofallelespergroup
mean(Na)
Arnica_L18@loc.n.all
A18 <- read.genalex("AGAIN.csv", sep = ";")
AR <- as.matrix(allelic.richness(Arnica_S18))
summary(AM32_subset_s18)
Arnica_L18$loci
allelic.richness(Arnica_D)
boxplot(t(allelic.richness(Arnica_S18)$Ar[1:18,]),
         main = "Allelreichtum", xlab = "Loci", ylab = "AR-Werte")
legend(2, 9, c("Arm01"))

#Heat-map und genetische Distanz---------------------------------------------------


fst <- genet.dist(Arnica_D, diploid=TRUE, method="WC84")
fst

fst_boot <- boot.ppfst(Arnica_S18, nboot=1000, quant=c(0.025,0.975), diploid = TRUE)
fst_boot$ul
fst_boot$ll

fst_matrix <- as.matrix(fst,)
view(fst_matrix)
fst_matrix1 <- matrix(c(0.42,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																															
                        0.54,	0.71,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																															
                        0.32,	0.49,	0.56,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																														
                        0.50,	0.63,	0.47,	0.53,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																													
                        0.31,	0.43,	0.49,	0.29,	0.45,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																												
                        0.31,	0.48,	0.56,	0.32,	0.50,	0.27,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																											
                        0.44,	0.60,	0.70,	0.32,	0.64,	0.42,	0.41,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																										
                        0.30,	0.45,	0.64,	0.25,	0.51,	0.19,	0.33,	0.44,	NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA,																									
                        0.18,	0.49,	0.64,	0.01,	0.56,	0.23,	0.25,	0.41,	0.15,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																								
                        0.29,	0.57,	0.71,	0.25,	0.54,	0.29,	0.34,	0.35,	0.15,	0.14,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,	NA,																					
                        0.35,	0.56,	0.72,	0.23,	0.61,	0.38,	0.36,	0.46,	0.35,	0.19,	0.25,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																						
                        0.32,	0.52,	0.62,	0.23,	0.53,	0.32,	0.41,	0.39,	0.34,	0.22,	0.13,	0.18,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																					
                        0.17,	0.39,	0.60,	0.23,	0.49,	0.17,	0.24,	0.41,	0.12,	0.10,	0.18,	0.26,	0.28,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,																				
                        0.35, 0.46,	0.51,	0.30,	0.47,	0.24,	0.30,	0.40,	0.21,	0.25,	0.23,	0.40,	0.36,	0.16,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA,																			
                        0.23,	0.40,	0.49,	0.18,	0.38,	0.18,	0.25,	0.33,	0.08,	0.10,	0.04,	0.23,	0.20,	0.05,	0.06,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA,
                        0.39,	0.54,	0.62,	0.34,	0.59,	0.35,	0.37,	0.46,	0.34,	0.36,	0.41,	0.46,	0.45,	0.32,	0.37, 0.28,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,	NA,  																	
                        0.35,	0.56,	0.70,	0.26,	0.57,	0.35,	0.37,	0.40,	0.37,	0.27,	0.24,	0.22,	0.23,	0.32,	0.33,	0.18,	0.41,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA,																
                        0.21,	0.53,	0.66,	0.30,	0.56,	0.34,	0.26,	0.41,	0.29,	0.20,	0.24,	0.22,	0.30,	0.21,	0.32,	0.20,	0.40,	0.17,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  	NA,														
                        0.39,	0.61,	0.71,	0.42,	0.62,	0.38,	0.41,	0.37,	0.46,	0.42,	0.45,	0.42,	0.42,	0.38,	0.42,	0.33,	0.47,	0.42,	0.41, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA,
                        0.45,	0.54,	0.65,	0.38,	0.59,	0.36,	0.38,	0.52,	0.37,	0.40,	0.43,	0.45,	0.47,	0.35,	0.37,	0.30,	0.37,	0.44,	0.43, 0.49,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,	NA, 												
                        0.40,	0.60,	0.68,	0.36,	0.60,	0.37,	0.41,	0.35,	0.44,	0.38,	0.40,	0.33,	0.33,	0.35,	0.38,	0.28,	0.45,	0.35,	0.41, 0.37,	0.49,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,	 NA,												
                        0.28,	0.49,	0.59,	0.28,	0.51,	0.33,	0.34,	0.33,	0.31,	0.18,	0.21,	0.17,	0.12,	0.21,	0.34,	0.18,	0.42,	0.12,	0.09, 0.35,	0.45,	0.30,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,		NA,											
                        0.41,	0.57,	0.65,	0.43,	0.57,	0.38,	0.40,	0.49,	0.44,	0.39,	0.43,	0.46,	0.43,	0.36,	0.42,	0.28,	0.43,	0.41,	0.42, 0.46,	0.52,	0.44,	0.35,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,	 	NA,									
                        0.17,	0.51,	0.70,	0.13,	0.57,	0.27,	0.24,	0.30,	0.22,	0.00,	0.10,	0.26,	0.27,	0.11,	0.30,	0.09,	0.36,	0.22,	0.15, 0.45,	0.41,	0.37,	0.16,	0.33,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,		NA,								
                        0.37,	0.55,	0.61,	0.34,	0.53,	0.33,	0.36,	0.33,	0.37,	0.26,	0.25,	0.23,	0.25,	0.30,	0.36,	0.18,	0.40,	0.15,	0.31, 0.35,	0.45,	0.21,	0.14,	0.32,	0.26,	NA, NA, NA, NA, NA, NA, NA, NA, NA,	 	NA,							
                        0.22,	0.52,	0.68,	0.29,	0.58,	0.29,	0.34,	0.35,	0.29,	0.03,	0.28,	0.18,	0.33,	0.12,	0.34,	0.16,	0.41,	0.18,	0.18, 0.36,	0.44,	0.29,	0.23,	0.40,	0.15,	0.27,	NA, NA, NA, NA, NA, NA, NA, NA,		NA,						
                        0.27,	0.43,	0.51,	0.17,	0.46,	0.24,	0.21,	0.27,	0.22,	0.13,	0.19,	0.18,	0.24,	0.19,	0.28,	0.16,	0.27,	0.28,	0.25, 0.31,	0.33,	0.30,	0.21,	0.26,	0.08,	0.24,	0.13,	NA, NA, NA, NA, NA, NA, NA,	 	NA,					
                        0.25,	0.37,	0.45,	0.19,	0.38,	0.12,	0.18,	0.30,	0.14,	0.01,	0.12,	0.22,	0.23,	0.04,	0.19,	0.06,	0.28,	0.23,	0.21, 0.32,	0.29,	0.30,	0.23,	0.32,	0.13,	0.28,	0.20,	0.11,	NA, NA, NA, NA, NA, NA,		NA,				
                        0.22,	0.44,	0.52,	0.17,	0.47,	0.28,	0.30,	0.33,	0.23,	0.07,	0.21,	0.20,	0.19,	0.17,	0.30,	0.15,	0.34,	0.18,	0.26, 0.36,	0.39,	0.30,	0.20,	0.32,	0.06,	0.25,	0.16,	0.18,	0.17,	NA, NA, NA, NA, NA,	 	NA,			
                        0.35,	0.38,	0.54,	0.34,	0.49,	0.27,	0.41,	0.47,	0.16,	0.29,	0.37,	0.42,	0.39,	0.17,	0.34,	0.21,	0.42,	0.43,	0.41, 0.50,	0.36,	0.49,	0.40,	0.46,	0.34,	0.43,	0.39,	0.28,	0.22,	0.34,	NA, NA, NA, NA,	 	NA,		
                        0.33,	0.50,	0.58,	0.20,	0.52,	0.33,	0.40,	0.34,	0.31,	0.17,	0.27,	0.24,	0.21,	0.29,	0.34,	0.16,	0.39,	0.25,	0.30, 0.40,	0.40,	0.36,	0.16,	0.42,	0.29,	0.24,	0.32,	0.23,	0.23,	0.23,	0.38,	NA, NA, NA,		NA,
                        0.27,	0.47,	0.56,	0.22,	0.51,	0.31,	0.33,	0.38,	0.30,	0.06,	0.22,	0.19,	0.22,	0.20,	0.33, 0.15,	0.35,	0.26,	0.30, 0.37,	0.38,	0.29,	0.25,	0.29,	0.12,	0.24,	0.20,	0.16,	0.20,	0.12,	0.37,	0.23,	NA, NA,		 NA,
                        0.21,	0.44,	0.62,	0.11,	0.54,	0.32,	0.29,	0.36,	0.31,	0.00, 0.15,	0.17,	0.22,	0.10,	0.28,	0.09,	0.35,	0.20,	0.10, 0.38,	0.38,	0.28,	0.11,	0.38,	0.00,	0.22,	0.12,	0.10,	0.13,	0.09,	0.32,	0.22,	0.08,	NA,	NA,
                        0.31,	0.59,	0.78,	0.19,	0.61,	0.28,	0.34,	0.31,	0.15,	0.06,	0.08,	0.30,	0.26,	0.07,	0.26, 0.06, 0.37, 0.21, 0.23, 0.43,	0.42,	0.22,	0.15,	0.47,	0.17,	0.16,	0.05,	0.15,	0.17,	0.21,	0.37,	0.26,	0.22,	0.03, NA,), nrow = 35, byrow = TRUE, ncol = 35)

fst_matrix2 <- matrix(c(0.28,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                        0.18,	0.41,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                        0.20,	0.40,	0.17,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                        0.33,	0.47,	0.42,	0.41, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                        0.30,	0.37,	0.44,	0.43, 0.49,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 												
                        0.28,	0.45,	0.35,	0.41, 0.37,	0.49,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 												
                        0.18,	0.42,	0.12,	0.09, 0.35,	0.45,	0.30,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 											
                        0.28,	0.43,	0.41,	0.42, 0.46,	0.52,	0.44,	0.35,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 										
                        0.09,	0.36,	0.22,	0.15, 0.45,	0.41,	0.37,	0.16,	0.33,	NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 									
                        0.18,	0.40,	0.15,	0.31, 0.35,	0.45,	0.21,	0.14,	0.32,	0.26,	NA, NA, NA, NA, NA, NA, NA, NA, NA, 								
                        0.16,	0.41,	0.18,	0.18, 0.36,	0.44,	0.29,	0.23,	0.40,	0.15,	0.27,	NA, NA, NA, NA, NA, NA, NA, NA, 							
                        0.16,	0.27,	0.28,	0.25, 0.31,	0.33,	0.30,	0.21,	0.26,	0.08,	0.24,	0.13,	NA, NA, NA, NA, NA, NA, NA, 						
                        0.06,	0.28,	0.23,	0.21, 0.32,	0.29,	0.30,	0.23,	0.32,	0.13,	0.28,	0.20,	0.11,	NA, NA, NA, NA, NA, NA, 					
                        0.15,	0.34,	0.18,	0.26, 0.36,	0.39,	0.30,	0.20,	0.32,	0.06,	0.25,	0.16,	0.18,	0.17,	NA, NA, NA, NA, NA, 				
                        0.21,	0.42,	0.43,	0.41, 0.50,	0.36,	0.49,	0.40,	0.46,	0.34,	0.43,	0.39,	0.28,	0.22,	0.34,	NA, NA, NA, NA, 			
                        0.16,	0.39,	0.25,	0.30, 0.40,	0.40,	0.36,	0.16,	0.42,	0.29,	0.24,	0.32,	0.23,	0.23,	0.23,	0.38,	NA, NA, NA, 	
                        0.15,	0.35,	0.26,	0.30, 0.37,	0.38,	0.29,	0.25,	0.29,	0.12,	0.24,	0.20,	0.16,	0.20,	0.12,	0.37,	0.23,	NA, NA, 
                        0.09,	0.35,	0.20,	0.10, 0.38,	0.38,	0.28,	0.11,	0.38,	0.00,	0.22,	0.12,	0.10,	0.13,	0.09,	0.32,	0.22,	0.08,	NA,
                        0.06, 0.37, 0.21, 0.23, 0.43,	0.42,	0.22,	0.15,	0.47,	0.17,	0.16,	0.05,	0.15,	0.17,	0.21,	0.37,	0.26,	0.22,	0.03), nrow = 19, byrow = TRUE, ncol = 19)
colors(distinct = FALSE)
heatmap.2(fst_matrix1, 
          Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          main = "FST-Werte",
          margins = c(5, 5),
          cellnote = fst_matrix1,
          notecol = "black", 
          notecex = 0.9,
          labRow = c("PIND",	"POZO",	"RIPO",	"SALD",	"AUS01",	"CZE01",	"DEN03",	"FRA01",	"FRA02",	"GER03",	"GER06",	"GER07",	"ITA01",	"ITA05",	"ITA08",	"LIT02",	"NET01",	"NOR01",	"NOR03",	"ROM01",	"SWE01",	"AM10",	"AM15",	"AM21",	"AM30",	"AM32",	"AM37",	"AM40",	"AM42",	"AM43",	"AM46",	"CAR",	"ELS",	"LAG"),
          labCol = c("COUT",	"PIND",	"POZO",	"RIPO",	"SALD",	"AUS01",	"CZE01",	"DEN03",	"FRA01",	"FRA02",	"GER03",	"GER06",	"GER07",	"ITA01",	"ITA05",	"ITA08",	"LIT02",	"NET01",	"NOR01",	"NOR03",	"ROM01",	"SWE01",	"AM10",	"AM15",	"AM21",	"AM30",	"AM32",	"AM37",	"AM40",	"AM42",	"AM43",	"AM46",	"CAR",	"ELS",	"LAG"),
          cexRow = 0.8,
          cexCol = 0.9,
          col = colorRampPalette(c("#c8e4fa", "#9ac3e6", "#4b6296"))(100),
          trace = "none",
          key = TRUE,
          keysize = 1,
          density.info = "none",
          symkey = FALSE,
          key.title = "FST")

heatmap.2(fst_matrix2, 
          Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          margins = c(5, 5),
          cellnote = fst_matrix2,
          notecol = "black", 
          notecex = 0.9,
          labRow = c("LIT02",	"NET01",	"NOR01",	"NOR03",	"ROM01",	"SWE01",	"AM10",	"AM15",	"AM21",	"AM30",	"AM32",	"AM37",	"AM40",	"AM42",	"AM43",	"AM46",	"CAR",	"ELS",	"LAG"),
          labCol = c("ITA08",	"LIT02",	"NET01",	"NOR01",	"NOR03",	"ROM01",	"SWE01",	"AM10",	"AM15",	"AM21",	"AM30",	"AM32",	"AM37",	"AM40",	"AM42",	"AM43",	"AM46",	"CAR",	"ELS",	"LAG"),
          cexRow = 0.9,
          cexCol = 0.9,
          col = colorRampPalette(c("#c8e4fa", "#9ac3e6", "#4b6296"))(100),
          trace = "none",
          key = TRUE,
          keysize = 1.0,
          density.info = "none",
          symkey = FALSE,
          key.title = "FST")

heatmap.2(fst_matrix, 
          Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          main = "FST Values",
          margins = c(5, 5),
          cellnote = fst_matrix,
          notecol = "black",
          labRow = c("SB-COUT", "SB-PIND", "SB-POZO", "SB-RIPO", "SB-SALD", ""),
          labCol = c("AM102", "AM202", "AM205", "AM30"),
          cexRow = 1,
          cexCol = 1,
          col = colorRampPalette(c("#c8e4fa", "#9ac3e6", "#4b6296"))(100),
          trace = "none",
          key = TRUE,
          keysize = 0.5,
          density.info = "none",
          symkey = FALSE,
          key.title = "FST")


arnica_genind <- genclone2genind(Arnica_L18)

arnica_gendist <- prevosti.dist(as.matrix(arnica_genind$tab))
arnica_gendist <- prevosti.dist(arnica_genind)

repeat_lengths <- c(4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 6, 2, 2, 2, 2, 3, 2, 2)
arnica_gendist <- bruvo.dist(arnica_genind, replen = repeat_lengths)

# PCoA --------------------------------------------------------------------

# Prevosti's distance can be used for non size-based markers
# (e.g. sequence-based microsatellites)
# Bruvo's can be used for size based-analysis (stepwise mutation model)

# convert genclone to genind object
arnica_genind <- genclone2genind(Arnica_D)
#arnica_genind_clean <- na.omit(arnica_genind), this doesn't seem to be necessary

# create distance matrix (PREVOSTI)

#arnica_gendist <- dist(as.matrix(arnica_genind$tab))
arnica_gendist <- prevosti.dist(as.matrix(arnica_genind$tab))
arnica_gendist <- prevosti.dist(arnica_genind)

# create distance matrix (BRUVO)
# needs definition of repeat length

repeat_lengths <- c(4, 4, 4, 4, 4, 4, 2)
arnica_gendist <- bruvo.dist(arnica_genind, replen = repeat_lengths)

####
# can also use other distances
# https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html
####

arnica_gendist

pcoa_result <- cmdscale(arnica_gendist, k = 2)
plot(pcoa_result, main = "PCoA Plot")

pcoa_result2 <- cmdscale(arnica_gendist, k = 68, eig = T)

pcoa_result2$eig[1]/sum(pcoa_result2$eig) * 100

pcoa_result2$eig[2]/sum(pcoa_result2$eig) * 100


pcoa_df <- data.frame(pc_1 = pcoa_result[, 1], pc_2 = pcoa_result[, 2])
pcoa_df$Population <- as.factor(pop(arnica_genind))

ggplot(pcoa_df, position = jitter, aes(x = pc_1, y = pc_2, color = Population, shape = Population)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "red1", "red2", "red3", "red4", "navy", "blue", "skyblue", "cyan", "cyan2", "steelblue", "steelblue1", "steelblue2", "turquoise", "turquoise2", "turquoise3", "steelblue4", "cyan4", "royalblue", "royalblue1", "royalblue4", "skyblue3", "magenta1", "maroon1", "maroon2", "maroon3", "maroon4", "magenta2", "magenta3", "magenta4", "orchid4", "orchid3", "orange1", "orange2", "orange3")) +
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 18, 12, 13, 14, 15, 16, 17, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3)) +
  ggtitle("PCoA - Alter Datensatz") + 
  xlab("Axis 1 (20.37%)") +
  ylab("Axis 2 (11.80%)") +
  theme_bw()

colours()
shape
# type in the values you get above for the xlab/ylab

# screeplot eigenvalues
"#18b6b2","#82e748","#edae49","#ED5564","#6958b5", "cyan", "coral", "olivedrab", "brown", "red", "blue", "bisque4", "gray", "gold", "deeppink", "plum4", "lightpink", "khaki", "turquoise", "maroon", "orange", "green", "orchid", "sienna", "peachpuff2", "powderblue", "tan4", "purple3", "tomato", "tan", "salmon", "black", "seagreen2", "thistle", "seagreen"
pcoa_result2$eig
explained <- pcoa_result2$eig / sum(pcoa_result2$eig)

tibble(pe = explained, 
       axis = 1:length(explained)) %>%
  ggplot(aes(x = axis, y=pe)) +
  geom_line() +
  ggtitle("Screeplot Sequence-based data") +
  ylab("Percentage explained") +
  xlab("Axis") +
  theme_bw()

view(allelic.richness(Arnica_S18))


total_hw <- hw.test(Arnica_S, B = 0)
pop_hw <- seppop(Arnica_S) %>% lapply(hw.test, B = 0)
pop_hw_mat <- sapply(pop_hw, "[", i = TRUE, j = 3)

x <- as.loci(Arnica_S)
hw.test(x)
alpha <- 0.05
pop_hw_new <- pop_hw_mat
pop_hw_new[pop_hw_new > alpha] <- 1

library("lattice")
levelplot(t(pop_hw_new), xlab = "Populations", ylab = "Loci",
          main = "Deviations from Hardy- Weinberg")

view(allelic.richness(SWE01_subset_s18)$Ar)
mean(AR1)




Arnica_C <- read.genalex("BELGIEN.csv", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = ";", recode = FALSE)

AUS01 <- popsub(Arnica_C, sublist = "SC-AUS01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
CZE01 <- popsub(Arnica_C, sublist = "SC-CZE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
DEN03 <- popsub(Arnica_C, sublist = "SC-DEN03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA01 <- popsub(Arnica_C, sublist = "SC-FRA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
FRA02 <- popsub(Arnica_C, sublist = "SC-FRA02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER03 <- popsub(Arnica_C, sublist = "SC-GER03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER06 <- popsub(Arnica_C, sublist = "SC-GER06", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
GER07 <- popsub(Arnica_C, sublist = "SC-GER07", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA01 <- popsub(Arnica_C, sublist = "SC-ITA01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA05 <- popsub(Arnica_C, sublist = "SC-ITA05", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ITA08 <- popsub(Arnica_C, sublist = "SC-ITA08", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
LIT02 <- popsub(Arnica_C, sublist = "SC-LIT02", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NET01 <- popsub(Arnica_C, sublist = "SC-NET01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR01 <- popsub(Arnica_C, sublist = "SC-NOR01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
NOR03 <- popsub(Arnica_C, sublist = "SC-NOR03", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
ROM01 <- popsub(Arnica_C, sublist = "SC-ROM01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)
SWE01 <- popsub(Arnica_C, sublist = "SC-SWE01", exclude = NULL, blacklist = NULL, mat = NULL, drop = TRUE)


poppr(Arnica_C)

basic.stats(SWE01)
summary(AUS01)

#Isolation by Distanz---------------------------------------------------------------------------------------------------
#bzw Korrelationen herausfinden----------------------------------------------------------------------------------------

fst.tran <- fst_matrix / (1 - fst_matrix)

library(hierfstat)   # für Fst
library(geosphere)   # für geographische Distanz
library(ggplot2)
geopop <- read.csv("LocationPop11.csv", header = TRUE, sep = ";")
dist(geopop, method = "maximum")
view(pop.distgeo)
install.packages("geosphere")


coords <- data.frame(
          Pop = 1:35,
          Lon = c(-7.098171422, -9.110212004, -6.832622226, 2.190319275, -3.386262099, 12.40865,
                  16.35439, 8.361139, 2.20942, 7.2174, 9.471056, 8.807, 9.8161, 9.172926,
                  10.92444, 8.311028, 24.35278, 6.22, 6.742692, 9.565, 22.835, 13.903, 12.703056,
                  13.772317, 8.020739, 11.643639, 8.429342, 12.2915, 12.865022, 11.223953, 6.419611,
                  12.13525, 1.973496, 6.221667, 5.722778),
          Lat = c(42.60765393, 42.88454738, 42.8820907, 42.20063349, 43.14786888, 46.84454,
                  50.37576, 56.482, 42.63571, 44.05773, 50.62905, 50.945, 54.0458, 44.58704,
                  46.09739, 46.35628, 54.14945, 52.80181, 58.0738, 59.531, 46.526, 55.808,
                  54.392778, 51.51006, 47.868375, 53.305294, 50.249122, 50.238436, 47.642931,
                  48.110811, 44.486361, 53.3592, 47.98458, 50.47528, 49.64639)
)
geo.mat <- distm(coords[,c("Lon","Lat")], fun = distHaversine)
geo.mat.km <- geo.mat / 1000

getPairs <- function(mat, pops) {
  combs <- combn(pops, 2)
  data.frame(
    Pop1 = combs[1,],
    Pop2 = combs[2,],
    Value = mat[lower.tri(mat)]
  )
}
fst.df <- getPairs(as.matrix(fst.tran), coords$Pop)
geo.df <- getPairs(as.matrix(geo.mat.km), coords$Pop)

ibd.df <- merge(fst.df, geo.df, by = c("Pop1","Pop2"))
names(ibd.df) <- c("Pop1","Pop2","GeneticDist","GeoDist")

ggplot(ibd.df, aes(x = GeoDist, y = GeneticDist)) +  # km statt Meter
  geom_point(shape = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Distanz (km)",
       y = expression(F[ST]/(1-F[ST])),
       title = "Isolation durch Distanz (Alte Daten)") +
  theme_bw(base_size = 14)+
  annotate("text", x = max(ibd.df$GeoDist)*0.7, y = max(ibd.df$GeneticDist)*0.9,
           label = label_text, size = 5, hjust = 0)

cor.test(ibd.df$GeoDist, ibd.df$GeneticDist)

# 1. Korrelationsanalyse
test <- cor.test(ibd.df$GeoDist, ibd.df$GeneticDist)

r_value <- round(test$estimate, 2)  # gerundeter r-Wert
p_value <- test$p.value

# 2. Signifikanz-Sternchen
stars <- ifelse(p_value < 0.001, "***",
                ifelse(p_value < 0.01, "**",
                       ifelse(p_value < 0.05, "*", "ns")))

label_text <- paste0("r = ", r_value, stars)


# Beispiel: zwei Distanzmatrizen (müssen "dist"-Objekte sein)
gen.dist <- as.dist(fst.tran)   # genetische Distanzen
geo.dist <- as.dist(geo.mat.km)     # geographische Distanzen

mantel.result <- mantel(gen.dist, geo.dist, method = "pearson", permutations = 9999)
print(mantel.result)

poppr.amova(Arnica_D)
pop(Arnica_S18) <- c(rep("Pop1", 5),
                     rep("Pop2", 5),
                     rep("Pop3", 5),
                     rep("Pop4", 5),
                     rep("Pop5", 6),
                     rep("Pop6", 6),
                     rep("Pop7", 5),
                     rep("Pop8", 5),
                     rep("Pop9", 2),
                     rep("Pop10", 2),
                     rep("Pop11", 2),
                     rep("Pop12", 5),
                     rep("Pop13", 4),
                     rep("Pop14", 2),
                     rep("Pop15", 4),
                     rep("Pop16", 4),
                     rep("Pop17", 5),
                     rep("Pop18", 5),
                     rep("Pop19", 3),
                     rep("Pop20", 5),
                     rep("Pop21", 5),
                     rep("Pop22", 5),
                     rep("Pop23", 5),
                     rep("Pop24", 5),
                     rep("Pop25", 2),
                     rep("Pop26", 6),
                     rep("Pop27", 4),
                     rep("Pop28", 4),
                     rep("Pop29", 4),
                     rep("Pop30", 5),
                     rep("Pop31", 5),
                     rep("Pop32", 6),
                     rep("Pop33", 6),
                     rep("Pop34", 2),
                     rep("Pop35", 1))


# Welche Populationen gibt es?
levels(pop(Arnica_S18))

# Wie viele Individuen pro Population?
table(pop(Arnica_S18))

library(poppr)

# AMOVA mit vorhandener Populationshierarchie
amova.result <- poppr.amova(Arnica_S18, ~Pop/Subpop)

# Ergebnisse ansehen
amova.result

# Signifikanz mit Permutationstest
amova.test <- randtest(amova.result, nrepet = 9999)
amova.test


data(Arnica_S18)


strata(Arnica_S18) <- data.frame(
  Population = pop(Arnica_S18)
)
gc_pruned <- clonecorrect(Arnica_S18, strata = ~Population)
amova.result <- 
poppr.amova(gc_pruned, ~Population, missing = "mean")


# Jetzt kannst du explizit mit ~Population rechnen
amova.result <- poppr.amova(Arnica_S18, ~Population, missing = "ignore")
amova.result



strata(Arnica_S18) <- data.frame(Population = pop(Arnica_S18))
amova.result <- poppr.amova(Arnica_S18, ~Population, missing = "mean")



class(Arnica_S18)
popNames(Arnica_S18)      # zeigt dir die Populationsnamen
nPop(Arnica_S18) 
amova.result <- poppr.amova(Arnica_S18, ~Population, missing = "mean")
amova.result

poppr.amova(
  Arnica_S18,
  ~Population,
  clonecorrect = TRUE,
  within = TRUE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  sep = "_",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 1L,
  missing = "loci",
  cutoff = 0.5,
  quiet = FALSE,
  method = c("pegas"),
  nperm = 0
)
#Korrelationen ----------------------------------------------------------------------------------
# Pakete laden
library(ggplot2)
library(dplyr)

# Beispiel-Daten (bitte durch deine ersetzen)
# ------------------------------------------------
# Angenommen, du hast ein Dataframe namens df:
# Spalten: Latitude, uHe, Ar, Fis
set.seed(123)
df <- data.frame(
  Latitude = c(42.60765393, 42.88454738, 42.8820907, 42.20063349, 43.14786888, 46.84454,
                   50.37576, 56.482, 42.63571, 44.05773, 50.62905, 50.945, 54.0458, 44.58704,
                   46.09739, 46.35628, 54.14945, 52.80181, 58.0738, 59.531, 46.526, 55.808,
                   54.392778, 51.51006, 47.868375, 53.305294, 50.249122, 50.238436, 47.642931,
                   48.110811, 44.486361, 53.3592, 47.98458, 50.47528),
  Ho = c(0.460, 0.321, 0.068, 0.505, 0.258, 0.602, 0.511, 0.341, 0.500, 0.562, 0.441, 0.387, 0.555, 0.785, 0.525, 0.684, 0.338, 0.491, 0.324, 0.385, 0.403, 0.363, 0.434,
        0.449, 0.558, 0.456, 0.324, 0.597, 0.632, 0.611, 0.580, 0.611, 0.511, 0.392 ),
  Ar = c(1.4, 1.1, 1.0, 3.1, 1.1, 2.2,
         1.5, 1.8, 1.0, 1.4, 1.3, 1.4, 2.1, 1.4, 1.6, 1.3, 1.9, 1.4, 1.4, 1.4, 1.5, 1.4, 2.6,
         2.0, 1.5, 2.5, 1.4, 2.5, 2.6, 2.4, 1.6, 2.4, 2.8, 1.2 ),
  Fis = c(0.166, -0.045, 0.225, 0.087, 0.019, 0.062,
              0.013, 0.046, 0.263, 0.250, 0.227, 0.106, -0.095, -0.047, 0.114,
              0.108, 0.262, -0.096,  0.395, -0.103, 0.090, 0.062, 0.187,
              -0.114, 0.022,  0.038,  0.405,  0.042, 0.165, 0.018, -0.051, -0.174,
              0.077, 0.292),
  Elevation = c(1330, 535, 1495, 700, 1002, 2478, 
                    766, 24, 2094, 1468, 352, 342, 10, 1222, 1122, 2334, 113, 8, 11, 543, 1343, 181, 0,
                    97, 1424, 35, 538, 585, 1425, 569, 1900, 62, 131, 570),
  Allelanzahl = c(45, 29, 19, 59, 28, 61,  47, 37, 23, 30, 30, 33, 49, 29, 53,
                  45, 43, 34, 37, 34, 42, 39, 53, 41, 31, 49, 32,  69, 64, 60, 56, 49,55,
                  25),
  private_Allele = c(12,24, 9, 12, 30, 29, 1, 0, 3, 0, 0, 0, 2, 5, 26, 8, 11, 2, 1, 21, 24, 4, 2, 0, 0, 2, 0, 4, 12, 9, 13, 10, 8, 1)
)

# ------------------------------------------------
# Funktion zum Erstellen einzelner Plots
plot_corr <- function(data, x, y, y_label) {
  # lineares Modell
  model <- lm(data[[y]] ~ data[[x]], data = data)
  r2 <- summary(model)$r.squared
  pval <- summary(model)$coefficients[2,4]
  
  # Signifikanzcode
  stars <- ifelse(pval <= 0.001, "***",
                  ifelse(pval <= 0.01, "**",
                         ifelse(pval <= 0.05, "*", "ns")))
  
  ggplot(data, aes_string(x = x, y = y)) +
    geom_point(shape = 21, fill = "white", color = "black", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    annotate("text",
             x = max(data[[x]]) - 2,
             y = max(data[[y]]),
             label = paste0("R² = ", round(r2, 2), stars),
             size = 5,
             hjust = 1) +
    labs(x = "Elevation", y = y_label) +
    theme_bw(base_size = 14)
}

# ------------------------------------------------
# Drei Diagramme generieren
p1 <- plot_corr(df, "Latitude", "Ho", "Heterozygotie (Ho)")
p2 <- plot_corr(df, "Latitude", "Ar", "Allelenreichtum (Ar)")
p3 <- plot_corr(df, "Latitude", "Fis", "Inzucht (Fis)")
p7 <- plot_corr(df, "Latitude", "private_Allele", "Anzahl privater Allele")
p8 <- plot_corr(df, "Latitude", "Allelanzahl", "Allelanzahl")
p4 <- plot_corr(df, "Elevation", "Ho", "Heterozygotie (Ho)")
p5 <- plot_corr(df, "Elevation", "Ar", "Allelenreichtum (Ar)")
p6 <- plot_corr(df, "Elevation", "Fis", "Inzucht (Fis)")
p9 <- plot_corr(df, "Elevation", "private_Allele", "Anzahl privater Allele")
p10 <- plot_corr(df, "Elevation", "Allelanzahl", "Allelanzahl")
# ------------------------------------------------
# Alles untereinander anordnen
install.packages("patchwork")
library(patchwork)
final_plot <- p10 / p9 
final_plot2 <- p8 / p7 
final_plot
final_plot2

cite.pa
cite.pack("hierfstat")



median(c(0.434,
         0.449,
         0.558,
         0.456,
         0.324,
         0.597,
         0.632,
         0.611,
         0.580,
         0.611,
         0.602,
         0.511,
         0.341,
         0.500,
         0.562,
         0.441,
         0.387,
         0.555,
         0.785,
         0.525,
         0.684,
         0.338,
         0.491,
         0.324,
         0.385,
         0.403,
         0.363,
         0.460,
         0.321,
         0.068,
         0.505,
         0.258,
         0.511,
         0.392
))
median(c(0.534,
         0.402,
         0.571,
         0.474,
         0.545,
         0.623,
         0.757,
         0.622,
         0.552,
         0.521,
         0.643,
         0.518,
         0.357,
         0.678,
         0.750,
         0.571,
         0.433,
         0.506,
         0.750,
         0.593,
         0.767,
         0.458,
         0.448,
         0.535,
         0.349,
         0.443,
         0.389,
         0.553,
         0.307,
         0.088,
         0.553,
         0.263,
         0.554,
         0.555
))
median(c(0.187,
         -0.114,
         0.022,
         0.038,
         0.405,
         0.042,
         0.165,
         0.018,
         -0.051,
         -0.174,
         0.062,
         0.013,
         0.046,
         0.263,
         0.250,
         0.227,
         0.106,
         -0.095,
         -0.047,
         0.114,
         0.108,
         0.262,
         -0.096,
         0.395,
         -0.103,
         0.090,
         0.062,
         0.166,
         -0.045,
         0.225,
         0.087,
         0.019,
         0.077,
         0.292
))
median(c(53,
         41,
         31,
         49,
         32,
         69,
         64,
         60,
         56,
         49,
         61,
         47,
         37,
         23,
         30,
         30,
         33,
         49,
         29,
         53,
         45,
         43,
         34,
         37,
         34,
         42,
         39,
         45,
         29,
         19,
         59,
         28,
         55,
         25
))
median(c(2.6,
         2.0,
         1.5,
         2.5,
         1.4,
         2.5,
         2.6,
         2.4,
         1.6,
         2.4,
         2.2,
         1.5,
         1.8,
         1.0,
         1.4,
         1.3,
         1.4,
         2.1,
         1.4,
         1.6,
         1.3,
         1.9,
         1.4,
         1.4,
         1.4,
         1.5,
         1.4,
         1.4,
         1.1,
         1.0,
         3.1,
         1.1,
         2.8,
         1.2
))
