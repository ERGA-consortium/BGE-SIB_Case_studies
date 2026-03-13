## Script to make chromosome feature drawings ##

install.packages("chromoMap")
library("chromoMap")
setwd("~/Desktop/ArnicaChrSSR") ## replace text between "" by own path

## prepare data files 
# in Excel / programmatically
# see here for format 
# https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html#Getting_Started

## set data file names 
chr <- "ChromosomeFile.csv"
elm <- "ElementFile.csv"

## test format of files
#head(read.table(elm,sep = ";"))
# if separator not "\t", read in files 

chr <- list(read.table(chr,sep = ";"))
#elm2 <- read.table(elm,sep = ";")
#elm1 <- list(elm2[elm2$V5 != "gap"])
elm <- list(read.table(elm,sep = ";"))


## test
#chromoMap(chr,elm)

### adapt the graphical output
wfac <- 4
dcol <- list(c("orange","darkgrey","turquoise", "yellow"))
#acol <- list(c(rgb(1,0,0,0), "orange","turquoise", "yellow"))

chromoMap(chr,elm, chr_color = c("lightgrey"), #c("#2144FF"),
          #legend = T, lg_x = 100, lg_y = 250, 
          #labels=T, anno_col = acol,
          win.summary.display=T,
          n_win.factor = wfac,
          #fixed.window = T, window.size = wfac)
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = dcol, 
          # y_chr_scale = "vertical",#ch_gap = 10,
          export.options=T)
