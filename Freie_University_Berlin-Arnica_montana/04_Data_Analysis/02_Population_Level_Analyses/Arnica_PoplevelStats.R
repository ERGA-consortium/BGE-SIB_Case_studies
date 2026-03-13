setwd("~/Documents/Arbeit/Projekte/ArnicaMontana/ArnicaSNP-BGE/06_Genotyping/")

#install.packages("poppr")
library("poppr")

popdata <- read.genalex("Genotypes2_all.csv", sep=";")

View(popdata)
popdata
summary(popdata)
poppr(popdata)

#splitStrata(popdata) <- ~Pop

setwd("~/Documents/Arbeit/Projekte/ArnicaMontana/ArnicaSNP-BGE/10_Analysis/")

## get lists of the population and locus names
poplist <- unique(popdata$pop)
loclist <- unique(popdata$loc.fac)

## print out information on alleles per population (& locus)
for (i in poplist) {
  locdat <- locus_table(popdata, pop=i)
  print(i)
  print(locdat)
}

## get per-pop allele presence/absence matrix
alperpop <- rowsum(popdata@tab,popdata@pop, na.rm=TRUE)
alperpop2 <- alperpop>0

## pairwise shared alleles between populations
sharedalleles <- matrix(nrow=length(poplist), ncol=length(poplist))

for (i in seq(length(poplist))) {
  for (j in seq(length(poplist))) {
    sharedalleles[i,j] <- sum((alperpop2[i,]+alperpop2[j,])==2)
  }
}


## look at the output of "summary" to find out how to access the Heterozygosity values
#View(summary(popdata[popdata$pop=="AM10"]))

## calculate FIS per locus & population
Fistable <- data.frame(row.names=loclist)
for (i in seq(length(poplist))) {
  print(poplist[i])
  sumdat <- summary(popdata[popdata$pop==poplist[i]])
  print(sumdat)
  Fis = 1- sumdat$Hobs/sumdat$Hexp
  print(Fis)
  Fistable[[poplist[i]]] <- Fis
}
colnames(Fistable) <- poplist

## calculate variability ("Streuung") of FIS across loci per population
streu <- rep(0,length(Fistable))
for (i in seq(length(Fistable))) {
  streu[i] <- abs(min(Fistable[,i], na.rm =TRUE)) + max(Fistable[,i], na.rm =TRUE)
  print(streu[i])
  }

# Private alleles
pa <- rowSums(private_alleles(popdata))

# Allelic diversity
nalleles <- summary(popdata)$pop.n.all / summary(popdata)$n.by.pop #length(loclist)


## calculate genetic distances & save to file
# between individuals (Prevosti distance)
indist <- prevosti.dist(popdata)
indist_m <- as.matrix(indist)
write.table(indist_m, "AM_Indist.csv", sep=";", row.names = indNames(popdata), col.names=indNames(popdata))

# between populations (pairwise FST)
#install.packages("hierfstat")
library("hierfstat")

popdist <- genet.dist(popdata, diploid=TRUE, method="Fst")
popdist_m <- as.matrix(popdist)
write.table(popdist_m, "AM_Popdist.csv", sep=";")

# allelic richness
ar <- allelic.richness(popdata)
#allele.count(popdata)[1][allele.count(popdata)[1]>0]

#####
library(ape)
library(ggplot2)

#setwd("~/Documents/Arbeit/Projekte/ArnicaMontana/ArnicaSNP-BGE/02_Com/Presentations/Portugal")
dat <- read.table("AM_Popdist.csv", sep=";", dec=".", header=TRUE, na.strings = " ")

pops <- as.factor(rownames(dat))
Fdist <- as.dist(dat[,2:length(dat)])
Fpcoa <- pcoa(Fdist)

poplist <- sub("..-","",rownames(dat))
studylist <- as.factor(sub("-.*","",rownames(dat)))

Fpcoa$pops <- poplist

eigenplot <- ggplot(Fpcoa$values, aes(x=as.factor(seq(1,dim(Fpcoa$values)[1])), y=Rel_corr_eig)) + 
  geom_bar(stat="identity", fill=c(rep("#000000",2), rep("#999999",dim(Fpcoa$values)[1]-2))) + 
  labs(x= "Axis", y="Explained variance") + theme_light() 

Fpcoa2 <- as.data.frame(Fpcoa$vectors)
Fpcoa2$pops <- factor(poplist, levels=poplist)

Fpcoa2$study <- studylist
#Fpcoa2$symbols <- symbollist
cc <- c(scales::seq_gradient_pal("lightgreen", "forestgreen", "Lab")(seq(0,1,length.out=sum(studylist=="SB"))), 
        scales::seq_gradient_pal("skyblue", "blue", "Lab")(seq(0,1,length.out=sum(studylist=="SC")/2-0.5)),
        scales::seq_gradient_pal("turquoise", "darkblue", "Lab")(seq(0,1,length.out=sum(studylist=="SC")/2)), 
        scales::seq_gradient_pal("yellow", "darkorange", "Lab")(seq(0,1,length.out=sum(studylist=="SD"))), 
        scales::seq_gradient_pal("red", "darkred", "Lab")(seq(0,1,length.out=sum(studylist=="SR"))))

pcoaplot <- ggplot(Fpcoa2, aes(x=Axis.1, y=Axis.2, color=pops, 
                               shape=pops, stroke=2, cex=1.5)) + 
  labs(x= paste("Axis 1:", round(Fpcoa$values$Rel_corr_eig[1]*100, 2), "%"), 
       y=paste("Axis 2:", round(Fpcoa$values$Rel_corr_eig[2]*100, 2), "%")) + 
  geom_point() + theme_light() + #+ theme(legend.position = "none") +
  scale_shape_manual(values = rep(c(21, 22, 23, 24, 25),7)) +
  scale_color_manual(values=cc)
  #geom_text(aes(label = pops), nudge_x=0.01, nudge_y=0.01) 

pcoaplot

png("AM_FST_PcoA.png", width=10, height=7, units="in", res=600)
print(pcoaplot)
dev.off()

png("AM_FST_PcoA_eig.png", width=10, height=2, units="in", res=600)
print(eigenplot)
dev.off()

###
datout <- data.frame(poplist)
datout$A1 <- Fpcoa2$Axis.1
datout$A2 <- Fpcoa2$Axis.2
datout$A3 <- Fpcoa2$Axis.3
datout$rFis2 <- streu/2
datout$palleles <- pa
datout$nalleles <- nalleles
datout <- cbind(datout, alperpop)
write.table(datout, "AM_PopPCoA.csv", sep=";", row.names=FALSE)

###
library(igraph)
net <- graph_from_adjacency_matrix(as.matrix(Fdist), 
                                       weighted=TRUE, mode="undirected", diag=F)
net2 <- graph_from_adjacency_matrix(sharedalleles/max(sharedalleles), 
                                       weighted=TRUE, mode="undirected", diag=F)

png("AM_FST_heatmap.png", width=10, height=10, units="in", res=600)
heatmap(as.matrix(Fdist), labRow=poplist, labCol=poplist, symm=TRUE, main="Pairwise FST")
dev.off()

png("AM_lowFST.png", width=10, height=10, units="in", res=600)
plot(net, layout=layout.circle, main=expression("Highest pairwise similarity 1-F"["ST"]), # layout_with_mds(network)
     vertex.label.cex=0.7, edge.width=(1-edge_attr(net,"weight"))^16*50,
     vertex.label.color="black", edge.color="red", vertex.label=poplist,
     vertex.color=c("limegreen", "skyblue", "yellow", "magenta")[as.factor(studylist)]
) 
dev.off()

png("AM_highFST.png", width=10, height=10, units="in", res=600)
plot(net, layout=layout.circle, main=expression("Highest pairwise dissimilarity F"["ST"]), # layout_with_mds(network)
     vertex.label.cex=0.7, edge.width=(edge_attr(net,"weight"))^16*50, 
     vertex.label.color="black", edge.color="grey", vertex.label=poplist,
     vertex.color=c("limegreen", "skyblue", "yellow", "magenta")[as.factor(studylist)]
) 
dev.off()

png("AM_sharedalleles.png", width=10, height=10, units="in", res=600)
plot(net2, layout=layout.circle, main="Number of shared alleles", # layout_with_mds(network)
     vertex.label.cex=0.7, edge.width=(edge_attr(net2,"weight")/max(edge_attr(net2,"weight")))^16*20,
     vertex.label.color="black", edge.color="red", vertex.label=poplist,
     vertex.color=c("limegreen", "skyblue", "yellow", "magenta")[as.factor(studylist)]
) 
dev.off()

png("AM_sharedalleles_heatmap.png", width=10, height=10, units="in", res=600)
heatmap(sharedalleles, labRow=poplist, labCol=poplist, symm=TRUE, main="Pairwise number of shared alleles")#, cexRow=0.5, cexCol=0.5)
dev.off()

######

oridata <- read.genalex("~/Documents/Arbeit/Projekte/ArnicaMontana/ArnicaSNP-BGE/09_Data/20241222_Arnica_SSRlenGenAlEx.csv", sep=";")

View(oridata)
oridata
summary(oridata)

## get lists of the population and locus names
opoplist <- unique(oridata$pop)
oloclist <- unique(oridata$loc.fac)

mlg(oridata)
oristat <- poppr(oridata)
er <- (oristat$MLG -1)/(oristat$N-1)

oFistable <- data.frame(row.names=oloclist)
for (i in seq(length(opoplist))) {
  print(opoplist[i])
  sumdat <- summary(oridata[oridata$pop==opoplist[i]])
  print(sumdat)
  Fis = 1- sumdat$Hobs/sumdat$Hexp
  print(Fis)
  oFistable[[opoplist[i]]] <- Fis
}
colnames(oFistable) <- opoplist

ostreu <- rep(0,length(oFistable))
for (i in seq(length(oFistable))) {
  ostreu[i] <- abs(min(oFistable[,i], na.rm =TRUE)) + max(oFistable[,i], na.rm =TRUE)
  print(ostreu[i])
}
ostreuout <- data.frame(row.names=opoplist)
ostreuout$streu <- ostreu/2
ostreuout$R <- er[0:length(opoplist)]
write.table(ostreuout, "Arnica_SSRlen_FisStreu2.csv", sep=";", col.names = NA)

