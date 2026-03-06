
plink.dist <- read.table(file = "...", header = FALSE, sep = "\t", dec = ".")

plink.id <- as.data.frame(read.table(file = "...", header = FALSE, sep = "\t", dec = "."))

metadata <- as.data.frame(read.csv(file = "..."))

plink.id$Species <- metadata$Species[match(plink.id$V1, metadata$Sample)]

str(plink.dist)

plink.matrix <- as.matrix(plink.dist)

dimnames(plink.matrix) <- list(plink.id[,3],plink.id[,3])

plink.matrix

plink.dist.df<-  as_tibble(plink.matrix, rownames="A") %>%
  pivot_longer(-A,names_to = "B", values_to = "distances") 

plink.dist.df %>%
  ggplot(aes(x=A, y=B, fill=distances)) +
  geom_tile()

plink.dist.df %>%
  ggplot(aes(x=A, y=B, fill=distances)) +
  geom_tile()+
  coord_equal()+
  scale_fill_gradient(low = "#FF0000", high = "#FFFFFF", name=NULL) +
  labs(x="Samples", y="Samples") +
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=3), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





NJ.tree <- ape::nj(plink.matrix)

plot(x=NJ.tree)

write.tree(NJ.tree, file = "NJ.tree.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)

tree <- read.tree("...")

tree_data <- NJ.tree$data

metadata <- read.csv("...")

ggplot(tree) + geom_tree() + theme_tree()+ geom_tiplab()

ggtree(tree)+ geom_text(aes(label=node), hjust=-.8) 

rerooted_tree <- root(tree, outgroup = "...", resolve.root = TRUE)

ggtree(rerooted_tree)+ geom_text(aes(label=node), hjust=-.8) 

ggplot(rerooted_tree) + geom_tree() + theme_tree()+ geom_tiplab()

ggtree(rerooted_tree) + 
  geom_tiplab() + 
  geom_cladelabel(node=21, label="...", 
                  color="...", offset=5000000, align=FALSE) + 
  geom_cladelabel(node=33, label="...", 
                  color="...", offset=8000000, align=FALSE) + 
  geom_cladelabel(node=19, label="...", 
                  color="...", offset=5000000, align=FALSE) + 
  geom_cladelabel(node=14, label="...", 
                  color="...", offset=5000000, align=FALSE) + 
  geom_cladelabel(node=44, label="...", 
                  color="...", offset=5000000, align=FALSE) + 
  theme_tree2() + 
  xlim(0, 60000000) + 
  theme_tree()



eigenvec <- as.data.frame(read_delim("...", delim = " ", col_names = FALSE))
eigenval <- read_delim("...", delim = " ", col_names = FALSE)

metadata <- as.data.frame(read.csv(file = "..."))

eigenvalue_percent <- round((eigenval$X1/sum(eigenval$X1))*100)

eigenvec$Species <- metadata$Species[match(eigenvec$X1, metadata$Sample)]


ggplot(data = eigenvec) +
  geom_point(mapping = aes(x = X3, y = X4, colour=Species), size = 3, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") + 
  labs(title = "PCA of Croatian sand gobies",
       x = paste0("Principal component 1 (",eigenvalue_percent[[1]]," %)"),
       y = paste0("Principal component 2 (",eigenvalue_percent[[2]]," %)"), 
       colour = "Species") +
  theme_minimal() 



samples<-read.table("...")[,1]
print(samples)


runs<-list()

for (i in 1:9){
  runs[[i]]<-read.table(paste0("....", i, ".Q"))
}

par(mar = c(6, 4, 2, 2) + 0.2) 
par(mfrow = c(2,3)) 
palt <- brewer.pal(12, "Set3")

for(i in 1:6){
  barplot(t(as.matrix(runs[[i]])),xlim = c(0,20), width=0.8,main = paste0("K=", i), col=palt[c(1:as.numeric(i+1))], ylab="Ancestry", border="black", names.arg = samples, cex.names = 0.4, las=1)
}

for(i in 1:9){
  runs[[i]]<-runs[[i]][c(2:10,12:14,11,19:26,29,17:18,27:28,30,15:16,1,31:38),]
}

samples_reorganized <- samples[c(2:10,12:14,11,19:26,29,17:18,27:28,30,15:16,1,31:38)]

for(i in 1:6){
  barplot(t(as.matrix(runs[[i]])),xlim = c(0,20), width=0.4,main = paste0("K=", i), col=palt[c(2:as.numeric(i+1))], ylab="Ancestry", border="black", names.arg = samples_reorganized, cex.names = 0.4, las=1)
}

for(i in 7:9){
  barplot(t(as.matrix(runs[[i]])),xlim = c(0,20), width=0.4,main = paste0("K=", i), col=palt[c(2:as.numeric(i+1))], ylab="Ancestry", border="black", names.arg = samples_reorganized, cex.names = 0.4, las=1)
}
