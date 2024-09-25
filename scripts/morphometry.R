#Collection of morphometric data and statistical analyses - Myrcia neoobscura species complex (Myrtaceae)
#Author: Paulo Henrique Gaem, PPG Biologia Vegetal, Unicamp (phgaem@gmail.com)

#EXTRACTING LEAF SHAPES AS X,Y COORDINATES
library(Momocs)
lf <- list.files('Leaves_Binary', full.names=TRUE) #getting images paths
coo <- import_jpg(lf) #import images
coo2<-Out(coo) #transform in coo-type objects
panel(coo2,col='yellow') #inspect shapes
stack(coo2) #taking a family picture of raw unscaled outlines 

#SETTING LEAF APICES AND BASES AS LANDMARKS
mlkds = NULL
coof<-coo2
for(i in 1:length(coof[[1]])) {
  print(i)
  x = coof[[1]][[i]]
  g1 = which(x[,1]==min(x[,1]))
  g2 = x[g1,2]
  d2 = abs(g2-mean(g2))
  idx1 = g1[d2==min(d2)][1]
  
  g1 = which(x[,1]==max(x[,1]))
  g2 = x[g1,2]
  d2 = abs(g2-mean(g2))
  idx2 = g1[d2==min(d2)][1]
  g1 = c(idx1,idx2)
  plot(x)
  points(x[g1,1],x[g1,2],pch=21,col='red',cex=3)
  mlkds[[i]] = g1
}
coof[['ldk']] = mlkds

#ALIGNING SHAPES
a<-coo_bookstein(coof) #setting a new baseline for the shapes based on the landmarks
panel(a,col='aquamarine') #inspect shapes again
stack(a) #taking a family picture of raw scaled smoothed outlines 

#TRANSFORMING SHAPES INTO MATHEMATICAL VARIABLES
b<-efourier(a,norm=F) #performing Elliptical Fourier Transforms to shapes
c<-PCA(b) #performing principal component analysis to Fourier-transformed shapes

#PLOTTING PCA (FIRST TWO COMPONENTS)
plot_PCA(x=c, axes=c(1,2))

#GETTING A PLOT DEPICTING THE CONTRIBUTIONS OF THE FIRST 10 COMPONENTS
pdf(file="pccontribs-leaves.pdf",width=8,height=8) #create a .pdf file
PCcontrib(c, nax=1:10) #plotting the contributions
dev.off()

#ATTENTION: this step needs inspection of the PCs truly contributing taxonomically to leaf shape

#EXTRACTING THE INFORMATIVE COMPONENTS INTO AN OBJECT
inform=c("PC1","PC2","PC3","PC5") # [REQUIRES EDITION] creating a vector containing the names of the taxonomically informative components
eigen_inform=c$x[,inform] #extracting the taxonomically informative components
colnames(eigen_inform)<-sub("PC","LEAF_PC",colnames(eigen_inform)) #renaming variables


#####
#GETTING LEAF AREAS
library(EBImage)
pics <- sapply(lf, readImage, simplify=F) #loading images

areas <- vector("list", length=length(pics)) #creating an empty list
names(areas) <- sub(".jpg", "", lf, fixed=T) #naming objects with paths of files

for (i in 1:length(pics)) {
  pics[[i]] -> pic
  pic.data <- imageData(channel(pic, mode = "blue")) # extract the 'blue' channel
  pic.data <- 1-pic.data  # reverse the image 
  pic.data[pic.data < .7] <- 0
  
  pic.lab <- bwlabel(pic.data)  # attempt to enclose the holes
  
  # calculate the area 
  pic.shape <- computeFeatures.shape(pic.lab)
  pic.shape[,1] -> pic.area
  sort(pic.area, decreasing=T) -> pic.area
  pic.area <- pic.area[1:2]
  names(pic.area) <- c("Leaf (px)", "Scale (px)")
  area.sq <- (round(pic.area[1]/pic.area[2],2))
  names(area.sq)<-c("Leaf (cm2)")
  
  c(pic.area, area.sq) -> area.out
  areas[[i]] <- area.out
}

do.call(rbind, areas) -> areas #getting area data into a table
areas <- as.data.frame(areas)

path<-rownames(areas)
areas<-cbind(path,areas)

library(tidyverse)
areas <- areas %>% 
  mutate(path = sub(".jpg","", lf, fixed=T)) %>% 
  separate(path, sep = "Leaves_Binary/", into = c("path_short", "names"))
areas$path_short=NULL
rownames(areas)<-areas$names; areas$names=NULL
detach("package:tidyverse", unload=TRUE)

area_leaves<-areas$'Leaf (cm2)'
names(area_leaves)<-rownames(areas)

inform_leaves=merge(eigen_inform,area_leaves,by=0,all=T) #combining leaf descriptors with leaf area
colnames(inform_leaves)[6]<-"LEAF_sqcm" #renaming leaf area variable
rownames(inform_leaves)<-inform_leaves$Row.names; inform_leaves$Row.names=NULL #setting row names

#EXTRACTING FLOWER BUD SHAPES AS X,Y COORDINATES
#####
lf2 <- list.files('Flowers_Binary', full.names=TRUE) #getting images' paths
coo <- import_jpg(lf2) #import images
coo2<-Out(coo) #transform in coo-type objects
panel(coo2,col='cyan') #inspect shapes
stack(coo2) #taking a family picture of raw unscaled outlines 

#SETTING FLOWER BUD APEX AND BASE AS LANDMARK
mlkds = NULL
coof<-coo2
for(i in 1:length(coof[[1]])) {
  print(i)
  x = coof[[1]][[i]]
  g1 = which(x[,1]==min(x[,1]))
  g2 = x[g1,2]
  d2 = abs(g2-mean(g2))
  idx1 = g1[d2==min(d2)][1]
  
  g1 = which(x[,1]==max(x[,1]))
  g2 = x[g1,2]
  d2 = abs(g2-mean(g2))
  idx2 = g1[d2==min(d2)][1]
  g1 = c(idx1,idx2)
  plot(x)
  points(x[g1,1],x[g1,2],pch=21,col='red',cex=3)
  #return(g1)
  mlkds[[i]] = g1
}
coof[['ldk']] = mlkds

#ALIGNING SHAPES
a<-coo_bookstein(coof) #setting a new baseline for the shapes based on the landmarks
a<-coo_smooth(a,10) #smoothing shapes
panel(a,col='green') #inspect shapes again
stack(a) #taking a family picture of raw procrustes-scaled smoothed outlines 

#TRANSFORMING SHAPES INTO MATHEMATICAL VARIABLES
b<-efourier(a,norm=F) #performing Elliptical Fourier Transforms to shapes
c<-PCA(b) #performing principal component analysis to Fourier-transformed shapes

#PLOTTING PCA (FIRST TWO COMPONENTS)
plot_PCA(x=c, axes=c(1,2))

#GETTING A PLOT DECIPTING THE CONTRIBUTIONS OF THE FIRST 10 COMPONENTS
pdf(file="pccontribs-buds.pdf",width=8,height=8) #create a .pdf file
PCcontrib(c, nax=1:10) #plotting the contributions
dev.off()

#ATTENTION: this step needs inspection of the PCs truly contributing taxonomically to flower bud shape

#EXTRACTING THE INFORMATIVE COMPONENTS INTO AN OBJECT
inform=c("PC1","PC2","PC4","PC7") # [REQUIRES EDITION] creating a vector containing the names of the taxonomically informative components
eigen_inform2=c$x[,inform] #extracting the taxonomically informative components
colnames(eigen_inform2)<-sub("PC","BUD_PC",colnames(eigen_inform2)) #renaming variables

#GETTING FLOWER BUD AREAS
#####
pics <- sapply(lf2, readImage, simplify=F) #loading images

areas <- vector("list", length=length(pics)) #creating an empty list
names(areas) <- sub(".jpg", "", lf2, fixed=T) #naming objects with paths of files

for (i in 1:length(pics)) {
  pics[[i]] -> pic
  pic.data <- imageData(channel(pic, mode = "blue")) # extract the 'blue' channel
  pic.data <- 1-pic.data  # reverse the image 
  pic.data[pic.data < .7] <- 0
  
  pic.lab <- bwlabel(pic.data)  # attempt to enclose the holes
  
  # calculate the area 
  pic.shape <- computeFeatures.shape(pic.lab)
  pic.shape[,1] -> pic.area
  sort(pic.area, decreasing=T) -> pic.area
  pic.area <- pic.area[1:2]
  names(pic.area) <- c("Bud (px)", "Scale (px)")
  area.sq <- (round(pic.area[1]/pic.area[2],2))
  names(area.sq)<-c("Bud (mm2)")
  
  c(pic.area, area.sq) -> area.out
  areas[[i]] <- area.out
}

do.call(rbind, areas) -> areas #getting area data into a table
areas <- as.data.frame(areas)

path<-rownames(areas)
areas<-cbind(path,areas)

library(tidyverse)
areas <- areas %>% 
  mutate(path = sub(".jpg","", lf2, fixed=T)) %>% 
  separate(path, sep = "Flowers_Binary/", into = c("path_short", "names"))
areas$path_short=NULL
rownames(areas)<-areas$names; areas$names=NULL
detach("package:tidyverse", unload=TRUE)

area_buds<-areas$'Bud (mm2)'
names(area_buds)<-rownames(areas)

inform_buds=merge(eigen_inform2,area_buds,by=0,all=T) #combining leaf descriptors with leaf area
colnames(inform_buds)[6]<-"BUD_sqcm" #renaming leaf area variable
rownames(inform_buds)<-inform_buds$Row.names; inform_buds$Row.names=NULL #setting row names

#COMBINING LINEAR AND GEOMETRIC MORPHOMETRIC DATA IN A SINGLE OBJECT
#####
morph<-data.table::fread('data.txt',header=T,na.strings="") #calling morphological matrix
morph<-as.data.frame(morph)
rownames(morph)<-morph$ID; morph$ID=NULL #setting row names

morph2=merge(morph,inform_leaves,by=0,all=T) #combining morphological data with leaf descriptors
morph2<-arrange(morph2, id) #reorganising rows by a single identifier
rownames(morph2)<-morph2$Row.names; morph2$Row.names=NULL #setting row names

morph3=merge(morph2,inform_buds,by=0,all.x=T) #combining morphological data+leaf descriptors with bud descriptors
morph3<-arrange(morph3, id) #reorganising rows by a single identifier
rownames(morph3)<-morph3$Row.names; morph3$Row.names=NULL #setting row names

write.csv(morph3,file='morph_total.csv')


#####
#LOADING DATA
data=read.csv('morph_total.csv')
rownames(data)<-data$X; data$X=NULL #setting row names

#FILTERING SPECIES NAMES AND MORPH DATA
final.names<-data[,2]
data<-cbind(final.names,data[,12:ncol(data)])

#EXPLORING VARIATION OF RAW MORPHOMETRIC DATA

str(data) #some columns are vectors of character-type variables

chr_cols<-unlist(lapply(data, is.character)) #identifying what columns have character-type variables
data[chr_cols]<-lapply(data[chr_cols], as.factor) #turning character-type variables into factor-type variables

str(data) #now we only have numeric and factor-type variables

#calculating variables that are ratios
data$leaf.marg<-data$leaf.marg.vein0/data$leaf.half.width0
data$infl.ped.prop<-data$infl.ped0/data$infl.length1
data$infl.lat.prop<-data$infl.lat0/data$infl.length1
data$infl.bran.prop<-data$infl.branches0/data$infl.length1

#erasing variables that are not to be used
data$leaf.marg.vein0=NULL
data$leaf.half.width0=NULL
data$infl.ped0=NULL
data$infl.lat0=NULL
data$infl.branches0=NULL

str(data)

#removing insufficiently varying variables
data$fruit.diam=NULL
data$fruit.calyx.persist=NULL

#transforming some short-range varying integer into factorial variables
data$infl.quant<-as.factor(data$infl.quant)
data$infl.veg<-as.factor(data$infl.veg)
data$infl.compos<-as.factor(data$infl.compos)
data$hair.golden<-as.factor(data$hair.golden)
data$hair.rufous<-as.factor(data$hair.rufous)
data$hair.ochre<-as.factor(data$hair.ochre)
data$hair.white<-as.factor(data$hair.white)

str(data)

#PLOTTING NUMERIC VARIABLES IN BOXPLOTS

column.names<-c('Final names', 'Plant height (m)','Branch exfoliation','Internode length (mm)',
                'Branch indumentum','Midvein position','Number of lateral veins',
                'Leaf margin revolution','Leaf gland display','Basal bracts display',
                'Inflorescence quantity','Inflorescence agglomeration','Number of vegetative nodes above flowering node',
                'Inflorescence branching pattern', 'Infl. composition','Inflorescence length (mm)',
                'Bract length (mm)','Bract persistence','Number of flowers on basal infl. branch',
                'Infl. glands display','Infl. hair density','Infl. hair length', 'Infl. hair golden',
                'Infl. hair rufous','Infl. hair ochraceous','Infl. hair white','Bracteole shape',
                'Bracteole length (mm)','Bracteole persistence','Pedicel length (mm)',
                'Hypanthium indum.','Calyx indum.','LEAF PC1','LEAF PC2', 'LEAF PC3', 'LEAF PC5',
                'Leaf area (cm²)','Fl. bud PC1','Fl. bud PC2','Fl. bud PC4', 'Fl. bud PC7',
                'Flower bud area (mm²)','Marginal vein dist. from margin (% leaf width)','Peduncle length (% infl. main axis length)',
                'Lower lat. branch lth./infl. lth. (widewness index)','Number of lat. bran./infl. lth. (compaction index)'
                )

op <- par(family='serif')

for (i in 1:ncol(data)){
if (is.factor(data[[i]])!=T){
  tiff(paste(names(data[i]), "boxplot.tiff", sep="_"), units="mm", width=250, height=100,
       res=300)
  boxplot(data[[i]] ~ data$final.names, outpch=20,
          names=c('ARE', 'EXC'),
          xlab="Taxa", 
          ylab=column.names[i]
  )
  dev.off()
  }
}

par(op) 
 
#PREPARING DATASET FOR MORPHOMETRIC ANALYSES 
#SCALING NUMERIC VARS

data <- data[,2:ncol(data)] #removing column containing names

escala <- function(x){
  desvio <- sd(x, na.rm = TRUE)
  media <- mean(x, na.rm=TRUE)
  
  return((x-media)/desvio)
}

library(magrittr)
library(dplyr)

numeric.vars=c("height","branch.intern","leaf.lat.veins","infl.length1",
               "infl.br.length","infl.flowers","flower.br.length","flower.pedic",
               "LEAF_PC1","LEAF_PC2","LEAF_PC3", "LEAF_PC5","LEAF_sqcm",
               "BUD_PC1", "BUD_PC2","BUD_PC4", "BUD_PC7", "BUD_sqcm",
               "leaf.marg","infl.ped.prop","infl.lat.prop","infl.bran.prop")

datamut <- data %>%
  mutate_at(numeric.vars, escala)

glimpse(datamut)

#GETTING A DENDROGRAM
library(cluster)
dd<-daisy(datamut, metric='gower') #building a distance matrix using Gower's index
hc<-hclust(dd, method="ward.D2") #building a dendrogram using Ward's method
coph<-cophenetic(hc); cor(dd,coph) #calculating the cophenetic correlation coefficient

library(factoextra)
pdf(file=paste(Sys.Date(), 'dendrogram.pdf',sep="-"), width=13, height=13) #creating an empty pdf
fviz_dend(hc, cex = 0.5,
          type='circular') #getting a circular dendrogram with 12 automatic coloured groups
dev.off()

# Calculate silhouette width for many k using hclust
#####
sil_width <- c(NA)
for(i in 2:9){
  agnes_clust <- hcut(dd, k=i, isdiss=T, hc_func="hclust",
                      hc_method="ward.D2", hc_metric=NULL)
  
  sil_width[i] <- agnes_clust$silinfo$avg.width
}

# Plot sihouette width (higher is better)
tiff(filename=paste(Sys.Date(), 'silhouette_widths.tiff', sep='-'), height=300, width=600)
plot(1:9, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width",
     las=2,
     xaxt='n')
lines(1:9, sil_width)
axis(1, at=2:9)
dev.off()

#CLUSTERBOOT
#####
library(fpc)

clusters <- 3
nboots <- 1000
clus.boot <- clusterboot(dd, 
                         B=nboots, # Number of bootstrap resamples
                         distances=F, #input data contains a dissimilarity matrix
                         clustermethod=hclustCBI, # for hierarchical clustering 
                         method="ward.D2", # use what we used in "hclust"
                         k=clusters, 
                         count=T) #show progress on screen?


AvgJaccard <- clus.boot$bootmean
Instability <- clus.boot$bootbrd/nboots
Clusters <- c(1:clusters)
Eval <- cbind(Clusters, AvgJaccard, Instability)
Eval
write.csv(Eval,file='Eval.csv')

groups <- cutree(hc, k=clusters) # Assign observations to groups

hierarchical_clusters <- as.factor(groups) # First, make it a factor so that it is not confused with numeric data
species <- row.names(data)   # Extract row names from the working data. If we removed 
# any rows from the original data due to missingness, then
# we will need this to merge the clusters with the original data.
# Combine the two into a data frame
HC_results <- cbind(species, hierarchical_clusters) # Make cluster data

write.csv(HC_results,file='HC_results.csv')