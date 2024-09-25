#Phenological patterns of the species of the Myrcia neoobscura species complex (Myrtaceae)
# Paulo Henrique Gaem - PPG Biologial Vegetal - University of Campinas (phgaem@gmail.com)

# Loading packages
library(monographaR)

# loading data
data<-data.table::fread('data.txt',header=T,na.strings="")

#transforming data for analysis
pheno<-as.data.frame(cbind(data$`Final name`,data$`Coll. month`,data$Phenology))
pheno$Species<-as.factor(pheno$V1); pheno$Month<-as.numeric(pheno$V2); pheno$Phenology<-as.factor(pheno$V3)
pheno$V1=NULL; pheno$V2=NULL; pheno$V3=NULL

str(pheno) #checking data

#running the command
phenoHist(pheno, mfrow=c(2,2), shrink=1.2, axis.cex=1.5, title.cex=1.5, 
          pdf=TRUE)