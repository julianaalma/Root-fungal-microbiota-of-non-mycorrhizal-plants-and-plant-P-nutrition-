####Code for Figure 5 and 6 
###Heatmap for OTUs associated to higher plant P accumulation (Fig. 5), 
####dendogram of these OTUs, and venn diagram  (Fig. 6)

###libraries needed####
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(tibble)
library(vegan)
library(ade4)
library(devtools)
library(pairwiseAdonis)
library(RSQLite)    
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(car)
library(psych)
library(corrplot)
library(MASS)
library(agricolae)
library(stats)
library(RColorBrewer)
library(corrplot)
library(psych)
library(cowplot)
library(stringr)
library(adegraphics)
library(DescTools)
library(metagMisc)
library(ggpubr)
library(forcats)
library(ggord)
library(randomForest)
library(paletteer)
library(ggnewscale)
library(ape)
library(ggvenn)
library(ggVennDiagram)

#set wd
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Machine_Learning_Paccum_prediction_Fig6")

###import data###
otu_mat<- read_excel("tableitsx_phyloseq_corrected.xlsx", sheet = "Table OTU plants")
tax_mat<- read_excel("tableitsx_phyloseq_corrected.xlsx", sheet = "Taxonomy")
samples_df <- read_excel("tableitsx_phyloseq_corrected.xlsx", sheet = "Metadata plants")
###format file###
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

phy_obj_plant<- phyloseq(OTU, TAX, samples)
phy_obj_plant<- prune_species(speciesSums(phy_obj_plant) > 0, phy_obj_plant) ### to discard OTUs with rel. abu. in roots = 0 (only present in soil samples)


###transform data in log10(RA+1)######
##relative abundance for each OTU abundOTU/abundtotOTU
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
phy_obj_plant_ra<-otu_table(phy_obj_plant_ra)+1
phy_obj_plant_log<-log10(otu_table(phy_obj_plant_ra))
###reassemble otu_table with taxonomy and metadata
phy_obj_plant_log<- phyloseq(otu_table(phy_obj_plant_log), TAX, samples)



############ Analyses: RF + Spearman + Capscale to identify best predictors of plant P accumulation #####

#####Make files#####
##### 2% occurence file: keep only most OTUs present in more than 2% of the samples
phy_obj_plant_log0.02<-phyloseq_filter_prevalence(
  phy_obj_plant_log,
  prev.trh = 0.02,
  abund.trh = NULL) #4748 OTUs

#get the mean relative abundance (RA) for each OTU
tax.mean <- taxa_sums(phy_obj_plant_log0.02)/nsamples(phy_obj_plant_log0.02)

#prune to keep OTUs with mean relative abundance (RA) > 0.0001 <=> log10(RA+1) > 0.0000434
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02) #1132 OTUs
class(sites.prune)
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)

#####Make Alpine and Lowland datasets for AM mycorrhizal (myc) and non-AM non-mycorrhizal (non-myc) plants
#remove "0" only taxa
phy_obj_plant_log0.02alpine<-subset_samples(phy_obj_plant_log0.02,Type_site=="Alpine")
phy_obj_plant_log0.02alpineno0<- prune_species(speciesSums(phy_obj_plant_log0.02alpine) > 0, phy_obj_plant_log0.02alpine)

phy_obj_plant_log0.02alpinenonmyc<-subset_samples(phy_obj_plant_log0.02alpine,Status=="Non-mycorhizal")
phy_obj_plant_log0.02alpinenonmycno0<- prune_species(speciesSums(phy_obj_plant_log0.02alpinenonmyc) > 0, phy_obj_plant_log0.02alpinenonmyc)

phy_obj_plant_log0.02alpinemyc<-subset_samples(phy_obj_plant_log0.02alpine,Status=="Mycorhizal")
phy_obj_plant_log0.02alpinemycno0<- prune_species(speciesSums(phy_obj_plant_log0.02alpinemyc) > 0, phy_obj_plant_log0.02alpinemyc)


phy_obj_plant_log0.02plaine<-subset_samples(phy_obj_plant_log0.02,Type_site=="Plaine")
phy_obj_plant_log0.02plaineno0<- prune_species(speciesSums(phy_obj_plant_log0.02plaine) > 0, phy_obj_plant_log0.02plaine)


phy_obj_plant_log0.02plainenonmyc<-subset_samples(phy_obj_plant_log0.02plaine,Status=="Non-mycorhizal")
phy_obj_plant_log0.02plainenonmycno0<- prune_species(speciesSums(phy_obj_plant_log0.02plainenonmyc) > 0, phy_obj_plant_log0.02plainenonmyc)

phy_obj_plant_log0.02plainemyc<-subset_samples(phy_obj_plant_log0.02plaine,Status=="Mycorhizal")
phy_obj_plant_log0.02plainemycno0<- prune_species(speciesSums(phy_obj_plant_log0.02plainemyc) > 0, phy_obj_plant_log0.02plainemyc)


#######################################################

#####Analysis for non-myc plants on alpine sites ######

### Sub Select only taxa for which mean relative abundance (RA) > 0.0001 <=> log10(RA+1) > 0.0000434
tax.mean <- taxa_sums(phy_obj_plant_log0.02alpinenonmycno0)/nsamples(phy_obj_plant_log0.02alpinenonmycno0)
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02alpinenonmycno0)

### Make a dataframe
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)
tab0.02<-data.frame(tab0.02)
s<-sample_data(sites.prune)
total<-cbind(tab0.02,s)


#### Random Forest analysis per Biome and Status ####

###For numeric vector Pplant/soil###
view(s)
y <- as.numeric(s$Pplant.soil)
x <- data.matrix(tab0.02)
id <- seq(1:length(y))
x2<-colnames(tab0.02)
numsets <- 10

#actual y representing 10 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)

#predictions are probabilities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)

#create data.frame for importance of each otu
rf.importance<-data.frame()
rf.imp.matrix <- matrix(nrow=length(x2),ncol=1)

#Metrics for model accuracy
rf.var.explained.matrix<-matrix(nrow=5,ncol=numsets)

###Model
par(new=T)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  subset <- matrix(0,nrow=72,ncol=5)
  subset[1:14,1] <- sort(sample(id,14))
  subset[1:14,2] <- sort(sample(which(!(id %in% subset)),14))
  subset[1:14,3] <- sort(sample(which(!(id %in% subset)),14))
  subset[1:15,4] <- sort(sample(which(!(id %in% subset)),15))
  subset[1:15,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  rf.importance<-NULL
  rf.var.explained<-NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x[-subset[,i],]
    #remove OTUs that have all zeroes in training set
    x.test <- x[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=y.train,ntree=100,importance=T)
    rf.pred <- c(rf.pred, predict(rf.train,x.test))
    #predictions are added as a queue with all samples being predicted one in each the 5 CV loops (here 72 predictions per loop). At the end we obtain 10 columns with 1 col for each of the iterations
    rf.importance<-cbind(rf.importance,rf.train$importance[,1])
    rf.var.explained<-c(rf.var.explained,(rf.train$rsq[100]))
    plot(rf.train)
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  x.sample <- x[y.sample,]
  rf.imp.matrix<-cbind(rf.imp.matrix,rf.importance)
  row.names(rf.imp.matrix)<-row.names(rf.importance)
  rf.var.explained.matrix[,n]<-rf.var.explained
}

##### Assess algorithm performance by comparing real vs predicted values

#Test model performance by comparing predicted vs real data that he has never seen in the test sets, for each independent iteration round (with 10 iteration rounds)
#calculate "average variance explained" from each iteration round (72 predictions x 10 independent iterations = 10 columns)
#cor(rf.pred.matrix[,1], y.test.matrix[,1])^2 # R Squared for each pred vs observed column
#cor(rf.pred.matrix[,2], y.test.matrix[,2])^2 # R Squared for each pred vs observed column, and so on...
#get the diagonal
NM_Alp_mean_R2_test<-mean(diag(cor(rf.pred.matrix,y.test.matrix)^2))
NM_Alp_mean_R2_test
#~68.6% of the variability in plant shoot P accumulation could be explained by the model (prediction test)

#visualize the correlation between predictions and actual values, from the 10 iterations
plot(x=rf.pred.matrix,y=y.test.matrix, xlab="Predicted shoot P accumulation", ylab="Actual shoot P accumulation in test set")
#add line if perfect prediction
abline(a=0, b=1, col="red", lwd = 3)


###Create values dataframes
#Each variable (OTU) has 50 values for importance (5 CV loops X 10 independent iterations)
rf.imp.df<-as.data.frame(rf.imp.matrix)
rf.imp.df$mean<-rowMeans(rf.imp.df[,2:51])
rf.imp.df$sd <- apply(rf.imp.df[,2:51], 1, sd, na.rm=TRUE)
rf.imp.df<-rf.imp.df[order(rf.imp.df$mean,decreasing=TRUE),]

###Threshold of 100 most important OTUs to keep based on Mean Decrase Accuracy
plot(rf.imp.df$mean)
grid(nx=10,ny=10)
abline(v=100,col="red")
rf.imp.df.100<-rf.imp.df[1:100,]

###Add taxonomical information to OTUs of interest
tax<-tax_table(phy_obj_plant_log0.02alpinenonmycno0)
tax<-data.frame(tax)
RF_non_myc_alpine_taxo<-merge(rf.imp.df.100,tax,by="row.names")
RF_non_myc_alpine_taxo$name<-paste(RF_non_myc_alpine_taxo$Row.names,RF_non_myc_alpine_taxo$Genus,sep="_")


###### Spearman correlation: Paccum ~ RA OTUs ##### 

##Make correlations with BH corrections for multiple testing
corspearman<-corr.test(tab0.02,total$Pplant.soil,method="spearman",adjust="BH")
corspearman1<-cbind(corspearman$p.adj,corspearman$r)
corspearman1<-data.frame(corspearman1)
corspearman1<-corspearman1 %>%
  rename(
    p.adj=X1,
    rho=X2
  )
cp2<-data.frame()

### Keep only significant OTUs with positive correlations
cp2<-corspearman1[corspearman1$p.adj<0.05,]
cp3<-cp2[cp2$rho>0,]
head(cp3[1:10,])
view(cp3)
nrow(cp3)
#89 OTUS


###### CAPscale analysis ######

###Capscale function from vegan used to perfrom db-rda. Makes an NMDS on Bray-curtis dissimilarities then constrained analysis using P accumulation
capscaletry<- capscale(tab0.02 ~ Pplant.soil, data=total,dist="bray")
anova(capscaletry)
summary(capscaletry)

###extract sample scores (capscaletry) for plotting
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))
total$CAP1 <- x$CAP1
total$MDS1 <- x$MDS1

###scores OTUS
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))

###scores Paccum added
scoreotus<-rbind(scoresotus,scores(capscaletry,display = "bp"))

###plot data
ggord(x,scoreotus,grp_in=total$Status,ellipse=TRUE, ellipse_pro=0.95,xlims=c(-3, 3),
      ylim=c(-3, 3))+geom_vline(xintercept = -0.8220, linetype="dotted")+
  geom_vline(xintercept = 0.8155, linetype="dotted")+
  geom_vline(xintercept = 0.1, linetype="dotted")+
  geom_vline(xintercept = -0.1, linetype="dotted")

###Extract OTUs and cap1 scores
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))
#get negative CAP1 scores = correlation with Pplant.soil (negative vector)
arrangedscorsotus<-scoresotus %>% filter(CAP1 <0)
arrangedscorsotus<-arrangedscorsotus %>% arrange(CAP1)
nrow(arrangedscorsotus)
#238 OTUs

# add taxonomy to selected OTUS
taxspeciesscores1<-merge(tax,arrangedscorsotus,by="row.names")

# #plot distribution and top 200
# taxspeciesscores1<-taxspeciesscores1 %>% arrange(CAP1)
# plot(abs(taxspeciesscores1$CAP1))
# grid(nx=10,ny=10)
# abline(v=200,col="red")

######Merge results######
#assemble tables
tabtoprintrf<-merge(RF_non_myc_alpine_taxo,taxspeciesscores1,by="Row.names")
cp3$Row.names<-rownames(cp3)
tabtoprintrf<-merge(tabtoprintrf,cp3,by="Row.names")
nrow(tabtoprintrf)
#57 OTUS

#export CSV
tabtoprintrf<-tabtoprintrf %>%
  rename(OTU=Row.names)
tabtoprintrf$name1<-paste(tabtoprintrf$OTU,tabtoprintrf$Species.x,sep="_")  
write.csv2(tabtoprintrf,"Capwilx_RF_Spearman_nonmyc_alpinegoodone_JA.csv")

##########################################################################

#####Analysis for myc plants on alpine sites ######

### Sub Select only taxa for which mean RA is > than 0.0001 -> 0.0000434 = log10(0.0001+1)
tax.mean <- taxa_sums(phy_obj_plant_log0.02alpinemycno0)/nsamples(phy_obj_plant_log0.02alpinemycno0)
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02alpinemycno0)

### Make a dataframe
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)
tab0.02<-data.frame(tab0.02)
s<-sample_data(sites.prune)
total<-cbind(tab0.02,s)

#### Random Forest analysis per Biome and Status ####

###For numeric vector Pplant/soil###
view(s)
y <- as.numeric(s$Pplant.soil)
x <- data.matrix(tab0.02)
id <- seq(1:length(y))
x2<-colnames(tab0.02)
numsets <- 10

#actual y representing 10 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)

#predictions are probabilities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)

#create data.frame for importance of each otu
rf.importance<-data.frame()
rf.imp.matrix <- matrix(nrow=length(x2),ncol=1)

#Metrics for model accuracy
rf.var.explained.matrix<-matrix(nrow=5,ncol=numsets)

###Model
par(new=T)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  subset <- matrix(0,nrow=72,ncol=5)
  subset[1:14,1] <- sort(sample(id,14))
  subset[1:14,2] <- sort(sample(which(!(id %in% subset)),14))
  subset[1:14,3] <- sort(sample(which(!(id %in% subset)),14))
  subset[1:15,4] <- sort(sample(which(!(id %in% subset)),15))
  subset[1:15,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  rf.importance<-NULL
  rf.var.explained<-NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x[-subset[,i],]
    #remove OTUs that have all zeroes in training set
    x.test <- x[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=y.train,ntree=100,importance=T)
    rf.pred <- c(rf.pred, predict(rf.train,x.test))
    rf.importance<-cbind(rf.importance,rf.train$importance[,1])
    rf.var.explained<-c(rf.var.explained,(rf.train$rsq[100]))
    plot(rf.train)
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  x.sample <- x[y.sample,]
  rf.imp.matrix<-cbind(rf.imp.matrix,rf.importance)
  row.names(rf.imp.matrix)<-row.names(rf.importance)
  rf.var.explained.matrix[,n]<-rf.var.explained
}

##### Asses algorithm performance by comparing real vs predicted values

#Test model performance by comparing predicted vs real data that he has never seen in the test sets, for each independent iteration round (with 10 iteration rounds)
#calculate "average variance explained" from each iteration round (72 predictions x 10 independent iterations = 10 columns)
#cor(rf.pred.matrix[,1], y.test.matrix[,1])^2 # R Squared for each pred vs observed column
#cor(rf.pred.matrix[,2], y.test.matrix[,2])^2 # R Squared for each pred vs observed column, and so on...
#get the diagonal
Myc_Alp_mean_R2_test<-mean(diag(cor(rf.pred.matrix,y.test.matrix)^2))
Myc_Alp_mean_R2_test
#~78.8% of the variability in plant shoot P accumulation could be explained by the model (prediction test)

#visualize all predictions from 10 iterations
plot(x=rf.pred.matrix,y=y.test.matrix, xlab="Predicted shoot P accumulation", ylab="Actual shoot P accumulation in test set")
#add line if perfect prediction
abline(a=0, b=1, col="red", lwd = 3)


###Create values dataframes
rf.imp.df<-as.data.frame(rf.imp.matrix)
rf.imp.df$mean<-rowMeans(rf.imp.df[,2:51])
rf.imp.df$sd <- apply(rf.imp.df[,2:51], 1, sd, na.rm=TRUE)
rf.imp.df<-rf.imp.df[order(rf.imp.df$mean,decreasing=TRUE),]

###Threshold of OTUs to keep
plot(rf.imp.df$mean)
grid(nx=10,ny=10)
abline(v=100,col="red")
rf.imp.df.100<-rf.imp.df[1:100,]

###Add taxonomical information to OTUs of interest
tax<-tax_table(phy_obj_plant_log0.02alpinemycno0)
tax<-data.frame(tax)
RF_myc_alpine_taxo<-merge(rf.imp.df.100,tax,by="row.names")
RF_myc_alpine_taxo$name<-paste(RF_myc_alpine_taxo$Row.names,RF_myc_alpine_taxo$Genus,sep="_")


###### Spearman correlation: Paccum ~ RA OTUs##### 

##Make correlations with BH corrections for multiple testing
corspearman<-corr.test(tab0.02,total$Pplant.soil,method="spearman",adjust="BH")
corspearman1<-cbind(corspearman$p.adj,corspearman$r)
corspearman1<-data.frame(corspearman1)
corspearman1<-corspearman1 %>%
  rename(
    p.adj=X1,
    rho=X2
  )
cp2<-data.frame()

###Keep only significant OTUs with positive correlations
cp2<-corspearman1[corspearman1$p.adj<0.05,]
cp3<-cp2[cp2$rho>0,]
head(cp3[1:10,])
view(cp3)
nrow(cp3)
#162 OTUs


###### CAPscale analysis ######

###Capscale function from vegan used to perfrom db-rda. Makes an NMDS on Bray-curtis dissimilarities then constrained analysis using P accumulation
capscaletry<- capscale(tab0.02 ~ Pplant.soil, data=total,dist="bray")
anova(capscaletry)
summary(capscaletry)

###extract sample scores (capscaletry) for plotting
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))
total$CAP1 <- x$CAP1
total$MDS1 <- x$MDS1

###scores OTUS
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))

###scores Paccum added
scoreotus<-rbind(scoresotus,scores(capscaletry,display = "bp"))

###plot data
ggord(x,scoreotus,grp_in=total$Status,ellipse=TRUE, ellipse_pro=0.95,xlims=c(-3, 3),
      ylim=c(-3, 3))+geom_vline(xintercept = -0.8220, linetype="dotted")+
  geom_vline(xintercept = 0.8155, linetype="dotted")+
  geom_vline(xintercept = 0.1, linetype="dotted")+
  geom_vline(xintercept = -0.1, linetype="dotted")

###Extract OTUs and cap1 scores
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))
#get positive CAP1 scores = correlation with Pplant.soil (positive vector)
arrangedscorsotus<-scoresotus %>% filter(CAP1 >0)
arrangedscorsotus<-arrangedscorsotus %>% arrange(desc(CAP1))
nrow(arrangedscorsotus)
#258 OTUs

# add taxonomy to selected OTUS
taxspeciesscores1<-merge(tax,arrangedscorsotus,by="row.names")

#plot distribution and top 200
taxspeciesscores1<-taxspeciesscores1 %>% arrange(desc(CAP1))
plot(abs(taxspeciesscores1$CAP1))
grid(nx=10,ny=10)
abline(v=200,col="red")

######Merge results######
#assemble tables
tabtoprintrf<-merge(RF_myc_alpine_taxo,taxspeciesscores1,by="Row.names")
cp3$Row.names<-rownames(cp3)
tabtoprintrf<-merge(tabtoprintrf,cp3,by="Row.names")
nrow(tabtoprintrf)
#70 OTUs

#export CSV
tabtoprintrf$name1<-paste(tabtoprintrf$OTU,tabtoprintrf$Species.x,sep="_")  
write.csv2(tabtoprintrf,"Capwilx_RF_Spearman_myc_alpinegoodone_JA.csv")

#############################################################################

#####Analysis for non-myc plants on lowland sites ######

### Sub Select only taxa for which mean RA is > than 0.0001 -> 0.0000434 = log10(0.0001+1)
tax.mean <- taxa_sums(phy_obj_plant_log0.02plainenonmycno0)/nsamples(phy_obj_plant_log0.02plainenonmycno0)
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02plainenonmycno0)#LOG10+1 de 0.0001 

### Make a dataframe
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)
tab0.02<-data.frame(tab0.02)
s<-sample_data(sites.prune)
total<-cbind(tab0.02,s)


#### Random Forest analysis per Biome and Status ####

###For numeric vector Pplant/soil###
view(s)
y <- as.numeric(s$Pplant.soil)
x <- data.matrix(tab0.02)
id <- seq(1:length(y))
x2<-colnames(tab0.02)
numsets <- 10

#actual y representing 10 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)

#predictions are probabilities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)

#create data.frame for importance of each otu
rf.importance<-data.frame()
rf.imp.matrix <- matrix(nrow=length(x2),ncol=1)

#Metrics for model accuracy
rf.var.explained.matrix<-matrix(nrow=5,ncol=numsets)

###Model
par(new=T)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  subset <- matrix(0,nrow=55,ncol=5)
  #correction: nrow=55 not 54
  subset[1:11,1] <- sort(sample(id,11))
  subset[1:11,2] <- sort(sample(which(!(id %in% subset)),11))
  subset[1:11,3] <- sort(sample(which(!(id %in% subset)),11))
  subset[1:11,4] <- sort(sample(which(!(id %in% subset)),11))
  subset[1:11,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  rf.importance<-NULL
  rf.var.explained<-NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x[-subset[,i],]
    #remove OTUs that have all zeroes in training set
    x.test <- x[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=y.train,ntree=100,importance=T)
    rf.pred <- c(rf.pred, predict(rf.train,x.test))
    rf.importance<-cbind(rf.importance,rf.train$importance[,1])
    rf.var.explained<-c(rf.var.explained,(rf.train$rsq[100]))
    plot(rf.train)
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  x.sample <- x[y.sample,]
  rf.imp.matrix<-cbind(rf.imp.matrix,rf.importance)
  row.names(rf.imp.matrix)<-row.names(rf.importance)
  rf.var.explained.matrix[,n]<-rf.var.explained
}

##### Asses algorithm performance by comparing real vs predicted values

#Test model performance by comparing predicted vs real data that he has never seen in the test sets, for each independent iteration round (with 10 iteration rounds)
#calculate "average variance explained" from each iteration round (72 predictions x 10 independent iterations = 10 columns)
#cor(rf.pred.matrix[,1], y.test.matrix[,1])^2 # R Squared for each pred vs observed column
#cor(rf.pred.matrix[,2], y.test.matrix[,2])^2 # R Squared for each pred vs observed column, and so on...
#get the diagonal
NM_Lowland_mean_R2_test<-mean(diag(cor(rf.pred.matrix,y.test.matrix)^2))
NM_Lowland_mean_R2_test
#~56.6% of the variability in plant shoot P accumulation could be explained by the model (prediction test)

#visualize the correlation between predictions and actual values, from the 10 iterations
plot(x=rf.pred.matrix,y=y.test.matrix, xlab="Predicted shoot P accumulation", ylab="Actual shoot P accumulation in test set")
#add line if perfect prediction
abline(a=0, b=1, col="red", lwd = 3)

###Create values dataframes
#Each variable (OTU) has 50 values for importance (5 CV loops X 10 independent iterations)
rf.imp.df<-as.data.frame(rf.imp.matrix)
rf.imp.df$mean<-rowMeans(rf.imp.df[,2:51])
rf.imp.df$sd <- apply(rf.imp.df[,2:51], 1, sd, na.rm=TRUE)
rf.imp.df<-rf.imp.df[order(rf.imp.df$mean,decreasing=TRUE),]

###Threshold of OTUs to keep
plot(rf.imp.df$mean)
grid(nx=10,ny=10)
abline(v=100,col="red")
rf.imp.df.100<-rf.imp.df[1:100,]

###Add taxonomical information to OTUs of interest
tax<-tax_table(phy_obj_plant_log0.02plainenonmycno0)
tax<-data.frame(tax)
RF_nonmyc_plaine_taxo<-merge(rf.imp.df.100,tax,by="row.names")
RF_nonmyc_plaine_taxo$name<-paste(RF_nonmyc_plaine_taxo$Row.names,RF_nonmyc_plaine_taxo$Genus,sep="_")


###### Spearman correlation: Paccum ~ RA OTUs##### 

##Make correlations with BH corrections for multiple testing
corspearman<-corr.test(tab0.02,total$Pplant.soil,method="spearman",adjust="BH")
corspearman1<-cbind(corspearman$p.adj,corspearman$r)
corspearman1<-data.frame(corspearman1)
corspearman1<-corspearman1 %>%
  rename(
    p.adj=X1,
    rho=X2
  )
cp2<-data.frame()

###Keep only significant OTUs with positive correlations
cp2<-corspearman1[corspearman1$p.adj<0.05,]
cp3<-cp2[cp2$rho>0,]
head(cp3[1:10,])
view(cp3)
nrow(cp3)
#73 OTUs


###### CAPscale analysis ######

###Capscale function from vegan used to perfrom db-rda. Makes an NMDS on Bray-curtis dissimilarities then constrained analysis using P accumulation
capscaletry<- capscale(tab0.02 ~ Pplant.soil, data=total,dist="bray")
anova(capscaletry) 
summary(capscaletry)

###extract sample scores (capscaletry) for plotting
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))
total$CAP1 <- x$CAP1
total$MDS1 <- x$MDS1

###scores OTUS
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))

###scores Paccum added
scoreotus<-rbind(scoresotus,scores(capscaletry,display = "bp"))

###plot data
ggord(x,scoreotus,grp_in=total$Status,ellipse=TRUE, ellipse_pro=0.95,xlims=c(-3, 3),
      ylim=c(-3, 3))+geom_vline(xintercept = -0.8220, linetype="dotted")+
  geom_vline(xintercept = 0.8155, linetype="dotted")+
  geom_vline(xintercept = 0.1, linetype="dotted")+
  geom_vline(xintercept = -0.1, linetype="dotted")

###Extract OTUs and cap1 scores
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))
#get positive CAP1 scores = correlation with Pplant.soil (positive vector)
arrangedscorsotus<-scoresotus %>% filter(CAP1 >0)
arrangedscorsotus<-arrangedscorsotus %>% arrange(desc(CAP1))
nrow(taxspeciesscores1)
#292 OTUs

# add taxonomy to selected OTUS
taxspeciesscores1<-merge(tax,arrangedscorsotus,by="row.names")

#plot distribution and top 200
taxspeciesscores1<-taxspeciesscores1 %>% arrange(desc(CAP1))
plot(abs(taxspeciesscores1$CAP1))
grid(nx=10,ny=10)
abline(v=200,col="red")

######Merge results######
#assemble tables
tabtoprintrf<-merge(RF_nonmyc_plaine_taxo,taxspeciesscores1,by="Row.names")
cp3$Row.names<-rownames(cp3)
tabtoprintrf<-merge(tabtoprintrf,cp3,by="Row.names")
nrow(tabtoprintrf)
#32 OTUS

#export CSV
tabtoprintrf$name1<-paste(tabtoprintrf$OTU,tabtoprintrf$Species.x,sep="_")  
write.csv2(tabtoprintrf,"Capwilx_RF_Spearman_nonmyc_plainegoodone_JA.csv")

##########################################################################

#####Analysis for myc plants on lowland sites ######

### Sub Select only taxa for which mean RA is > than 0.0001 -> 0.0000434 = log10(0.0001+1)
tax.mean <- taxa_sums(phy_obj_plant_log0.02plainemycno0)/nsamples(phy_obj_plant_log0.02plainemycno0)
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02plainemycno0)

### Make a dataframe
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)
tab0.02<-data.frame(tab0.02)
s<-sample_data(sites.prune)
total<-cbind(tab0.02,s)


#### Random Forest analysis per Biome and Status ####

###For numeric vector Pplant/soil###
view(s)
y <- as.numeric(s$Pplant.soil)
x <- data.matrix(tab0.02)
id <- seq(1:length(y))
x2<-colnames(tab0.02)
numsets <- 10

#actual y representing 10 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)

#predictions are probabilities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)

#create data.frame for importance of each otu
rf.importance<-data.frame()
rf.imp.matrix <- matrix(nrow=length(x2),ncol=1)

#Metrics for model accuracy
rf.var.explained.matrix<-matrix(nrow=5,ncol=numsets)

###Model
par(new=T)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  subset <- matrix(0,nrow=54,ncol=5)
  subset[1:11,1] <- sort(sample(id,11))
  subset[1:11,2] <- sort(sample(which(!(id %in% subset)),11))
  subset[1:11,3] <- sort(sample(which(!(id %in% subset)),11))
  subset[1:11,4] <- sort(sample(which(!(id %in% subset)),11))
  subset[1:10,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  rf.importance<-NULL
  rf.var.explained<-NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x[-subset[,i],]
    #remove OTUs that have all zeroes in training set
    x.test <- x[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=y.train,ntree=100,importance=T)
    rf.pred <- c(rf.pred, predict(rf.train,x.test))
    #predictions are added as a queue with all samples being predicted one in each the 5 CV loops (here 72 predictions per loop). At the end we obtain 10 columns with 1 col for each of the iterations
    rf.importance<-cbind(rf.importance,rf.train$importance[,1])
    rf.var.explained<-c(rf.var.explained,(rf.train$rsq[100]))
    plot(rf.train)
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  x.sample <- x[y.sample,]
  rf.imp.matrix<-cbind(rf.imp.matrix,rf.importance)
  row.names(rf.imp.matrix)<-row.names(rf.importance)
  rf.var.explained.matrix[,n]<-rf.var.explained
}

##### Asses algorithm performance by comparing real vs predicted values

#Test model performance by comparing predicted vs real data that he has never seen in the test sets, for each independent iteration round (with 10 iteration rounds)
#calculate "average variance explained" from each iteration round (72 predictions x 10 independent iterations = 10 columns)
#cor(rf.pred.matrix[,1], y.test.matrix[,1])^2 # R Squared for each pred vs observed column
#cor(rf.pred.matrix[,2], y.test.matrix[,2])^2 # R Squared for each pred vs observed column, and so on...
#get the diagonal
Myc_Lowland_mean_R2_test<-mean(diag(cor(rf.pred.matrix,y.test.matrix)^2))
Myc_Lowland_mean_R2_test
#~57.6% of the variability in plant shoot P accumulation could be explained by the model (prediction test)

#visualize the correlation between predictions and actual values, from the 10 iterations
plot(x=rf.pred.matrix,y=y.test.matrix, xlab="Predicted shoot P accumulation", ylab="Actual shoot P accumulation in test set")
#add line if perfect prediction
abline(a=0, b=1, col="red", lwd = 3)


###Create values dataframes
#Each variable (OTU) has 50 values for importance (5 CV loops X 10 independent iterations)
rf.imp.df<-as.data.frame(rf.imp.matrix)
rf.imp.df$mean<-rowMeans(rf.imp.df[,2:51])
rf.imp.df$sd <- apply(rf.imp.df[,2:51], 1, sd, na.rm=TRUE)
rf.imp.df<-rf.imp.df[order(rf.imp.df$mean,decreasing=TRUE),]

###Threshold of OTUs to keep
plot(rf.imp.df$mean)
grid(nx=10,ny=10)
abline(v=100,col="red")
rf.imp.df.100<-rf.imp.df[1:100,]

###Add taxonomical information to TOUs of interest
tax<-tax_table(phy_obj_plant_log0.02plainemycno0)
tax<-data.frame(tax)
RF_myc_plaine_taxo<-merge(rf.imp.df.100,tax,by="row.names")
RF_myc_plaine_taxo$name<-paste(RF_myc_plaine_taxo$Row.names,RF_myc_plaine_taxo$Genus,sep="_")


###### Spearman correlation: Paccum ~ RA OTUs##### 

##Make correlations with BH corrections for multiple testing
corspearman<-corr.test(tab0.02,total$Pplant.soil,method="spearman",adjust="BH")
corspearman1<-cbind(corspearman$p.adj,corspearman$r)
corspearman1<-data.frame(corspearman1)
corspearman1<-corspearman1 %>%
  rename(
    p.adj=X1,
    rho=X2
  )
cp2<-data.frame()

###Keep only significant OTUs with positive correlations
cp2<-corspearman1[corspearman1$p.adj<0.05,]
cp3<-cp2[cp2$rho>0,]
head(cp3[1:10,])
view(cp3)
nrow(cp3)
#98 OTUs


###### CAPscale analysis ######

###Capscale function from vegan used to perfrom db-rda. Makes an NMDS on Bray-curtis dissimilarities then constrained analysis using P accumulation
capscaletry<- capscale(tab0.02 ~ Pplant.soil, data=total,dist="bray")
anova(capscaletry)
summary(capscaletry)

###extract sample scores (capscaletry) for plotting
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))
total$CAP1 <- x$CAP1
total$MDS1 <- x$MDS1

###scores OTUS
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))

###scores Paccum added
scoreotus<-rbind(scoresotus,scores(capscaletry,display = "bp"))

###plot data
ggord(x,scoreotus,grp_in=total$Status,ellipse=TRUE, ellipse_pro=0.95,xlims=c(-3, 3),
      ylim=c(-3, 3))+geom_vline(xintercept = -0.8220, linetype="dotted")+
  geom_vline(xintercept = 0.8155, linetype="dotted")+
  geom_vline(xintercept = 0.1, linetype="dotted")+
  geom_vline(xintercept = -0.1, linetype="dotted")

###Extract OTUs and cap1 scores
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))
#get positive CAP1 scores = correlation with Pplant.soil (positive vector)
arrangedscorsotus<-scoresotus %>% filter(CAP1 >0)
arrangedscorsotus<-arrangedscorsotus %>% arrange(desc(CAP1))
nrow(taxspeciesscores1)
#315 OTUs

# add taxonomy to selected OTUS
taxspeciesscores1<-merge(tax,arrangedscorsotus,by="row.names")

#plot distribution and top 200
taxspeciesscores1<-taxspeciesscores1 %>% arrange(desc(CAP1))
plot(abs(taxspeciesscores1$CAP1))
grid(nx=10,ny=10)
abline(v=200,col="red")


######Merge results######
#assemble tables
tabtoprintrf<-merge(RF_myc_plaine_taxo,taxspeciesscores1,by="Row.names")
cp3$Row.names<-rownames(cp3)
tabtoprintrf<-merge(tabtoprintrf,cp3,by="Row.names")
nrow(tabtoprintrf)
#41 OTUS

#export CSV
tabtoprintrf$name1<-paste(tabtoprintrf$OTU,tabtoprintrf$Species.x,sep="_")  
write.csv2(tabtoprintrf,"Capwilx_RF_Spearman_myc_plainegoodone_JA.csv")


##Figure 6B: Dendograms and venn diagram######
alpinenmotu<-read.csv2("Capwilx_RF_Spearman_nonmyc_alpinegoodone.csv",header=TRUE)
alpinemycotu<-read.csv2("Capwilx_RF_Spearman_myc_alpinegoodone.csv",header=TRUE)
plainemycotu<-read.csv2("Capwilx_RF_Spearman_myc_plainegoodone.csv",header=TRUE)
plainenmotu<-read.csv2("Capwilx_Spearman_RF_nonmyc_plainegoodone.csv",header=TRUE)


####Merged and curated files

Paccum_OTUs<-read_excel("heatmap_annotation_correlP_tot.xlsx",sheet = "heatmap_annotation_correlP_tot")
Paccum_OTUs<-Paccum_OTUs[-c(1:9),]

Paccum_alpine_NM<-filter(Paccum_OTUs,sites=="alpinenonmyc")
Cap_RF_spearman_alpine_NM<-read.csv2("Capwilx_RF_Spearman_nonmyc_alpinegoodone.csv")
Paccum_alpine_NM_sd<-merge(Paccum_alpine_NM,Cap_RF_spearman_alpine_NM,by="OTU")

Paccum_plaine_NM<-filter(Paccum_OTUs,sites=="plainenonmyc")
Cap_RF_spearman_plaine_NM<-read.csv2("Capwilx_Spearman_RF_nonmyc_plainegoodone.csv")
Paccum_plaine_NM_sd<-merge(Paccum_plaine_NM,Cap_RF_spearman_plaine_NM,by="OTU")

Paccum_nonmyc<-rbind(Paccum_plaine_NM_sd,Paccum_alpine_NM_sd)

Paccum_alpine_Myc<-filter(Paccum_OTUs,sites=="alpinemyc")
Cap_RF_spearman_alpine_Myc<-read.csv2("Capwilx_RF_Spearman_myc_alpinegoodone.csv")
Paccum_alpine_Myc_sd<-merge(Paccum_alpine_Myc,Cap_RF_spearman_alpine_Myc,by="OTU")

Paccum_plaine_Myc<-filter(Paccum_OTUs,sites=="plainemyc")
Cap_RF_spearman_plaine_Myc<-read.csv2("Capwilx_RF_Spearman_myc_plainegoodone.csv")
Paccum_plaine_Myc_sd<-merge(Paccum_plaine_Myc,Cap_RF_spearman_plaine_Myc,by="OTU")

Paccum_myc<-rbind(Paccum_alpine_Myc_sd,Paccum_plaine_Myc_sd)

Paccum<-rbind(Paccum_myc,Paccum_nonmyc)
Paccum$OTU<-as.factor(Paccum$OTU)
Paccum$OTU_sites<-paste(Paccum$OTU,Paccum$sites)

###Importance MSE variable



paletteer_d("ggsci::default_igv")



palette <- c("o__Magnaporthales"="#C74F26",
             "o__Helotiales"="#00ffffff",
             "o__Hypocreales"="#827A8C",
             "o__Glomerales"="#68165cff",
             "o__Glomerellales"="#68D768",
             "o__Geoglossales"="#74C475",
             "o__Coniochaetales"="#5CB0DC",
             "o__Chaetothyriales"="#B86135",
             "o__Cantharellales"="#F0E685",
             "p__Ascomycota_unclassified"="#b5cf6bff",
             "o__Agaricales"="#739957",
             "k__Fungi_unclassified"="#E1D5DC",
             "o__Minutisphaerales"="#ffcb3d",
             "o__Mycosphaerellales"= "#7863A3",
             "o__Myrmecridiales" ="#E3AE68",
             "o__Pezizales"="#CCDEB5",
             "o__Microascales"="#D60047",
             "o__Mortierellales"="#D48F5C",
             "o__Pleosporales"="#AE1E61",
             "o__Tremellales"="#99CC01",
             "o__Sebacinales"="#66007B",
             "o__Sordariales" ="#E6C76E",
             "o__Thelebolales" ="#58635C",
             "o__Tubeufiales"="#A8A8A8",
             "o__Xylariales"="#CC9900",
             "o__Capnodiales"= "#456882ff",
             "o__Holtermanniales"="#914721ff",
             "o__Filobasidiales"="#7f2166ff",
             "c__Sordariomycetes_unclassified"="#cc3a30ff",
             "o__Pezizomycotina_ord_Incertae_sedis"="#612878ff",
             "o__Mytilinidales"="#3a1951ff",
             "c__Dothideomycetes_unclassified"="#4f4fffff",
             "o__Thelephorales"="#a898ccff"
)




Paccum$meanMSE<-as.numeric(Paccum$meanMSE)
Paccum$meanMSE2<-Paccum$meanMSE
Paccum$meanMSE<-log10(Paccum$meanMSE)
Paccum$meanMSE[Paccum$sites %in% c("alpinemyc", "plainemyc")] <- -Paccum$meanMSE[Paccum$sites %in% c("alpinemyc", "plainemyc")]


P1<-Paccum%>%
  mutate(name = fct_reorder(OTU, desc(...1))) %>%
  ggplot(aes(x=name, y=meanMSE, fill=Order.x,colour=Order.x))+
  geom_segment(aes(x = name, xend = name, y = 0, yend = meanMSE,
                   color =Order.x),lwd=1)+
  scale_color_manual(values=palette)+
  #c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF","#5DB1DDFF","#802268FF","#6BD76BFF","#D595A7FF","#924822FF","#837B8DFF","#C75127FF","#D58F5CFF","#7A65A5FF","#E4AF69FF","#3B1B53FF","#CDDEB7FF","#612A79FF","#AE1F63FF","#E7C76FFF","#5A655EFF","#CC9900FF","#99CC00FF","#A9A9A9FF","#CC9900FF","#99CC00FF","#33CC00FF","#00CC33FF","#00CC99FF","#0099CCFF","#0A47FFFF","#4775FFFF","#FFC20AFF","#FFD147FF","#990033FF"))+
  geom_point(aes(),size=2,stroke=2)+ scale_shape_manual(values=c(16))+
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
  # position=position_dodge(0.05))+facet_wrap(~sites,scales="free")+
  theme_classic()+
  theme(axis.text.x = element_text(NULL),
        axis.title.x = element_blank()) + 
  theme(axis.text.x.bottom = element_text(angle=90,size=14,hjust = 0.95))+ 
  #theme(legend.position="none")+
  theme(legend.text = element_text(size = 8),legend.position = "left")+
  theme(axis.title = element_text(size=24))+ 
  theme(axis.text.y = element_text(size=8))
P1<-P1+coord_flip()
P1

####################################

##Fig 6

###Make dendogram for AM plants

otu_arbre<-Paccum_myc[c(1,8:14,16)]
rownames(otu_arbre)<-otu_arbre$namesites
for(i in 1:ncol(otu_arbre)){
  otu_arbre[,i]<-as.factor(otu_arbre[,i])}

frm <- ~Kingdom.x/Phylum.x/Class.x/Order.x/Family.x/Genus.x/Species.x/namesites
tr <- as.phylo(frm, data = otu_arbre, collapse=FALSE)
tr$edge.length <- rep(1, nrow(tr$edge))
plot(tr, show.node.label=TRUE)
Nnode(tr)
write.tree(phy=tr, file="dendo_tree_Pcorrel_myc_plants130924_PDcalc_species_otusname.newick")

##Make dendogram for non-AM plants
otu_arbre<-Paccum_nonmyc[c(1,8:14,16)]
rownames(otu_arbre)<-otu_arbre$namesites
for(i in 1:ncol(otu_arbre)){
  otu_arbre[,i]<-as.factor(otu_arbre[,i])}

frm <- ~Kingdom.x/Phylum.x/Class.x/Order.x/Family.x/Genus.x/Species.x/namesites
tr <- as.phylo(frm, data = otu_arbre, collapse=FALSE)
tr$edge.length <- rep(1, nrow(tr$edge))
plot(tr, show.node.label=TRUE)
Nnode(tr)
write.tree(phy=tr, file="dendo_tree_Pcorrel_nonmyc_plants130924_species_otusname.newick")

## these files were then used to produce final figures in itol


#####Make venn diagram

heat<-read.csv2("heatmap_annotation_correlP_tot.csv",stringsAsFactors = T)

heat1<-heat[10:209,]
for(i in 16:ncol(heat1)){
  heat1[,i]<-as.numeric(heat1[,i])}
heat1$meanabund<-rowMeans(heat1[,c(16:268)])
heat1<-separate(heat1,sites,into=c("Biome","Status"),sep="ne")
heatforlist<-heat1[,c(1,15)]

# Get unique occurrences from the second column
unique_occurrences <- unique(heatforlist$Status)

# Create an empty list to store results
result_list <- list()

# Loop through unique occurrences and extract corresponding data
for (occurrence in unique_occurrences) {
  subset_data <- heatforlist[heatforlist$Status == occurrence, ]
  result_list[[occurrence]] <- subset_data
}
x<-result_list

x1<-list(myc=x$myc[,1],
         nonmyc=x$nonmyc[,1])


# Default plot
names(x1)<-c("AM plants","Non-AM plants")
ggvenn(x1,fill_color = c("#f8766bff", "#00bdc2ff"),
       stroke_size = 0.5, set_name_size = 5)
