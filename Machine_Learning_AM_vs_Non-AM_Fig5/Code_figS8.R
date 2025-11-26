###Code for Figure S8
###Heatmap for OTUs enriched in AM or non-AM (non-Myc) plants

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
library(randomForest)
library(stats)

#set wd
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Machine_Learning_Fig5_Fig6")

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
##relative abundance for each OTU
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
phy_obj_plant_ra<-otu_table(phy_obj_plant_ra)+1
phy_obj_plant_log<-log10(otu_table(phy_obj_plant_ra))
###reassemble otu_table with taxonomy and metadata
phy_obj_plant_log<- phyloseq(otu_table(phy_obj_plant_log), TAX, samples)



############ Analyses: RF + Wilcoxon's test + Capscale to identify best predictors of AM vs non-AM differences #####

#####Make files#####
##### 2% occurence file: keep only OTUs present in more than 2% of the samples
phy_obj_plant_log0.02<-phyloseq_filter_prevalence(
  phy_obj_plant_log,
  prev.trh = 0.02,
  abund.trh = NULL) #4748 OTUs

#get the mean log10(RA+1) for each OTU
tax.mean <- taxa_sums(phy_obj_plant_log0.02)/nsamples(phy_obj_plant_log0.02)

#prune to keep OTUs with mean relative abundance (RA) > 0.0001 <=> log10(RA+1) > 0.0000434
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02) #1132 OTUs
class(sites.prune)
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)

#####Make Alpine and Lowland datasets
#remove "0" only taxa
phy_obj_plant_log0.02alpine<-subset_samples(phy_obj_plant_log0.02,Type_site=="Alpine")
phy_obj_plant_log0.02alpineno0<- prune_species(speciesSums(phy_obj_plant_log0.02alpine) > 0, phy_obj_plant_log0.02alpine)

phy_obj_plant_log0.02plaine<-subset_samples(phy_obj_plant_log0.02,Type_site=="Plaine")
phy_obj_plant_log0.02plaineno0<- prune_species(speciesSums(phy_obj_plant_log0.02plaine) > 0, phy_obj_plant_log0.02plaine)

#######################################################

#####Alpine sites######

### Sub Select again only taxa for which mean relative abundance (RA) > 0.0001 <=> log10(RA+1) > 0.0000434
tax.mean <- taxa_sums(phy_obj_plant_log0.02alpineno0)/nsamples(phy_obj_plant_log0.02alpineno0)
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02alpineno0) ####822 OTUs

### Make a dataframe
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)
tab0.02<-data.frame(tab0.02)
s<-sample_data(sites.prune)
total<-cbind(tab0.02,s)


#### Random Forest analysis between AM vs non-AM plants ####

####vector of discrete trait coded 1=Myc and 0=Non-myc
view(s)
y <- as.numeric(s$Status=="Mycorhizal")
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
rf.errorrate.matrix<-matrix(nrow=5,ncol=numsets)
par(new=T)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  subset <- matrix(0,nrow=144,ncol=5)
  subset[1:29,1] <- sort(sample(id,29))
  subset[1:29,2] <- sort(sample(which(!(id %in% subset)),29))
  subset[1:29,3] <- sort(sample(which(!(id %in% subset)),29))
  subset[1:29,4] <- sort(sample(which(!(id %in% subset)),29))
  subset[1:28,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  rf.importance<-NULL
  rf.errorrate<-NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x[-subset[,i],]
    x.test <- x[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=as.factor(y.train),importance=TRUE,ntree=300)
    rf.pred <- c(rf.pred, predict(rf.train,x.test,type="prob")[,2])
    #Get feature importance = Mean Decrease Accuracy
    rf.importance<-cbind(rf.importance,rf.train$importance[,3])
    rf.errorrate<-c(rf.errorrate,mean(rf.train$confusion[,3]))
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
  rf.errorrate.matrix[,n]<-rf.errorrate
}

save(y.sample.matrix,y.test.matrix,rf.pred.matrix,file="test_predict_alpine.Rdata")

##### Asses classifier performance with ROC curve
library(pROC)
load("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Machine_Learning_Fig5_Fig6/test_predict_alpine.Rdata")
#load("C:/Users/pbruyant/Documents/git_Microbiota/codes-article-principal/Microbiome analyses/test_predict.Rdata")
#calculate average predicted y by sample, 
#then generate ROC between predicted and actual y
for (n in 1:numsets) {
  y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
  y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
  y.test.matrix[,n] <- y.test.temp[,2]
  rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
  rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
  rf.pred.matrix[,n] <- rf.pred.temp[,2]
}
identical(y,y.test.matrix[,1])
rf.pred.avg <- rowMeans(rf.pred.matrix)
rf <- roc(y,rf.pred.avg)
rf$auc #Area under the curve: 0.9414
plot(rf$specificities,rf$sensitivities,type="n",xlim=c(1,0),xlab="Specificity",ylab="Sensitivity",main="Myc vs Non-myc all")
lines(rf$specificities,rf$sensitivities,col="black")
dev.off()

####error
rf.error<-as.data.frame(rf.errorrate.matrix)
write.csv2(rf.error,"rf.error_class_alpine.csv")


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
tax<-tax_table(phy_obj_plant_log0.02alpineno0)
tax<-data.frame(tax)
RF_alpine_taxo<-merge(rf.imp.df.100,tax,by="row.names")
RF_alpine_taxo$name<-paste(RF_alpine_taxo$Row.names,RF_alpine_taxo$Genus,sep="_")


###### CAPscale analysis ######

######Capscale function from vegan used to perform db-rda. Makes an NMDS on Bray-curtis dissimilarities then constrained analysis by "Status"
capscaletry<- capscale(tab0.02~ as.factor(Status), data=total,dist="bray")
anova(capscaletry)
summary(capscaletry)

###extract sample scores (capscaletry) for plotting
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$MDS1

#plot separation
p3<-ggplot(total, aes(x= CAP1, y=MDS1,color = Status)) + 
  stat_ellipse(aes(fill=factor(total$Status)), geom = "polygon",alpha = 0,show.legend=FALSE,lwd=1) +
  geom_point(size=2) + theme_classic()+
  scale_color_manual(values=c("#f8756bff","#00bdc2ff"))+
  scale_fill_manual(values=c("#f8756bff","#00bdc2ff"))+
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  theme(plot.title = element_text(hjust = 0.5),legend.text = element_text(size=18,family="arial"),
        axis.text = element_text(size=18,family="arial"),legend.title=element_text(size=18,family="arial"),axis.title=element_text(size=18,family="arial"))+ 
  labs(color="Status",x="CAP1 (1.6%)",y="NMDS1 (6%)")+theme(legend.position="bottom")+
  coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 4.5))
p3


###Extract CAP1 scores for all OTUs
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))
arrangedscorsotus<-scoresotus %>% arrange(desc(abs(CAP1)))
taxo<-tax_table(phy_obj_plant_log)
taxo<-data.frame(taxo)
taxspeciesscores1<-merge(taxo,arrangedscorsotus,by="row.names")
taxspeciesscores1<-taxspeciesscores1 %>% arrange(desc(abs(CAP1)))
#822 OTUs

#add taxonomy info
rownames(taxspeciesscores1)<-taxspeciesscores1$Row.names
taxspeciesscores1<-t(taxspeciesscores1)
taxspeciesscores1<-data.frame(taxspeciesscores1)
taxspeciesscores2<-total %>% dplyr::select(which(colnames(total) %in% names(taxspeciesscores1)))


###### Wilcoxon's tests Myc vs non-Myc on all OTUs #####
library(stats)

#Add plant Status info
taxspeciesscores2$Status<-total$Status
taxspeciesscores2<-taxspeciesscores2[order(taxspeciesscores2$Status),]
taxspeciesscores3<-taxspeciesscores2[,1:ncol(taxspeciesscores2)] #822 OTUs = 822 columns

#Run Wilcoxon test comparing OTU RA between Myc vs Non-Myc plants, get only OTUs showing significant differences
estimates_P.value<-data.frame()
for (i in 1:822){
  mod <- wilcox.test(taxspeciesscores2[1:72,i],taxspeciesscores2[73:144,i])
  estimates_P.value[i,1] <- mod$p.value
  rownames(estimates_P.value)[i]<-colnames(taxspeciesscores3)[i]
}
estimates_P.value$Row.names<-rownames(estimates_P.value)
estimatessigni<-estimates_P.value[estimates_P.value$V1 <= 0.05,]
estimatessigni$Row.names<-rownames(estimatessigni)
#105 OTUS

######Merge results####
#OTUs kept if in top 100 in RF AND Significantly differential abundant based of Wilcoxon's test => 57 OTUs
tabtoprint<-merge(estimatessigni,t(taxspeciesscores1),by="Row.names")
tabtoprint<-tabtoprint%>% rename("p-value"=V1)
tabtoprintrf_alp<-merge(RF_alpine_taxo,tabtoprint,by="Row.names")
write.csv2(tabtoprintrf_alp,"Capwilx_RF_myc_vs_nonmyc_alpine_57otus.csv")


#####################################################################

#####Lowland sites######

### Sub Select only taxa for which mean RA is > than 0.0001 -> 0.0000434 = log10(0.0001+1)
tax.mean <- taxa_sums(phy_obj_plant_log0.02plaineno0)/nsamples(phy_obj_plant_log0.02plaineno0)
sites.prune <- prune_taxa(tax.mean > 0.0000434, phy_obj_plant_log0.02plaineno0)#850 OTUs

### Make a dataframe
tab0.02<-otu_table(sites.prune)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)
tab0.02<-data.frame(tab0.02)
s<-sample_data(sites.prune)
total<-cbind(tab0.02,s)


#### Random Forest analysis between AM vs non-AM plants ####

####vector of discrete trait coded 1=Myc and 0=Non-myc
view(s)
y <- as.numeric(s$Status=="Mycorhizal")
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
rf.errorrate.matrix<-matrix(nrow=5,ncol=numsets)
par(new=T)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  subset <- matrix(0,nrow=109,ncol=5)
  subset[1:21,1] <- sort(sample(id,21))
  subset[1:22,2] <- sort(sample(which(!(id %in% subset)),22))
  subset[1:22,3] <- sort(sample(which(!(id %in% subset)),22))
  subset[1:22,4] <- sort(sample(which(!(id %in% subset)),22))
  subset[1:22,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  rf.importance<-NULL
  rf.errorrate<-NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x[-subset[,i],]
    x.test <- x[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=as.factor(y.train),importance=TRUE,ntree=300)
    rf.pred <- c(rf.pred, predict(rf.train,x.test,type="prob")[,2])
    #Get feature importance = Mean Decrease Accuracy
    rf.importance<-cbind(rf.importance,rf.train$importance[,3])
    rf.errorrate<-c(rf.errorrate,mean(rf.train$confusion[,3]))
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
  rf.errorrate.matrix[,n]<-rf.errorrate
}

save(y.sample.matrix,y.test.matrix,rf.pred.matrix,file="test_predict_lowland.Rdata")


##### Asses classifier performance with ROC curve

library(pROC)
load("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Machine_Learning_Fig5_Fig6/test_predict_lowland.Rdata")
#load("C:/Users/pbruyant/Documents/git_Microbiota/codes-article-principal/Microbiome analyses/test_predict.Rdata")
#calculate average predicted y by sample, 
#then generate ROC between predicted and actual y
for (n in 1:numsets) {
  y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
  y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
  y.test.matrix[,n] <- y.test.temp[,2]
  rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
  rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
  rf.pred.matrix[,n] <- rf.pred.temp[,2]
}
identical(y,y.test.matrix[,1])
rf.pred.avg <- rowMeans(rf.pred.matrix)
rf <- roc(y,rf.pred.avg)
rf$auc #Area under the curve: 0.865
plot(rf$specificities,rf$sensitivities,type="n",xlim=c(1,0),xlab="Specificity",ylab="Sensitivity",main="Myc vs Non-myc all")
lines(rf$specificities,rf$sensitivities,col="black")
dev.off()


####error
rf.error<-as.data.frame(rf.errorrate.matrix)
write.csv2(rf.error,"rf.error_class_lowland.csv")


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
tax<-tax_table(phy_obj_plant_log0.02plaineno0)
tax<-data.frame(tax)
RF_plaine_taxo<-merge(rf.imp.df.100,tax,by="row.names")
RF_plaine_taxo$name<-paste(RF_plaine_taxo$Row.names,RF_plaine_taxo$Genus,sep="_")



###### CAPscale analysis ######

######Capscale function from vegan used to perform db-rda. Makes an NMDS on Bray-curtis dissimilarities then constrained analysis by "Status"
capscaletry<- capscale(tab0.02~ as.factor(Status), data=total,dist="bray")
anova(capscaletry)
summary(capscaletry)

###extract sample scores (capscaletry) for plotting
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$MDS1

#plot separation
p3<-ggplot(total, aes(x= CAP1, y=MDS1,color = Status)) + 
  stat_ellipse(aes(fill=factor(total$Status)), geom = "polygon",alpha = 0,show.legend=FALSE,lwd=1) +
  geom_point(size=2) + theme_classic()+
  scale_color_manual(values=c("#f8756bff","#00bdc2ff"))+
  scale_fill_manual(values=c("#f8756bff","#00bdc2ff"))+
  #scale_shape_manual(values=c(1,16))+
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  theme(plot.title = element_text(hjust = 0.5),legend.text = element_text(size=18,family="arial"),
        axis.text = element_text(size=18,family="arial"),legend.title=element_text(size=18,family="arial"),axis.title=element_text(size=18,family="arial"))+ 
  labs(color="Status",x="CAP1 (1.6%)",y="NMDS1 (10.4%)")+theme(legend.position="bottom")+
  coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 4.5))
p3


###Extract CAP1 scores for all OTUs
scoresotus <- as.data.frame(scores(capscaletry, display = "species"))
arrangedscorsotus<-scoresotus %>% arrange(desc(abs(CAP1)))
taxo<-tax_table(phy_obj_plant_log)
taxo<-data.frame(taxo)
taxspeciesscores1<-merge(taxo,arrangedscorsotus,by="row.names")
taxspeciesscores1<-taxspeciesscores1 %>% arrange(desc(abs(CAP1)))
#850 OTUs


#add taxonomy info
rownames(taxspeciesscores1)<-taxspeciesscores1$Row.names
taxspeciesscores1<-t(taxspeciesscores1)
taxspeciesscores1<-data.frame(taxspeciesscores1)
taxspeciesscores2<-total %>% dplyr::select(which(colnames(total) %in% names(taxspeciesscores1)))


###### Wilcoxon tests #####
library(stats)

#Add plant Status info
taxspeciesscores2$Status<-total$Status
taxspeciesscores2<-taxspeciesscores2[order(taxspeciesscores2$Status),]
taxspeciesscores3<-taxspeciesscores2[,1:ncol(taxspeciesscores2)] #850 OTUs

estimates_P.value<-data.frame()
for (i in 1:850){
  mod <- wilcox.test(taxspeciesscores2[1:54,i],taxspeciesscores2[55:109,i])
  estimates_P.value[i,1] <- mod$p.value
  rownames(estimates_P.value)[i]<-colnames(taxspeciesscores3)[i]
}
estimates_P.value$Row.names<-rownames(estimates_P.value)
estimatessigni<-estimates_P.value[estimates_P.value$V1 <= 0.05,]
estimatessigni$Row.names<-rownames(estimatessigni)
#82 OTUs

######Merge results####
#OTUs kept if in top 100 in RF AND Significantly differential abundant based of Wilcoxon's test => 47 OTUs
tabtoprint<-merge(estimatessigni,t(taxspeciesscores1),by="Row.names")
tabtoprint<-tabtoprint%>% rename("p-value"=V1)
tabtoprintrf_plaine<-merge(RF_plaine_taxo,tabtoprint,by="Row.names")
write.csv2(tabtoprintrf_plaine,"Capwilx_RF_myc_vs_nonmyc_plaine_47otus.csv")


#####################################################################################################

####Make dataframe for visualizing as heatmap in Morpheus
###Calculate mean log10+1(RA) per family x sites for each OTU found within the analyses (57 Otus Alpine + 47 Otus Lowland)

library(dplyr)

##Indicate in which dataset OTUs are enriched###
tabtoprintrf_plaine$CAP1<-as.numeric(tabtoprintrf_plaine$CAP1)
tabtoprintrf_plaine <- tabtoprintrf_plaine %>% 
  dplyr::mutate(Biome_Status = dplyr::case_when(.$CAP1 < 0 ~ "Lowland_Myc",
                        .$CAP1 > 0 ~ "Lowland_Non_Myc"))

tabtoprintrf_alp$CAP1<-as.numeric(tabtoprintrf_alp$CAP1)
tabtoprintrf_alp <- tabtoprintrf_alp %>% 
  mutate(Biome_Status = case_when(CAP1 < 0 ~ "Alpine_Myc",
                                  CAP1 > 0 ~ "Alpine_Non_Myc"))

###merge files for lowland and alpine biome
tab_to_print_all<-rbind(tabtoprintrf_plaine,tabtoprintrf_alp)

###Any non_unique OTUs ? 
nrow(tab_to_print_all) #104 rows
length(unique(tab_to_print_all$Row.names)) #102 rows
#2 OTUs identified in both biomes
tab_to_print_all[duplicated(tab_to_print_all$Row.names), ]
##Otu000135 in both biomes
##Otu000761 in both biomes

###Get RA values  
tab0.02<-otu_table(phy_obj_plant_log0.02)
tab0.02<-data.frame(tab0.02)
tab0.02<-t(tab0.02)
tab0.02<-as.data.frame(tab0.02)

results_with_RA<-tab0.02 %>% dplyr::select(which(colnames(tab0.02) %in% tab_to_print_all$Row.names))
results_with_RA<-t(results_with_RA)

###Calculate OTU averages per treatment: Family x Biome 
ast_columns <- grep("Per1Ast|Com1Ast|Doua1Ast", colnames(results_with_RA), value = TRUE)
poa_columns <- grep("Per1Poa|Com1Poa|Doua1Poa", colnames(results_with_RA), value = TRUE)
ger_columns <- grep("Per1Ger|Com1Ger|Doua1Ger", colnames(results_with_RA), value = TRUE)
bra_columns <- grep("Per1Bra|Com1Bra|Doua1Bra", colnames(results_with_RA), value = TRUE)
car_columns <- grep("Per1Car|Com1Car|Doua1Car", colnames(results_with_RA), value = TRUE)
cyp_columns <- grep("Per1Cyp|Com1Cyp|Doua1Cyp", colnames(results_with_RA), value = TRUE)

ast_col_alp <- grep("Gal1Ast|Lau1Ast|Cha1Ast|Cla1Ast", colnames(results_with_RA), value = TRUE)
poa_col_alp <- grep("Gal1Poa|Lau1Poa|Cha1Poa|Cla1Poa", colnames(results_with_RA), value = TRUE)
ger_col_alp <- grep("Gal1Ger|Lau1Ger|Cha1Ger|Cla1Ger", colnames(results_with_RA), value = TRUE)
ren_col_alp <- grep("Gal1Ren|Lau1Ren|Cha1Ren|Cla1Ren", colnames(results_with_RA), value = TRUE)
bra_col_alp <- grep("Gal1Bra|Lau1Bra|Cha1Bra|Cla1Bra", colnames(results_with_RA), value = TRUE)
car_col_alp <- grep("Gal1Car|Lau1Car|Cha1Car|Cla1Car", colnames(results_with_RA), value = TRUE)
cyp_col_alp <- grep("Gal1Cyp|Lau1Cyp|Cha1Cyp|Cla1Cyp", colnames(results_with_RA), value = TRUE)

# Calculate mean for each row 
ast_lowland_means <- rowMeans(results_with_RA[, ast_columns, drop = FALSE])
poa_lowland_means <- rowMeans(results_with_RA[, poa_columns, drop = FALSE])
ger_lowland_means <- rowMeans(results_with_RA[, ger_columns, drop = FALSE])
bra_lowland_means <- rowMeans(results_with_RA[, bra_columns, drop = FALSE])
car_lowland_means <- rowMeans(results_with_RA[, car_columns, drop = FALSE])
cyp_lowland_means <- rowMeans(results_with_RA[, cyp_columns, drop = FALSE])

ast_alp_means <- rowMeans(results_with_RA[, ast_col_alp, drop = FALSE])
poa_alp_means <- rowMeans(results_with_RA[, poa_col_alp, drop = FALSE])
ger_alp_means <- rowMeans(results_with_RA[, ger_col_alp, drop = FALSE])
ren_alp_means <- rowMeans(results_with_RA[, ren_col_alp, drop = FALSE])
bra_alp_means <- rowMeans(results_with_RA[, bra_col_alp, drop = FALSE])
car_alp_means <- rowMeans(results_with_RA[, car_col_alp, drop = FALSE])
cyp_alp_means <- rowMeans(results_with_RA[, cyp_col_alp, drop = FALSE])

# Combine the means into a data frame
means_df <- data.frame(ast_lowland_means, poa_lowland_means, ger_lowland_means, bra_lowland_means, car_lowland_means, cyp_lowland_means,
                       ast_alp_means, poa_alp_means,ger_alp_means,ren_alp_means,bra_alp_means, car_alp_means,cyp_alp_means)

###Combine RA with taxonomy and RF results 
means_df<-merge(tab_to_print_all,means_df,by.x="Row.names",by.y="row.names",all.x = TRUE)
###Keep column of interest
means_df<-means_df[,c(1,53,54,63:71,73:86)]

##Change the enrichment of the two OTUs found in both biomes
means_df$Biome_Status[means_df$Row.names =="Otu000135"] <- "Lowland and Alpine Myc"
means_df$Biome_Status[means_df$Row.names =="Otu000761"] <- "Lowland and Alpine Myc"

#write table
write.csv2(means_df,"means_abund_by_fam_biomes_OTUs_myc_vs_non_myc.csv")

###remove both rows duplicated for heatmap visualization
means_df<-means_df[-c(20,66),]

###Calculate overall mean RA for each OTU
means_df$mean_RA<-rowMeans(means_df[,14:26])

##reorder so all are not on the left-side of the df
means_df<-means_df[,c(27,1:26)]

write.csv2(means_df,"means_abund_by_fam_biomes_OTUs_myc_vs_non_myc_for_heatmap.csv")
#Use this file for visualization in Morpheus