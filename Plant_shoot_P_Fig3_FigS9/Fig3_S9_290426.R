### Code for figures 3 and S9 on plant shoot P

###########
#Figure 3
### P shoots
### Mind that in the code the different biomes are called "type_sites" and Lowland is called "Plaine"
### libraries needed####
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

###phyloseq files (only plant data)####
#setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/1_Plant_shoot_P_Fig3")
setwd("C:/Users/pbruyant/Documents/Article thèse/Nvelle version/New_code_fig3_S9/")

###import data###
otu_mat<- read_excel("tableitsx_phyloseq.xlsx", sheet = "Table_OTU")
tax_mat<- read_excel("tableitsx_phyloseq.xlsx", sheet = "Taxonomy")
samples_df <- read_excel("tableitsx_phyloseq.xlsx", sheet = "Metadata")
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

#####Shoot P ####
pplants<-samples_df
pplants$Family<-as.factor(pplants$Family)
pplants$Site<-as.factor(pplants$Site)

###order data
pplants<-pplants%>% rstatix::reorder_levels("Family",order=c("Asteraceae","Geraniaceae","Ranunculaceae","Poaceae","Caryophyllaceae","Brassicaceae","Cyperaceae"))
pplants<-pplants%>% rstatix::reorder_levels("Site",order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","La Doua","Commelle"))
pplants$abbrev<-pplants$Family
pplants$abbrev<-dplyr::recode_factor(pplants$abbrev,Asteraceae="Ast",Ranunculaceae="Ran",
                                     Geraniaceae="Ger",Poaceae="Poa", Caryophyllaceae= "Car",
                                     Brassicaceae="Bra",Cyperaceae="Cyp")


p1<-ggboxplot(pplants,x="abbrev",y="Pplant",add="jitter",color="Status",
              fill="Status",size=0.1, width=0.8,ylab = "Shoot [P] (mg/g)",alpha=0.1)+
  theme_classic()+theme(axis.line.x = element_line())+ 
  theme(axis.text.y= element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text.x= element_text(size=16,angle=90),axis.title.x = element_blank())+ 
  theme(strip.text.x = element_text(size=16),strip.background = element_blank())+ 
  facet_wrap(~Site, scales="free",ncol=7)+ ylim(0,8)+
  theme(legend.position="bottom",legend.text = element_text(size=16),
        legend.title=element_text(size=18))
p1



###Statistical analysis####

#Test normality of shoot P data
shapiro.test(pplants$Pplant)
#not normal

#test Homoscedasticity between groups
car::leveneTest(Pplant ~Site/Status/Family, data=pplants)
#different variances

#Normalize data for ANOVA using box-cox transformation
library(bestNormalize)
obj_bestnormP<-bestNormalize(pplants$Pplant)
obj_bestnormP
pplants$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data

#verify normality of transformed data
shapiro.test(pplants$bestnormP)
#OK

#verify Homoscedasticity of transformed data
car::leveneTest(bestnormP ~Site/Status/Family, data=pplants)
#OK


###### ANOVA with nested model shoot P ~ Site/Status/Family ####
res.aov<-aov(bestnormP ~ Site/Status/Family, data=pplants)
summary(res.aov)

# Df Sum Sq Mean Sq F value Pr(>F)    
# Site                 6  76.38  12.730   80.52 <2e-16 ***
#   Site:Status          7  17.97   2.567   16.23 <2e-16 ***
#   Site:Status:Family  28 124.29   4.439   28.07 <2e-16 ***
#   Residuals          211  33.36   0.158                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Get R2 values

df<-summary(res.aov)
sum(df[[1]][,2])
sqR<-df[[1]][,2]/sum(df[[1]][,2])
sqR
# [1] 0.30309739 0.07130223 0.49321577 0.13238462





###within sites  

##### ANOVA between Status for each Site ####
library(stats)
anova_results <- list()
leveneTest_results <-list()
shapiro_results<-list()
anova_results_bestnorm <- list()
leveneTest_results_bestnorm  <-list()
shapiro_results_bestnorm <-list()
results_wilx<-list()
# Get unique sites in the dataframe
sites <- unique(pplants$Site)

# Iterate over each site
for (Site in sites) {
  # Filter the dataframe for the current site
  site_df <- pplants[pplants$Site == Site, ]
  # Perform ANOVA
  res.aov <- aov(Pplant ~ factor(Status), data =site_df)
  # Store ANOVA result
  anova_results[[as.character(Site)]] <- summary(res.aov)
  #perform levene_test
  res.levene<-leveneTest(Pplant ~ factor(Status), data =site_df)
  leveneTest_results[[Site]] <-res.levene
  #shapiro test
  res.shap<-shapiro.test(site_df$Pplant)
  shapiro_results[[Site]] <- res.shap
  #bestnorm
  obj_bestnorm<-boxcox(site_df$Pplant)
  site_df$bestnorm<-obj_bestnorm$x.t
  res.aov <- aov(bestnorm ~ factor(Status), data =site_df)
  anova_results_bestnorm[[Site]] <- summary(res.aov)
  #perform levene_test
  res.levene<-leveneTest(bestnorm ~ factor(Status), data =site_df)
  leveneTest_results_bestnorm[[Site]] <- res.levene
  #shapiro test
  res.shap<-shapiro.test(site_df$bestnorm)
  shapiro_results_bestnorm[[Site]] <-res.shap
  ##for_sites_when_condition not respected
  res.wilx<-rstatix::wilcox_test(site_df,Pplant~ Status,p.adjust.method = "none")
  results_wilx[[Site]]<-res.wilx
}

# Print ANOVA results for each site
anova_results
leveneTest_results
shapiro_results
anova_results_bestnorm
leveneTest_results_bestnorm
shapiro_results_bestnorm
results_wilx
####Only significant differences in Shoot P for Galibier et Lautaret

###See mean for each site
mean<- pplants %>% group_by(Site,Status) %>% summarize(., mean_ShootP=mean(`Pplant`))
mean

mean_myc<-filter(mean,Status=="Mycorhizal")
mean_nonmyc<-filter(mean,Status=="Non-mycorhizal")

# Calculate the percentage increase or decrease of total P of non-myc plants compared to myc plants by sites
mean_myc$percentage_change <- (mean_nonmyc$mean_ShootP - mean_myc$mean_ShootP) / mean_myc$mean_ShootP* 100
mean_myc$mean_nm<-mean_nonmyc$mean_ShootP
mean_myc

##### ANOVA between Families for each site#####

library(stats)
library(car)
library(rstatix)

anova_results <- list()
leveneTest_results <- list()
shapiro_results <- list()
anova_results_bestnorm <- list()
leveneTest_results_bestnorm <- list()
shapiro_results_bestnorm <- list()
tukey_results <- list()
tukey_results_bn <- list()
pairwise_wilcox_results <- list()

sites <- unique(pplants$Site)

for (Site in sites) {
  
  site_df <- pplants[pplants$Site == Site, ]
  #Anova
  res.aov <- aov(Pplant~ factor(Family), data = site_df)
  anova_results[[Site]] <- summary(res.aov)
  #Variance check
  res.levene <- car::leveneTest(Pplant ~ factor(Family), data = site_df)
  leveneTest_results[[Site]] <- res.levene
  #Shapiro on residues
  res.shap <- shapiro.test(residuals(res.aov))  
  shapiro_results[[Site]] <- res.shap
  # Transform using boxcox
  obj_bestnorm <- bestNormalize::boxcox(site_df$Pplant)
  site_df$bestnorm<-obj_bestnorm$x.t
  # Anova on Transformed data
  res.aov_bn <- aov(bestnorm ~ factor(Family), data = site_df)
  anova_results_bestnorm[[Site]] <- summary(res.aov_bn)
  #Verify assumptions
  res.levene_bn <- car::leveneTest(bestnorm ~ factor(Family), data = site_df)
  leveneTest_results_bestnorm[[Site]] <- res.levene_bn
  #Veirfy normality
  res.shap_bn <- shapiro.test(residuals(res.aov_bn))
  shapiro_results_bestnorm[[Site]] <- res.shap_bn
  #Pairwise comparisons
  if(res.shap_bn$p.value > 0.04 & res.levene_bn$`Pr(>F)`[1] > 0.05){
    
    # Tukey on BN
    tukey_bn <- TukeyHSD(res.aov_bn)
    tukey_results_bn[[Site]] <- tukey_bn
    
  } else {
    
    # airwise Wilcoxon
    pw <- pairwise.wilcox.test(site_df$Pplant,
                               site_df$Family,
                               p.adjust.method = "fdr")
    
    pairwise_wilcox_results[[Site]] <- pw
  }
}

anova_results
leveneTest_results
shapiro_results

anova_results_bestnorm
leveneTest_results_bestnorm
shapiro_results_bestnorm

tukey_results
tukey_results_bn
pairwise_wilcox_results

library(multcompView)

tukey_to_letters <- function(tukey_obj) {
  
  # extract Tukey p-values
  tuk <- tukey_obj$`factor(Family)`
  
  # get p-values
  pvals <- tuk[, "p adj"]
  
  # names must match comparisons
  names(pvals) <- rownames(tuk)
  
  # convert to compact letters
  letters <- multcompLetters(pvals)$Letters
  
  return(letters)
}

letters_results <- lapply(tukey_results_bn, tukey_to_letters)
letters_results


#####################################################
#Fig. S9
#Shoot P accumulation

#####Shoot P ####
pplants<-samples_df
pplants$Family<-as.factor(pplants$Family)
pplants$Site<-as.factor(pplants$Site)
pplants$`Pplant/soil`<-as.factor(pplants$Site)

###order data
pplants<-pplants%>% rstatix::reorder_levels("Family",order=c("Asteraceae","Geraniaceae","Ranunculaceae","Poaceae","Caryophyllaceae","Brassicaceae","Cyperaceae"))
pplants<-pplants%>% rstatix::reorder_levels("Site",order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","La Doua","Commelle"))
pplants$abbrev<-pplants$Family
pplants$abbrev<-dplyr::recode_factor(pplants$abbrev,Asteraceae="Ast",Ranunculaceae="Ran",
                                     Geraniaceae="Ger",Poaceae="Poa", Caryophyllaceae= "Car",
                                     Brassicaceae="Bra",Cyperaceae="Cyp")


p2<-ggboxplot(pplants,x="abbrev",y="`Pplant/soil`",add="jitter",color="Status",
              fill="Status",size=0.1, width=0.8,ylab = "Shoot P accumulation",alpha=0.1)+
  theme_classic()+theme(axis.line.x = element_line())+ 
  theme(axis.text.y= element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text.x= element_text(size=16,angle=90),axis.title.x = element_blank())+ 
  theme(strip.text.x = element_text(size=16),strip.background = element_blank())+ 
  facet_wrap(~Site, scales="free",ncol=7)+
  theme(legend.position="bottom",legend.text = element_text(size=16),
        legend.title=element_text(size=18))
p2


###Statistical analysis P accumulation ####

#Test normality of shoot P data
shapiro.test(pplants$`Pplant/soil`)
#not normal

#test Homoscedasticity between groups
car::leveneTest(`Pplant/soil` ~Site/Status/Family, data=pplants)
#different variances

#Normalize data for ANOVA using box-cox transformation
library(bestNormalize)
obj_bestnormP<-bestNormalize(pplants$`Pplant/soil`)
obj_bestnormP
pplants$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data

#verify normality of transformed data
shapiro.test(pplants$bestnormP)
#OK

#verify Homoscedasticity of transformed data
car::leveneTest(bestnormP ~Site/Status/Family, data=pplants)
#OK


###### ANOVA with nested model shoot P ~ Site/Status/Family ####
res.aov<-aov(bestnormP ~ Site/Status/Family, data=pplants)
summary(res.aov)

# Df Sum Sq Mean Sq F value   Pr(>F)    
# Site                 6 184.45  30.741  509.02  < 2e-16 ***
#   Site:Status          7   6.23   0.890   14.74 1.34e-15 ***
#   Site:Status:Family  28  48.58   1.735   28.73  < 2e-16 ***
#   Residuals          211  12.74   0.060                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Get R2 values

df<-summary(res.aov)
sum(df[[1]][,2])
sqR<-df[[1]][,2]/sum(df[[1]][,2])
sqR

#[1] 0.73193619 0.02472710 0.19276973 0.05056698


###within sites  

##### ANOVA between Status for each Site ####
library(stats)
anova_results <- list()
leveneTest_results <-list()
shapiro_results<-list()
anova_results_bestnorm <- list()
leveneTest_results_bestnorm  <-list()
shapiro_results_bestnorm <-list()
results_wilx<-list()
# Get unique sites in the dataframe
sites <- unique(pplants$Site)

# Iterate over each site
for (Site in sites) {
  # Filter the dataframe for the current site
  site_df <- pplants[pplants$Site == Site, ]
  # Perform ANOVA
  res.aov <- aov(`Pplant/soil` ~ factor(Status), data =site_df)
  # Store ANOVA result
  anova_results[[as.character(Site)]] <- summary(res.aov)
  #perform levene_test
  res.levene<-leveneTest( `Pplant/soil` ~ factor(Status), data =site_df)
  leveneTest_results[[Site]] <-res.levene
  #shapiro test
  res.shap<-shapiro.test(site_df$`Pplant/soil`)
  shapiro_results[[Site]] <- res.shap
  #bestnorm
  obj_bestnorm<-boxcox(site_df$`Pplant/soil`)
  site_df$bestnorm<-obj_bestnorm$x.t
  res.aov <- aov(bestnorm ~ factor(Status), data =site_df)
  anova_results_bestnorm[[Site]] <- summary(res.aov)
  #perform levene_test
  res.levene<-leveneTest(bestnorm ~ factor(Status), data =site_df)
  leveneTest_results_bestnorm[[Site]] <- res.levene
  #shapiro test
  res.shap<-shapiro.test(site_df$bestnorm)
  shapiro_results_bestnorm[[Site]] <-res.shap
  ##for_sites_when_condition not respected
  res.wilx<-rstatix::wilcox_test(site_df,`Pplant/soil` ~ Status,p.adjust.method = "none")
  results_wilx[[Site]]<-res.wilx
}

# Print ANOVA results for each site
anova_results
leveneTest_results
shapiro_results
anova_results_bestnorm
leveneTest_results_bestnorm
shapiro_results_bestnorm
results_wilx
####Only significant differences in Shoot P for Galibier et Lautaret

##### ANOVA between Families for each site#####

library(stats)
library(car)
library(rstatix)

anova_results <- list()
leveneTest_results <- list()
shapiro_results <- list()
anova_results_bestnorm <- list()
leveneTest_results_bestnorm <- list()
shapiro_results_bestnorm <- list()
tukey_results <- list()
tukey_results_bn <- list()
pairwise_wilcox_results <- list()

sites <- unique(pplants$Site)

for (Site in sites) {
  
  site_df <- pplants[pplants$Site == Site, ]
  #Anova
  res.aov <- aov(`Pplant/soil`~ factor(Family), data = site_df)
  anova_results[[Site]] <- summary(res.aov)
  #Variance check
  res.levene <- car::leveneTest(`Pplant/soil` ~ factor(Family), data = site_df)
  leveneTest_results[[Site]] <- res.levene
  #Shapiro on residues
  res.shap <- shapiro.test(residuals(res.aov))  
  shapiro_results[[Site]] <- res.shap
  # Transform using boxcox
  obj_bestnorm <- bestNormalize::boxcox(site_df$`Pplant/soil`)
  site_df$bestnorm<-obj_bestnorm$x.t
  # Anova on Transformed data
  res.aov_bn <- aov(bestnorm ~ factor(Family), data = site_df)
  anova_results_bestnorm[[Site]] <- summary(res.aov_bn)
  #Verify assumptions
  res.levene_bn <- car::leveneTest(bestnorm ~ factor(Family), data = site_df)
  leveneTest_results_bestnorm[[Site]] <- res.levene_bn
  #Veirfy normality
  res.shap_bn <- shapiro.test(residuals(res.aov_bn))
  shapiro_results_bestnorm[[Site]] <- res.shap_bn
  #Pairwise comparisons
  if(res.shap_bn$p.value > 0.04 & res.levene_bn$`Pr(>F)`[1] > 0.05){
    
    # Tukey on BN
    tukey_bn <- TukeyHSD(res.aov_bn)
    tukey_results_bn[[Site]] <- tukey_bn
    
  } else {
    
    # airwise Wilcoxon
    pw <- pairwise.wilcox.test(site_df$`Pplant/soil`,
                               site_df$Family,
                               p.adjust.method = "fdr")
    
    pairwise_wilcox_results[[Site]] <- pw
  }
}

anova_results
leveneTest_results
shapiro_results

anova_results_bestnorm
leveneTest_results_bestnorm
shapiro_results_bestnorm

tukey_results
tukey_results_bn
pairwise_wilcox_results

library(multcompView)

tukey_to_letters <- function(tukey_obj) {
  
  # extract Tukey p-values
  tuk <- tukey_obj$`factor(Family)`
  
  # get p-values
  pvals <- tuk[, "p adj"]
  
  # names must match comparisons
  names(pvals) <- rownames(tuk)
  
  # convert to compact letters
  letters <- multcompLetters(pvals)$Letters
  
  return(letters)
}

letters_results <- lapply(tukey_results_bn, tukey_to_letters)
letters_results

########## Phylo-ANOVA analysis for within site P level differences between AM vs non-AM plants, give phylogeny

#do a subtree from plant megatree with selected species
#The megatree is ‘GBOTB.extended.WP.tre’ reported by Jin and Qian (2022).
devtools::install_github("jinyizju/U.PhyloMaker")
library("U.PhyloMaker")
#setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/codes-article-principal - Copie/Shoot P values")
sp.list <- read.csv("Species_list.csv", sep=';')
megatree <- read.tree('https://raw.githubusercontent.com/megatrees/plant_20221117/main/plant_megatree.tre')
gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/plant_20221117/main/plant_genus_list.csv', sep=",")
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type = 1, scenario = 3)
result
#tree as newick to open in itol
write.tree(result$phylo, "output_tree.tre")
#write table with plant families
write.table(result$sp.list, file="sp_list_fam.txt", quote=FALSE, sep="\t")

library(readxl)
samples_df <- read_excel("tableitsx_phyloseq.xlsx", sheet = "Metadata")
as.data.frame(samples_df)

#get mean Shoot P concentration per species
shootP_means_sp<-aggregate(samples_df[, 10], list(samples_df$Espèce), mean)
colnames(shootP_means_sp)<-c("Species", "avg_shootP")
#Reorder the values
sp.list<-sp.list[order(sp.list$Species),]
#adjust names: no spaces
sp.list$Species<-sub(" ", "_", sp.list$Species)
row.names(sp.list)<-sp.list$Species
status<-as.factor(as.vector(sp.list$Status))
names(status)=sp.list$Species

shootP_means_sp<-shootP_means_sp[order(shootP_means_sp$Species),]
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)=names(status)


#install.packages('geiger')
library(geiger)


#Site by site comparissons give the phylogeny
#GAL
#get mean Shoot P per species
gal<-subset(samples_df, Site=="Galibier")
shootP_means_sp<-aggregate(gal[, 10], list(gal$Espèce, gal$Status), mean)
colnames(shootP_means_sp)<-c("Species", "Status", "avg_shootP")
shootP_means_sp$Species<-sub(" ", "_", shootP_means_sp$Species)
gp<-as.factor(as.vector(shootP_means_sp$Status))
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)<-shootP_means_sp$Species
names(gp)<-shootP_means_sp$Species

aov.phylo(shootP~gp, result$phylo, nsim = 50)

# Analysis of Variance Table
# 
# Response: dat
# Df Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1 0.9322 0.93220  1.3493 0.30999           0.1961
# Residuals  4 2.7635 0.69088                                 
# 
# Call:
#   lm(formula = dat ~ group)
# 
# Coefficients:
#   (Intercept)  groupNon-mycorhizal  
# 2.4806              -0.7883  

#Lau
#get mean Shoot P per species
lau<-subset(samples_df, Site=="Lautaret")
shootP_means_sp<-aggregate(lau[, 10], list(lau$Espèce, lau$Status), mean)
colnames(shootP_means_sp)<-c("Species", "Status", "avg_shootP")
shootP_means_sp$Species<-sub(" ", "_", shootP_means_sp$Species)
gp<-as.factor(as.vector(shootP_means_sp$Status))
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)<-shootP_means_sp$Species
names(gp)<-shootP_means_sp$Species

aov.phylo(shootP~gp, result$phylo, nsim = 50)

# Analysis of Variance Table
# 
# Response: dat
# Df Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1 1.4243 1.42431  1.8115 0.24956           0.1569
# Residuals  4 3.1451 0.78629                                 
# 
# Call:
#   lm(formula = dat ~ group)
# 
# Coefficients:
#   (Intercept)  groupNon-mycorhizal  
# 2.6583              -0.9744  

#Cla
#get mean Shoot P per species
cla<-subset(samples_df, Site=="Clarée")
shootP_means_sp<-aggregate(cla[, 10], list(cla$Espèce, cla$Status), mean)
colnames(shootP_means_sp)<-c("Species", "Status", "avg_shootP")
shootP_means_sp$Species<-sub(" ", "_", shootP_means_sp$Species)
gp<-as.factor(as.vector(shootP_means_sp$Status))
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)<-shootP_means_sp$Species
names(gp)<-shootP_means_sp$Species

aov.phylo(shootP~gp, result$phylo, nsim = 50)

# Analysis of Variance Table
# 
# Response: dat
# Df  Sum-Sq  Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1 0.03604 0.036037  0.1403 0.72699           0.6863
# Residuals  4 1.02741 0.256854                                 
# 
# Call:
#   lm(formula = dat ~ group)
# 
# Coefficients:
#   (Intercept)  groupNon-mycorhizal  
# 1.288                0.155  


#Cham
#get mean Shoot P per species
cha<-subset(samples_df, Site=="Chamrousse")
shootP_means_sp<-aggregate(cha[, 10], list(cha$Espèce, cha$Status), mean)
colnames(shootP_means_sp)<-c("Species", "Status", "avg_shootP")
shootP_means_sp$Species<-sub(" ", "_", shootP_means_sp$Species)
gp<-as.factor(as.vector(shootP_means_sp$Status))
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)<-shootP_means_sp$Species
names(gp)<-shootP_means_sp$Species

aov.phylo(shootP~gp, result$phylo, nsim = 50)

# Analysis of Variance Table
# 
# Response: dat
# Df  Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1 0.15952 0.15952 0.27266 0.62914           0.6863
# Residuals  4 2.34027 0.58507                                 
# 
# Call:
#   lm(formula = dat ~ group)
# 
# Coefficients:
#   (Intercept)  groupNon-mycorhizal  
# 1.2472               0.3261  

#Com
#get mean Shoot P per species
com<-subset(samples_df, Site=="Commelle")
shootP_means_sp<-aggregate(com[, 10], list(com$Espèce, com$Status), mean)
colnames(shootP_means_sp)<-c("Species", "Status", "avg_shootP")
shootP_means_sp$Species<-sub(" ", "_", shootP_means_sp$Species)
gp<-as.factor(as.vector(shootP_means_sp$Status))
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)<-shootP_means_sp$Species
names(gp)<-shootP_means_sp$Species

aov.phylo(shootP~gp, result$phylo, nsim = 50)

# Analysis of Variance Table
# 
# Response: dat
# Df  Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1  1.9684  1.9685 0.77186 0.42925            0.451
# Residuals  4 10.2011  2.5503                                 
# 
# Call:
#   lm(formula = dat ~ group)
# 
# Coefficients:
#   (Intercept)  groupNon-mycorhizal  
# 3.626               -1.146  

#Per
#get mean Shoot P per species
per<-subset(samples_df, Site=="Perouges")
shootP_means_sp<-aggregate(per[, 10], list(per$Espèce, per$Status), mean)
colnames(shootP_means_sp)<-c("Species", "Status", "avg_shootP")
shootP_means_sp$Species<-sub(" ", "_", shootP_means_sp$Species)
gp<-as.factor(as.vector(shootP_means_sp$Status))
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)<-shootP_means_sp$Species
names(gp)<-shootP_means_sp$Species

aov.phylo(shootP~gp, result$phylo, nsim = 50)

# Analysis of Variance Table
# 
# Response: dat
# Df  Sum-Sq Mean-Sq  F-value  Pr(>F) Pr(>F) given phy
# group      1 0.02492 0.02492 0.042678 0.84642           0.9412
# Residuals  4 2.33551 0.58388
# 
# Call:
#   lm(formula = dat ~ group)
# 
# Coefficients:
#   (Intercept)  groupNon-mycorhizal
# 1.9450              -0.1289

#Doua
#get mean Shoot P per species
dou<-subset(samples_df, Site=="La Doua")
shootP_means_sp<-aggregate(dou[, 10], list(dou$Espèce, dou$Status), mean)
colnames(shootP_means_sp)<-c("Species", "Status", "avg_shootP")
shootP_means_sp$Species<-sub(" ", "_", shootP_means_sp$Species)
gp<-as.factor(as.vector(shootP_means_sp$Status))
shootP<-as.vector(shootP_means_sp$avg_shootP)
names(shootP)<-shootP_means_sp$Species
names(gp)<-shootP_means_sp$Species

aov.phylo(shootP~gp, result$phylo, nsim = 50)

# Analysis of Variance Table
# 
# Response: dat
# Df  Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1  3.2178  3.2178  1.4244 0.28622           0.2941
# Residuals  5 11.2955  2.2591                                 
# 
# Call:
#   lm(formula = dat ~ group)
# 
# Coefficients:
#   (Intercept)  groupNon-mycorhizal  
# 4.081               -1.370  