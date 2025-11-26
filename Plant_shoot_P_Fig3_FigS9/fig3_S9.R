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
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/1_Plant_shoot_P_Fig3")
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
