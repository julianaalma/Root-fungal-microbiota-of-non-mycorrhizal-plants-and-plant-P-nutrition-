###Fig. S7A and C: Glomeromycota Relative abundance and Mycorrhization levels

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

setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/AM_Glomeromycota_FigS5")

###Import OTU data ####
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

#create phyloseq object
phy_obj_plant<- phyloseq(OTU, TAX, samples)
phy_obj_plant<- prune_species(speciesSums(phy_obj_plant) > 0, phy_obj_plant) ### to keep OTUs only present in plants (not only in soils)

#calculate relative abundances and transform
phylum<-tax_glom(phy_obj_plant,taxrank = rank_names(phy_obj_plant)[2]) #agglomerate by phylum : tax. rank =2
phylum_ra<-transform_sample_counts(phylum, function(x) x/sum(x) ) 
phylum_log<-otu_table(phylum_ra)+1
phylum_log<-log10(otu_table(phylum_log))
phylum2 <- phyloseq(otu_table(phylum_log), tax_table(phylum), sample_data(phylum))

#select Glomeromycota data 
glom_log<-subset_taxa(phylum2,Phylum=="p__Glomeromycota")


#### Fig S7 C: Rel. Abu. Glomeromycota

#Rename factors
glom_log_tab<-psmelt(glom_log)
glom_log_tab$Abundance<-as.numeric(glom_log_tab$Abundance)
glom_log_tab$abbrev<-dplyr::recode_factor(glom_log_tab$sample_Family,Asteraceae="Ast",Ranunculaceae="Ran",
                                           Geraniaceae="Ger",Poaceae="Poa", Caryophyllaceae= "Car",
                                           Brassicaceae="Bra",Cyperaceae="Cyp")
glom_log_tab<-rstatix::reorder_levels(glom_log_tab,sample_Family,c("Asteraceae","Geraniaceae","Ranunculaceae","Poaceae","Caryophyllaceae","Brassicaceae","Cyperaceae"))
glom_log_tab<-glom_log_tab%>% reorder_levels("Site",order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","La Doua","Commelle"))

#Plot Rel. Abu. Glomeromycota

p1<-ggboxplot(glom_log_tab,x="abbrev",y="Abundance",add="jitter",color="Status",
              fill="Status",size=0.3, width=0.8, ylab = "Glomeromycota relative abundance (log10(x+1))",alpha=0.1)+
  stat_summary(fun.y=mean, geom="point", shape=20, color="black")+
  theme_classic()+theme(axis.line.x = element_line())+ 
  theme(axis.text.y= element_text(size=10),axis.title.y = element_text(size=10))+ 
  theme(axis.text.x= element_text(size=10,angle=90),axis.title.x = element_blank())+ 
  theme(strip.text.x = element_text(size=10),strip.background = element_blank())+ 
  facet_wrap(~Site, scales="free",ncol=7)+ ylim(0,0.15)+
  scale_color_manual(values=c("#f8756bff","#00bdc2ff"))+
  scale_fill_manual(values=c("#f8756bff","#00bdc2ff"))+
  theme(legend.position="bottom",legend.text = element_text(size=10), legend.title=element_text(size=10))
p1


#Statistics
#comparing Glomeromycota Relative abundance between Myc. and non-Myc. plants per site

library(car)
library(bestNormalize)

#List all analyses
anova_results <- list()
leveneTest_results <-list()
shapiro_results<-list()
anova_results_bestnorm <- list()
leveneTest_results_bestnorm  <-list()
shapiro_results_bestnorm <-list()
results_wilx<-list()

# Get unique sites in the dataframe
sites <- unique(glom_log_tab$Site)

# Iterate all analyses over each site
for (Site in sites) {
  # Filter the dataframe for the current site
  site_df <- glom_log_tab[glom_log_tab$Site == Site, ]
  # Perform ANOVA
  res.aov <- aov(Abundance ~ factor(Status), data =site_df)
  # Store ANOVA result
  anova_results[[as.character(Site)]] <- summary(res.aov)
  #perform levene_test
  res.levene<-leveneTest(Abundance ~ factor(Status), data =site_df)
  leveneTest_results[[Site]] <-res.levene
  #shapiro test
  res.shap<-shapiro.test(site_df$Abundance)
  shapiro_results[[Site]] <- res.shap
  #Normalize data
  obj_bestnorm<-bestNormalize(site_df$Abundance)
  site_df$bestnorm<-obj_bestnorm$x.t
  res.aov <- aov(bestnorm ~ factor(Status), data =site_df)
  anova_results_bestnorm[[Site]] <- summary(res.aov)
  #perform levene_test on normalized data
  res.levene<-leveneTest(bestnorm ~ factor(Status), data =site_df)
  leveneTest_results_bestnorm[[Site]] <- res.levene
  #perform levene_test on normalized data
  res.shap<-shapiro.test(site_df$bestnorm)
  shapiro_results_bestnorm[[Site]] <-res.shap
  #Non-parametric test for non-normal data 
  res.wilx<-rstatix::wilcox_test(site_df,Abundance~ Status,p.adjust.method = "none")
  results_wilx[[Site]]<-res.wilx
}


# Print results for each site
anova_results
leveneTest_results
shapiro_results
anova_results_bestnorm
leveneTest_results_bestnorm
shapiro_results_bestnorm
results_wilx

#No difference between Myc and Non-Myc in sites Clarée and Chamrousse, Significant differences in all other sites 


#####################################

#### Fig S7 A: Mycorrhization levels per treatment (percentage of root pieces mycorrhized, per treatment)
#values are not available for individual plants. 
#From each treatment, approx. 10 fine root pieces from each individual plant were taken and pooled into a 
#composite sample to asses mycorrhization levels within the treatment. 
#20-30 root pieces were analysed for each treatment.

#Plot distribution of mycorrhization degree

myc_score<-read_excel("Score_mycorrhization_JA.xlsx", sheet="Myc_Scores")

score_summary <- myc_score %>%
  group_by(SiteFamily, Mycorrhization_Score) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n / sum(n)) * 100)

score_summary$Site <- sapply(score_summary$SiteFamily, function(x) {
  if (grepl("Cham", x, ignore.case = TRUE)) {
    return("Chamrousse")
  } else if (grepl("Cla", x, ignore.case = TRUE)) {
    return("Clarée")
  } else if (grepl("Com", x, ignore.case = TRUE)) {
    return("Commelle")
  } else if (grepl("Lau", x, ignore.case = TRUE)) {
    return("Lautaret")
  } else if (grepl("Gal", x, ignore.case = TRUE)) {
    return("Galibier")
  } else if (grepl("Dou", x, ignore.case = TRUE)) {
    return("Doua")
  } else if (grepl("Per", x, ignore.case = TRUE)) {
    return("Perouges")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

score_summary$sample_Family <- sapply(score_summary$SiteFamily, function(x) {
  if (grepl("Ast", x, ignore.case = TRUE)) {
    return("Asteraceae")
  } else if (grepl("Bra", x, ignore.case = TRUE)) {
    return("Brassicaceae")
  } else if (grepl("Car", x, ignore.case = TRUE)) {
    return("Caryophyllaceae")
  } else if (grepl("Cyp", x, ignore.case = TRUE)) {
    return("Cyperaceae")
  } else if (grepl("Ger", x, ignore.case = TRUE)) {
    return("Geraniaceae")
  } else if (grepl("Poa", x, ignore.case = TRUE)) {
    return("Poaceae")
  } else if (grepl("Ran", x, ignore.case = TRUE)) {
    return("Ranunculaceae")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

#Plot distributions
score_summary$good_sample_name<-paste(score_summary$Site,score_summary$sample_Family,sep="_")
score_summary$Mycorrhization_Score<-as.factor(score_summary$Mycorrhization_Score)
score_summary$abbrev<-dplyr::recode_factor(score_summary$sample_Family,Asteraceae="Ast",Ranunculaceae="Ran",
                                           Geraniaceae="Ger",Poaceae="Poa", Caryophyllaceae= "Car",
                                           Brassicaceae="Bra",Cyperaceae="Cyp")
score_summary<-rstatix::reorder_levels(score_summary,sample_Family,c("Asteraceae","Geraniaceae","Ranunculaceae","Poaceae","Caryophyllaceae","Brassicaceae","Cyperaceae"))
score_summary<-score_summary%>% reorder_levels("Site",order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","Doua","Commelle"))
score_summary$Percentage_of_mycorrhization<-dplyr::recode_factor(score_summary$Mycorrhization_Score,"0"="0","1"="< 1 %",
                                                                 "2"="< 10 %","3" ="< 50 %", "4" ="> 50 %",
                                                                 "5" ="> 90 %")

colors<-c("#ffffff","#ffcdfe","#fe66cb","#ca51afff","#6a2d65ff","#40415aff")

p2<-ggplot(score_summary,aes(abbrev,percentage,color=Percentage_of_mycorrhization,fill=Percentage_of_mycorrhization))+
  geom_bar(position="fill", stat="identity",width = 0.8)+scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+facet_wrap(~Site,nrow=1,scales="free_y")+theme_minimal()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(),
        axis.text.x = element_text(size=10, angle=90),axis.title.x = element_blank())+
  theme(axis.text.y= element_text(size=10), axis.title.y = element_text(size=10), axis.line = element_line(), axis.ticks = element_line())+
  theme(strip.text.x = element_text(size=10),strip.background = element_blank())+ scale_y_continuous ()+
  ylab("Percentage of mycorrhized root pieces (presence of arbuscules and/or vesicules)")

p2



#### Correlation mycorrhization level X Glomeromycota rel. abu.

#Average Rel. Abu. Glomeromycota per treatment
glom_log_tab<-psmelt(glom_log)
glom_log_avg<-aggregate(glom_log_tab$Abundance, list(glom_log_tab$SiteFamily), FUN=mean)
names(glom_log_avg)<-c("SiteFamily","Abundance")


#Attribute Site names
glom_log_avg$Site <- sapply(glom_log_avg$SiteFamily, function(x) {
  if (grepl("Cham", x, ignore.case = TRUE)) {
    return("Chamrousse")
  } else if (grepl("Cla", x, ignore.case = TRUE)) {
    return("Claree")
  } else if (grepl("Com", x, ignore.case = TRUE)) {
    return("Commelle")
  } else if (grepl("Lau", x, ignore.case = TRUE)) {
    return("Lautaret")
  } else if (grepl("Gal", x, ignore.case = TRUE)) {
    return("Galibier")
  } else if (grepl("Dou", x, ignore.case = TRUE)) {
    return("Doua")
  } else if (grepl("Per", x, ignore.case = TRUE)) {
    return("Perouges")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

#Attribute Family names
glom_log_avg$sample_Family <- sapply(glom_log_avg$SiteFamily, function(x) {
  if (grepl("Ast", x, ignore.case = TRUE)) {
    return("Asteraceae")
  } else if (grepl("Bra", x, ignore.case = TRUE)) {
    return("Brassicaceae")
  } else if (grepl("Car", x, ignore.case = TRUE)) {
    return("Caryophyllaceae")
  } else if (grepl("Cyp", x, ignore.case = TRUE)) {
    return("Cyperaceae")
  } else if (grepl("Ger", x, ignore.case = TRUE)) {
    return("Geraniaceae")
  } else if (grepl("Poa", x, ignore.case = TRUE)) {
    return("Poaceae")
  } else if (grepl("Ran", x, ignore.case = TRUE)) {
    return("Ranunculaceae")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

glom_log_avg$good_sample_name<-paste(glom_log_avg$Site,glom_log_avg$sample_Family,sep="_")


##### Mycorrhization data: percentage of root pieces colonized

#Get percentages of root pieces with 0 colonization (score=0)
score_summary2 <- score_summary %>% filter(Mycorrhization_Score == 0)
score_summary3 <- score_summary2 %>% mutate(percentage_colonized = 100 - percentage)

### Adjust treatment names
score_summary3$good_sample_name<-paste(score_summary3$Site,score_summary3$sample_Family,sep="_")

###merge both data files
glom_score_RA<-merge(score_summary3,glom_log_avg,by="good_sample_name")

## Add Biome information
glom_score_RA$Biome<- sapply(glom_score_RA$SiteFamily.x, function(x) {
  if (grepl("Cham", x, ignore.case = TRUE)) {
    return("Alpine")
  } else if (grepl("Cla", x, ignore.case = TRUE)) {
    return("Alpine")
  } else if (grepl("Com", x, ignore.case = TRUE)) {
    return("Lowland")
  } else if (grepl("Lau", x, ignore.case = TRUE)) {
    return("Alpine")
  } else if (grepl("Gal", x, ignore.case = TRUE)) {
    return("Alpine")
  } else if (grepl("Dou", x, ignore.case = TRUE)) {
    return("Lowland")
  } else if (grepl("Per", x, ignore.case = TRUE)) {
    return("Lowland")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

## Add Myc. Status information
glom_score_RA$Status <- sapply(glom_score_RA$SiteFamily.x, function(x) {
  if (grepl("Ast", x, ignore.case = TRUE)) {
    return("Myc")
  } else if (grepl("Bra", x, ignore.case = TRUE)) {
    return("Non-myc")
  } else if (grepl("Car", x, ignore.case = TRUE)) {
    return("Non-myc")
  } else if (grepl("Cyp", x, ignore.case = TRUE)) {
    return("Non-myc")
  } else if (grepl("Ger", x, ignore.case = TRUE)) {
    return("Myc")
  } else if (grepl("Poa", x, ignore.case = TRUE)) {
    return("Myc")
  } else if (grepl("Ran", x, ignore.case = TRUE)) {
    return("Myc")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

## Plot overall correlation across treatments (42 points)
p3<-ggplot(glom_score_RA, aes(x = Abundance, y = percentage_colonized)) +
  geom_point(aes(color = Status, shape=Biome),size=3) +
  geom_smooth(method="lm",se=TRUE,level=0.95, color="black") +
  theme_classic2()+ylab("Percentage of mycorrhized root pieces (presence of arbuscules and/or vesicules)") +xlab("Mean Glomeromycota relative abundance (log10(x+1)")
p3

#Average mycorrhization rate in AM plants
Myc_AM_avg<-aggregate(glom_score_RA$percentage_colonized, list(glom_score_RA$Status), FUN=mean)
Myc_AM_avg

# Group.1        x
# 1     Myc 69.78011
# 2 Non-myc 14.09742

#Average Glomeromycota RA in AM plants
phylum3 <- phyloseq(otu_table(phylum_ra), tax_table(phylum), sample_data(phylum))
glom_RA<-subset_taxa(phylum3,Phylum=="p__Glomeromycota")
glom_RA2<-psmelt(glom_RA)
glom_RA_myc<-aggregate(glom_RA2$Abundance, list(glom_RA2$Status), FUN=mean)
glom_RA_myc

# Group.1           x
# 1     Mycorhizal 0.040235165
# 2 Non-mycorhizal 0.005992171

#Correlation statistics
cor.test(glom_score_RA$Abundance,glom_score_RA$percentage_colonized,method="pearson")

# Pearson's product-moment correlation
# 
# data:  glom_score_RA$Abundance and glom_score_RA$percentage_colonized
# t = 2.3686, df = 34, p-value = 0.02368
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.05455022 0.62731484
# sample estimates:
#       cor 
# 0.3763413 



cor.test(glom_score_RA$Abundance,glom_score_RA$percentage_colonized,method="spearman")

# Spearman's rank correlation rho
# 
# data:  glom_score_RA$Abundance and glom_score_RA$percentage_colonized
# S = 4967.5, p-value = 0.03069
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3606794 


############################################
