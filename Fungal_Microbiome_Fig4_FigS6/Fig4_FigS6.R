###Fig 4 and S6
### Comparison of fungal communities between treatments and pant families, with and without AM Glomeromycota taxa 
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


#set wd
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/PB_old_codes-article-principal - Copie/Microbiome analyses")

###import data######
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

phy_obj_plant<- phyloseq(OTU, TAX, samples)
phy_obj_plant<- prune_species(speciesSums(phy_obj_plant) > 0, phy_obj_plant) ### to discard OTUs with rel. abu. in roots = 0 (only present in soil samples)

#transform data in log(RA+1)
##relative abundance for each OTU abundOTU/abundtotOTU
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
phy_obj_plant_ra<-otu_table(phy_obj_plant_ra)+1
phy_obj_plant_log<-log10(otu_table(phy_obj_plant_ra))
###reassemble otu_table with taxonomy and metadata
phy_obj_plant_log<- phyloseq(otu_table(phy_obj_plant_log), TAX, samples)


#####Permanova on whole dataset#####

#####character to factor to conduct pcoA
tab_otu<-otu_table(phy_obj_plant_log)
tab_otu<-data.frame(tab_otu)
tab_otu<-t(tab_otu)
total<-cbind(tab_otu,samples_df)
total$Status[total$Status=="Nonmyc"] <- "Non-mycorhizal"
total$Site<-as.factor(total$Site)
total$Status<-as.factor(total$Status)
total$Family<-as.factor(total$Family)
total$Type_site<-as.factor(total$Type_site)


#####Permanova with fully nested model, and permutations constrained by "site" using strata
permanova<-adonis2(formula = tab_otu ~ Site/Status/Family, data = total, permutations = 99999, method = "bray", by = "terms", strata=total$Site)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = tab_otu ~ Site/Status/Family, data = total, permutations = 9999, method = "bray", by = "terms", strata = total$Site)
# Df SumOfSqs      R2      F Pr(>F)    
# Site                 6   18.684 0.16297 9.8305  1e-04 ***
#   Site:Status          7    4.817 0.04201 2.1723  1e-04 ***
#   Site:Status:Family  28   24.310 0.21204 2.7408  1e-04 ***
#   Residual           211   66.838 0.58298                  
# Total              252  114.649 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Calculate exact P values from pseudo-F results of Permutations
#exact P = proportion of Pseudo F > Observed F, for each term

##Site
#Extract the observed F-value
observed_F<- permanova$F[1]
# Extract all the permutation F-values
perm_F <- attr(permanova,"F.perm")[,1]
# Calculate the exact p-value
exact_p <- sum(perm_F >= observed_F) / length(perm_F)
# Print the exact p-value
print(exact_p)
#equals 0

##Status
#Extract the observed F-value
observed_F<- permanova$F[2]
# Extract all the permutation F-values
perm_F <- attr(permanova,"F.perm")[,2]
# Calculate the exact p-value
exact_p <- sum(perm_F >= observed_F) / length(perm_F)
# Print the exact p-value
print(exact_p)
#equals 0

##Family
#Extract the observed F-value
observed_F<- permanova$F[3]
# Extract all the permutation F-values
perm_F <- attr(permanova,"F.perm")[,3]
# Calculate the exact p-value
exact_p <- sum(perm_F >= observed_F) / length(perm_F)
# Print the exact p-value
print(exact_p)
#equals 0


#####Figure 4 A, B and C####

######Panel A: Capscale by site###############
capscaletry<- capscale(tab_otu ~ Site, data=total,dist="bray")
summary(capscaletry)
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$CAP2
p1<-ggplot(total, aes(x= CAP1, y=MDS1,color = Site)) + 
  stat_ellipse(aes(fill=factor(total$Site)), geom = "polygon",alpha = 0,show.legend=FALSE,lwd=1) +
  geom_point(size=3) + theme_classic()+
  scale_color_manual(values=c("#84eb7e","#83a205","#e95c3e","#01d0b5","#c39741","#8eaaff","#8245b8"))+
  scale_fill_manual(values=c( "#84eb7e","#83a205","#e95c3e","#01d0b5","#c39741","#8eaaff","#8245b8"))+
  coord_cartesian(xlim=c(-2, 2), ylim=c(-3, 3)) +
  #scale_shape_manual(values=c(1,16))+
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),legend.text = element_text(size=24),
        axis.text = element_text(size=24),legend.title=element_text(size=24),axis.title=element_text(size=18))+ 
  labs(color="Sampling sites",x="CAP1 (6%)",y="CAP2 (3.2%)")+theme(legend.position="bottom")
p1

#######Panel B: Capescale by family######
capscaletry<- capscale(tab_otu ~ Family, data=total,dist="bray")
summary(capscaletry)
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$CAP2
p2<-ggplot(total, aes(x= CAP1, y=MDS1,color = Family)) + 
  stat_ellipse(aes(fill=factor(total$Family)), geom = "polygon",alpha = 0,lwd=1,show.legend=FALSE) +
  geom_point(size=3) + theme_classic()+
  scale_color_manual(values=c("#8c7209ff","#dfafd6ff","#c2f2efff","#e3cec5ff","#ad4498ff","#3ac42dff","#492359ff"))+
  scale_fill_manual(values=c("#8c7209ff","#dfafd6ff","#c2f2efff","#e3cec5ff","#ad4498ff","#3ac42dff","#492359ff"))+
  coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 4.5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  #scale_shape_manual(values=c(1,16))+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),legend.text = element_text(size=24),
        axis.text = element_text(size=24),legend.title=element_text(size=24),axis.title=element_text(size=24))+
  labs(color="Plant Families",x="CAP1 (2.2%)",y="CAP2 (1.4%)")+
  theme(legend.position="bottom")
p2

######Panel C: Capscale by biomes and status####
total$Biome_Status<-paste(total$Type_site,total$Status)
capscaletry<- capscale(tab_otu ~ Type_site*Status, data=total,dist="bray")
summary(capscaletry)
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$CAP2

###make capscale fig
p3<-ggplot(total, aes(x= CAP1, y=MDS1,color = Type_site)) + 
  stat_ellipse(aes(fill=factor(total$Status)), geom = "polygon", alpha = 0,show.legend=FALSE,
               lwd=1) +
  geom_point(size=4,aes(shape=Status))+ theme_classic()+
  scale_color_manual(labels=c("Alpine (low P)","Lowland (high P)"),
                     values=c("#1b9e77ff","#FFC107"))+
  scale_fill_manual(values=c("#1b9e77ff","#FFC107"))+
  scale_shape_manual(values=c(1,16))+
  coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 4.5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),legend.text = element_text(size=24),
        axis.text = element_text(size=24),legend.title = element_text(size=24),axis.title = element_text(size=24))+ 
  labs(color="Biomes",x="CAP1 (5.9%)",y="CAP2 (0.8%)")+theme(legend.position="bottom")
p3


### put together figure
ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1, legend="none")


####Figure S6 : the same analyses without AM Glomeromycota fungi ####

######PerMANOVA WITHOUT GLOMEROMYCOTA#####

##recalculate RA log
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

phy_obj_plant<- phyloseq(OTU, TAX, samples)
phy_obj_plant<- prune_species(speciesSums(phy_obj_plant) > 0, phy_obj_plant) ### to discard OTUs with rel. abu. in roots = 0 (only present in soil samples) 
no_glomero<-subset_taxa(phy_obj_plant,!Phylum=="p__Glomeromycota") #exclude Glomeromycota OTUs

####transform data in log(RA+1)#
##relative abundance for each OTU abundOTU/abundtotOTU
no_glomero_ra<-transform_sample_counts(no_glomero, function(x) x/sum(x) )
no_glomero_ra<-otu_table(no_glomero_ra)+1
no_glomero_log<-log10(otu_table(no_glomero_ra))

###reassemble otu_table with taxonomy and metadata
no_glomero_log<- phyloseq(otu_table(no_glomero_log), no_glomero@tax_table, samples)

#####character to factor to conduct pcoA
tab_otu<-otu_table(no_glomero_log)
tab_otu<-data.frame(tab_otu)
tab_otu<-t(tab_otu)
total<-cbind(tab_otu,samples_df)
total<-cbind(tab_otu,no_glomero_log@sam_data)
total$Status[total$Status=="Nonmyc"] <- "Non-mycorhizal"
total$Site<-as.factor(total$Site)
total$Status<-as.factor(total$Status)
total$Family<-as.factor(total$Family)
total$Type_site<-as.factor(total$Type_site)


##Permanova##
permanova<-adonis2(formula = tab_otu ~ Site/Status/Family, data = total, permutations = 99999, method = "bray", by = "terms", strata=total$Site)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 99999
# 
# adonis2(formula = tab_otu ~ Site/Status/Family, data = total, permutations = 99999, method = "bray", by = "terms", strata = total$Site)
# Df   SumOfSqs         R2       F Pr(>F)    
# Site                 6  18.869534 0.16489162 9.98577  1e-05 ***
#   Site:Status          7   4.723794 0.04127893 2.14271  1e-05 ***
#   Site:Status:Family  28  24.390194 0.21313397 2.76585  1e-05 ***
#   Residual           211  66.452456 0.58069548                   
# Total              252 114.435979 1.00000000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#####Figure S8 A, B and C####

######Capscale for datsets without Glomeromycota####

###Panel A: capscale by sites
capscaletry<- capscale(tab_otu ~ Site, data=total,dist="bray")
summary(capscaletry)
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$CAP2
p1<-ggplot(total, aes(x= CAP1, y=MDS1,color = Site)) + 
  stat_ellipse(aes(fill=factor(total$Site)), geom = "polygon",alpha = 0,show.legend=FALSE,lwd=1) +
  geom_point(size=3) + theme_classic()+
  scale_color_manual(values=c("#84eb7e","#83a205","#e95c3e","#01d0b5","#c39741","#8eaaff","#8245b8"))+
  scale_fill_manual(values=c( "#84eb7e","#83a205","#e95c3e","#01d0b5","#c39741","#8eaaff","#8245b8"))+
  coord_cartesian(xlim=c(-2, 2), ylim=c(-3, 3)) +
  #scale_shape_manual(values=c(1,16))+
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),legend.text = element_text(size=24),
        axis.text = element_text(size=24),legend.title=element_text(size=24),axis.title=element_text(size=18))+ 
  labs(color="Sampling sites",x="CAP1 (6%)",y="CAP2 (3.2%)")+theme(legend.position="bottom")
p1

###Panel B: capscale by family###
capscaletry<- capscale(tab_otu ~ Family, data=total,dist="bray")
summary(capscaletry)
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$CAP2
p2<-ggplot(total, aes(x= CAP1, y=MDS1,color = Family)) + 
  stat_ellipse(aes(fill=factor(total$Family)), geom = "polygon",alpha = 0,lwd=1,show.legend=FALSE) +
  geom_point(size=3) + theme_classic()+
  scale_color_manual(values=c("#8c7209ff","#dfafd6ff","#c2f2efff","#e3cec5ff","#ad4498ff","#3ac42dff","#492359ff"))+
  scale_fill_manual(values=c("#8c7209ff","#dfafd6ff","#c2f2efff","#e3cec5ff","#ad4498ff","#3ac42dff","#492359ff"))+
  coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 4.5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  #scale_shape_manual(values=c(1,16))+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),legend.text = element_text(size=24),
        axis.text = element_text(size=24),legend.title=element_text(size=24),axis.title=element_text(size=24))+
  labs(color="Plant Families",x="CAP1 (2.2%)",y="CAP2 (1.4%)")+
  theme(legend.position="bottom")
p2

###Panel C: Capscale by biomes and status
total$Biome_Status<-paste(total$Type_site,total$Status)
capscaletry<- capscale(tab_otu ~ Type_site*Status, data=total,dist="bray")
summary(capscaletry)
x <- as.data.frame(scores(capscaletry, display = "sites", choices=c(1,2)))###extract coordinates
total$CAP1 <- x$CAP1
total$MDS1 <- x$CAP2
###make capscale fig
p3<-ggplot(total, aes(x= CAP1, y=MDS1,color = Type_site)) + 
  stat_ellipse(aes(fill=factor(total$Status)), geom = "polygon", alpha = 0,show.legend=FALSE,
               lwd=1) +
  geom_point(size=4,aes(shape=Status))+ theme_classic()+
  scale_color_manual(labels=c("Alpine (low P)","Lowland (high P)"),
                     values=c("#1b9e77ff","#FFC107"))+
  scale_fill_manual(values=c("#1b9e77ff","#FFC107"))+
  scale_shape_manual(values=c(1,16))+
  coord_cartesian(xlim=c(-3, 3), ylim=c(-3, 4.5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted")+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),legend.text = element_text(size=24),
        axis.text = element_text(size=24),legend.title = element_text(size=24),axis.title = element_text(size=24))+ 
  labs(color="Biomes",x="CAP1 (5.9%)",y="CAP2 (0.8%)")+theme(legend.position="bottom")
p3

### put together figure
ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1, legend="none")

