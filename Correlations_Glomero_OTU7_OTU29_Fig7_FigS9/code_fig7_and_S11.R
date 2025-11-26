####Code for Figure 7 and S11
### Correlations for each group of interest, presented in figures 7 and S9
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
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Correlations_Glomero_OTU7_OTU29_Fig7_FigS9")

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


#### Figure 7 ####

####transform data in log(RA+1)##
##relative abundance for each OTU abundOTU/abundtotOTU
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
phy_obj_plant_ra<-otu_table(phy_obj_plant_ra)+1
phy_obj_plant_log<-log10(otu_table(phy_obj_plant_ra))
###reassemble otu_table with taxonomy and metadata
phy_obj_plant_log<- phyloseq(otu_table(phy_obj_plant_log), TAX, samples)


##extract abundance OTU7 and OTU 29 #####
y4<- psmelt(phy_obj_plant_log)
y4$Abundance<-as.numeric(y4$Abundance)

y5<-filter(y4,OTU=="Otu000029"| OTU=="Otu000007")
myc<-filter(y5,Status=="Mycorhizal")
NM<-filter(y5,Status=="Non-mycorhizal")


#####Graphics OTU 7 and 29 in NM plants ####

p1<-ggplot(NM, aes(x = Abundance, y = Pplant.soil)) +
  geom_point(aes(color = Type_site,shape=sample_Family),size=3) + xscale("log10",.format=TRUE)+ 
  scale_color_manual(values=c("#004D40","#FFC107"))+
  scale_shape_manual(values=c(15,19,17))+
  geom_smooth(aes(color=Type_site),method="lm",se=TRUE,level=0.95) + facet_wrap(~OTU,nrow=1,axes = "all",axis.labels = "all")+
  theme_classic2()+ylab("P accumulation")+xlab("log(relative abundance + 1)")+ylim(0,1270)
p1


#Statistics Spearman correlations

#OTU29
#Alpine
NM_OTU29_Alpine<-filter(NM,OTU =="Otu000029" & Type_site == "Alpine")
cor.test(NM_OTU29_Alpine$Abundance,NM_OTU29_Alpine$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho

# data:  NM_OTU29_Alpine$Abundance and NM_OTU29_Alpine$Pplant.soil
# S = 33780, p-value = 5.482e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4568706 

#Brassicaceae only
NM_OTU29_Alpine_Brass<-filter(NM,OTU =="Otu000029" & Type_site == "Alpine" & sample_Family == "Brassicaceae")
cor.test(NM_OTU29_Alpine_Brass$Abundance,NM_OTU29_Alpine_Brass$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU29_Alpine_Brass$Abundance and NM_OTU29_Alpine_Brass$Pplant.soil
# S = 554.39, p-value = 1.714e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7589592 

#Cyperaceae only
NM_OTU29_Alpine_Cyp<-filter(NM,OTU =="Otu000029" & Type_site == "Alpine" & sample_Family == "Cyperaceae")
cor.test(NM_OTU29_Alpine_Cyp$Abundance,NM_OTU29_Alpine_Cyp$Pplant.soil,method="spearman",adjust="none")
plot(rank(NM_OTU29_Alpine_Cyp$Abundance),rank(NM_OTU29_Alpine_Cyp$Pplant.soil))

#Spearman's rank correlation rho

# data:  NM_OTU29_Alpine_Cyp$Abundance and NM_OTU29_Alpine_Cyp$Pplant.soil
# S = 1045.6, p-value = 0.005845
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5453856 


#Caryophyllaceae only
NM_OTU29_Alpine_Car<-filter(NM,OTU =="Otu000029" & Type_site == "Alpine" & sample_Family == "Caryophyllaceae")
cor.test(NM_OTU29_Alpine_Car$Abundance,NM_OTU29_Alpine_Car$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU29_Alpine_Car$Abundance and NM_OTU29_Alpine_Car$Pplant.soil
# S = 2180.1, p-value = 0.8089
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.05212058 



#OTU29
#Lowland (Plaine)
NM_OTU29_Lowland<-filter(NM,OTU =="Otu000029" & Type_site == "Plaine")
cor.test(NM_OTU29_Lowland$Abundance,NM_OTU29_Lowland$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU29_Lowland$Abundance and NM_OTU29_Lowland$Pplant.soil
# S = 13455, p-value = 5.831e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5146207 

#Brassicaceae only
NM_OTU29_Lowland_Bra<-filter(NM,OTU =="Otu000029" & Type_site == "Plaine" & sample_Family == "Brassicaceae")
cor.test(NM_OTU29_Lowland_Bra$Abundance,NM_OTU29_Lowland_Bra$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU29_Lowland_Bra$Abundance and NM_OTU29_Lowland_Bra$Pplant.soil
# S = 457.74, p-value = 0.02443
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5276201 

#Cyperaceae only
NM_OTU29_Lowland_Cyp<-filter(NM,OTU =="Otu000029" & Type_site == "Plaine" & sample_Family == "Cyperaceae")
cor.test(NM_OTU29_Lowland_Cyp$Abundance,NM_OTU29_Lowland_Cyp$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU29_Lowland_Cyp$Abundance and NM_OTU29_Lowland_Cyp$Pplant.soil
# S = 308, p-value = 0.000575
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7298246 

#Caryophyllaceae only
NM_OTU29_Lowland_Car<-filter(NM,OTU =="Otu000029" & Type_site == "Plaine" & sample_Family == "Caryophyllaceae")
cor.test(NM_OTU29_Lowland_Car$Abundance,NM_OTU29_Lowland_Car$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU29_Lowland_Car$Abundance and NM_OTU29_Lowland_Car$Pplant.soil
# S = 489.01, p-value = 0.03659
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4953475 



#OTU7
#Alpine
NM_OTU7_Alpine<-filter(NM,OTU =="Otu000007" & Type_site == "Alpine")
cor.test(NM_OTU7_Alpine$Abundance,NM_OTU7_Alpine$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Alpine$Abundance and NM_OTU7_Alpine$Pplant.soil
# S = 24261, p-value = 1.292e-08
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6099208 

#Brassicaceae only
NM_OTU7_Alpine_Bra<-filter(NM,OTU =="Otu000007" & Type_site == "Alpine" & sample_Family == "Brassicaceae")
cor.test(NM_OTU7_Alpine_Bra$Abundance,NM_OTU7_Alpine_Bra$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Alpine_Bra$Abundance and NM_OTU7_Alpine_Bra$Pplant.soil
# S = 716.13, p-value = 0.0001988
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6886391 

#Cyperaceae only
NM_OTU7_Alpine_Cyp<-filter(NM,OTU =="Otu000007" & Type_site == "Alpine" & sample_Family == "Cyperaceae")
cor.test(NM_OTU7_Alpine_Cyp$Abundance,NM_OTU7_Alpine_Cyp$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Alpine_Cyp$Abundance and NM_OTU7_Alpine_Cyp$Pplant.soil
# S = 765.57, p-value = 0.0003696
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6671432 

#Caryophyllaceae only
NM_OTU7_Alpine_Car<-filter(NM,OTU =="Otu000007" & Type_site == "Alpine" & sample_Family == "Caryophyllaceae")
cor.test(NM_OTU7_Alpine_Car$Abundance,NM_OTU7_Alpine_Car$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Alpine_Car$Abundance and NM_OTU7_Alpine_Car$Pplant.soil
# S = 792.07, p-value = 0.0005053
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6556221 


#OTU7
#Lowland (PLaine)
NM_OTU7_Lowland<-filter(NM,OTU =="Otu000007" & Type_site == "Plaine")
cor.test(NM_OTU7_Lowland$Abundance,NM_OTU7_Lowland$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Lowland$Abundance and NM_OTU7_Lowland$Pplant.soil
# S = 27338, p-value = 0.9205
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.01377536 

#Brassicaceae only
NM_OTU7_Lowland_Bra<-filter(NM,OTU =="Otu000007" & Type_site == "Plaine" & sample_Family == "Brassicaceae")
cor.test(NM_OTU7_Lowland_Bra$Abundance,NM_OTU7_Lowland_Bra$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Lowland_Bra$Abundance and NM_OTU7_Lowland_Bra$Pplant.soil
# S = 1087.6, p-value = 0.6286
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1223746 

#Cyperaceae only
NM_OTU7_Lowland_Cyp<-filter(NM,OTU =="Otu000007" & Type_site == "Plaine" & sample_Family == "Cyperaceae")
cor.test(NM_OTU7_Lowland_Cyp$Abundance,NM_OTU7_Lowland_Cyp$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Lowland_Cyp$Abundance and NM_OTU7_Lowland_Cyp$Pplant.soil
# S = 1293.8, p-value = 0.5818
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1349239 

#Caryophyllaceae only
NM_OTU7_Lowland_Car<-filter(NM,OTU =="Otu000007" & Type_site == "Plaine" & sample_Family == "Caryophyllaceae")
cor.test(NM_OTU7_Lowland_Car$Abundance,NM_OTU7_Lowland_Car$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_OTU7_Lowland_Car$Abundance and NM_OTU7_Lowland_Car$Pplant.soil
# S = 912.71, p-value = 0.8189
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.05809522


#####Graphics Glomeromycota in NM plants ####

## Get Glomeromycota relative abundances
phylum<-tax_glom(
  phy_obj_plant,
  taxrank = rank_names(phy_obj_plant)[2])

phylum_ra<-transform_sample_counts(phylum, function(x) x/sum(x) )
phylum_log<-otu_table(phylum_ra)+1
phylum_log<-log10(otu_table(phylum_log))
glom_log <- phyloseq(otu_table(phylum_log), tax_table(phylum), samples)

glomero<-subset_taxa(glom_log,Phylum=="p__Glomeromycota")
glomero_abund<- psmelt(glomero)

## Make Graph
glomero_abund<-filter(glomero_abund,Status=="Mycorhizal")

p2<-ggplot(glomero_abund, aes(x = Abundance, y =Pplant.soil)) +
  geom_point(aes(color = Type_site,shape=sample_Family),size=3) + xscale("log10",.format=TRUE)+ 
  scale_color_manual(values=c("#004D40","#FFC107"))+
  scale_shape_manual(values=c(1,0,2,7))+
  geom_smooth(aes(color=Type_site),method="lm",se=TRUE,level=0.95) + 
  theme_classic2()+ylab("P accumulation")+xlab("log(relative abundance + 1)")+ylim(0,1270)

p2



#Statistics Spearman correlations

#Glomeromycota
#Alpine
Myc_Glom_Alpine<-filter(glomero_abund, Type_site == "Alpine")
cor.test(Myc_Glom_Alpine$Abundance,Myc_Glom_Alpine$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Alpine$Abundance and Myc_Glom_Alpine$Pplant.soil
# S = 18617, p-value = 7.292e-12
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7006793 

#only Asteraceae
Myc_Glom_Alpine_Ast<-filter(glomero_abund, Type_site == "Alpine" & sample_Family=="Asteraceae")
cor.test(Myc_Glom_Alpine_Ast$Abundance,Myc_Glom_Alpine_Ast$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Alpine_Ast$Abundance and Myc_Glom_Alpine_Ast$Pplant.soil
# S = 632.55, p-value = 6.138e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7249785 


#only Poaceae
Myc_Glom_Alpine_Poa<-filter(glomero_abund, Type_site == "Alpine" & sample_Family=="Poaceae")
cor.test(Myc_Glom_Alpine_Poa$Abundance,Myc_Glom_Alpine_Poa$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Alpine_Poa$Abundance and Myc_Glom_Alpine_Poa$Pplant.soil
# S = 449.6, p-value = 2.152e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8045227 

#only Geraniaceae or Ranunculaceae
Myc_Glom_Alpine_Ger<-filter(glomero_abund, Type_site == "Alpine" & (sample_Family=="Geraniaceae"| sample_Family=="Ranunculaceae"))
cor.test(Myc_Glom_Alpine_Ger$Abundance,Myc_Glom_Alpine_Ger$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Alpine_Ger$Abundance and Myc_Glom_Alpine_Ger$Pplant.soil
# S = 832, p-value = 0.001031
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6382609 

#Glomeromycota
#Lowland
Myc_Glom_Lowland<-filter(glomero_abund, Type_site == "Plaine")
cor.test(Myc_Glom_Lowland$Abundance,Myc_Glom_Lowland$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Lowland$Abundance and Myc_Glom_Lowland$Pplant.soil
# S = 17286, p-value = 0.01159
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3411157 

#only Asteraceae
Myc_Glom_Lowland_Ast<-filter(glomero_abund, Type_site == "Plaine" & sample_Family=="Asteraceae")
cor.test(Myc_Glom_Lowland_Ast$Abundance,Myc_Glom_Lowland_Ast$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Lowland_Ast$Abundance and Myc_Glom_Lowland_Ast$Pplant.soil
# S = 676, p-value = 0.222
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3023736 

#only Poaceae
Myc_Glom_Lowland_Poa<-filter(glomero_abund, Type_site == "Plaine" & sample_Family=="Poaceae")
cor.test(Myc_Glom_Lowland_Poa$Abundance,Myc_Glom_Lowland_Poa$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Lowland_Poa$Abundance and Myc_Glom_Lowland_Poa$Pplant.soil
# S = 1028, p-value = 0.8114
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.06088751 

#only Geraniaceae or Ranunculaceae
Myc_Glom_Lowland_Ger<-filter(glomero_abund, Type_site == "Plaine" & (sample_Family=="Geraniaceae"| sample_Family=="Ranunculaceae"))
cor.test(Myc_Glom_Lowland_Ger$Abundance,Myc_Glom_Lowland_Ger$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  Myc_Glom_Lowland_Ger$Abundance and Myc_Glom_Lowland_Ger$Pplant.soil
# S = 1030, p-value = 0.805
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.0629515 

plot_grid(p2,p1,rel_widths = c(0.40,0.60))



####Figure S11####


#####Graphics OTU7 and OTU29 in AM-Mycorrhizal plants####

y5<-filter(y4,OTU=="Otu000029"| OTU=="Otu000007")

myc<-filter(y5,Status=="Mycorhizal")

p4<-ggplot(myc, aes(x = Abundance, y =Pplant.soil)) +
  geom_point(aes(color = Type_site,shape=sample_Family),size=3) + xscale("log10",.format=TRUE)+ 
  scale_color_manual(values=c("#004D40","#FFC107"))+
  scale_shape_manual(values=c(1,0,2,7))+
  geom_smooth(aes(color=Type_site),method="lm",se=TRUE,level=0.95) + facet_wrap(~OTU,nrow=1,axes = "all",axis.labels = "all")+
  theme_classic2()+ylab("P accumulation")+xlab("log(relative abundance + 1)")+ylim(0,1270)

p4

#Statistics Spearman correlations

#OTU29
#Alpine
myc_OTU29_Alpine<-filter(myc,OTU =="Otu000029" & Type_site == "Alpine")
cor.test(myc_OTU29_Alpine$Abundance,myc_OTU29_Alpine$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Alpine$Abundance and myc_OTU29_Alpine$Pplant.soil
# S = 15010, p-value = 1.164e-14
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7586741 

#only Asteraceae
myc_OTU29_Alpine_Ast<-filter(myc,OTU =="Otu000029" & Type_site == "Alpine" & sample_Family == "Asteraceae")
cor.test(myc_OTU29_Alpine_Ast$Abundance,myc_OTU29_Alpine_Ast$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Alpine_Ast$Abundance and myc_OTU29_Alpine_Ast$Pplant.soil
# S = 406.33, p-value = 7.766e-07
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8233328 

#only Poaceae
myc_OTU29_Alpine_Poa<-filter(myc,OTU =="Otu000029" & Type_site == "Alpine" & sample_Family == "Poaceae")
cor.test(myc_OTU29_Alpine_Poa$Abundance,myc_OTU29_Alpine_Poa$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Alpine_Poa$Abundance and myc_OTU29_Alpine_Poa$Pplant.soil
# S = 680.28, p-value = 0.0001226
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7042264

#only Geraniaceae or Ranunculaceae
myc_OTU29_Alpine_Ger<-filter(myc,OTU =="Otu000029" & Type_site == "Alpine" & (sample_Family == "Geraniaceae"|sample_Family == "Ranunculaceae"))
cor.test(myc_OTU29_Alpine_Ger$Abundance,myc_OTU29_Alpine_Ger$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Alpine_Ger$Abundance and myc_OTU29_Alpine_Ger$Pplant.soil
# S = 711.08, p-value = 0.000186
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6908363


#OTU29
#Lowland
myc_OTU29_Lowland<-filter(myc,OTU =="Otu000029" & Type_site == "Plaine")
cor.test(myc_OTU29_Lowland$Abundance,myc_OTU29_Lowland$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Lowland$Abundance and myc_OTU29_Lowland$Pplant.soil
# S = 18253, p-value = 0.02531
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3042409

#only Asteraceae
myc_OTU29_Lowland_Ast<-filter(myc,OTU =="Otu000029" & Type_site == "Plaine" & sample_Family == "Asteraceae")
cor.test(myc_OTU29_Lowland_Ast$Abundance,myc_OTU29_Lowland_Ast$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Lowland_Ast$Abundance and myc_OTU29_Lowland_Ast$Pplant.soil
# S = 579.8, p-value = 0.0985
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4016521 

#only Poaceae
myc_OTU29_Lowland_Poa<-filter(myc,OTU =="Otu000029" & Type_site == "Plaine" & sample_Family == "Poaceae")
cor.test(myc_OTU29_Lowland_Poa$Abundance,myc_OTU29_Lowland_Poa$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Lowland_Poa$Abundance and myc_OTU29_Lowland_Poa$Pplant.soil
# S = 390, p-value = 0.01019
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5975232 

#only Geraniaceae or Ranunculaceae
myc_OTU29_Lowland_Ger<-filter(myc,OTU =="Otu000029" & Type_site == "Plaine" & (sample_Family == "Geraniaceae" |sample_Family == "Ranunculaceae"))
cor.test(myc_OTU29_Lowland_Ger$Abundance,myc_OTU29_Lowland_Ger$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU29_Lowland_Ger$Abundance and myc_OTU29_Lowland_Ger$Pplant.soil
# S = 1020, p-value = 0.837
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.05263158


#OTU7
#Alpine
myc_OTU7_Alpine<-filter(myc,OTU =="Otu000007" & Type_site == "Alpine")
cor.test(myc_OTU7_Alpine$Abundance,myc_OTU7_Alpine$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Alpine$Abundance and myc_OTU7_Alpine$Pplant.soil
# S = 17591, p-value = 1.375e-12
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7171745 

#only Asteraceae
myc_OTU7_Alpine_Ast<-filter(myc,OTU =="Otu000007" & Type_site == "Alpine" & sample_Family == "Asteraceae")
cor.test(myc_OTU7_Alpine_Ast$Abundance,myc_OTU7_Alpine_Ast$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Alpine_Ast$Abundance and myc_OTU7_Alpine_Ast$Pplant.soil
# S = 599.82, p-value = 3.683e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.7392069 

#only Poaceae
myc_OTU7_Alpine_Poa<-filter(myc,OTU =="Otu000007" & Type_site == "Alpine" & sample_Family == "Poaceae")
cor.test(myc_OTU7_Alpine_Poa$Abundance,myc_OTU7_Alpine_Poa$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Alpine_Poa$Abundance and myc_OTU7_Alpine_Poa$Pplant.soil
# S = 428.39, p-value = 1.325e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8137424 

#only Geraniaceae or Ranunculaceae
myc_OTU7_Alpine_Ger<-filter(myc,OTU =="Otu000007" & Type_site == "Alpine" & (sample_Family == "Geraniaceae"|sample_Family == "Ranunculaceae"))
cor.test(myc_OTU7_Alpine_Ger$Abundance,myc_OTU7_Alpine_Ger$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Alpine_Ger$Abundance and myc_OTU7_Alpine_Ger$Pplant.soil
# S = 816.79, p-value = 0.0006688
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6448742 

#OTU7
#Lowland
myc_OTU7_Lowland<-filter(myc,OTU =="Otu000007" & Type_site == "Plaine")
cor.test(myc_OTU7_Lowland$Abundance,myc_OTU7_Lowland$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Lowland$Abundance and myc_OTU7_Lowland$Pplant.soil
# S = 22280, p-value = 0.2765
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1507639 

#only Asteraceae
myc_OTU7_Lowland_Ast<-filter(myc,OTU =="Otu000007" & Type_site == "Plaine" & sample_Family == "Asteraceae")
cor.test(myc_OTU7_Lowland_Ast$Abundance,myc_OTU7_Lowland_Ast$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Lowland_Ast$Abundance and myc_OTU7_Lowland_Ast$Pplant.soil
# S = 629.86, p-value = 0.1545
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3499883

#only Poaceae
myc_OTU7_Lowland_Poa<-filter(myc,OTU =="Otu000007" & Type_site == "Plaine" & sample_Family == "Poaceae")
cor.test(myc_OTU7_Lowland_Poa$Abundance,myc_OTU7_Lowland_Poa$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Lowland_Poa$Abundance and myc_OTU7_Lowland_Poa$Pplant.soil
# S = 595.66, p-value = 0.1143
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3852853

#only Geraniaceae or Ranunculaceae
myc_OTU7_Lowland_Ger<-filter(myc,OTU =="Otu000007" & Type_site == "Plaine" & (sample_Family == "Geraniaceae"|sample_Family == "Ranunculaceae"))
cor.test(myc_OTU7_Lowland_Ger$Abundance,myc_OTU7_Lowland_Ger$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  myc_OTU7_Lowland_Ger$Abundance and myc_OTU7_Lowland_Ger$Pplant.soil
# S = 957.61, p-value = 0.9631
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.01175447



#####Graphics Glomeromycota in Non-AM plants####

glomero_abund<- psmelt(glomero)
glomero_abund_NM<-filter(glomero_abund,Status=="Non-mycorhizal")

p3<-ggplot(glomero_abund_NM, aes(x = Abundance, y =Pplant.soil)) +
  geom_point(aes(color = Type_site,shape=sample_Family),size=3) + xscale("log10",.format=TRUE)+ 
  scale_color_manual(values=c("#004D40","#FFC107"))+
  scale_shape_manual(values=c(15,19,17))+
  geom_smooth(aes(color=Type_site),method="lm",se=TRUE,level=0.95) + 
  theme_classic2()+ylab("P accumulation")+xlab("log(relative abundance + 1)")+ylim(0,1270)+
  theme(legend.position = "right")
p3

#Statistics Spearman correlations

#Glomeromycota
#Alpine
NM_Glom_Alpine<-filter(glomero_abund_NM, Type_site == "Alpine")
cor.test(NM_Glom_Alpine$Abundance,NM_Glom_Alpine$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Alpine$Abundance and NM_Glom_Alpine$Pplant.soil
# S = 50287, p-value = 0.1071
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#     rho 
# 0.19148

#only Brassicaeae
NM_Glom_Alpine_Bra<-filter(glomero_abund_NM, Type_site == "Alpine" & sample_Family == "Brassicaceae")
cor.test(NM_Glom_Alpine_Bra$Abundance,NM_Glom_Alpine_Bra$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Alpine_Bra$Abundance and NM_Glom_Alpine_Bra$Pplant.soil
# S = 1570, p-value = 0.1308
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3173913 

#only Cyperaceae
NM_Glom_Alpine_Cyp<-filter(glomero_abund_NM, Type_site == "Alpine" & sample_Family == "Cyperaceae")
cor.test(NM_Glom_Alpine_Cyp$Abundance,NM_Glom_Alpine_Cyp$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Alpine_Cyp$Abundance and NM_Glom_Alpine_Cyp$Pplant.soil
# S = 1923.9, p-value = 0.4452
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1635138 

#only Caryophyllaceae
NM_Glom_Alpine_Car<-filter(glomero_abund_NM, Type_site == "Alpine" & sample_Family == "Caryophyllaceae")
cor.test(NM_Glom_Alpine_Car$Abundance,NM_Glom_Alpine_Car$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Alpine_Car$Abundance and NM_Glom_Alpine_Car$Pplant.soil
# S = 1756, p-value = 0.2646
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2365217 


#Glomeromycota
#Lowland
NM_Glom_Lowland<-filter(glomero_abund_NM, Type_site == "Plaine")
cor.test(NM_Glom_Lowland$Abundance,NM_Glom_Lowland$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Lowland$Abundance and NM_Glom_Lowland$Pplant.soil
# S = 33693, p-value = 0.1141
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2154801 

#only Brassicaeae
NM_Glom_Lowland_Bra<-filter(glomero_abund_NM, Type_site == "Plaine" & sample_Family == "Brassicaceae")
cor.test(NM_Glom_Lowland_Bra$Abundance,NM_Glom_Lowland_Bra$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Lowland_Bra$Abundance and NM_Glom_Lowland_Bra$Pplant.soil
# S = 1093.1, p-value = 0.6127
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1280331 

#only Cyperaceae
NM_Glom_Lowland_Cyp<-filter(glomero_abund_NM, Type_site == "Plaine" & sample_Family == "Cyperaceae")
cor.test(NM_Glom_Lowland_Cyp$Abundance,NM_Glom_Lowland_Cyp$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Lowland_Cyp$Abundance and NM_Glom_Lowland_Cyp$Pplant.soil
# S = 1272, p-value = 0.6361
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1157895 

#only Caryophyllaceae
NM_Glom_Lowland_Car<-filter(glomero_abund_NM, Type_site == "Plaine" & sample_Family == "Caryophyllaceae")
cor.test(NM_Glom_Lowland_Car$Abundance,NM_Glom_Lowland_Car$Pplant.soil,method="spearman",adjust="none")

# Spearman's rank correlation rho
# 
# data:  NM_Glom_Lowland_Car$Abundance and NM_Glom_Lowland_Car$Pplant.soil
# S = 1184, p-value = 0.3747
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2218782 

plot_grid(p3,p4,rel_widths = c(0.40,0.60))
