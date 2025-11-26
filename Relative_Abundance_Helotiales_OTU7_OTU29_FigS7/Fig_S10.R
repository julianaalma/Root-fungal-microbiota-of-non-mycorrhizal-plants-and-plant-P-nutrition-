###Fig S10
###Relative abundance Helotiales OTU 7 and 29
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

###phyloseq files (only plant data)####
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Relative_Abundance_Helotiales_OTU7_OTU29_FigS7")
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

####transform data in log(RA+1)##
##relative abundance for each OTU abundOTU/abundtotOTU
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
phy_obj_plant_ra<-otu_table(phy_obj_plant_ra)+1
phy_obj_plant_log<-log10(otu_table(phy_obj_plant_ra))
###reassemble otu_table with taxonomy and metadata
phy_obj_plant_log<- phyloseq(otu_table(phy_obj_plant_log), TAX, samples)


##extract abundance OTU7 and OTU 29 #####
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
y4<- psmelt(phy_obj_plant_ra)
y4$Abundance<-as.numeric(y4$Abundance)

y5<-filter(y4,OTU=="Otu000029"| OTU=="Otu000007")
y5$Site<-as.factor(y5$Site)
y5<-y5%>% rstatix::reorder_levels("Site",order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","La Doua","Commelle"))
y5$abbrev<-y5$sample_Family
y5$abbrev<-dplyr::recode_factor(y5$abbrev,Asteraceae="Ast",Ranunculaceae="Ran",
                                Geraniaceae="Ger",Poaceae="Poa", Caryophyllaceae= "Car",
                                Brassicaceae="Bra",Cyperaceae="Cyp")
#Find smallest value in abundance for each OTU
min_abundance_per_otu <- y5 %>%
  filter(Abundance > 0) %>%
  group_by(OTU) %>%
  summarise(Smallest_Non_Zero_Abundance = min(Abundance)) %>%
  mutate(Half_Smallest_Abundance = Smallest_Non_Zero_Abundance / 2)

# Replace 0 values in the original dataframe with half the smallest abundance for each OTU
half_smallest_otu29 <- min_abundance_per_otu %>%
  filter(OTU == "Otu000029") %>%
  pull(Half_Smallest_Abundance)

# Replace zeros with smallest value based on OTU
y5_replaced_zero <- y5 %>%
  mutate(Abundance = ifelse(OTU == "Otu000029" & Abundance == 0, 
                            half_smallest_otu29,
                            ifelse(Abundance == 0 & OTU != "Otu000029",
                                   min_abundance_per_otu$Half_Smallest_Abundance[match(OTU, min_abundance_per_otu$OTU)],
                                   Abundance)))


y5_replaced_zero$Abundance_log<-log10(y5_replaced_zero$Abundance)

p1<-ggplot(y5_replaced_zero,aes(x=Sample,y=Abundance_log,color=Status,fill=Status))+
  geom_point(alpha=0.5)+ geom_vline(xintercept=-6.0)+geom_segment(aes(y = case_when(
    OTU == "Otu000029" ~ -5.5,
    OTU == "Otu000007" ~ -5.3,
    TRUE ~ -5.4  # Default value if OTU doesn't match specific cases
  ), yend=Abundance_log), size=1,alpha=0.5)+
  facet_wrap(~OTU+Site+abbrev,nrow=2,scales="free_x")+
  theme_classic()+theme(axis.line.x = element_line())+
  theme(axis.text.y= element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text.x= element_blank(),axis.title.x = element_blank())+ 
  theme(strip.text.x = element_text(size=16),strip.background = element_blank())+ 
  theme(legend.position="bottom",legend.text = element_text(size=16),
        legend.title=element_text(size=18))+facet_wrap(~OTU+Site+abbrev,nrow=2,scales="free_x")+
  labs(y="log10 (Relative abundance)")
p1


