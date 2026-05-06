###Figs S5A 
###Diversity analysis 
###libraries needed####
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)        # filter and reformat data frames
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
setwd("C:/Users/pbruyant/Documents/Article thèse/Nvelle version/Git_new_clean-20251123T112115Z-1-001/Git_new_clean/Plant_shoot_P_Fig3")
###import data###
otu_mat<- read_excel("tableitsx_phyloseq.xlsx", sheet = "Table_OTU")
tax_mat<- read_excel("tableitsx_phyloseq.xlsx", sheet = "Taxonomy")
samples_df <- read_excel("tableitsx_phyloseq.xlsx", sheet = "Metadata")
###formate file###
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
phy_obj_plant<- prune_species(speciesSums(phy_obj_plant) > 0, phy_obj_plant) ### to keep OTUs only prensent in plants (not only in soils) 


####Rarefaction curves####
###Do not forget to turn the otu_mat from the phyloseq object !!! not the same orientation as everyone else of course
OTU_data<-as.data.frame(phy_obj_plant@otu_table)
samples_df<-as.data.frame(phy_obj_plant@sam_data)
OTU_data<-t(OTU_data)

###table for ggplot2
rarecurve_data <- vegan::rarecurve(as.matrix(OTU_data), 
                                   col =  as.factor(samples_df$Family),
                                   step = 100, label = TRUE,tidy=TRUE)

rarecurve_data$Sites <- sapply(rarecurve_data$Site, function(x) {
  if (grepl("Cha", x, ignore.case = TRUE)) {
    return("Chamrousse")
  } else if (grepl("Cla", x, ignore.case = TRUE)) {
    return("Clarée")
  } else if (grepl("Gal", x, ignore.case = TRUE)) {
    return("Galibier")
  } else if (grepl("Lau", x, ignore.case = TRUE)) {
    return("Lautaret")
  } else if (grepl("Per", x, ignore.case = TRUE)) {
    return("Perouges")
  } else if (grepl("Doua", x, ignore.case = TRUE)) {
    return("La Doua")
  } else if (grepl("Com", x, ignore.case = TRUE)) {
    return("Commelle")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

rarecurve_data$Families <- sapply(rarecurve_data$Site, function(x) {
  if (grepl("Car", x, ignore.case = TRUE)) {
    return("Caryophyllaceae")
  } else if (grepl("Cyp", x, ignore.case = TRUE)) {
    return("Cyperaceae")
  } else if (grepl("Bra", x, ignore.case = TRUE)) {
    return("Brassicaceae")
  } else if (grepl("Poa", x, ignore.case = TRUE)) {
    return("Poaceae")
  } else if (grepl("Ren", x, ignore.case = TRUE)) {
    return("Ranunculaceae")
  } else if (grepl("Ger", x, ignore.case = TRUE)) {
    return("Geraniaceae")
  } else if (grepl("Ast", x, ignore.case = TRUE)) {
    return("Asteraceae")
  } else {
    return("unknown")  # In case none of the above words are found
  }
})

rarecurve_data$Families<-as.factor(rarecurve_data$Families)
col_compartment<-c("Asteraceae"="#8c7209ff",
                   "Brassicaceae"="#dfafd6ff",
                   "Caryophyllaceae"="#c2f2efff",
                   "Cyperaceae"="#e3cec5ff",
                   "Geraniaceae"="#ad4498ff",
                   "Poaceae"="#3ac42dff",
                   "Ranunculaceae"="#492359ff")


rarecurve_data$Sites<-as.factor(rarecurve_data$Sites)
rarecurve_data<-rarecurve_data%>% rstatix::reorder_levels("Sites",order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","La Doua","Commelle"))

###plot rarecurve with ggplot 
rarecurve_ggplot <- ggplot(rarecurve_data, aes(x = Sample, y = Species, color = Families, group = Site)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~Sites, scales = "free", nrow = 2) +
  labs(title = "Rarefaction Curves") +
  scale_color_manual(values = col_compartment) +
  scale_x_continuous(labels = scales::label_number()) +
  theme_classic() +
  theme(
    axis.line.x = element_line(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 0),
    axis.title.x = element_blank(),
    strip.text.x = element_text(size = 16),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

rarecurve_ggplot

