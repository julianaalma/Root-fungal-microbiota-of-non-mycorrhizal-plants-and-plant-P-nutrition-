###Figs S5: Fungal communities' alpha diversity and composition (order level)  
###Diversity analysis 
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
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Fungal_Microbiome_Fig4_FigS6_FigS8")
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

####transform data in log(RA+1)######
##relative abundance for each OTU abundOTU/abundtotOTU
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
phy_obj_plant_ra<-otu_table(phy_obj_plant_ra)+1
phy_obj_plant_log<-log10(otu_table(phy_obj_plant_ra))
###reassemble otu_table with taxonomy and metadata
phy_obj_plant_log<- phyloseq(otu_table(phy_obj_plant_log), TAX, samples)

#####alpha diversity analyses#####
alpha_index<-estimate_richness(phy_obj_plant,split=TRUE,measures = NULL)###provide raw data
alpha_index_df_sample<-cbind(alpha_index,samples_df)

alpha_index_df_sample<-rstatix::reorder_levels(alpha_index_df_sample,Family,c("Asteraceae","Geraniaceae","Ranunculaceae","Poaceae","Caryophyllaceae","Brassicaceae","Cyperaceae"))
alpha_index_df_sample<-alpha_index_df_sample%>% reorder_levels("Site",order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","La Doua","Commelle"))
alpha_index_df_sample$abbrev<- alpha_index_df_sample$Family
alpha_index_df_sample$abbrev<-dplyr::recode_factor(alpha_index_df_sample$abbrev,Asteraceae="Ast",Ranunculaceae="Ran",
                                                   Geraniaceae="Ger",Poaceae="Poa", Caryophyllaceae= "Car",
                                                   Brassicaceae="Bra",Cyperaceae="Cyp")

#Plot Shannon's H
p1<-ggboxplot(alpha_index_df_sample,x="abbrev",y="Shannon",add="jitter",
              color = "Status",fill="Status",size=0.5, width=0.3,alpha=0.5,
              xlab="Plant families",ylab="Shannon index")+ 
  theme_classic()+theme(axis.text.x = element_text(size=16,angle=90))+
  theme(axis.text.y = element_text(size=16))+facet_wrap(~Site,ncol=7,scale="free_x")+
  theme(strip.text.x = element_text(size=16),strip.background = element_blank())+
  theme(axis.text.y = element_text(size=16),axis.title.y = element_text(size=16),legend.position = "none")
p1


p2<-ggboxplot(alpha_index_df_sample,x="abbrev",y="Chao1",add="jitter",
              color = "Status",fill="Status",size=0.5, width=0.3,alpha=0.5,
              xlab="Plant families",ylab="Chao1 index")+ 
  theme_classic()+theme(axis.text.x = element_text(size=16,angle=90))+
  theme(axis.text.y = element_text(size=16))+facet_wrap(~Site,ncol=7,scale="free_x")+
  theme(strip.text.x = element_text(size=16),strip.background = element_blank())+
  theme(axis.text.y = element_text(size=16),axis.title.y = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))
p2

library(cowplot)
plot_alpha<-plot_grid(p1,p2,ncol=2,labels=c("A","B"),hjust=-0.5,align="h",rel_widths =c(0.9,1.1),rel_heights = c(0.25,0.75))
plot_alpha

#Statistics
####### Anova on alpha diversity#######

#Chao1
alpha_index_df_sample$VARIABLE<-alpha_index_df_sample$Chao1

#test Homoscedasticity between groups
leveneTest(VARIABLE~Site/Status/Family, data=alpha_index_df_sample)
#normality test (must be >0.05)
shapiro.test(alpha_index_df_sample$VARIABLE)
# -> Non-normal data

#Normalize data
library(bestNormalize)
obj_bestnorm<-bestNormalize(alpha_index_df_sample$VARIABLE)
obj_bestnorm
alpha_index_df_sample$bestnorm<-obj_bestnorm$x.t
res.aov <- aov(bestnorm ~ Site/Status/Family, data=alpha_index_df_sample)
summary(res.aov)

# 
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Site                 6  71.85  11.974  23.842  < 2e-16 ***
#   Site:Status          7  10.55   1.507   3.001  0.00504 ** 
#   Site:Status:Family  28  63.63   2.272   4.524 7.63e-11 ***
#   Residuals          211 105.98   0.502                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Shannon's H
alpha_index_df_sample$VARIABLE<-alpha_index_df_sample$Shannon

#test Homoscedasticity between groups
leveneTest(VARIABLE~Site/Status/Family, data=alpha_index_df_sample)
#normality test (must be >0.05)
shapiro.test(alpha_index_df_sample$VARIABLE)
# -> Non-normal data

#Normalize data
library(bestNormalize)
obj_bestnorm<-bestNormalize(alpha_index_df_sample$VARIABLE)
obj_bestnorm
alpha_index_df_sample$bestnorm<-obj_bestnorm$x.t
res.aov <- aov(bestnorm ~ Site/Status/Family, data=alpha_index_df_sample)
summary(res.aov)
# 
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Site                 6  59.38   9.896  17.669  < 2e-16 ***
#   Site:Status          7   6.19   0.885   1.579    0.143    
# Site:Status:Family  28  68.26   2.438   4.353 2.55e-10 ***
#   Residuals          211 118.17   0.560                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 




#### Beta diversity plotting #####
##### Stackplot, relative abundance taxa at "order" level ########
y <- tax_glom(phy_obj_plant, taxrank = "Order") # agglomerate taxa
y_ra<-transform_sample_counts(y, function(x) x/sum(x) )
y_ra <- phyloseq(otu_table(y_ra), tax_table(y), samples)
y4 <- psmelt(y_ra) # create dataframe from phyloseq object
y4$order <- as.character(y4$Order) #convert to character
y4$Abundance<-as.numeric(y4$Abundance)
y4$order[y4$Abundance < 0.05] <- "x_Other" #rename as "Other" orders with < 5% abundance

y4<-rstatix::reorder_levels(y4,sample_Family,order=c("Asteraceae","Geraniaceae","Ranunculaceae","Poaceae","Caryophyllaceae","Brassicaceae","Cyperaceae"))
y4<-rstatix::reorder_levels(y4,Site,order=c("Galibier","Lautaret","Chamrousse","Clarée","Perouges","La Doua","Commelle"))
y4$abbrev<-y4$sample_Family
y4$abbrev<-dplyr::recode_factor(y4$abbrev,Asteraceae="Ast",Ranunculaceae="Ran",
                                Geraniaceae="Ger",Poaceae="Poa", Caryophyllaceae= "Car",
                                Brassicaceae="Bra",Cyperaceae="Cyp")

palette <- c(
  "o__Helotiales"                         = "#00ffffff",
  "o__Mortierellales"                     = "#00A676",
  "o__Glomerellales"                      = "#008A8B",
  "o__Thelebolales"                       = "#005353",
  "o__Platygloeales"                      = "#0277bc",
  "o__Sebacinales"                        = "#1A6233",
  "o__Mytilinidales"                      = "#2C6B49",
  "o__Capnodiales"                        = "#013013",
  "o__Gomphales"                          = "#000E05",
  "o__Phomatosporales"                    = "#5A9A8F",
  "o__Chaetosphaeriales"                 = "#A9DAD1",
  "o__Cantharellales"                     = "#84DEBB",
  "o__Savoryellales"                      = "#86AEB0",
  "o__Chaetothyriales"                    = "#CAFF70",
  "o__Verrucariales"                      = "#e5ffa1",
  "o__Minutisphaerales"                   = "#8A6E21",
  "o__Russulales"                         = "#F5B700",
  "o__Amphisphaeriales"                   = "#FFE8A3",
  "o__Atheliales"                         = "#BFD200",
  "p__Ascomycota_unclassified"            = "#8c8888",
  "o__Archaeosporales"                    = "#006E90",
  "o__Sclerococcales"                     = "#457B9D",
  "c__Leotiomycetes_unclassified"         = "#08519C",
  "o__Boliniales"                         = "#1B4662",
  "o__Pezizales"                          = "#881AA9",
  "o__Xylariales"                         = "#bed8ff",
  "o__Mycosphaerellales"                  = "#4B66CC",
  "o__Trechisporales"                     = "#8D9EF2",
  "o__Auriculariales"                     = "#D7C2D7",
  "o__Pleosporales"                       = "#ccabff",
  "o__Sordariomycetes_ord_Incertae_sedis"= "#B95EFF",
  "o__Olpidiales"                         = "#DC5EFF",
  "o__Hypocreales"                        = "#9996C3",
  "o__Leotiales"                          = "#3C1642",
  "o__Rhytismatales"                      = "#2A1F00",
  "c__Dothideomycetes_unclassified"       = "#5C4033",
  "o__Corticiales"                        = "#6F1A07",
  "o__Pezizomycotina_ord_Incertae_sedis" = "#988200",
  "o__Branch06"                           = "#AA6F73",
  "o__Sakaguchiales"                      = "#C57B57",
  "o__Polyporales"                        = "#CDAA7D",
  "k__Fungi_unclassified"                 = "#AFABAB",
  "x_Other"                                 = "grey85",
  "o__Geoglossales"                       = "#74C475",
  "o__Myrmecridiales"                     = "#7da800",
  "o__Glomerales"                         = "#68165cff",
  "c__Agaricomycetes_unclassified"         = "#AD1800",
  "o__Microascales"                       = "#DB3C50",
  "o__Acarosporales"                      = "#CC5500",
  "o__Thelephorales"                      = "#ff6700",
  "o__Sordariales"                        = "#ff728c",
  "o__Tremellales"                        = "#ff83cf",
  "o__Agaricales"                         = "#FF56D6",
  "o__Tubeufiales"                        = "#860042",
  "o__Coniochaetales"                     = "#F3BEC5",
  "o__Helicobasidiales"                   = "#4C1F00",
  "o__Magnaporthales"                     = "#994454",
  "o__Orbiliales"                         = "#FF9CAF",
  "o__Amylocorticiales"                   = "#803400"
)







p <- ggplot(y4, aes(fill=order, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity",width = 0.8)+facet_wrap(~Site+y4$abbrev, nrow=1, scale="free_x")+
  theme_classic()
p<-p+theme(legend.title = element_text(size = 14), 
           legend.text = element_text(size = 10),panel.spacing = unit(0.05, "lines"),
           axis.text.x = element_text(angle=90),legend.key.size = unit(0.5, "lines"),
           legend.position="bottom",strip.text.x = element_text(size=10),
           strip.background = element_blank())+
  scale_fill_manual(values=palette) + 
  labs(y="Relative abundance of main fungal orders",fill="Fungal order")
p

