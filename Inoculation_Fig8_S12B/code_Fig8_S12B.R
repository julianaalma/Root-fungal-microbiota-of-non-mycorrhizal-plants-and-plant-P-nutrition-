###Code for Figure 8 and S12B
###Fungal inoculation on non-mycorrhizal plants
###Results from two independent experiments

###libraries needed####
library(biostat)
library(readxl)
library(cowplot)

##set wd
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Inoculation_Fig8")

###get data####
exp1<-read_excel("Analyses_P_tourbe_exp1_270524.xlsx")
exp2<-read_excel("Analyses_P_Tourbe_2_220724.xlsx")
total<-rbind(exp1,exp2)
total$Experiment<-as.factor(total$Experiment)
total$Fresh_biomass_mg<-as.numeric(total$Fresh_biomass_mg)
mean(total$Dry_biomass_mg)
sd(total$Dry_biomass_mg)
library(dplyr)
library(ggplot2)
total<-total%>% rstatix::reorder_levels("Species",order=c("Arabis alpina","Minuartia verna","Carex sempervirens"))
total$PO4_µg_mgFW<-as.numeric(total$PO4_µg_mgFW)
total$PO4_µg_gFW<-total$PO4_µg_mgFW/1000

install.packages("ggbreak")
library(ggbreak) 
library(patchwork)


#### shoot fresh biomass (mg)#####

ShootFW<-total %>% 
  rstatix::reorder_levels(.,Inoculation,c("NI","Isolate_1-39","Isolate_1-40")) %>%
  rstatix::reorder_levels(.,Species,c("Arabis alpina","Minuartia verna","Carex sempervirens")) %>%
  ggplot(.,aes(y=(Fresh_biomass_mg), x= Inoculation,color=Inoculation))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Inoculation),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Shoot biomass (mg FW)"))+xlab("Inoculation")+
  stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=total)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Species,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
  scale_fill_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
ShootFW


#### Shoot P concentration (µg/mg FW) ##### 

Pplants_concentration<-total %>% 
  rstatix::reorder_levels(.,Inoculation,c("NI","Isolate_1-39","Isolate_1-40")) %>%
  rstatix::reorder_levels(.,Species,c("Arabis alpina","Minuartia verna","Carex sempervirens")) %>%
  ggplot(.,aes(y=(PO4_µg_mgFW), x= Inoculation,color=Inoculation))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Inoculation),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Shoot P conc. (µg/mg FW)"))+xlab("Inoculation")+
  stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=total)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Species,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
  scale_fill_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))

Pplants_concentration

#### Total shoot P (mg)#####

ShootP<-total %>% 
  rstatix::reorder_levels(.,Inoculation,c("NI","Isolate_1-39","Isolate_1-40")) %>%
  rstatix::reorder_levels(.,Species,c("Arabis alpina","Minuartia verna","Carex sempervirens")) %>%
  ggplot(.,aes(y=(PO4_µg_total), x= Inoculation,color=Inoculation))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Inoculation),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Total shoot P (µg)"))+xlab("Inoculation")+
  stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=total)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Species,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
  scale_fill_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.3)))
  ShootP
  
  #all plots together
  plot_grid(ShootFW, Pplants_concentration, ShootP, ncol = 1, nrow = 3)
  
  
  #statistics in comparison to NI per plant, per variable
  Arabis<-filter(total, Species=="Arabis alpina")
  Minuartia<-filter(total, Species=="Minuartia verna")
  Carex<-filter(total,Species=="Carex sempervirens")

  pairwise.wilcox.test(Arabis$Fresh_biomass_mg, Arabis$Inoculation, p.adjust.method="none")
  pairwise.wilcox.test(Arabis$PO4_µg_mgFW, Arabis$Inoculation, p.adjust.method="none")
  pairwise.wilcox.test(Arabis$PO4_µg_total, Arabis$Inoculation, p.adjust.method="none")
  
  pairwise.wilcox.test(Minuartia$Fresh_biomass_mg, Minuartia$Inoculation, p.adjust.method="none")
  pairwise.wilcox.test(Minuartia$PO4_µg_mgFW, Minuartia$Inoculation, p.adjust.method="none")
  pairwise.wilcox.test(Minuartia$PO4_µg_total, Minuartia$Inoculation, p.adjust.method="none")
  
  pairwise.wilcox.test(Carex$Fresh_biomass_mg, Carex$Inoculation, p.adjust.method="none")
  pairwise.wilcox.test(Carex$PO4_µg_mgFW, Carex$Inoculation, p.adjust.method="none")
  pairwise.wilcox.test(Carex$PO4_µg_total, Carex$Inoculation, p.adjust.method="none")
  
  
  
  ########################
  
  ###Fig S12B
  ##### Comparison of fungal root colonization levels among treatments
  
  colonisation<-read.table("Data_colonisation.txt", sep="\t", h=T)

  colonisation2<-colonisation %>% 
    rstatix::reorder_levels(.,Isolate,c("1_39","1_40")) %>%
    rstatix::reorder_levels(.,Plant,c("Arabis_alpina","Minuartia_verna","Carex_sempervirens")) %>%
    ggplot(.,aes(y=(Proportion_intersections_hyphae), x= Isolate, color=Isolate))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(aes(color=Isolate),width = 0.2, alpha=0.3, size=1) + 
    ylab(expression("Root colonisation"))+xlab("Inoculation")+
    stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=colonisation)+
    theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_blank(),
                          axis.title.x= element_blank(),axis.title.y= element_text(size=12),legend.text = element_text(size=12),
                          legend.title = element_blank(),axis.ticks = element_line(),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                          panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
    facet_wrap(~Plant,nrow = 1,scales = "free_y")+
    scale_color_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
    scale_fill_manual(values=c("#bebebeff", "#5fbcd3ff","#2c89a0ff"))+
    scale_y_continuous(limits = c(0.00,1.2))
  colonisation2
  
  
#Statistics
#non-normal distribution
shapiro.test(colonisation$Proportion_intersections_hyphae)
#pairwise wilcoxon test (non-parametric)
colonisation$treatment<-paste(colonisation$Isolate, colonisation$Plant)
test<-pairwise.wilcox.test(colonisation$Proportion_intersections_hyphae, colonisation$treatment, p.adjust.method="fdr")
make_cld(test)
