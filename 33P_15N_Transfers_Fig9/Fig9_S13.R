#code for Figure 9B and S13
#fungus to plant 33P and 15N transfer upon inoculation

###mind that in this code, poireau means "leek (A. ampeloprasum), Champignon means "fungi"
####import data + librairies######
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

#set wd
setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/33P_15N_Transfers_Fig9")

##### data #####
res_P<- read_excel("2024-PaulineB-Lyon_final.xlsx")

########## Fig. S12

### shoot dry weight (mg)#####

ShootDW <- res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=(`Pois Sec PA (mg)`), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Shoot biomass (mg DW)"))+xlab("Inoculation")+
  
  theme_minimal()+theme(axis.text.y = element_text(size=5),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=5),legend.text = element_text(size=5),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=5),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
ShootDW

####Statistics shoot DW ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$`Pois Sec PA (mg)`,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$`Pois Sec PA (mg)`,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$`Pois Sec PA (mg)`,Poireau$Champignon,p.adjust.method="none")

########### Phosphorus

#### Shoot P concentration (µg/mg DW) ##### 

Pplants_concentration<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=(`Quantité de P (mg/g PA)`), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Shoot P conc. (µg/mg DW)"))+xlab("Inoculation")+
  
  theme_minimal()+theme(axis.text.y = element_text(size=5),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=5),legend.text = element_text(size=5),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=5),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))

Pplants_concentration

####Statistics shoot P concentration ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$`Quantité de P (mg/g PA)`,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$`Quantité de P (mg/g PA)`,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$`Quantité de P (mg/g PA)`,Poireau$Champignon,p.adjust.method="none")


#### Total shoot P (mg)#####

ShootP<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=(`Quantité de P (mg P/PA)`), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Total shoot P (mg)"))+xlab("Inoculation")+
  
  theme_minimal()+theme(axis.text.y = element_text(size=5),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=5),legend.text = element_text(size=5),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=5),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
ShootP

####Statistics total shoot P ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$`Quantité de P (mg P/PA)`,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$`Quantité de P (mg P/PA)`,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$`Quantité de P (mg P/PA)`,Poireau$Champignon,p.adjust.method="none")


##33P per shoot sample ####
Shoot33P<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=(`33P (KBq dans PA)`), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("33P (KBq in shoots)"))+xlab("Inoculation")+
  
  theme_minimal()+theme(axis.text.y = element_text(size=5),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=5),legend.text = element_text(size=5),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=5),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
Shoot33P

####Statistics 33P per shoot sample ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$`33P (KBq dans PA)`,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$`33P (KBq dans PA)`,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$`33P (KBq dans PA)`,Poireau$Champignon,p.adjust.method="none")


########### Nitrogen

##N concentration µg/mg DW ####

res_P$Nconc<-res_P$`Quantité de N (mg N/PA)`/(res_P$`Pois Sec PA (mg)`/1000)

ShootN_conc<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=Nconc, x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Shoot N conc. (µg/mg DW)"))+xlab("Inoculation")+
  
  theme_minimal()+theme(axis.text.y = element_text(size=5),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=5),legend.text = element_text(size=5),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=5),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
ShootN_conc

####Statistics shoot N conc. ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$Nconc,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$Nconc,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$Nconc,Poireau$Champignon,p.adjust.method="none")


## Total shoot N (mg) ####

ShootN_tot<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=(`Quantité de N (mg N/PA)`), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("Total shoot N (mg)"))+xlab("Inoculation")+
  
  theme_minimal()+theme(axis.text.y = element_text(size=5),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=5),legend.text = element_text(size=5),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=5),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
ShootN_tot

####Statistics total shoot ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$`Quantité de N (mg N/PA)`,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$`Quantité de N (mg N/PA)`,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$`Quantité de N (mg N/PA)`,Poireau$Champignon,p.adjust.method="none")


## 15N per shoot sample (mg) ####

res_P$N15tot<-res_P$`Quantité de 15N (µg N/PA)`/1000

Shoot15N<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=N15tot, x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=1) + 
  ylab(expression("15N in shoots (mg)"))+xlab("Inoculation")+
  
  theme_minimal()+theme(axis.text.y = element_text(size=5),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=5),legend.text = element_text(size=5),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=5),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1,scales = "free_y")+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
Shoot15N

####Statistics shoot 15N ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$N15tot,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$N15tot,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$N15tot,Poireau$Champignon,p.adjust.method="none")

#assemble plot
plot_grid(ShootDW, Pplants_concentration, ShootP, Shoot33P, ShootN_conc, ShootN_tot, Shoot15N, ncol = 1, nrow = 7)



######################### Main Fig. 9B


##33P per g DW in shoots####
res_P$P33PA_per_PA<-(res_P$`33P (KBq dans PA)`/res_P$`Pois Sec PA (mg)`)*1000

res_P %>%
  group_by(Plante,Champignon) %>%
  summarise(.,mean(P33PA_per_PA))

Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")


ShootP_perPA<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=(P33PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("33P KBq/g shoot DW)"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1)+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
ShootP_perPA


####Statistics shoot 33P KBq/g shoot DW ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$P33PA_per_PA,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$P33PA_per_PA,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$P33PA_per_PA,Poireau$Champignon,p.adjust.method="none")


## adjust scale for each plant
ShootP_perPA_Ara<-Arabis %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  ggplot(.,aes(y=(P33PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("33P KBq/g shoot DW"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12, angle=30, hjust=1),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12,),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"), legend.position = "none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)),limits = c(0,6))+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  ggtitle("A. alpina")
  
ShootP_perPA_Ara


ShootP_perPA_Min<-Minuartia %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  ggplot(.,aes(y=(P33PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("33P KBq/g shoot DW"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12, angle=30, hjust=1),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12,),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"), legend.position = "none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)),limits = c(0,6))+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  ggtitle("M. verna")

ShootP_perPA_Min

ShootP_perPA_Poi<-Poireau %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  ggplot(.,aes(y=(P33PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("33P KBq/g shoot DW"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12, angle=30, hjust=1),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12,),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"), legend.position = "none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)),limits = c(0,16))+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  ggtitle("A. ampeloprasum")

ShootP_perPA_Poi


##### 15N

##15N per g DW in shoots (mg/g DW)####
res_P$N15PA_per_PA<-(res_P$`Quantité de 15N (µg N/PA)`/res_P$`Pois Sec PA (mg)`)

res_P %>%
  group_by(Plante,Champignon) %>%
  summarise(.,mean(N15PA_per_PA))

Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")


Shoot15N_perPA<-res_P %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  rstatix::reorder_levels(.,Plante,c("A. alpina","M verna","Poireau")) %>%
  ggplot(.,aes(y=(N15PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("15N mg/g shoot DW"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_blank(),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"))+
  facet_wrap(~Plante,nrow = 1)+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))
Shoot15N_perPA

####Statistics shoot 33P KBq/g shoot DW ####
Arabis<-filter(res_P, Plante=="A. alpina")
Minuartia<-filter(res_P, Plante=="M verna")
Poireau<-filter(res_P,Plante=="Poireau")

pairwise.wilcox.test(Arabis$N15PA_per_PA,Arabis$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Minuartia$N15PA_per_PA,Minuartia$Champignon,p.adjust.method="none")
pairwise.wilcox.test(Poireau$N15PA_per_PA,Poireau$Champignon,p.adjust.method="none")



## adjust scale for each plant
ShootN_perPA_Ara<-Arabis %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  ggplot(.,aes(y=(N15PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("15N mg/g shoot DW"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12, angle=30, hjust=1),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12,),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"), legend.position = "none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)),limits = c(0,15))+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  ggtitle("A. alpina")

ShootN_perPA_Ara


ShootN_perPA_Min<-Minuartia %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  ggplot(.,aes(y=(N15PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("15N mg/g shoot DW"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12, angle=30, hjust=1),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12,),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"), legend.position = "none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)),limits = c(0,15))+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  ggtitle("M. verna")

ShootN_perPA_Min

ShootN_perPA_Poi<-Poireau %>% 
  rstatix::reorder_levels(.,Champignon,c("Non-inoculé","1-39","1-40","Rhizophagus irregularis")) %>%
  ggplot(.,aes(y=(N15PA_per_PA), x= Champignon,color=Champignon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=Champignon),width = 0.2, alpha=0.3, size=2) + 
  ylab(expression("15N mg/g shoot DW"))+xlab("Inoculation")+
  #stat_summary(fun="mean", geom="point", shape=20, size=4, color="grey50", fill="grey50",data=res_P)+
  theme_minimal()+theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12, angle=30, hjust=1),
                        axis.title.x= element_blank(),axis.title.y= element_text(size=12,),legend.text = element_text(size=12),
                        legend.title = element_blank(),axis.ticks = element_line(),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size=12),
                        panel.border = element_rect(color = "black",fill = "NA"),axis.line = element_line(colour = "black"), legend.position = "none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)),limits = c(0,15))+
  scale_color_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  scale_fill_manual(values=c("grey", "#5fbcd3ff","#2c89a0ff","#C07F80"))+
  ggtitle("A. ampeloprasum")

ShootN_perPA_Poi


#Assemble plots
plot_grid(ShootP_perPA_Ara, ShootP_perPA_Min, ShootP_perPA_Poi, ShootN_perPA_Ara, ShootN_perPA_Min, ShootN_perPA_Poi, ncol = 3, nrow = 2)

