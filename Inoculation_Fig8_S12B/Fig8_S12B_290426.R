###Code for Figure 8 and S12B
###Fungal inoculation on non-mycorrhizal plants
###Results from two independent experiments

###libraries needed####
library(biostat)
library(readxl)
library(cowplot)

##set wd
#setwd("~/Documents/Labo/MS/4_Pauline_xp/V4_resub_2025/Git_new/Inoculation_Fig8")
setwd("C:/Users/pbruyant/Documents/Article thèse/Nvelle version/New_code_Fig8/")

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
  
##############################
  ####Add pot effect :####
  
  ###test for Arabis####
  Arabis$effect_exp_pot<-paste(Arabis$Experiment,Arabis$Inoc_Pot,sep="_")
  
  #Test normality of shoot P data
  shapiro.test(Arabis$Fresh_biomass_mg)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(Fresh_biomass_mg ~ Inoculation*effect_exp_pot, data=Arabis)
  #Ok
  
  #Normalize data for ANOVA using box-cox transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Arabis$Fresh_biomass_mg)
  obj_bestnormP
  Arabis$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data yeo-johnson
  
  #verify normality of transformed data
  shapiro.test(Arabis$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Arabis)
  #Not ok 0.04 so okayish
  
  ###### ANOVA with effect of exp_pot
  res.aov<-aov(bestnormP ~Inoculation*effect_exp_pot, data=Arabis)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value   Pr(>F)    
  #Inoculation      2  13.54   6.768   9.125 0.000229 ***
  #  effect_exp_pot  24  38.29   1.596   2.151 0.004479 ** 
  #  Residuals      100  74.17   0.742 
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  #Get R2 values
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  # [1] 0.1074336 0.3039172 0.5886492
  ####Petit effet significatif de l'experience X pot
  
  ###Effect of pot so either mean by pot or mixed model with a nested pot as random effect. Did both 
  
  ###Arabis with pot acounted for####
  
  ###for Biomass####
  ####Sum by pot and divide by number of plant per pot####
  Arabis$effect_exp_pot<-paste(Arabis$Experiment,Arabis$Inoculation,Arabis$Pot,sep=";")
  
  
  Arabis_mean<-Arabis %>% 
    group_by(effect_exp_pot) %>%
    summarize(
      across(c(Fresh_biomass_mg,PO4_µg_mgFW,PO4_µg_total),
             ~mean(.x,na.rm=TRUE), .names="mean_{.col}")) 
  Arabis_mean<-separate(Arabis_mean, effect_exp_pot, c("Exp","Inoculation","Pot"),sep=";")
  
  ###Biomass
  #Test normality of shoot P data
  shapiro.test(Arabis_mean$mean_Fresh_biomass_mg)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_Fresh_biomass_mg ~ Inoculation, data=Arabis_mean)
  #Ok
  
  # ANOVA
  res.aov<-aov(mean_Fresh_biomass_mg ~Inoculation, data=Arabis_mean)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value Pr(>F)  
  #Inoculation  2  808.9   404.5   3.461 0.0478 *
  #  Residuals   24 2804.7   116.9                 
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.2238594 0.7761406
  
  TukeyHSD(res.aov)
  
  #Fit: aov(formula = mean_Fresh_biomass_mg ~ Inoculation, data = Arabis_mean)
  
  #Inoculation
  #diff        lwr        upr     p adj
  #Isolate_1-40-Isolate_1-39   3.585635  -9.140485 16.3117552 0.7636840
  #NI-Isolate_1-39            -9.395569 -22.121689  3.3305515 0.1771097
  #NI-Isolate_1-40           -12.981204 -25.707324 -0.2550834 0.0449693
  
  
  ###model with pot as random effect####
  library(lme4)
  library(lmerTest)
  ###use bestnorm data as not normal 
  
  Arabis$Inoculation <- relevel(as.factor(Arabis$Inoculation), ref = "NI")
  model <- lmer(bestnormP ~ Inoculation + (1|effect_exp_pot), data = Arabis)
  summary(model)
  
  #Fixed effects:
  #  Estimate Std. Error     df t value Pr(>|t|)    
  #(Intercept)              -0.4283     0.1973 24.3616  -2.171  0.03989 * 
  #InoculationIsolate_1-39   0.4538     0.2754 23.1138   1.648  0.11292   
  #InoculationIsolate_1-40   0.8139     0.2771 23.8144   2.937  0.00724 **
  
  
  ###P concentration#####
  ###mean####
  #Test normality of shoot P data
  shapiro.test(Arabis_mean$mean_PO4_µg_mgFW)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_PO4_µg_mgFW ~ Inoculation, data=Arabis_mean)
  #Ok
  
  #ANOVA 
  res.aov<-aov(mean_PO4_µg_mgFW ~Inoculation, data=Arabis_mean)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value Pr(>F)  
  #Inoculation  2 0.2072 0.10358   19.15 1.07e-05 ***
  #  Residuals   24 0.1298 0.00541                     
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.6147219 0.3852781
  
  TukeyHSD(res.aov)
  #Fit: aov(formula = mean_PO4_µg_mgFW ~ Inoculation, data = Arabis_mean)
  
  #$Inoculation
  #diff        lwr          upr     p adj
  #Isolate_1-40-Isolate_1-39 -0.12150650 -0.2080958 -0.034917236 0.0049956
  #NI-Isolate_1-39           -0.21390303 -0.3004923 -0.127313773 0.0000066
  #NI-Isolate_1-40           -0.09239654 -0.1789858 -0.005807277 0.0349162
  
  
  ###Mixed model#####
  
  shapiro.test(Arabis$PO4_µg_mgFW)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(PO4_µg_mgFW ~ Inoculation*effect_exp_pot, data=Arabis)
  #Ok
  
  #Normalize data for ANOVA using box-cox transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Arabis$PO4_µg_mgFW)
  obj_bestnormP
  Arabis$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data boxcox
  
  #verify normality of transformed data
  shapiro.test(Arabis$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Arabis)
  #Ok
  
  Arabis$Inoculation <- relevel(as.factor(Arabis$Inoculation), ref = "NI")
  model <- lmer(bestnormP ~ Inoculation + (1|effect_exp_pot), data = Arabis)
  summary(model)
  
  #Fixed effects:
  #  Estimate Std. Error       df t value Pr(>|t|)    
  #(Intercept)              -0.7673     0.1733 26.9014  -4.428 0.000143 ***
  #InoculationIsolate_1-39   1.4579     0.2422 25.7125   6.018 2.46e-06 ***
  #InoculationIsolate_1-40   0.8162     0.2436 26.3600   3.351 0.002442 **
  
  ###P total#####
  ###mean####
  #Test normality of shoot P data
  shapiro.test(Arabis_mean$mean_PO4_µg_total)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_PO4_µg_total ~ Inoculation, data=Arabis_mean)
  #Ok
  
  ###### ANOVA 
  res.aov<-aov(mean_PO4_µg_total ~Inoculation, data=Arabis_mean)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value Pr(>F)  
  #Inoculation  2  216.1  108.05   9.266 0.00104 **
  #Residuals   24  279.9   11.66                       
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.4357101 0.5642899
  
  TukeyHSD(res.aov)
  #Fit: aov(formula = mean_PO4_µg_total ~ Inoculation, data = Arabis_mean)
  
  #Inoculation
  #diff        lwr       upr     p adj
  #Isolate_1-40-Isolate_1-39 -0.8404183  -4.860504  3.179667 0.8613586
  #NI-Isolate_1-39           -6.3772862 -10.397372 -2.357201 0.0016231
  #NI-Isolate_1-40           -5.5368679  -9.556953 -1.516782 0.0058424
  
  ###Mixed-model#####
  library(lme4)
  library(lmerTest)
  
  shapiro.test(Arabis$PO4_µg_total)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(PO4_µg_total ~ Inoculation*effect_exp_pot, data=Arabis)
  #Ok
  
  #Normalize data for ANOVA using box-cox transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Arabis$PO4_µg_total)
  obj_bestnormP
  Arabis$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data boxcox
  
  #verify normality of transformed data
  shapiro.test(Arabis$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Arabis)
  #Ok
  
  
  Arabis$Inoculation <- relevel(as.factor(Arabis$Inoculation), ref = "NI")
  model <- lmer(bestnormP ~ Inoculation + (1|effect_exp_pot), data = Arabis)
  summary(model)
  #Fixed effects:
  #  Estimate Std. Error      df t value Pr(>|t|)    
  #(Intercept)              -0.7648     0.1933 25.6362  -3.957 0.000534 ***
  #  InoculationIsolate_1-39   1.1327     0.2708 24.6968   4.184 0.000315 ***
  #  InoculationIsolate_1-40   1.1345     0.2720 25.1851   4.172 0.000315 ***
  
  
  ###test for Minuartia#####
  Minuartia$effect_exp_pot<-paste(Minuartia$Experiment,Minuartia$Inoc_Pot,sep="_")
  
  #Test normality of shoot P data
  shapiro.test(Minuartia$Fresh_biomass_mg)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(Fresh_biomass_mg ~ Inoculation*effect_exp_pot, data=Minuartia)
  #Ok
  
  #Normalize data for ANOVA using asinh transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Minuartia$Fresh_biomass_mg)
  obj_bestnormP
  Minuartia$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data asinh
  
  #verify normality of transformed data
  shapiro.test(Minuartia$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Minuartia)
  #Ok
  
  ###### ANOVA
  res.aov<-aov(bestnormP ~Inoculation*effect_exp_pot, data=Minuartia)
  summary(res.aov)
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  # [1] 0.04162281 0.24827741 0.22527398 0.48482580
  
  
  ####Sum by pot and divide by number of plant per pot for Biomass
  Minuartia$effect_exp_pot<-paste(Minuartia$Experiment,Minuartia$Inoculation,Minuartia$Pot,sep=";")
  
  
  Minuartia_mean<-Minuartia %>% 
    group_by(effect_exp_pot) %>%
    summarize(
      across(c(Fresh_biomass_mg,PO4_µg_mgFW,PO4_µg_total),
             ~mean(.x,na.rm=TRUE), .names="mean_{.col}")) 
  Minuartia_mean<-separate(Minuartia_mean, effect_exp_pot, c("Exp","Inoculation","Pot"),sep=";")
  
  ###Biomass####
  ###Mean####
  #Test normality of shoot P data
  shapiro.test(Minuartia_mean$mean_Fresh_biomass_mg)
  #not ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_Fresh_biomass_mg ~ Inoculation, data=Minuartia_mean)
  #Ok
  
  #Normalize data for ANOVA using asinh
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Minuartia_mean$mean_Fresh_biomass_mg)
  obj_bestnormP
  Minuartia_mean$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data asinh
  
  #verify normality of transformed data
  shapiro.test(Minuartia_mean$bestnormP)
  #OK
  
  #test Homoscedasticity between groups
  car::leveneTest(bestnormP ~ Inoculation, data=Minuartia_mean)
  #Ok
  
  ###### ANOVA with bestnorm
  res.aov<-aov(bestnormP ~Inoculation, data=Minuartia_mean)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value Pr(>F)
  #Inoculation  2   2.79   1.396   1.428  0.253
  #Residuals   36  35.21   0.978  
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.07349711 0.92650289
  TukeyHSD(res.aov)
  
  #$Inoculation
  #diff        lwr       upr     p adj
  #Isolate_1-40-Isolate_1-39  0.2265899 -0.7215257 1.1747054 0.8294253
  #NI-Isolate_1-39           -0.4193863 -1.3675018 0.5287292 0.5316923
  #NI-Isolate_1-40           -0.6459762 -1.5940917 0.3021393 0.2322748
  
  
  
  ##Mixed model####
  Minuartia$Inoculation <- relevel(as.factor(Minuartia$Inoculation), ref = "NI")
  model <- lmer(bestnormP ~ Inoculation + (1|effect_exp_pot), data = Minuartia)
  summary(model)
  #Fixed effects:
  #  Estimate Std. Error      df t value Pr(>|t|)
  #(Intercept)              -0.2586     0.1970 35.6581  -1.313    0.198
  #InoculationIsolate_1-39   0.2866     0.2798 36.2649   1.024    0.312
  #InoculationIsolate_1-40   0.4558     0.2788 35.8039   1.634    0.111
  
  ###P concentration#####
  ###Mean####
  #Test normality of shoot P data
  shapiro.test(Minuartia_mean$mean_PO4_µg_mgFW)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_PO4_µg_mgFW ~ Inoculation, data=Minuartia_mean)
  #Ok
  
  ###### ANOVA 
  res.aov<-aov(mean_PO4_µg_mgFW ~Inoculation, data=Minuartia_mean)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value Pr(>F)  
  #Inoculation  2  1.171  0.5854   19.46 1.87e-06 ***
  #Residuals   36  1.083  0.0301                     
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.5194581 0.4805419
  TukeyHSD(res.aov)
  
  #$Inoculation
  #diff        lwr         upr     p adj
  #Isolate_1-40-Isolate_1-39 -0.09111231 -0.2574088  0.07518417 0.3832296
  #NI-Isolate_1-39           -0.40454015 -0.5708366 -0.23824367 0.0000024
  #NI-Isolate_1-40           -0.31342784 -0.4797243 -0.14713136 0.0001440
  
  ##Mixed model####
  #Test normality of shoot P data
  shapiro.test(Minuartia$PO4_µg_mgFW)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(PO4_µg_mgFW ~ Inoculation*effect_exp_pot, data=Minuartia)
  #Ok
  
  #Normalize data for ANOVA using yeo-johnson transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Minuartia$PO4_µg_mgFW)
  obj_bestnormP
  Minuartia$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data yeo-Johnson
  
  #verify normality of transformed data
  shapiro.test(Minuartia$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Minuartia)
  #Ok
  
  Minuartia$Inoculation <- relevel(as.factor(Minuartia$Inoculation), ref = "NI")
  model <- lmer(bestnormP~ Inoculation + (1|effect_exp_pot), data = Minuartia)
  summary(model)
  
  #Fixed effects:
  #  Estimate Std. Error      df t value Pr(>|t|)    
  #(Intercept)              -0.8488     0.1510 35.1870  -5.623 2.39e-06 ***
  #  InoculationIsolate_1-39   1.4557     0.2146 35.8757   6.783 6.45e-08 ***
  # InoculationIsolate_1-40   1.1584     0.2138 35.3531   5.419 4.37e-06 ***
  
  
  
  
  
  ###P total#####
  ###mean####
  #Test normality of shoot P data
  shapiro.test(Minuartia_mean$mean_PO4_µg_total)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_PO4_µg_total ~ Inoculation, data=Minuartia_mean)
  #Ok
  
  ###### ANOVA 
  res.aov<-aov(mean_PO4_µg_total ~Inoculation, data=Minuartia_mean)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value Pr(>F)  
  #            Df Sum Sq Mean Sq F value   Pr(>F)    
  #Inoculation  2  187.0   93.51   14.67 2.19e-05 ***
  #  Residuals   36  229.4    6.37                      
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #0.4490789 0.5509211
  
  TukeyHSD(res.aov)
  
  #Df Sum Sq Mean Sq F value   Pr(>F)    
  #                                 diff       lwr       upr     p adj
  #Isolate_1-40-Isolate_1-39 -0.06109351 -2.481460  2.359273 0.9979037
  #NI-Isolate_1-39           -4.67567382 -7.096041 -2.255307 0.0001019
  #NI-Isolate_1-40           -4.61458032 -7.034947 -2.194213 0.0001227
  
  ##Mixed model####
  
  #Test normality of shoot P data
  shapiro.test(Minuartia$PO4_µg_total)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(PO4_µg_total~ Inoculation, data=Minuartia)
  #pas ok
  
  #Normalize data for ANOVA using box transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Minuartia$PO4_µg_total)
  obj_bestnormP
  Minuartia$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data boxcox
  
  #verify normality of transformed data
  shapiro.test(Minuartia$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Minuartia)
  #Ok
  
  
  Minuartia$Inoculation <- relevel(as.factor(Minuartia$Inoculation), ref = "NI")
  model <- lmer(bestnormP~ Inoculation + (1|effect_exp_pot), data = Minuartia)
  summary(model)
  
  #Fixed effects:
  #  Estimate Std. Error      df t value Pr(>|t|)    
  #(Intercept)              -0.8149     0.1503 36.1179  -5.421 4.10e-06 ***
  #InoculationIsolate_1-39   1.2489     0.2138 36.8833   5.841 1.04e-06 ***
  #InoculationIsolate_1-40   1.2088     0.2129 36.3029   5.678 1.82e-06 ***
  
  
  ###test for Carex#####
  Carex$effect_exp_pot<-paste(Carex$Experiment,Carex$Inoc_Pot,sep="_")
  
  #Test normality of shoot P data
  shapiro.test(Carex$Fresh_biomass_mg)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(Fresh_biomass_mg ~ Inoculation*effect_exp_pot, data=Carex)
  #Ok
  
  #Normalize data for ANOVA using sqrt transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Carex$Fresh_biomass_mg)
  obj_bestnormP
  Carex$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data sqrt
  
  #verify normality of transformed data
  shapiro.test(Carex$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Carex)
  #Ok
  
  ###### ANOVA
  res.aov<-aov(bestnormP ~Inoculation*effect_exp_pot, data=Carex)
  summary(res.aov)
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  # [1] 0.04162281 0.24827741 0.22527398 0.48482580
  
  
  ####Sum by pot and divide by number of plant per pot for Biomass
  Carex$effect_exp_pot<-paste(Carex$Experiment,Carex$Inoculation,Carex$Pot,sep=";")
  
  
  Carex_mean<-Carex %>% 
    group_by(effect_exp_pot) %>%
    summarize(
      across(c(Fresh_biomass_mg,PO4_µg_mgFW,PO4_µg_total),
             ~mean(.x,na.rm=TRUE), .names="mean_{.col}")) 
  Carex_mean<-separate(Carex_mean, effect_exp_pot, c("Exp","Inoculation","Pot"),sep=";")
  
  ###Biomass####
  ###Mean####
  #Test normality of shoot P data
  shapiro.test(Carex_mean$mean_Fresh_biomass_mg)
  #not ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_Fresh_biomass_mg ~ Inoculation, data=Carex_mean)
  #Ok
  
  ###### ANOVA 
  res.aov<-aov(mean_Fresh_biomass_mg ~Inoculation, data=Carex_mean)
  summary(res.aov)
  
  #Df Sum Sq Mean Sq F value Pr(>F)
  #Inoculation  2   2.86   1.430   0.282  0.757
  #Residuals   22 111.62   5.073 
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.07349711 0.92650289
  TukeyHSD(res.aov)
  
  #$Inoculation
  #diff        lwr       upr     p adj
  #Isolate_1-40-Isolate_1-39 0.1172338 -2.632195 2.866663 0.9936965
  #NI-Isolate_1-39           0.7797917 -2.049348 3.608931 0.7703488
  #NI-Isolate_1-40           0.6625579 -2.086871 3.411987 0.8186965
  
  
  
  ##Mixed model####
  Carex$Inoculation <- relevel(as.factor(Carex$Inoculation), ref = "NI")
  model <- lmer(bestnormP~ Inoculation + (1|effect_exp_pot), data = Carex)
  summary(model)
  #Fixed effects:
  #  Estimate Std. Error      df t value Pr(>|t|)
  #(Intercept)              0.06988    0.20775 19.28770   0.336    0.740
  #InoculationIsolate_1-39 -0.13960    0.29625 19.61348  -0.471    0.643
  #InoculationIsolate_1-40 -0.09815    0.28302 18.68055  -0.347    0.733
  
  ###P concentration#####
  ###Mean####
  #Test normality of shoot P data
  shapiro.test(Carex_mean$mean_PO4_µg_mgFW)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_PO4_µg_mgFW ~ Inoculation, data=Carex_mean)
  #Ok
  
  ###### ANOVA 
  res.aov<-aov(mean_PO4_µg_mgFW ~Inoculation, data=Carex_mean)
  summary(res.aov)
  
  
  #            Df Sum Sq Mean Sq F value Pr(>F)
  #Inoculation  2 0.5234  0.2617   2.119  0.144
  #Residuals   22 2.7177  0.1235 
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.1614949 0.8385051
  TukeyHSD(res.aov)
  
  #$Inoculation
  #diff        lwr         upr     p adj
  #Isolate_1-40-Isolate_1-39  0.1938120 -0.2352099 0.62283385 0.5034736
  #NI-Isolate_1-39           -0.1561867 -0.5976465 0.28527320 0.6528787
  #NI-Isolate_1-40           -0.3499986 -0.7790205 0.07902322 0.1240793
  
  ##Mixed model####
  #Test normality of shoot P data
  shapiro.test(Carex$PO4_µg_mgFW)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(PO4_µg_mgFW ~ Inoculation*effect_exp_pot, data=Carex)
  #Ok
  
  
  Carex$Inoculation <- relevel(as.factor(Carex$Inoculation), ref = "NI")
  model <- lmer(PO4_µg_mgFW~ Inoculation + (1|effect_exp_pot), data = Carex)
  summary(model)
  
  #Fixed effects:
  #  Fixed effects:
  #Estimate Std. Error      df t value Pr(>|t|)    
  #(Intercept)               1.1049     0.1118 18.4331   9.881 8.55e-09 ***
  #  InoculationIsolate_1-39   0.1643     0.1593 18.7632   1.031   0.3155    
  #InoculationIsolate_1-40   0.3393     0.1525 17.9157   2.225   0.0392 *  
  
  
  
  ###P total#####
  ###mean####
  #Test normality of shoot P data
  shapiro.test(Carex_mean$mean_PO4_µg_total)
  #ok
  
  #test Homoscedasticity between groups
  car::leveneTest(mean_PO4_µg_total ~ Inoculation, data=Carex_mean)
  #Ok
  
  ###### ANOVA 
  res.aov<-aov(mean_PO4_µg_total ~Inoculation, data=Carex_mean)
  summary(res.aov)
  
  #            Df Sum Sq Mean Sq F value   Pr(>F)    
  #Inoculation  2    7.5   3.751   0.174  0.842
  #Residuals   22  474.9  21.586                      
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  df<-summary(res.aov)
  sum(df[[1]][,2])
  sqR<-df[[1]][,2]/sum(df[[1]][,2])
  sqR
  #[1] 0.01555358 0.98444642
  
  TukeyHSD(res.aov)
  
  #$Inoculation
  #diff       lwr      upr     p adj
  #Isolate_1-40-Isolate_1-39  1.3282208 -4.342930 6.999371 0.8277391
  #NI-Isolate_1-39            0.6276135 -5.207952 6.463179 0.9606330
  #NI-Isolate_1-40           -0.7006072 -6.371758 4.970543 0.9484192
  
  ##Mixed model####
  
  #Test normality of shoot P data
  shapiro.test(Carex$PO4_µg_total)
  #not normal
  
  #test Homoscedasticity between groups
  car::leveneTest(PO4_µg_total~ Inoculation, data=Carex)
  #pas ok
  
  #Normalize data for ANOVA using box transformation
  library(bestNormalize)
  obj_bestnormP<-bestNormalize(Carex$PO4_µg_total)
  obj_bestnormP
  Carex$bestnormP<-obj_bestnormP$x.t #contains the transformed (normalized) data yeo-Johnson
  
  #verify normality of transformed data
  shapiro.test(Carex$bestnormP)
  #OK
  
  #verify Homoscedasticity of transformed data
  car::leveneTest(bestnormP ~Inoculation*effect_exp_pot, data=Carex)
  #Ok
  
  
  Carex$Inoculation <- relevel(as.factor(Carex$Inoculation), ref = "NI")
  model <- lmer(bestnormP~ Inoculation + (1|effect_exp_pot), data = Carex)
  summary(model)
  
  #Fixed effects:
  #  Estimate Std. Error      df t value Pr(>|t|)    
  #(Intercept)             -0.12180    0.23279 19.01676  -0.523    0.607
  #InoculationIsolate_1-39  0.03402    0.33148 19.35604   0.103    0.919
  #InoculationIsolate_1-40  0.22992    0.31787 18.55168   0.723    0.479
  
  
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
