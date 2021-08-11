rm(list=ls())

library(survival)
library(survminer)
library(dplyr)
library(gridExtra)

library(moonBook)
library(readxl)
library(ggplot2)

PDAC_data_r <- read_excel("~/PDAC_data_r_revision.xlsx", 
                          col_types = c("numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric","numeric",
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "date", "date", "numeric", 
                                        "numeric", "date", "numeric", "text", 
                                        "numeric", "numeric","date", "numeric", "numeric", 
                                        "text"))

table1<-mytable(PDAC_data_r)

## factor conversion

#age
PDAC_data_r<-PDAC_data_r %>% mutate(age_group=ifelse(Age>=65,"Old","Young"))
PDAC_data_r$age_group<-factor(PDAC_data_r$age_group)

#sex
PDAC_data_r$Sex <-factor(PDAC_data_r$Sex,
                          levels = c("0","1"),
                          labels = c("Male", "Female"))
#cyst
PDAC_data_r$Intratumoral_cystic_type <-factor(PDAC_data_r$Intratumoral_cystic_type,
                                            levels = c("0","1","2"),
                                            labels = c("None", "Neoplastic mucin cyst", "Necrosis"))
#Diffusion
PDAC_data_r$Diffusion_restriction <-factor(PDAC_data_r$Diffusion_restriction,
                        levels = c("0","1"),
                        labels = c("Neg", "Pos"))
#CA19-9
PDAC_data_r$CA19_9_cut_off <-factor(PDAC_data_r$CA19_9_cut_off,
                                  levels = c("0","1"),
                                  labels = c("<34", ">=34"))

#Location
PDAC_data_r$Location <-factor(PDAC_data_r$Location,
                                    levels = c("0","1","2"),
                                    labels = c("Head", "Body","Tail"))

#T stage
PDAC_data_r$T_stage <-factor(PDAC_data_r$T_stage,
                              levels = c("0","1","2","3"),
                              labels = c("T1", "T2","T3","T4"))
#N stage
PDAC_data_r$N_stage <-factor(PDAC_data_r$N_stage,
                              levels = c("0","1"),
                              labels = c("N0", "N1"))

#N stage 8th
PDAC_data_r$N_stage8th <-factor(PDAC_data_r$N_stage8th,
                             levels = c("0","1","2"),
                             labels = c("N0", "N1", "N2"))

#UE image
PDAC_data_r$UET1WI <-factor(PDAC_data_r$UET1WI,
                             levels = c("0","1"),
                             labels = c("Hyper", "Hypo"))

#AP
PDAC_data_r$AP <-factor(PDAC_data_r$AP,
                            levels = c("0","1"),
                            labels = c("Hyper", "Hypo"))

#PVP
PDAC_data_r$PVP <-factor(PDAC_data_r$PVP,
                            levels = c("0","1"),
                            labels = c("Hyper", "Hypo"))
#DP
PDAC_data_r$DP <-factor(PDAC_data_r$DP,
                         levels = c("0","1"),
                         labels = c("Hyper", "Hypo"))

#secondary signs
PDAC_data_r$Secondary_signs <-factor(PDAC_data_r$Secondary_signs,
                         levels = c("0","1"),
                         labels = c("Neg", "Pos"))

#ADC
PDAC_data_r<-PDAC_data_r %>% mutate(ADC_group=ifelse(ADC_value>=1400,"High","Low"))
PDAC_data_r$ADC_group<-factor(PDAC_data_r$ADC_group)







### KM curve

surv_object<-Surv(time=PDAC_data_r$RFS,event = PDAC_data_r$Recurrence_st2)
surv_object
fit1<-survfit(surv_object~Intratumoral_cystic_type, data = PDAC_data_r)
summary(fit1)

pairwise_survdiff(Surv(RFS,Recurrence_st2) ~ Intratumoral_cystic_type, PDAC_data_r)

SURVPLOT1<-ggsurvplot(fit1, data=PDAC_data_r, 
                     risk.table = TRUE,
                     risk.table.height = 0.2,
                     risk.table.col = "strata",
                     risk.table.fontsize = 7,
                     risk.table.y.text.col = TRUE,
                     tables.theme = theme_survminer(font.main = 12),
                     #legend.title = "",
                     legend.labs = c("Ordinary PDAC", "PDAC with \nneoplastic mucin cyst", "PDAC with \nimaging necrosis"),
                     legend=c(0.825,0.73),
                     conf.int = FALSE, 
                     conf.int.style = "step",
                     conf.int.alpha = 0.1,
                     cumevents = FALSE,
                     palette = "aaas",
                     pval= FALSE, 
                     pval.method = FALSE,
                     pval.coord = c(0.0,0.0),
                     pval.method.coord = c(0.05,0.05),
                     pval.size = 4,
                     pval.method.size = 4,
                     # ncensor.plot = TRUE
                     )

SURVPLOT1$plot = SURVPLOT1$plot + 
  ylab("Recurrence-Free Survival Rate (%)")+
  geom_text(aes(x=18.5,y=0.98,label= "italic('p')"),size = 5.5, parse = TRUE)+
  geom_text(aes(x=49.5,y=0.98,label= " = .002 for ordinary PDAC vs. PDAC with neoplastic mucin cyst"),size = 5.5)+
  
  geom_text(aes(x=18.5,y=0.93,label= "italic('p')"),size = 5.5, parse = TRUE)+
  geom_text(aes(x=47.3,y=0.93,label= " = .001 for ordinary PDAC vs. PDAC with imaging necrosis"),size = 5.5)+
  scale_y_continuous(breaks = seq(0,1,by = 0.1), labels = paste0(seq(0,100,by = 10),"%")) + 
  xlab("Time (months)")+
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18, vjust = -15),
        #legend.background = element_rect(linetype = "solid", colour = "black")
        )+
  theme(legend.text = element_text(size = 13, color = "black"))+
  theme(legend.title=element_blank())

SURVPLOT1$table = SURVPLOT1$table + 
  ylab(NULL)+xlab(NULL)+
  #ggtitle("PDAC")+
  scale_x_continuous(breaks = c(0,20,40,60,80), labels = c("","","","",""))
# +
#   theme(axis.text.x = element_text(size = 10), 
#         axis.text.y = element_text(size = 12, face = "bold"), 
#         axis.title.x = element_text(size = 10), 
#         axis.title.y = element_text(size = 10)
#         )

SURVPLOT1

ggsave("RFS_revision.tiff", plot = print(SURVPLOT1), dpi=200, dev='tiff', height=8.27, width=11.69, units="in")

fit2<-survfit(surv_object~N_stage, data = PDAC_data_r)
summary(fit2)

SURVPLOT2<-ggsurvplot(fit2, data=PDAC_data_r, 
                      risk.table = TRUE, 
                      risk.table.height = 0.2,
                      legend.title = "N stage",
                      legend.labs = c("N0", "N1 or N2"),
                      legend=c(0.875,0.085),
                      conf.int = FALSE, 
                      conf.int.style = "step",
                      conf.int.alpha = 0.1,
                      cumevents = FALSE,
                      palette = "lancet",
                      # tables.theme = theme_cleantable(),
                      pval= FALSE, 
                      pval.method = FALSE,
                      pval.coord = c(0.0,0.0),
                      pval.method.coord = c(0.05,0.05),
                      pval.size = 4,
                      pval.method.size = 4,
                      # ncensor.plot = TRUE
)

SURVPLOT2$plot = SURVPLOT2$plot + 
  ylab("Recurrence-Free Survival probability (%)")+
  geom_text(aes(x=65,y=0.9,label= "p = 0.016"),size = 6)+
  scale_y_continuous(breaks = seq(0,1,by = 0.1), labels = paste0(seq(0,100,by = 10),"%")) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.background = element_rect(linetype = "solid", colour = "black")
  )

SURVPLOT2$table = SURVPLOT2$table + 
  ylab(NULL)+
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)
  )

SURVPLOT2
ggpar(SURVPLOT2, font.legend = c(15,"bold"))


fit3<-survfit(surv_object~DP, data = PDAC_data_r)
summary(fit3)

SURVPLOT3<-ggsurvplot(fit3, data=PDAC_data_r, 
                      risk.table = TRUE, 
                      risk.table.height = 0.2,
                      legend.title = "Delayed Phase",
                      legend.labs = c("Iso or Hyperintensity","Hypointensity"),
                      legend=c(0.875,0.085),
                      conf.int = FALSE, 
                      conf.int.style = "step",
                      conf.int.alpha = 0.1,
                      cumevents = FALSE,
                      palette = "lancet",
                      # tables.theme = theme_cleantable(),
                      pval= FALSE, 
                      pval.method = FALSE,
                      pval.coord = c(0.0,0.0),
                      pval.method.coord = c(0.05,0.05),
                      pval.size = 4,
                      pval.method.size = 4,
                      # ncensor.plot = TRUE
)

SURVPLOT3$plot = SURVPLOT3$plot + 
  ylab("Recurrence-Free Survival probability (%)")+
  geom_text(aes(x=65,y=0.9,label= "p = 0.027"),size = 6)+
  scale_y_continuous(breaks = seq(0,1,by = 0.1), labels = paste0(seq(0,100,by = 10),"%")) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.background = element_rect(linetype = "solid", colour = "black")
  )

SURVPLOT3$table = SURVPLOT3$table + 
  ylab(NULL)+
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)
  )

SURVPLOT3
ggpar(SURVPLOT3, font.legend = c(15,"bold"))


fit4<-survfit(surv_object~Diffusion_restriction, data = PDAC_data_r)
summary(fit4)

SURVPLOT4<-ggsurvplot(fit4, data=PDAC_data_r, 
                      risk.table = TRUE, 
                      risk.table.height = 0.2,
                      legend.title = "Diffusion restriction",
                      legend.labs = c("Negative","Positive"),
                      legend=c(0.875,0.085),
                      conf.int = FALSE, 
                      conf.int.style = "step",
                      conf.int.alpha = 0.1,
                      cumevents = FALSE,
                      palette = "lancet",
                      # tables.theme = theme_cleantable(),
                      pval= TRUE, 
                      pval.method = FALSE,
                      pval.coord = c(0.0,0.0),
                      pval.method.coord = c(0.05,0.05),
                      pval.size = 4,
                      pval.method.size = 4,
                      # ncensor.plot = TRUE
)

SURVPLOT4$plot = SURVPLOT4$plot + 
  ylab("Recurrence-Free Survival probability (%)")+
  geom_text(aes(x=65,y=0.9,label= "p = 0.014"),size = 6)+
  scale_y_continuous(breaks = seq(0,1,by = 0.1), labels = paste0(seq(0,100,by = 10),"%")) + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        legend.background = element_rect(linetype = "solid", colour = "black")
  )

SURVPLOT4$table = SURVPLOT4$table + 
  ylab(NULL)+
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)
  )

SURVPLOT4
ggpar(SURVPLOT4, font.legend = c(15,"bold"))

splots<-list()
splots[[1]]<-SURVPLOT1
splots[[2]]<-SURVPLOT2
splots[[3]]<-SURVPLOT3
splots[[4]]<-SURVPLOT4
p<-arrange_ggsurvplots(splots, print = TRUE, ncol=2, nrow = 2)
#ggpar(p, font.legend = c(15,"bold"))
