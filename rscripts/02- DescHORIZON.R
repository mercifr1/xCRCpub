#' ####################################################
#'   
#'   2stgTGIOS
#'   Descriptive plots
#'   
#'   Q1-2022
#'   Francois Mercier 
#'   
#' ####################################################


#' sessionInfo()

os1<-readRDS("../HORIZONIII/data/HorizOSads.rds")
tk1<-readRDS("../HORIZONIII/data/HorizTGIads.rds")

#' ===============================================
#' Visualize SLD
#' ===============================================

#' Display BSLD distribution
xbreaks<-c(1, 3, 10, 30, 100, 300)
bslddf<-tk1 %>%
  filter(BL_VIS=="1") %>%
  mutate(BSLDmed=median(BSLD),
         BSLDq10=quantile(BSLD, probs=.1),
         BSLDq90=quantile(BSLD, probs=.9))
  
ggplot(bslddf, aes(x=BSLD))+
  geom_density(aes(y=..count..))+
  geom_vline(aes(xintercept=BSLDmed), lty=2)+
  geom_vline(aes(xintercept=BSLDq10), lty=2)+
  geom_vline(aes(xintercept=BSLDq90), lty=2)+
  scale_x_continuous("Baseline SLD (mm)", trans="log", breaks=xbreaks)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())

#' Display SLD spaghetti in Male vs Female
tk2<-left_join(tk1, os1, by="UID") 

ybreaks<-c(3, 30, 100, 300)
g0<-ggplot(tk2, aes(x=TKYEAR, y=SLD))+
  geom_line(aes(group=UID), colour="wheat4", alpha=0.2)+
  geom_point(colour="grey33", alpha=0.3, size=0.9)+
  facet_wrap(~SEX)+
  scale_x_continuous("Year", breaks=0.5*(0:5))+
  scale_y_continuous("SLD (mm)", breaks=ybreaks)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())
g0


#' ===============================================
#' Visualize OS
#' ===============================================

os2<-os1 %>%
  drop_na(LDH1_5)
os2 %>% count(LDH1_5)

#' Display OS KM in all
cox<-coxph(Surv(OSYEAR, 1-OSCEN)~1, data=os2)
os.kmest<-survfit(Surv(OSYEAR, 1-OSCEN)~1, data=os2)

g1<-survminer::ggsurvplot(os.kmest,
                          data=os2,
                          risk.table=T,
                          break.x.by=0.5, legend.title="",
                          xlab="Time (year)", ylab="Overall survival",
                          risk.table.fontsize=4, legend=c(0.8, 0.8))
g1


#' Display OS KM in LDH high vs LDH low
cox<-coxph(Surv(OSYEAR, 1-OSCEN)~LDH1_5, data=os2)
os.kmest<-survfit(Surv(OSYEAR, 1-OSCEN)~LDH1_5, data=os2)

g2<-survminer::ggsurvplot(os.kmest,
                          data=os2,
                          risk.table=T,
                          break.x.by=0.5, legend.title="",
                          xlab="Time (year)", ylab="Overall survival",
                          risk.table.fontsize=4, legend=c(0.8, 0.8))
g2

