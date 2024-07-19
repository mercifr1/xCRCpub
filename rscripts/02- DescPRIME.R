#' ####################################################
#'
#'   Descriptive plots
#'   PRIME study (from PDS)
#'
#'   Q3-2024
#'   Francois Mercier
#'
#' ####################################################


#' sessionInfo()

#' ===============================================
#' Visualize SLD
#' ===============================================

#' biom.df<-readRDS("./data/PRIME/PRIMETGIads.rds")
#' summary(biom.df)

#' Display SLD spaghetti
ybreaks<-c(3, 30, 100, 300)
mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])
ggplot(data=biom.df, aes(x=BIOMYR, y=BIOMVAL))+
    geom_point(colour="grey33", alpha=0.3, size=0.9)+
    geom_line(aes(group=SUBJID, colour=as.factor(ATRT)), alpha=0.6)+
    facet_wrap(~ATRT)+
    scale_x_continuous("Year", breaks=0.5*(0:5))+
    scale_y_continuous("SLD (mm)", breaks=ybreaks)+
    scale_colour_manual(values=mycols, guide="none")+
    theme_minimal()+
    theme(panel.grid.minor=element_blank())
ggsave("./data/PRIME/20240704 SpaghTGI.jpg", width=200, height=200/1.8, units="mm", dpi=600)


#' ===============================================
#' Visualize OS
#' ===============================================

#' event.df<-readRDS("./data/PRIME/PRIMEOSads.rds")
#' summary(event.df)

event.df |>
  group_by(ATRT) |>
  count(EVENTFL)

#' To extract the median survival in months (for prez to Basel Biostat Forum 2024-07-15)
temp<-event.df |> mutate(EVENTMO=EVENTYR*12)
survfit(formula = Surv(EVENTMO, EVENTFL)~ATRT, data = temp, type = "kaplan-meier")

#' CoxPH
cox<-coxph(Surv(EVENTYR, EVENTFL)~ATRT, data=event.df)
summary(cox)

os.kmest<-survfit(Surv(EVENTYR, EVENTFL)~ATRT, data=event.df)

#' Display OS KM
mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])
g1<-survminer::ggsurvplot(os.kmest,
  data=event.df, risk.table=T, break.x.by=.5, legend.title="",
  xlab="Year", ylab="Overall survival", palette = mycols,
  risk.table.fontsize=4, legend=c(0.8, 0.8))
g1
jpeg("./data/PRIME/20240704 KMOS.jpg", width=200, height=200, units="mm", res=600)
print(g1)
dev.off()





