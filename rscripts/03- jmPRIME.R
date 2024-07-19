#' ####################################################
#'
#'   Joint models
#'   PRIME study (from PDS)
#'
#'   Q3-2024
#'   Francois Mercier
#'
#' ####################################################


#' sessionInfo()


#' biom.df<-readRDS("./data/PRIME/PRIMETGIads.rds")
#' summary(biom.df)

#' event.df<-readRDS("./data/PRIME/PRIMEOSads.rds")
#' summary(event.df)


#' ===============================================
#' Model SLD
#' ===============================================

tgi.dat<-DataJoint(
  subject = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  longi   = DataLongitudinal(data=biom.df, threshold=3, formula= BIOMVAL~BIOMYR)
  )

tgi.in<-JointModel(longitudinal=LongitudinalSteinFojo(
    mu_bsld = prior_normal(log(70), .1),
    mu_ks = prior_normal(log(1.8), .1),
    mu_kg = prior_normal(log(0.15), .1),
    omega_bsld = prior_lognormal(log(0.1), .1),
    omega_ks = prior_lognormal(log(0.1), .5),
    omega_kg = prior_lognormal(log(0.1), .5),
    sigma = prior_lognormal(log(0.18), .5),
))

tgi.samples<-sampleStanModel(tgi.in, data=tgi.dat,
                             iter_sampling = 1000,
                             iter_warmup = 2000,
                             chains = 3, parallel_chains = 3)
tgi.out<-as.CmdStanMCMC(tgi.samples)
print(tgi.out, max_rows=500, digits=5)

#' Display profiles OBS vs IPRED for 10 random individuals
set.seed(1869)
selected_subjects<-sample(event.df$SUBJID, 10)
select_longquant_obs<-LongitudinalQuantities(tgi.samples, grid=GridObserved(subjects=selected_subjects))

mycols<-rep(ghibli::ghibli_palettes$MarnieMedium1, 2)
autoplot(select_longquant_obs)+
  labs(x="Time (years)", y="SLD (mm)")+
  theme_minimal()
ggsave("./data/PRIME/20240704 ObsPredTGI.jpg", width=300, height=300/1.8, units="mm", dpi=600)

#' extracting the individual parameters psi_xxx
#' for instance, ks_i
summary(tgi.out$draws("lm_sf_psi_ks"))
summary(tgi.out$draws("lm_sf_psi_kg"))
sigma<-summary(tgi.out$draws("lm_sf_sigma"))$median

#' IPRED vs. OBS
all_longquant_obs<-LongitudinalQuantities(tgi.samples, grid=GridObserved())
#' summary(all_longquant_obs)
all_longquant_obs_df<-tibble(summary(all_longquant_obs))
prepdf<-left_join(biom.df, all_longquant_obs_df, 
                  by=c("SUBJID"="group", "BIOMYR"="time"))
ggplot(prepdf, aes(x=median, y=BIOMVAL))+
  scale_x_continuous(trans="log", breaks=c(1, 3, 30, 300))+
  scale_y_continuous(trans="log", breaks=c(1, 3, 30, 300))+
  coord_equal(ratio=1)+
  labs(x="Indiv predicted SLD (mm)", y="Observed SLD (mm)")+
  geom_point(colour="grey22", alpha=0.3)+
  geom_abline(aes(slope=1, intercept=0), lty=2, lwd=1, colour="blue")+
  theme_minimal()
ggsave("./data/PRIME/20240704 ObsPredTGI-Identity.jpg", width=150, height=150, units="mm", dpi=600)

#' IWRES vs. IPRED
prepdf_iwres<-prepdf |>  mutate(IWRES=(BIOMVAL-median)/(sigma*median))
ggplot(prepdf_iwres, aes(x=median, y=IWRES))+
  labs(x="Indiv predicted SLD (mm)", y="IWRES")+
  geom_point(colour="grey22", alpha=0.3)+
  geom_abline(aes(slope=0, intercept=0), lty=2, lwd=1, colour="blue")+
  theme_minimal()
ggsave("./data/PRIME/20240704 IWRESPredTGI.jpg", width=150, height=50, units="mm", dpi=600)

prepdf_iwres<-prepdf |>  mutate(IWRES=(BIOMVAL-median)/(sigma*median))
ggplot(prepdf_iwres, aes(x=BIOMYR, y=IWRES))+
  labs(x="Time (year)", y="IWRES")+
  geom_point(colour="grey22", alpha=0.3)+
  geom_abline(aes(slope=0, intercept=0), lty=2, lwd=1, colour="blue")+
  theme_minimal()
ggsave("./data/PRIME/20240704 IWRESTimeTGI.jpg", width=150, height=50, units="mm", dpi=600)




#' ===============================================
#' Model OS
#' ===============================================

surv.dat<-DataJoint(
  subject  = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  survival = DataSurvival(data=event.df, formula=Surv(EVENTYR, EVENTFL)~ATRT)
)

surv.in<-JointModel(survival=SurvivalWeibullPH())

surv.samples<-sampleStanModel(surv.in, data=surv.dat,
                              iter_sampling = 1000,
                              iter_warmup = 2000,
                              chains = 3, parallel_chains = 3)
surv.out<-as.CmdStanMCMC(surv.samples)
print(surv.out, max_rows=500, digits=5)

#' LOOIC
surv.out$loo()


#' Display PRED vs OBS surv curves
surv.out.surv<-SurvivalQuantities(surv.samples, type="surv",
  grid=GridGrouped(times=seq(from=0, to=4, by=0.1),
  groups=split(event.df$SUBJID, event.df$ATRT)))

mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])

legendtitle<-"Treatment group"
autoplot(surv.out.surv, add_km=T, add_wrap=F)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_continuous(breaks=0:4)+
  labs(x="Time (years)", color=legendtitle, fill=legendtitle)+
  theme_minimal()
ggsave("./data/PRIME/20240704 ObsPredOS.jpg", width=200, height=200/1.8, units="mm", dpi=600)


#' ===============================================
#' Predicted HR by surv only model
#' ===============================================

surv.out.haz<-SurvivalQuantities(surv.samples, type="haz",
  grid=GridGrouped(times=seq(from=0, to=4, by=0.1),
  groups=split(event.df$SUBJID, event.df$ATRT)))

#' Extract samples as df
#' Derive the point wise mode and associated 95%HDCI, by group
surv.out.haz_df<-as.data.frame(surv.out.haz)

#' ExtractHR
surv.out.haz_fol<-surv.out.haz_df |> filter(group=="FOLFOX alone") |>
  rename(hazfol=values) |> select(time, hazfol) |> 
  mutate(iter=rep(1:3000, 41)) |>
  filter(time==2)
surv.out.haz_pan<-surv.out.haz_df |> filter(group=="Panitumumab + FOLFOX") |>
  rename(hazpan=values) |> select(time, hazpan) |> 
  mutate(iter=rep(1:3000, 41)) |>
  filter(time==2)

surv.out.HR_df<-left_join(surv.out.haz_fol, surv.out.haz_pan, by="iter") |>
  mutate(HRvalues=hazpan/hazfol)

surv.out.HR_df |>
  drop_na(HRvalues) |>
  summarize(median=median(HRvalues),
            p025=quantile(HRvalues, probs=0.025),
            p975=quantile(HRvalues, probs=0.975))


#' ===============================================
#' Model with no link = JM0.
#' ===============================================

jm0.in<-JointModel(
  longitudinal=LongitudinalSteinFojo(
    mu_bsld = prior_normal(log(70), .1),
    mu_ks = prior_normal(log(1.8), .1),
    mu_kg = prior_normal(log(0.15), .1),
    omega_bsld = prior_lognormal(log(0.1), .1),
    omega_ks = prior_lognormal(log(0.1), .5),
    omega_kg = prior_lognormal(log(0.1), .5),
    sigma = prior_lognormal(log(0.18), .5),
  ),
  survival = SurvivalWeibullPH(),
  link = Link(linkNone())
)


jm0.dat<-DataJoint(
  subject  = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  longitudinal = DataLongitudinal(data=biom.df, threshold=3, formula= BIOMVAL~BIOMYR),
  survival = DataSurvival(data=event.df, formula=Surv(EVENTYR, EVENTFL)~1)
)

jm0.samples<-sampleStanModel(jm0.in,
                             data = jm1.dat,
                             iter_sampling = 1000,
                             iter_warmup = 2000,
                             chains = 3, parallel_chains = 3)

jm0.out<-as.CmdStanMCMC(jm0.samples)

#' LOOIC
jm0.out$loo()


#' ===============================================
#' Model with no trt on OS = JM1.
#' ===============================================

jm.in<-JointModel(
  longitudinal=LongitudinalSteinFojo(
    mu_bsld = prior_normal(log(70), .1),
    mu_ks = prior_normal(log(1.8), .1),
    mu_kg = prior_normal(log(0.15), .1),
    omega_bsld = prior_lognormal(log(0.1), .1),
    omega_ks = prior_lognormal(log(0.1), .5),
    omega_kg = prior_lognormal(log(0.1), .5),
    sigma = prior_lognormal(log(0.18), .5),
  ),
  survival = SurvivalWeibullPH(),
  link = Link(linkTTG(prior_normal(0.01, 3)))
)

jm1.dat<-DataJoint(
  subject  = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  longitudinal = DataLongitudinal(data=biom.df, threshold=3, formula= BIOMVAL~BIOMYR),
  survival = DataSurvival(data=event.df, formula=Surv(EVENTYR, EVENTFL)~1)
)

jm1.samples<-sampleStanModel(jm.in,
                                 data = jm1.dat,
                                 iter_sampling = 1000,
                                 iter_warmup = 2000,
                                 chains = 3, parallel_chains = 3)

jm1.out<-as.CmdStanMCMC(jm1.samples)
print(jm1.out, max_rows=500, digits=5)

#' LOOIC
jm1.out$loo()

#' Display PRED vs OBS surv curves
jm1.out.surv<-SurvivalQuantities(jm1.samples, type="surv",
   grid=GridGrouped(times=seq(from=0, to=4, by=0.1),
   groups=split(event.df$SUBJID, event.df$ATRT)))

legendtitle<-"Treatment group"
mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])
autoplot(jm1.out.surv, add_km=T, add_wrap=F)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_continuous(breaks=0:4)+
  labs(x="Time (years)", color=legendtitle, fill=legendtitle)+
  theme_minimal()
ggsave("./data/PRIME/20240710 ObsPredOS-fromJM1.jpg", width=200, height=200/1.8, units="mm", dpi=600)


#' ===============================================
#' Predicted HR by JM1
#' ===============================================

jm1.out.haz<-SurvivalQuantities(jm1.samples, type="haz",
   grid=GridGrouped(times=seq(from=0, to=4, by=0.1),
   groups=split(event.df$SUBJID, event.df$ATRT)))

#' Extract samples as df
#' Derive the point wise mode and associated 95%HDCI, by group
jm1.out.haz_df<-as.data.frame(jm1.out.haz)

jm1.out.haz_ci<-jm1.out.haz_df |>
  group_by(group, time) |>
  summarize(median=median(values),
            p025=quantile(values, probs=0.025),
            p975=quantile(values, probs=0.975))

legendtitle<-"Treatment group"
mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])
ggplot(jm1.out.haz_ci, aes(time, median))+
  geom_line(aes(colour=as.factor(group)))+
  geom_ribbon(aes(ymin=p025, ymax=p975, fill=as.factor(group)), alpha=0.3)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_continuous(breaks=0:4)+
  labs(x="Time (years)", y="Hazard",
       color=legendtitle, fill=legendtitle)+
  theme_minimal()+
  theme(legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal')
ggsave("./data/PRIME/20240709 OShaz-fromJM1.jpg", width=200, height=200/1.8, units="mm", dpi=600)

#' ExtractHR
jm1.out.haz_fol<-jm1.out.haz_df |> filter(group=="FOLFOX alone") |>
  rename(hazfol=values) |> select(time, hazfol) |> 
  mutate(iter=rep(1:3000, 41)) |>
  filter(time==2)
jm1.out.haz_pan<-jm1.out.haz_df |> filter(group=="Panitumumab + FOLFOX") |>
  rename(hazpan=values) |> select(time, hazpan) |> 
  mutate(iter=rep(1:3000, 41)) |>
  filter(time==2)

jm1.out.HR_df<-left_join(jm1.out.haz_fol, jm1.out.haz_pan, by="iter") |>
  mutate(HRvalues=hazpan/hazfol)

jm1.out.HR_df |>
  drop_na(HRvalues) |>
  summarize(median=median(HRvalues),
            p025=quantile(HRvalues, probs=0.025),
            p975=quantile(HRvalues, probs=0.975))


#' ===============================================
#' Model incl trt on OS = JM2
#' note: jm.in is the same as the one used in JM1
#' ===============================================

jm2.dat<-DataJoint(
  subject  = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  longitudinal = DataLongitudinal(data=biom.df, threshold=3, formula= BIOMVAL~BIOMYR),
  survival = DataSurvival(data=event.df, formula=Surv(EVENTYR, EVENTFL)~ATRT)
)

jm2.samples<-sampleStanModel(jm.in,
  data = jm2.dat,
  iter_sampling = 1000,
  iter_warmup = 2000,
  chains = 3, parallel_chains = 3)

jm2.out<-as.CmdStanMCMC(jm2.samples)
print(jm2.out, max_rows=500, digits=5)

#' LOOIC
jm2.out$loo()

#' Display PRED vs OBS surv curves
jm2.out.surv<-SurvivalQuantities(jm2.samples, type="surv",
       grid=GridGrouped(times=seq(from=0, to=4, by=0.1),
       groups=split(event.df$SUBJID, event.df$ATRT)))

legendtitle<-"Treatment group"
autoplot(jm2.out.surv, add_km=T, add_wrap=F)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_continuous(breaks=0:4)+
  labs(x="Time (years)", color=legendtitle, fill=legendtitle)+
  theme_minimal()
ggsave("./data/PRIME/20240704 ObsPredOS-fromJM2.jpg", width=200, height=200/1.8, units="mm", dpi=600)


#' ===============================================
#' Predicted HR by JM2
#' ===============================================

jm2.out.haz<-SurvivalQuantities(jm2.samples, type="haz",
          grid=GridGrouped(times=seq(from=0, to=4, by=0.1),
          groups=split(event.df$SUBJID, event.df$ATRT)))

#' Extract samples as df
#' Derive the point wise mode and associated 95%HDCI, by group
jm2.out.haz_df<-as.data.frame(jm2.out.haz)

jm2.out.haz_ci<-jm2.out.haz_df |>
  group_by(group, time) |>
  summarize(median=median(values),
            p025=quantile(values, probs=0.025),
            p975=quantile(values, probs=0.975))

legendtitle<-"Treatment group"
mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])
ggplot(jm2.out.haz_ci, aes(time, median))+
  geom_line(aes(colour=as.factor(group)))+
  geom_ribbon(aes(ymin=p025, ymax=p975, fill=as.factor(group)), alpha=0.3)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_continuous(breaks=0:4)+
  labs(x="Time (years)", y="Hazard",
       color=legendtitle, fill=legendtitle)+
  theme_minimal()+
  theme(legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal')
ggsave("./data/PRIME/20240709 OShaz-fromJM2.jpg", width=150, height=200/1.8, units="mm", dpi=600)


#' ExtractHR

jm2.out.haz_fol<-jm2.out.haz_df |> filter(group=="FOLFOX alone") |>
  rename(hazfol=values) |> select(time, hazfol) |> 
  mutate(iter=rep(1:3000, 41)) |>
  filter(time==2)
jm2.out.haz_pan<-jm2.out.haz_df |> filter(group=="Panitumumab + FOLFOX") |>
  rename(hazpan=values) |> select(time, hazpan) |> 
  mutate(iter=rep(1:3000, 41)) |>
  filter(time==2)

jm2.out.HR_df<-left_join(jm2.out.haz_fol, jm2.out.haz_pan, by="iter") |>
  mutate(HRvalues=hazpan/hazfol)

jm2.out.HR_df |>
  drop_na(HRvalues) |>
  summarize(median=median(HRvalues),
            p025=quantile(HRvalues, probs=0.025),
            p975=quantile(HRvalues, probs=0.975))


#' ===============================================
#' Overlaying survival curves 
#' JM0=SurvOnly vs JM1 vs JM2
#' ===============================================

srv.over<-as.data.frame(surv.out.surv) |> mutate(model="SurvOnly")
jm1.over<-as.data.frame(jm1.out.surv) |> mutate(model="JM1")
jm2.over<-as.data.frame(jm2.out.surv) |> mutate(model="JM2")
surv.over<-rbind(srv.over, jm1.over, jm2.over)
surv.over_df<-surv.over |>
  group_by(model, group, time) |>
  summarize(median=median(values),
            p025=quantile(values, probs=0.025),
            p975=quantile(values, probs=0.975))


legendtitle<-"Model"
mycols<-c(rev(ghibli::ghibli_palettes$MononokeMedium)[c(1, 3, 5)])
ggplot(surv.over_df, aes(time, median))+
  facet_wrap(~group)+
  geom_line(aes(colour=as.factor(model)))+
  geom_ribbon(aes(ymin=p025, ymax=p975, fill=as.factor(model)), alpha=0.4)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_continuous(breaks=0:4)+
  labs(x="Time (years)", y="S(t)",
       color=legendtitle, fill=legendtitle)+
  theme_minimal()
ggsave("./data/PRIME/20240710 OSover.jpg", width=250, height=200/1.8, units="mm", dpi=600)



