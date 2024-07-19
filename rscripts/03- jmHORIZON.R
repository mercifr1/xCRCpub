#' ####################################################
#'
#'   Joint models
#'   HORIZONIII study (from PDS)
#'
#'   Q3-2024
#'   Francois Mercier
#'
#' ####################################################


#' sessionInfo()


#' biom.df<-readRDS("./data/PRIME/HorizTGIads.rds")
#' summary(biom.df)

#' event.df<-readRDS("./data/PRIME/HorizOSads.rds")
#' summary(event.df)


#' ===============================================
#' Model SLD
#' ===============================================

tgi.dat<-DataJoint(
  subject = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  longi   = DataLongitudinal(data=biom.df, threshold=3, formula= BIOMVAL~BIOMYR)
  )

tgi.in<-JointModel(longitudinal=LongitudinalSteinFojo(
  mu_bsld=prior_normal(log(70), .2),
  mu_ks=prior_normal(log(0.1), .3),
  mu_kg=prior_normal(log(0.5), .3),
  omega_bsld=prior_lognormal(log(0.5), 0.3),
  omega_ks=prior_lognormal(log(0.7), 0.3),
  omega_kg=prior_lognormal(log(0.7), 0.3),
  sigma=prior_lognormal(log(0.2), 0.1)))

tgi.samples<-sampleStanModel(tgi.in, data=tgi.dat, iter_warmup=2000,
  iter_sampling=1000, chains=3, parallel=3, refresh=500)
tgi.out<-as.CmdStanMCMC(tgi.samples)
print(tgi.out, max_rows=500, digits=5)

#' Display profiles OBS vs IPRED for 10 random individuals
selected_subjects<-sample(event.df$SUBJID, 10)
longquant_obs<-LongitudinalQuantities(tgi.samples, grid=GridObserved(subjects=selected_subjects))
g2<-autoplot(longquant_obs)+
    theme_minimal()+
    theme(panel.grid.minor=element_blank())
g2

#' om_BSLD: hist(rlnorm(1000, log(0.5), 0.3))
#' om_ks: hist(rlnorm(1000, log(0.7), 0.3))
#' om_sigma: hist(rlnorm(1000, log(0.2), 0.1))


#' ===============================================
#' Model OS
#' ===============================================

surv.dat<-DataJoint(
  subject  = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  survival = DataSurvival(data=event.df, formula=Surv(EVENTYR, EVENTFL)~1)
)

surv.in<-JointModel(survival=SurvivalWeibullPH())

surv.samples<-sampleStanModel(surv.in, data=surv.dat, iter_warmup=2000,
    iter_sampling=1000, chains=3, parallel=3, refresh=500)
surv.out<-as.CmdStanMCMC(surv.samples)
print(surv.out, max_rows=500, digits=5)

#' Display PRED vs OBS surv curves
expected.surv<-SurvivalQuantities(surv.samples, type="surv",
  grid=GridGrouped(times=seq(from=0, to=3, by=0.1),
                   groups=split(event.df$SUBJID, event.df$ATRT)))

mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])

legendtitle<-"Treatment group"
g3<-autoplot(expected.surv, add_km=T, add_wrap=F)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  labs(x="Time (years)", color=legendtitle, fill=legendtitle)+
  theme_minimal()
g3


#' ===============================================
#' Model JM
#' ===============================================

jm.dat<-DataJoint(
  subject  = DataSubject(data=event.df, subject="SUBJID", arm="ATRT", study="STUDY"),
  longitudinal = DataLongitudinal(data=biom.df, threshold=2, formula= BIOMVAL~BIOMYR),
  survival = DataSurvival(data=event.df, formula=Surv(EVENTYR, EVENTFL)~1)
)

jm.in<-JointModel(
  longitudinal=LongitudinalSteinFojo(
    mu_bsld=prior_normal(log(70), .2),
    mu_ks=prior_normal(log(0.1), .3),
    mu_kg=prior_normal(log(0.5), .3),
    omega_bsld=prior_lognormal(log(0.5), 0.3),
    omega_ks=prior_lognormal(log(0.7), 0.3),
    omega_kg=prior_lognormal(log(0.7), 0.3),
    sigma=prior_lognormal(log(0.2), 0.1)),
  survival=SurvivalWeibullPH(),
  link=Link(linkTTG(prior_normal(0.01, 3)))
)

jm.samples<-sampleStanModel(jm.in, data=jm.dat,
    iter_sampling=1000, iter_warmup=2000, chains=3, parallel_chains=3)

jm.out<-as.CmdStanMCMC(jm.samples)
print(jm.out, max_rows=500, digits=5)

#' Display PRED vs OBS surv curves
expected.surv<-SurvivalQuantities(jm.samples, type="surv",
    grid=GridGrouped(times=seq(from=0, to=3, by=0.1),
    groups=split(event.df$SUBJID, event.df$ATRT)))

mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])

legendtitle<-"Treatment group"
g4<-autoplot(expected.surv, add_km=T, add_wrap=F)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  labs(x="Time (years)", color=legendtitle, fill=legendtitle)+
  theme_minimal()
g4
