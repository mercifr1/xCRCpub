#' ####################################################
#'   
#'   2stgTGIOS
#'   Data preparation
#'   
#'   Q1-2022
#'   Francois Mercier 
#'   
#' ####################################################


#' sessionInfo()

#' ===============================================
#' Import and Select
#' ===============================================================

#' UID values
#' ---------------------------------
subj<-haven::read_sas("../HORIZONIII/data/rdpsubj.sas7bdat")
rcist<-haven::read_sas("../HORIZONIII/data/rdprcist.sas7bdat") 

keepin0<-c("UID", "LSVISDY", "SEX", "BEV_SDY", "BEV_EDY", "DEATDI", 
           "LDH1_5", "BL_VEGFN", "VEGF_STR")

#' - keep patients in PerProtocol set (PP_SET==1)
#' - keep rows from patients exposed to BEV for at least 3 weeks (i.e. at least 1st TA visit)
subj0<-subj %>% 
  filter(PP_SET==1, !is.na(BEV_SDY), !is.na(BEV_EDY), BEV_EDY>3*7) %>%
  mutate(UID=RANDCODE) %>%
  select(one_of(keepin0))
#' dim(subj0)
#' [1] 645  13

#' OS values
#' ---------------------------------
os0<-rcist %>% 
  mutate(UID=RANDCODE) %>%
  filter(UID %in% subj0$UID) %>%
  group_by(UID) %>%
    slice(1) %>%
  ungroup() %>%
  mutate(OSYEAR=OSTIM/365.25) %>%
  select(UID, DEATFLAG, OSTIM, OSCEN, OSYEAR)
#' length(os0$UID)
#' 645

#' SLD values
#' ---------------------------------
#' - remove rows where SLD is NA
#' - remove rows where patients have PBLCNT=NA i.e. remove patient with at least one post-baseline TA
tk0<-rcist %>% 
  mutate(UID=RANDCODE) %>%
  filter(UID %in% subj0$UID) %>%
  drop_na(STLDI) %>%
  filter(!is.na(PBLCNT)) %>%
  mutate(SLD=ifelse(STLDI==0, 2.5, STLDI*10), 
         BSLD=BL_STLDI*10,
         TKYEAR=ORDYTRT/365.25) 
#' length(unique(tk0$UID))
#' 640

#' ===============================================
#' Build analysis dataset
#' ===============================================================

os1<-left_join(subj0, os0, by="UID") %>%
  select(UID, SEX, OSTIM, OSCEN, OSYEAR, LDH1_5, BL_VEGFN, VEGF_STR)
#' length(unique(os1$UID))

tk1<-left_join(tk0, subj0, by="UID") %>%
  select(UID, ORDYTRT, TKYEAR, BL_VIS, BSLD, SLD, PBLCNT)
#' length(unique(tk1$UID))

#' saveRDS(os1, file="../HORIZONIII/data/HorizOSads.rds")
#' saveRDS(tk1, file="../HORIZONIII/data/HorizTGIads.rds")



