#'---
#' title: Case-crossover analysis - Example 2
#' author: ""
#' date: ""
#' output:
#'  github_document
#'---

#+ message=FALSE
pacman::p_load(
    rio,            # import and export files
    here,           # locate files 
    tidyverse,      # data management and visualization
    skimr,
    splines,
    gnm,            # Generalized Nonlinear Models
    Epi             # Statistical Analysis in Epidemiology
)

#' # Data
# data #-----------
ls_wdf <- readRDS(here("note_lehuynh/working_df.RDS"))

(dat.2city <- rbind(ls_wdf[[1]] %>% select(-allage1, -allage2, -ind),
                   ls_wdf[[2]] %>% select(-ind)))

names(dat.2city)

(dat.2city <- dat.2city %>% 
    mutate(city = as.factor(city),
           stratum = as.factor(city:year:month:dow)))

(df_ind <- dat.2city %>% 
    group_by(stratum) %>% 
    summarise(ind = sum(all)))

(wdf <- dat.2city %>% 
    left_join(df_ind,
              by = join_by(stratum)))

#' # Fit conditional Poisson with space-time-stratified strata
# space-time-stratified strata #-----------
model.cc.adj.city <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                         data = wdf, 
                         family = quasipoisson, 
                         subset = ind>0, 
                         eliminate = stratum)

summary(model.cc.adj.city)

Epi::ci.exp(model.cc.adj.city, subset = "l01pm10") 

#' # Stratified analysis by city
# stratified analysis by city #-----------------------
#' ## Valencia
## Valencia #-----------------------
model.cc.vlc <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data = wdf %>% filter(city == "Valencia"), 
                    family = quasipoisson, 
                    subset = ind>0, 
                    eliminate = stratum)

summary(model.cc.vlc)

#' ## London
## London #-------------
model.cc.ldn <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data = wdf %>% filter(city == "London"), 
                    family = quasipoisson, 
                    subset = ind>0, 
                    eliminate = stratum)

summary(model.cc.ldn)

#' ## Relative Risk
## relative risk #-----------------------
Epi::ci.exp(model.cc.vlc, subset = "l01pm10") 

Epi::ci.exp(model.cc.ldn, subset = "l01pm10") 

#' # Interaction analysis by location
# interaction by location #--------------
model.cc.city.int <- gnm(all ~ ns.tmean + ns.rh + city + l01pm10:city, 
                         data = wdf, 
                         family = quasipoisson, 
                         subset = ind>0, 
                         eliminate = stratum)

summary(model.cc.city.int)

Epi::ci.exp(model.cc.city.int, subset = "l01pm10") 

# Likelihood ratio test for effect modification.
anova(model.cc.adj.city, model.cc.city.int, test = "LRT")

