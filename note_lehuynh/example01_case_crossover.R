#'---
#' title: Case-crossover analysis - Example 1
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

ls_wdf

#' # Example 1
# example 1 #-----------
#' Exposure-outcome association and environmental time-varying confounders adjustment

ls_model1 <- ls_wdf %>% 
    map(\(data){
        # Fit fixed-effects conditional quasi-Poisson regression
        model <- gnm(all ~ ns.tmean + ns.rh + l01pm10,
                     data = data,
                     family = quasipoisson,
                     subset = ind > 0,
                     eliminate = stratum)
        
        return(model)
    })

#' ### Valencia
summary(ls_model1[[1]])

# Get Relative Risk for PM10.
Epi::ci.exp(ls_model1[[1]], subset = "l01pm10") 

#' ### London
summary(ls_model1[[2]])

# Get Relative Risk for PM10.
Epi::ci.exp(ls_model1[[2]], subset = "l01pm10") 

#' # Adjustment of sub-population time-invariant covariates
# sub-population #-----------------------
# reshape data
(df_model2 <- ls_wdf[[1]] %>% 
    pivot_longer(cols = c(allage1, allage2),
                 names_to = "age",
                 values_to = "allage"))

# Fit fixed-effects conditional quasi-Poisson regression adjusted by age
model.cc.adj.age <- gnm(allage ~ factor(age) + ns.tmean + ns.rh + l01pm10, 
                        data = df_model2,
                        family = quasipoisson,
                        eliminate = stratum)

summary(model.cc.adj.age)

Epi::ci.exp(model.cc.adj.age, subset = "l01pm10") 

#' Generate age-time-stratified strata
df_model2 %>% 
    mutate(age = as.factor(age),
           stratum4 = as.factor(year:month:dow:age)) %>% 
    select(year, month, dow, age, stratum4)

(df_model3 <- df_model2 %>% 
    mutate(age = as.factor(age),
           stratum4 = as.factor(year:month:dow:age)))

(df_ind4 <- df_model3 %>% 
        group_by(stratum4) %>% 
        summarise(ind4 = sum(all)))

(df_model3a <- df_model3 %>% 
    left_join(df_ind4,
              by = join_by(stratum4)))

# Fit conditional quasi-Poisson with age-time-stratified strata.
model.cc.str4 <- gnm(allage ~ ns.tmean + ns.rh + l01pm10, 
                     data = df_model3a,
                     family = quasipoisson,
                     subset = ind4>0,
                     eliminate = stratum4)

summary(model.cc.str4)

Epi::ci.exp(model.cc.str4, subset = "l01pm10") 
Epi::ci.exp(model.cc.adj.age, subset = "l01pm10") 
Epi::ci.exp(ls_model1[[1]], subset = "l01pm10") 

#' # Investigation of effect modification by age
# effect modification by age #---------------------
# age1 = <65 years
model.cc.age1 <- gnm(allage1 ~ ns.tmean + ns.rh + l01pm10, 
                     data = ls_wdf[["valencia"]], 
                     family = quasipoisson, 
                     subset = ind>0, 
                     eliminate = stratum)
summary(model.cc.age1)
Epi::ci.exp(model.cc.age1, subset = "l01pm10") 

# age2 = >65 years
model.cc.age2 <- gnm(allage2 ~ ns.tmean + ns.rh + l01pm10, 
                     data = ls_wdf[["valencia"]], 
                     family = quasipoisson, 
                     subset = ind>0, 
                     eliminate = stratum)
summary(model.cc.age2)
Epi::ci.exp(model.cc.age2, subset = "l01pm10") 

#' ## Interaction analysis by age 
model.cc.age.int <- gnm(allage ~ ns.tmean + ns.rh + factor(age) + l01pm10:factor(age), 
                        data = df_model3a, 
                        family = quasipoisson, 
                        subset = ind4>0, 
                        eliminate = stratum4)
summary(model.cc.age.int)
Epi::ci.exp(model.cc.age.int, subset = "l01pm10") 

# Likelihood ratio-test for effect for interaction.
anova(model.cc.str4, model.cc.age.int, test = "LRT")

