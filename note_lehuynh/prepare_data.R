#'---
#' title: Data preparation
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

#' # Check data
# check data #-----------
(df_valencia <- rio::import(here("valencia.csv")) %>% tibble())

skimr::skim(df_valencia)

(df_london <- rio::import(here("london.csv")) %>% tibble())

skimr::skim(df_london)

#' # Data preparation
# prepare data #-----------
#' ### Generate time-stratified strata
#' "1826 observations were split into 420 stratum sets defined by day-of-week 
#' within-month and year."
a1 <- df_valencia %>% 
    mutate(across(.cols = c(month, year, dow),
                  as.factor),
           stratum = as.factor(year:month:dow))

a1 %>% print(n = 20)

a1 %>% count(stratum)

#' ### Generate 'ind' = sum of events in each stratum
a1 %>%
    group_by(stratum) %>% 
    summarise(ind = sum(all))

#' ### Generate lagX for temperature, humidity, and PM10
df_valencia %>% 
    select(tmean, rh, pm10) %>% 
    mutate(across(.cols = c(tmean, rh, pm10),
                  ~ dplyr::lag(.x, n = 1),
                  .names = "{.col}_l1"),
           across(.cols = c(tmean, rh, pm10),
                  ~ dplyr::lag(.x, n = 3),
                  .names = "{.col}_l3"))

#' ### Generate splines for environmental factor adjustment
ns(df_valencia$tmean, df = 6)

df_valencia %>% 
    select(tmean) %>% 
    mutate(ns.tmean = ns(tmean, df = 6))

#' ### Generate working dataframes - 2 cities
wls <- list(df_valencia, df_london) %>% 
    set_names(c("valencia", "london")) %>% 
    map(\(data){
        
        wdf <- data %>% 
            mutate(across(.cols = c(month, year, dow),
                          as.factor),
                   # generate time-stratified strata
                   stratum = as.factor(year:month:dow),
                   # generate lag0 for temperature, humidity, and PM10
                   across(.cols = c(tmean, rh, pm10),
                          ~ dplyr::lag(.x, n = 0),
                          .names = "{.col}_l0"),
                   # generate lag1 for temperature, humidity, and PM10
                   across(.cols = c(tmean, rh, pm10),
                          ~ dplyr::lag(.x, n = 1),
                          .names = "{.col}_l1"),
                   # generate lag2 for temperature, humidity, and PM10
                   across(.cols = c(tmean, rh, pm10),
                          ~ dplyr::lag(.x, n = 2),
                          .names = "{.col}_l2"),
                   # generate lag3 for temperature, humidity, and PM10
                   across(.cols = c(tmean, rh, pm10),
                          ~ dplyr::lag(.x, n = 3),
                          .names = "{.col}_l3"),
                   l03tmean = (tmean_l0 + tmean_l1 + tmean_l2 + tmean_l3)/4,
                   l03rh = (rh_l0 + rh_l1 + rh_l2 + rh_l3)/4,
                   # generate averaged lag exposure for PM10 sale for a 10 ug/m3 increase
                   l01pm10 = ((pm10_l0 + pm10_l1)/2)/10,
                   # generate splines for temperature adjustment
                   ns.tmean = ns(l03tmean, df = 6),
                   # generate splines for humidity adjustment
                   ns.rh = ns(l03rh, df = 3))
        
        # Generate 'ind' = sum of events in each stratum
        df_ind <- data %>% 
                mutate(across(.cols = c(month, year, dow),
                              as.factor),
                       stratum = as.factor(year:month:dow)) %>%
                group_by(stratum) %>% 
                summarise(ind = sum(all))
        
        res <- wdf %>% 
            left_join(df_ind, by = join_by(stratum))
        
        return(res)
    })

#' Save working datasets
saveRDS(wls, here("note_lehuynh/working_df.RDS"))

