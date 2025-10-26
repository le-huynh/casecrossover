Supplementary_R_code.R
================

``` r
#######################################################################################
#######################################################################################
### Time-stratified case-crossover studies in environmental epidemiology:           ###
### a tutorial.                                                                     ###
### (2nd version submitted to IJE Educational Corner, 2023/04/01)                   ### 
#######################################################################################
#######################################################################################

# Install packages.
# install.packages("dplyr", "splines", "gnm", "Epi")
```

``` r
# Load packages.
library(dplyr); library(splines); library(gnm); library(Epi); library(here)

#######################################################################################
# Data management.
#######################################################################################

list.cities <- list()
list.cities[[1]]<-read.csv(here('valencia.csv'))
list.cities[[2]]<-read.csv(here('london.csv'))

list.cities[[1]] %>% tibble()
```

    ## # A tibble: 1,826 × 12
    ##    city     date       year month   day   dow   all allage1 allage2 tmean    rh  pm10
    ##    <chr>    <chr>     <int> <int> <int> <int> <int>   <int>   <int> <dbl> <dbl> <dbl>
    ##  1 Valencia 01jan2002  2002     1     1     2    31       5      26  8.90  85    17.8
    ##  2 Valencia 02jan2002  2002     1     2     3    29       3      26  8.5   89.4  21.7
    ##  3 Valencia 03jan2002  2002     1     3     4    23       4      19 10.8   86.8  43.2
    ##  4 Valencia 04jan2002  2002     1     4     5    17       2      15  9.60  88.1  19.8
    ##  5 Valencia 05jan2002  2002     1     5     6    19       6      13  9.20  88.8  15.9
    ##  6 Valencia 06jan2002  2002     1     6     0    26       4      22  9.70  87.2  20.5
    ##  7 Valencia 07jan2002  2002     1     7     1    25       6      19  8.80  72.4  48.9
    ##  8 Valencia 08jan2002  2002     1     8     2    21       0      21  9.40  78.5  26  
    ##  9 Valencia 09jan2002  2002     1     9     3    15       3      12  8.80  83.6  18.7
    ## 10 Valencia 10jan2002  2002     1    10     4    19       1      18  9.60  76.2  39.6
    ## # ℹ 1,816 more rows

``` r
list.cities[[2]] %>% tibble()
```

    ## # A tibble: 1,826 × 10
    ##    city   date       year month   day   dow   all   tmean    rh  pm10
    ##    <chr>  <chr>     <int> <int> <int> <int> <int>   <dbl> <dbl> <dbl>
    ##  1 London 01jan2002  2002     1     1     2   199 -0.225   75.7  71.7
    ##  2 London 02jan2002  2002     1     2     3   231  0.0875  77.5  40.2
    ##  3 London 03jan2002  2002     1     3     4   210  0.850   81.3  41.8
    ##  4 London 04jan2002  2002     1     4     5   203  0.538   85.4  50.4
    ##  5 London 05jan2002  2002     1     5     6   224  4.25    93.5  49.4
    ##  6 London 06jan2002  2002     1     6     0   198  7.07    96.4  31.1
    ##  7 London 07jan2002  2002     1     7     1   180  5.19    93.5  48.6
    ##  8 London 08jan2002  2002     1     8     2   188  3.51    81.5  48.6
    ##  9 London 09jan2002  2002     1     9     3   168  3.22    88.3  59.2
    ## 10 London 10jan2002  2002     1    10     4   194  5.32    85.4  48.7
    ## # ℹ 1,816 more rows

``` r
cities <- c('valencia', 'london')
for (city in 1:length(list.cities))
  {
  
  # Load raw datset.
  data <- list.cities[[city]]
  
  # Generate lags for temperature, humidity, and PM10.
  for (num in 0:3)
  {
    tmean <- paste("tmean_l", num, sep="")
    assign(tmean, dplyr::lag(data$tmean, num))
    rh <- paste("rh_l", num, sep="")
    assign(rh, dplyr::lag(data$rh, num))
    pm10 <- paste("pm10_l", num, sep="")
    assign(pm10, dplyr::lag(data$pm10, num))
  }
  
  # Generate splines for temperature adjustment.
  data$l03tmean <- (tmean_l0+tmean_l1+tmean_l2+tmean_l3)/4
  data$ns.tmean <- ns(data$l03tmean, df=6)
  
  # Generate splines for humidity adjustment.
  data$l03rh <- (rh_l0+rh_l1+rh_l2+rh_l3)/4
  data$ns.rh <- ns(data$l03rh, df=3)
  
  # Generate averaged lag exposure for PM10 sale for a 10 ug/m3 increase.
  data$l01pm10 <- (pm10_l0+pm10_l1)/2
  data$l01pm10 <- data$l01pm10/10

    # Save R data file.
    save(data, file = paste0(here(), "/note_original/", cities[city],".RData"))
}

#######################################################################################
# Example 1. Exposure-outcome association and environmental time-varying confounders 
# adjustment using the Valencia dataset.
#######################################################################################

# Load Valencia dataset.
load(here("note_original/valencia.RData"))

names(data)
```

    ##  [1] "city"     "date"     "year"     "month"    "day"      "dow"      "all"      "allage1" 
    ##  [9] "allage2"  "tmean"    "rh"       "pm10"     "l03tmean" "ns.tmean" "l03rh"    "ns.rh"   
    ## [17] "l01pm10"

``` r
# Generate time-stratified strata.
data$month <- as.factor(data$month)
data$year  <- as.factor(data$year)
data$dow   <- as.factor(data$dow)
data$stratum <- with(data, as.factor(year:month:dow))
data$ind <- tapply(data$all, data$stratum, sum)[data$stratum]

names(data)
```

    ##  [1] "city"     "date"     "year"     "month"    "day"      "dow"      "all"      "allage1" 
    ##  [9] "allage2"  "tmean"    "rh"       "pm10"     "l03tmean" "ns.tmean" "l03rh"    "ns.rh"   
    ## [17] "l01pm10"  "stratum"  "ind"

``` r
data %>% 
  tibble() %>% 
  select(date, year, month, dow, stratum, all, ind) %>% 
  print(n = 32)
```

    ## # A tibble: 1,826 × 7
    ##    date      year  month dow   stratum    all       ind
    ##    <chr>     <fct> <fct> <fct> <fct>    <int> <int[1d]>
    ##  1 01jan2002 2002  1     2     2002:1:2    31       116
    ##  2 02jan2002 2002  1     3     2002:1:3    29       113
    ##  3 03jan2002 2002  1     4     2002:1:4    23       109
    ##  4 04jan2002 2002  1     5     2002:1:5    17        82
    ##  5 05jan2002 2002  1     6     2002:1:6    19        93
    ##  6 06jan2002 2002  1     0     2002:1:0    26        97
    ##  7 07jan2002 2002  1     1     2002:1:1    25        80
    ##  8 08jan2002 2002  1     2     2002:1:2    21       116
    ##  9 09jan2002 2002  1     3     2002:1:3    15       113
    ## 10 10jan2002 2002  1     4     2002:1:4    19       109
    ## 11 11jan2002 2002  1     5     2002:1:5    24        82
    ## 12 12jan2002 2002  1     6     2002:1:6    21        93
    ## 13 13jan2002 2002  1     0     2002:1:0    24        97
    ## 14 14jan2002 2002  1     1     2002:1:1    17        80
    ## 15 15jan2002 2002  1     2     2002:1:2    29       116
    ## 16 16jan2002 2002  1     3     2002:1:3    16       113
    ## 17 17jan2002 2002  1     4     2002:1:4    18       109
    ## 18 18jan2002 2002  1     5     2002:1:5    18        82
    ## 19 19jan2002 2002  1     6     2002:1:6    28        93
    ## 20 20jan2002 2002  1     0     2002:1:0    24        97
    ## 21 21jan2002 2002  1     1     2002:1:1    15        80
    ## 22 22jan2002 2002  1     2     2002:1:2    16       116
    ## 23 23jan2002 2002  1     3     2002:1:3    28       113
    ## 24 24jan2002 2002  1     4     2002:1:4    26       109
    ## 25 25jan2002 2002  1     5     2002:1:5    23        82
    ## 26 26jan2002 2002  1     6     2002:1:6    25        93
    ## 27 27jan2002 2002  1     0     2002:1:0    23        97
    ## 28 28jan2002 2002  1     1     2002:1:1    23        80
    ## 29 29jan2002 2002  1     2     2002:1:2    19       116
    ## 30 30jan2002 2002  1     3     2002:1:3    25       113
    ## 31 31jan2002 2002  1     4     2002:1:4    23       109
    ## 32 01feb2002 2002  2     5     2002:2:5    20        83
    ## # ℹ 1,794 more rows

``` r
# Fit fixed-effects conditional quasi-Poisson regression.
model.cc <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                data=data, family=quasipoisson, subset=ind>0, eliminate=stratum)

summary(model.cc)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = data, subset = ind > 0)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -3.011741  -0.582200  -0.006527   0.554747   3.265350  
    ## 
    ## Coefficients of interest:
    ##            Estimate Std. Error t value Pr(>|t|)   
    ## ns.tmean1 -0.140549   0.062431  -2.251  0.02452 * 
    ## ns.tmean2 -0.233571   0.083941  -2.783  0.00547 **
    ## ns.tmean3 -0.091536   0.086605  -1.057  0.29072   
    ## ns.tmean4 -0.204639   0.080176  -2.552  0.01081 * 
    ## ns.tmean5 -0.194026   0.160232  -1.211  0.22614   
    ## ns.tmean6  0.195115   0.105515   1.849  0.06465 . 
    ## ns.rh1    -0.033250   0.032655  -1.018  0.30875   
    ## ns.rh2    -0.129567   0.119167  -1.087  0.27710   
    ## ns.rh3     0.001596   0.057949   0.028  0.97803   
    ## l01pm10    0.022684   0.010073   2.252  0.02448 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.9727476)
    ## 
    ## Residual deviance: 1376.1 on 1393 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
# Get Relative Risk for PM10.
Epi::ci.exp(model.cc, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002945 1.043341

``` r
#######################################################################################
# Adjustment of subpopulation time-invariant covariates.
#######################################################################################

# Reshape to long format.
long.all <- reshape2::melt(dplyr::select(data, date, allage1, allage2),
                           id.vars="date", variable.name="age", value.name="allage")

long.all %>% tibble()
```

    ## # A tibble: 3,652 × 3
    ##    date      age     allage
    ##    <chr>     <fct>    <int>
    ##  1 01jan2002 allage1      5
    ##  2 02jan2002 allage1      3
    ##  3 03jan2002 allage1      4
    ##  4 04jan2002 allage1      2
    ##  5 05jan2002 allage1      6
    ##  6 06jan2002 allage1      4
    ##  7 07jan2002 allage1      6
    ##  8 08jan2002 allage1      0
    ##  9 09jan2002 allage1      3
    ## 10 10jan2002 allage1      1
    ## # ℹ 3,642 more rows

``` r
long <- merge(long.all, dplyr::select(data, -c(allage1, allage2)), by="date", all.x=T)

long %>% names()
```

    ##  [1] "date"     "age"      "allage"   "city"     "year"     "month"    "day"      "dow"     
    ##  [9] "all"      "tmean"    "rh"       "pm10"     "l03tmean" "ns.tmean" "l03rh"    "ns.rh"   
    ## [17] "l01pm10"  "stratum"  "ind"

``` r
long %>% 
  tibble() %>% 
  select(date, age, allage, all, stratum, ind)
```

    ## # A tibble: 3,652 × 6
    ##    date      age     allage   all stratum        ind
    ##    <chr>     <fct>    <int> <int> <fct>    <int[1d]>
    ##  1 01apr2002 allage2     11    12 2002:4:1        61
    ##  2 01apr2002 allage1      1    12 2002:4:1        61
    ##  3 01apr2003 allage2     16    18 2003:4:2        86
    ##  4 01apr2003 allage1      2    18 2003:4:2        86
    ##  5 01apr2004 allage2     22    26 2004:4:4       106
    ##  6 01apr2004 allage1      4    26 2004:4:4       106
    ##  7 01apr2005 allage2     17    18 2005:4:5        93
    ##  8 01apr2005 allage1      1    18 2005:4:5        93
    ##  9 01apr2006 allage2      8     9 2006:4:6        83
    ## 10 01apr2006 allage1      1     9 2006:4:6        83
    ## # ℹ 3,642 more rows

``` r
# Fit fixed-effects conditional quasi-Poisson regression adjusted by age.
model.cc.adj.age <- gnm(allage ~ factor(age) + ns.tmean + ns.rh + l01pm10, 
                        data=long, family=quasipoisson, eliminate=stratum)

summary(model.cc.adj.age)
```

    ## 
    ## Call:
    ## gnm(formula = allage ~ factor(age) + ns.tmean + ns.rh + l01pm10, 
    ##     eliminate = stratum, family = quasipoisson, data = long)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.90129  -0.71047  -0.06278   0.56913   3.48587  
    ## 
    ## Coefficients of interest:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## factor(age)allage2  1.600814   0.015132 105.787  < 2e-16 ***
    ## ns.tmean1          -0.140549   0.062834  -2.237  0.02537 *  
    ## ns.tmean2          -0.233571   0.084483  -2.765  0.00573 ** 
    ## ns.tmean3          -0.091536   0.087164  -1.050  0.29372    
    ## ns.tmean4          -0.204639   0.080693  -2.536  0.01126 *  
    ## ns.tmean5          -0.194026   0.161266  -1.203  0.22901    
    ## ns.tmean6           0.195115   0.106196   1.837  0.06626 .  
    ## ns.rh1             -0.033250   0.032866  -1.012  0.31177    
    ## ns.rh2             -0.129567   0.119936  -1.080  0.28009    
    ## ns.rh3              0.001596   0.058323   0.027  0.97817    
    ## l01pm10             0.022684   0.010138   2.237  0.02533 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.9853419)
    ## 
    ## Residual deviance: 3330.6 on 3215 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 4

``` r
Epi::ci.exp(model.cc.adj.age, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002817 1.043473

``` r
# Generate age-time-stratified strata.
long$stratum4 <- with(long, as.factor(year:month:dow:age))
long$ind4 <- tapply(long$all, long$stratum4, sum)[long$stratum4]

long %>% tibble() %>% select(year, month, dow, stratum4, ind, ind4)
```

    ## # A tibble: 3,652 × 6
    ##    year  month dow   stratum4               ind      ind4
    ##    <fct> <fct> <fct> <fct>            <int[1d]> <int[1d]>
    ##  1 2002  4     1     2002:4:1:allage2        61        61
    ##  2 2002  4     1     2002:4:1:allage1        61        61
    ##  3 2003  4     2     2003:4:2:allage2        86        86
    ##  4 2003  4     2     2003:4:2:allage1        86        86
    ##  5 2004  4     4     2004:4:4:allage2       106       106
    ##  6 2004  4     4     2004:4:4:allage1       106       106
    ##  7 2005  4     5     2005:4:5:allage2        93        93
    ##  8 2005  4     5     2005:4:5:allage1        93        93
    ##  9 2006  4     6     2006:4:6:allage2        83        83
    ## 10 2006  4     6     2006:4:6:allage1        83        83
    ## # ℹ 3,642 more rows

``` r
# Fit conditional quasi-Poisson with age-time-stratified strata.
model.cc.str4 <- gnm(allage ~ ns.tmean + ns.rh + l01pm10, 
                     data=long, family=quasipoisson, subset=ind4>0, eliminate=stratum4)

summary(model.cc.str4)
```

    ## 
    ## Call:
    ## gnm(formula = allage ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum4, 
    ##     family = quasipoisson, data = long, subset = ind4 > 0)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.94222  -0.63605  -0.01766   0.55059   2.83927  
    ## 
    ## Coefficients of interest:
    ##            Estimate Std. Error t value Pr(>|t|)   
    ## ns.tmean1 -0.140549   0.061875  -2.272  0.02319 * 
    ## ns.tmean2 -0.233571   0.083193  -2.808  0.00503 **
    ## ns.tmean3 -0.091536   0.085834  -1.066  0.28632   
    ## ns.tmean4 -0.204639   0.079462  -2.575  0.01007 * 
    ## ns.tmean5 -0.194026   0.158805  -1.222  0.22189   
    ## ns.tmean6  0.195115   0.104575   1.866  0.06218 . 
    ## ns.rh1    -0.033250   0.032364  -1.027  0.30434   
    ## ns.rh2    -0.129567   0.118106  -1.097  0.27272   
    ## ns.rh3     0.001596   0.057433   0.028  0.97783   
    ## l01pm10    0.022684   0.009984   2.272  0.02315 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.955499)
    ## 
    ## Residual deviance: 2850.3 on 2796 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
Epi::ci.exp(model.cc.str4, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.003121 1.043157

``` r
#######################################################################################
# Investigation of effect modification by age.
#######################################################################################

# Stratified analysis by age. 
model.cc.age1 <- gnm(allage1 ~ ns.tmean + ns.rh + l01pm10, 
                     data=data, family=quasipoisson, subset=ind>0, eliminate=stratum)

summary(model.cc.age1)
```

    ## 
    ## Call:
    ## gnm(formula = allage1 ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = data, subset = ind > 0)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.88196  -0.65644  -0.02824   0.52184   2.85657  
    ## 
    ## Coefficients of interest:
    ##           Estimate Std. Error t value Pr(>|t|)  
    ## ns.tmean1 -0.08195    0.15652  -0.524   0.6006  
    ## ns.tmean2 -0.09791    0.20824  -0.470   0.6383  
    ## ns.tmean3 -0.05630    0.21130  -0.266   0.7900  
    ## ns.tmean4 -0.13930    0.19324  -0.721   0.4711  
    ## ns.tmean5 -0.19706    0.39791  -0.495   0.6205  
    ## ns.tmean6 -0.13649    0.25413  -0.537   0.5913  
    ## ns.rh1    -0.07536    0.07939  -0.949   0.3426  
    ## ns.rh2    -0.58427    0.29195  -2.001   0.0456 *
    ## ns.rh3    -0.18030    0.14458  -1.247   0.2126  
    ## l01pm10    0.01261    0.02460   0.513   0.6081  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.9525685)
    ## 
    ## Residual deviance: 1479.3 on 1393 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
model.cc.age2 <- gnm(allage2 ~ ns.tmean + ns.rh + l01pm10, 
                     data=data, family=quasipoisson, subset=ind>0, eliminate=stratum)

summary(model.cc.age2)
```

    ## 
    ## Call:
    ## gnm(formula = allage2 ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = data, subset = ind > 0)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.94292  -0.60608  -0.01156   0.57165   2.77253  
    ## 
    ## Coefficients of interest:
    ##           Estimate Std. Error t value Pr(>|t|)   
    ## ns.tmean1 -0.15083    0.06755  -2.233  0.02572 * 
    ## ns.tmean2 -0.25856    0.09102  -2.841  0.00456 **
    ## ns.tmean3 -0.09677    0.09424  -1.027  0.30466   
    ## ns.tmean4 -0.21551    0.08744  -2.465  0.01383 * 
    ## ns.tmean5 -0.18934    0.17372  -1.090  0.27594   
    ## ns.tmean6  0.26215    0.11507   2.278  0.02286 * 
    ## ns.rh1    -0.02548    0.03554  -0.717  0.47350   
    ## ns.rh2    -0.04253    0.12953  -0.328  0.74272   
    ## ns.rh3     0.03610    0.06276   0.575  0.56524   
    ## l01pm10    0.02474    0.01095   2.259  0.02406 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.9610814)
    ## 
    ## Residual deviance: 1365.3 on 1393 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
Epi::ci.exp(model.cc.age1, subset="l01pm10") 
```

    ##         exp(Est.)      2.5%    97.5%
    ## l01pm10  1.012693 0.9650332 1.062707

``` r
Epi::ci.exp(model.cc.age2, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.025051 1.003277 1.047297

``` r
# Interaction analysis by age. 
model.cc.age.int <- gnm(allage ~ ns.tmean + ns.rh + factor(age) + l01pm10:factor(age), 
                        data=long, family=quasipoisson, subset=ind4>0, eliminate=stratum4)

summary(model.cc.age.int)
```

    ## 
    ## Call:
    ## gnm(formula = allage ~ ns.tmean + ns.rh + factor(age) + l01pm10:factor(age), 
    ##     eliminate = stratum4, family = quasipoisson, data = long,     subset = ind4 > 0)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.9405  -0.6365  -0.0152   0.5525   2.8227  
    ## 
    ## Coefficients of interest:
    ##                             Estimate Std. Error t value Pr(>|t|)   
    ## ns.tmean1                  -0.140799   0.061886  -2.275  0.02297 * 
    ## ns.tmean2                  -0.233523   0.083204  -2.807  0.00504 **
    ## ns.tmean3                  -0.091537   0.085845  -1.066  0.28638   
    ## ns.tmean4                  -0.204634   0.079472  -2.575  0.01008 * 
    ## ns.tmean5                  -0.194203   0.158827  -1.223  0.22153   
    ## ns.tmean6                   0.195196   0.104589   1.866  0.06210 . 
    ## ns.rh1                     -0.033242   0.032368  -1.027  0.30451   
    ## ns.rh2                     -0.129200   0.118124  -1.094  0.27415   
    ## ns.rh3                      0.001885   0.057443   0.033  0.97383   
    ## factor(age)allage2          0.000000         NA      NA       NA   
    ## factor(age)allage1:l01pm10  0.012558   0.023618   0.532  0.59497   
    ## factor(age)allage2:l01pm10  0.024681   0.010837   2.277  0.02284 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.9557409)
    ## 
    ## Std. Error is NA where coefficient has been constrained or is unidentified
    ## 
    ## Residual deviance: 2850.1 on 2795 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
Epi::ci.exp(model.cc.age.int, subset="l01pm10") 
```

    ##                            exp(Est.)      2.5%    97.5%
    ## factor(age)allage1:l01pm10  1.012637 0.9668299 1.060614
    ## factor(age)allage2:l01pm10  1.024988 1.0034458 1.046992

``` r
# Likelihood ratio-test for effect for interaction.
anova(model.cc.str4, model.cc.age.int, test="LRT")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: allage ~ ns.tmean + ns.rh + l01pm10 - 1
    ## Model 2: allage ~ ns.tmean + ns.rh + factor(age) + factor(age):l01pm10 - 
    ##     1
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
    ## 1      2796     2850.3                     
    ## 2      2795     2850.1  1  0.21419   0.6359

``` r
#######################################################################################
# Example 2. Multi-location study.
#######################################################################################

# Load Valencia R data file.
load(here("note_original/valencia.RData"))
data.vlc <- dplyr::select(data, -allage1,-allage2)

# Load London R data file.
load(here("note_original/london.RData"))
data.ldn <- data

# Append datasets.
dat.2city <- rbind(data.vlc, data.ldn)

dat.2city %>% tibble()
```

    ## # A tibble: 3,652 × 15
    ##    city     date    year month   day   dow   all tmean    rh  pm10 l03tmean ns.tmean[,"1"] l03rh
    ##    <chr>    <chr>  <int> <int> <int> <int> <int> <dbl> <dbl> <dbl>    <dbl>          <dbl> <dbl>
    ##  1 Valencia 01jan…  2002     1     1     2    31  8.90  85    17.8    NA            NA      NA  
    ##  2 Valencia 02jan…  2002     1     2     3    29  8.5   89.4  21.7    NA            NA      NA  
    ##  3 Valencia 03jan…  2002     1     3     4    23 10.8   86.8  43.2    NA            NA      NA  
    ##  4 Valencia 04jan…  2002     1     4     5    17  9.60  88.1  19.8     9.45          0.177  87.3
    ##  5 Valencia 05jan…  2002     1     5     6    19  9.20  88.8  15.9     9.53          0.184  88.3
    ##  6 Valencia 06jan…  2002     1     6     0    26  9.70  87.2  20.5     9.82          0.215  87.7
    ##  7 Valencia 07jan…  2002     1     7     1    25  8.80  72.4  48.9     9.33          0.165  84.1
    ##  8 Valencia 08jan…  2002     1     8     2    21  9.40  78.5  26       9.27          0.161  81.7
    ##  9 Valencia 09jan…  2002     1     9     3    15  8.80  83.6  18.7     9.17          0.152  80.4
    ## 10 Valencia 10jan…  2002     1    10     4    19  9.60  76.2  39.6     9.15          0.150  77.7
    ## # ℹ 3,642 more rows
    ## # ℹ 3 more variables: ns.tmean[2:6] <dbl>, ns.rh <dbl[,3]>, l01pm10 <dbl>

``` r
# Generate space-time-stratified strata.
dat.2city$city  <- as.factor(dat.2city$city)
dat.2city$month <- as.factor(dat.2city$month)
dat.2city$year  <- as.factor(dat.2city$year)
dat.2city$dow   <- as.factor(dat.2city$dow)
dat.2city$stratum <- with(dat.2city, as.factor(city:year:month:dow))
dat.2city$ind <- tapply(dat.2city$all, dat.2city$stratum, sum)[dat.2city$stratum]

names(dat.2city)
```

    ##  [1] "city"     "date"     "year"     "month"    "day"      "dow"      "all"      "tmean"   
    ##  [9] "rh"       "pm10"     "l03tmean" "ns.tmean" "l03rh"    "ns.rh"    "l01pm10"  "stratum" 
    ## [17] "ind"

``` r
dat.2city %>% tibble() %>% select(city, year, month, dow, stratum, ind)
```

    ## # A tibble: 3,652 × 6
    ##    city     year  month dow   stratum                 ind
    ##    <fct>    <fct> <fct> <fct> <fct>             <int[1d]>
    ##  1 Valencia 2002  1     2     Valencia:2002:1:2       116
    ##  2 Valencia 2002  1     3     Valencia:2002:1:3       113
    ##  3 Valencia 2002  1     4     Valencia:2002:1:4       109
    ##  4 Valencia 2002  1     5     Valencia:2002:1:5        82
    ##  5 Valencia 2002  1     6     Valencia:2002:1:6        93
    ##  6 Valencia 2002  1     0     Valencia:2002:1:0        97
    ##  7 Valencia 2002  1     1     Valencia:2002:1:1        80
    ##  8 Valencia 2002  1     2     Valencia:2002:1:2       116
    ##  9 Valencia 2002  1     3     Valencia:2002:1:3       113
    ## 10 Valencia 2002  1     4     Valencia:2002:1:4       109
    ## # ℹ 3,642 more rows

``` r
dat.2city %>% tibble() %>% select(city, year, month, dow, stratum, ind) %>% tail(n = 10)
```

    ## # A tibble: 10 × 6
    ##    city   year  month dow   stratum                ind
    ##    <fct>  <fct> <fct> <fct> <fct>            <int[1d]>
    ##  1 London 2006  12    5     London:2006:12:5       746
    ##  2 London 2006  12    6     London:2006:12:6       696
    ##  3 London 2006  12    0     London:2006:12:0       582
    ##  4 London 2006  12    1     London:2006:12:1       565
    ##  5 London 2006  12    2     London:2006:12:2       593
    ##  6 London 2006  12    3     London:2006:12:3       570
    ##  7 London 2006  12    4     London:2006:12:4       549
    ##  8 London 2006  12    5     London:2006:12:5       746
    ##  9 London 2006  12    6     London:2006:12:6       696
    ## 10 London <NA>  <NA>  <NA>  <NA>                    NA

``` r
# Fit conditional Poisson with space-time-stratified strata.
model.cc.adj.city <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                         data=dat.2city, family=quasipoisson, subset=ind>0, eliminate=stratum)

summary(model.cc.adj.city)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = dat.2city, subset = ind > 0)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -3.08681  -0.61553  -0.01305   0.56604   3.59154  
    ## 
    ## Coefficients of interest:
    ##            Estimate Std. Error t value Pr(>|t|)    
    ## ns.tmean1 -0.075396   0.016057  -4.696 2.79e-06 ***
    ## ns.tmean2 -0.096497   0.022406  -4.307 1.71e-05 ***
    ## ns.tmean3 -0.091654   0.021974  -4.171 3.13e-05 ***
    ## ns.tmean4 -0.161615   0.021668  -7.459 1.16e-13 ***
    ## ns.tmean5  0.035013   0.043930   0.797 0.425517    
    ## ns.tmean6  0.337401   0.032908  10.253  < 2e-16 ***
    ## ns.rh1    -0.017364   0.013121  -1.323 0.185822    
    ## ns.rh2    -0.028924   0.042061  -0.688 0.491719    
    ## ns.rh3     0.018360   0.020235   0.907 0.364296    
    ## l01pm10    0.007139   0.001849   3.862 0.000115 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 1.057489)
    ## 
    ## Residual deviance: 2976.7 on 2795 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
Epi::ci.exp(model.cc.adj.city, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.007165 1.003522 1.010821

``` r
# Stratified analysis by city.
model.cc.vlc <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data=subset(dat.2city,city=="Valencia"), 
                    family=quasipoisson, subset=ind>0, eliminate=stratum)

summary(model.cc.vlc)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = subset(dat.2city, city == "Valencia"),     subset = ind > 0)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -3.011741  -0.582200  -0.006527   0.554747   3.265350  
    ## 
    ## Coefficients of interest:
    ##            Estimate Std. Error t value Pr(>|t|)   
    ## ns.tmean1 -0.140549   0.062431  -2.251  0.02452 * 
    ## ns.tmean2 -0.233571   0.083941  -2.783  0.00547 **
    ## ns.tmean3 -0.091536   0.086605  -1.057  0.29072   
    ## ns.tmean4 -0.204639   0.080176  -2.552  0.01081 * 
    ## ns.tmean5 -0.194026   0.160232  -1.211  0.22614   
    ## ns.tmean6  0.195115   0.105515   1.849  0.06465 . 
    ## ns.rh1    -0.033250   0.032655  -1.018  0.30875   
    ## ns.rh2    -0.129567   0.119167  -1.087  0.27710   
    ## ns.rh3     0.001596   0.057949   0.028  0.97803   
    ## l01pm10    0.022684   0.010073   2.252  0.02448 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.9727476)
    ## 
    ## Residual deviance: 1376.1 on 1393 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
model.cc.ldn <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data=subset(dat.2city,city=="London"), 
                    family=quasipoisson, subset=ind>0, eliminate=stratum)

summary(model.cc.ldn)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = subset(dat.2city, city == "London"),     subset = ind > 0)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -3.087601  -0.655932   0.003019   0.590973   3.606829  
    ## 
    ## Coefficients of interest:
    ##            Estimate Std. Error t value Pr(>|t|)    
    ## ns.tmean1 -0.071383   0.017236  -4.141 3.66e-05 ***
    ## ns.tmean2 -0.087034   0.024167  -3.601 0.000328 ***
    ## ns.tmean3 -0.091437   0.023591  -3.876 0.000111 ***
    ## ns.tmean4 -0.158039   0.023374  -6.761 2.00e-11 ***
    ## ns.tmean5  0.056185   0.047520   1.182 0.237270    
    ## ns.tmean6  0.357007   0.035927   9.937  < 2e-16 ***
    ## ns.rh1    -0.014247   0.015027  -0.948 0.343248    
    ## ns.rh2    -0.008091   0.046626  -0.174 0.862262    
    ## ns.rh3     0.023084   0.022533   1.024 0.305813    
    ## l01pm10    0.006390   0.001957   3.265 0.001120 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 1.139979)
    ## 
    ## Residual deviance: 1587.5 on 1392 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
Epi::ci.exp(model.cc.vlc, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002945 1.043341

``` r
Epi::ci.exp(model.cc.ldn, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10   1.00641 1.002558 1.010278

``` r
Epi::ci.exp(model.cc.adj.city, subset="l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.007165 1.003522 1.010821

``` r
# Interaction analysis by location.
model.cc.city.int <- gnm(all ~ ns.tmean + ns.rh + city + l01pm10:city, 
                        data=dat.2city, family=quasipoisson, subset=ind>0, eliminate=stratum)

summary(model.cc.city.int)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + city + l01pm10:city, eliminate = stratum, 
    ##     family = quasipoisson, data = dat.2city, subset = ind > 0)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -3.083227  -0.618923  -0.008919   0.565896   3.592595  
    ## 
    ## Coefficients of interest:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## ns.tmean1            -0.075573   0.016051  -4.708 2.62e-06 ***
    ## ns.tmean2            -0.096636   0.022399  -4.314 1.66e-05 ***
    ## ns.tmean3            -0.092180   0.021970  -4.196 2.81e-05 ***
    ## ns.tmean4            -0.162040   0.021664  -7.480 9.91e-14 ***
    ## ns.tmean5             0.034885   0.043917   0.794 0.427069    
    ## ns.tmean6             0.337616   0.032901  10.262  < 2e-16 ***
    ## ns.rh1               -0.017907   0.013122  -1.365 0.172448    
    ## ns.rh2               -0.028098   0.042052  -0.668 0.504083    
    ## ns.rh3                0.018794   0.020232   0.929 0.353000    
    ## cityValencia          0.000000         NA      NA       NA    
    ## cityLondon:l01pm10    0.006679   0.001872   3.567 0.000367 ***
    ## cityValencia:l01pm10  0.022315   0.010027   2.225 0.026135 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 1.056838)
    ## 
    ## Std. Error is NA where coefficient has been constrained or is unidentified
    ## 
    ## Residual deviance: 2974.2 on 2794 degrees of freedom
    ## AIC: NA
    ## 
    ## Number of iterations: 2

``` r
Epi::ci.exp(model.cc.city.int, subset="l01pm10") 
```

    ##                      exp(Est.)     2.5%    97.5%
    ## cityLondon:l01pm10    1.006701 1.003014 1.010402
    ## cityValencia:l01pm10  1.022566 1.002665 1.042861

``` r
# Likelihood ratio test for effect modification.
anova(model.cc.adj.city, model.cc.city.int, test="LRT")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: all ~ ns.tmean + ns.rh + l01pm10 - 1
    ## Model 2: all ~ ns.tmean + ns.rh + city + city:l01pm10 - 1
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
    ## 1      2795     2976.7                     
    ## 2      2794     2974.2  1   2.5022   0.1239

``` r
#######################################################################################
#######################################################################################
###                           End of script file                                    ###
#######################################################################################
#######################################################################################
```
