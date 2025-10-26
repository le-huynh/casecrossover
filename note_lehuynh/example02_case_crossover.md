Case-crossover analysis - Example 2
================

``` r
pacman::p_load(
    rio,            # import and export files
    here,           # locate files 
    tidyverse,      # data management and visualization
    skimr,
    splines,
    gnm,            # Generalized Nonlinear Models
    Epi             # Statistical Analysis in Epidemiology
)
```

# Data

``` r
# data #-----------
ls_wdf <- readRDS(here("note_lehuynh/working_df.RDS"))

(dat.2city <- rbind(ls_wdf[[1]] %>% select(-allage1, -allage2, -ind),
                   ls_wdf[[2]] %>% select(-ind)))
```

    ## # A tibble: 3,652 × 28
    ##    city     date  year  month   day dow     all tmean    rh  pm10 stratum tmean_l0 rh_l0 pm10_l0
    ##    <chr>    <chr> <fct> <fct> <int> <fct> <int> <dbl> <dbl> <dbl> <fct>      <dbl> <dbl>   <dbl>
    ##  1 Valencia 01ja… 2002  1         1 2        31  8.90  85    17.8 2002:1…     8.90  85      17.8
    ##  2 Valencia 02ja… 2002  1         2 3        29  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7
    ##  3 Valencia 03ja… 2002  1         3 4        23 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2
    ##  4 Valencia 04ja… 2002  1         4 5        17  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8
    ##  5 Valencia 05ja… 2002  1         5 6        19  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9
    ##  6 Valencia 06ja… 2002  1         6 0        26  9.70  87.2  20.5 2002:1…     9.70  87.2    20.5
    ##  7 Valencia 07ja… 2002  1         7 1        25  8.80  72.4  48.9 2002:1…     8.80  72.4    48.9
    ##  8 Valencia 08ja… 2002  1         8 2        21  9.40  78.5  26   2002:1…     9.40  78.5    26  
    ##  9 Valencia 09ja… 2002  1         9 3        15  8.80  83.6  18.7 2002:1…     8.80  83.6    18.7
    ## 10 Valencia 10ja… 2002  1        10 4        19  9.60  76.2  39.6 2002:1…     9.60  76.2    39.6
    ## # ℹ 3,642 more rows
    ## # ℹ 14 more variables: tmean_l1 <dbl>, rh_l1 <dbl>, pm10_l1 <dbl>, tmean_l2 <dbl>, rh_l2 <dbl>,
    ## #   pm10_l2 <dbl>, tmean_l3 <dbl>, rh_l3 <dbl>, pm10_l3 <dbl>, l03tmean <dbl>, l03rh <dbl>,
    ## #   l01pm10 <dbl>, ns.tmean <dbl[,6]>, ns.rh <dbl[,3]>

``` r
names(dat.2city)
```

    ##  [1] "city"     "date"     "year"     "month"    "day"      "dow"      "all"      "tmean"   
    ##  [9] "rh"       "pm10"     "stratum"  "tmean_l0" "rh_l0"    "pm10_l0"  "tmean_l1" "rh_l1"   
    ## [17] "pm10_l1"  "tmean_l2" "rh_l2"    "pm10_l2"  "tmean_l3" "rh_l3"    "pm10_l3"  "l03tmean"
    ## [25] "l03rh"    "l01pm10"  "ns.tmean" "ns.rh"

``` r
(dat.2city <- dat.2city %>% 
    mutate(city = as.factor(city),
           stratum = as.factor(city:year:month:dow)))
```

    ## # A tibble: 3,652 × 28
    ##    city     date  year  month   day dow     all tmean    rh  pm10 stratum tmean_l0 rh_l0 pm10_l0
    ##    <fct>    <chr> <fct> <fct> <int> <fct> <int> <dbl> <dbl> <dbl> <fct>      <dbl> <dbl>   <dbl>
    ##  1 Valencia 01ja… 2002  1         1 2        31  8.90  85    17.8 Valenc…     8.90  85      17.8
    ##  2 Valencia 02ja… 2002  1         2 3        29  8.5   89.4  21.7 Valenc…     8.5   89.4    21.7
    ##  3 Valencia 03ja… 2002  1         3 4        23 10.8   86.8  43.2 Valenc…    10.8   86.8    43.2
    ##  4 Valencia 04ja… 2002  1         4 5        17  9.60  88.1  19.8 Valenc…     9.60  88.1    19.8
    ##  5 Valencia 05ja… 2002  1         5 6        19  9.20  88.8  15.9 Valenc…     9.20  88.8    15.9
    ##  6 Valencia 06ja… 2002  1         6 0        26  9.70  87.2  20.5 Valenc…     9.70  87.2    20.5
    ##  7 Valencia 07ja… 2002  1         7 1        25  8.80  72.4  48.9 Valenc…     8.80  72.4    48.9
    ##  8 Valencia 08ja… 2002  1         8 2        21  9.40  78.5  26   Valenc…     9.40  78.5    26  
    ##  9 Valencia 09ja… 2002  1         9 3        15  8.80  83.6  18.7 Valenc…     8.80  83.6    18.7
    ## 10 Valencia 10ja… 2002  1        10 4        19  9.60  76.2  39.6 Valenc…     9.60  76.2    39.6
    ## # ℹ 3,642 more rows
    ## # ℹ 14 more variables: tmean_l1 <dbl>, rh_l1 <dbl>, pm10_l1 <dbl>, tmean_l2 <dbl>, rh_l2 <dbl>,
    ## #   pm10_l2 <dbl>, tmean_l3 <dbl>, rh_l3 <dbl>, pm10_l3 <dbl>, l03tmean <dbl>, l03rh <dbl>,
    ## #   l01pm10 <dbl>, ns.tmean <dbl[,6]>, ns.rh <dbl[,3]>

``` r
(df_ind <- dat.2city %>% 
    group_by(stratum) %>% 
    summarise(ind = sum(all)))
```

    ## # A tibble: 841 × 2
    ##    stratum           ind
    ##    <fct>           <int>
    ##  1 London:2002:1:0   738
    ##  2 London:2002:1:1   757
    ##  3 London:2002:1:2   925
    ##  4 London:2002:1:3   933
    ##  5 London:2002:1:4   948
    ##  6 London:2002:1:5   775
    ##  7 London:2002:1:6   766
    ##  8 London:2002:2:0   650
    ##  9 London:2002:2:1   638
    ## 10 London:2002:2:2   689
    ## # ℹ 831 more rows

``` r
(wdf <- dat.2city %>% 
    left_join(df_ind,
              by = join_by(stratum)))
```

    ## # A tibble: 3,652 × 29
    ##    city     date  year  month   day dow     all tmean    rh  pm10 stratum tmean_l0 rh_l0 pm10_l0
    ##    <fct>    <chr> <fct> <fct> <int> <fct> <int> <dbl> <dbl> <dbl> <fct>      <dbl> <dbl>   <dbl>
    ##  1 Valencia 01ja… 2002  1         1 2        31  8.90  85    17.8 Valenc…     8.90  85      17.8
    ##  2 Valencia 02ja… 2002  1         2 3        29  8.5   89.4  21.7 Valenc…     8.5   89.4    21.7
    ##  3 Valencia 03ja… 2002  1         3 4        23 10.8   86.8  43.2 Valenc…    10.8   86.8    43.2
    ##  4 Valencia 04ja… 2002  1         4 5        17  9.60  88.1  19.8 Valenc…     9.60  88.1    19.8
    ##  5 Valencia 05ja… 2002  1         5 6        19  9.20  88.8  15.9 Valenc…     9.20  88.8    15.9
    ##  6 Valencia 06ja… 2002  1         6 0        26  9.70  87.2  20.5 Valenc…     9.70  87.2    20.5
    ##  7 Valencia 07ja… 2002  1         7 1        25  8.80  72.4  48.9 Valenc…     8.80  72.4    48.9
    ##  8 Valencia 08ja… 2002  1         8 2        21  9.40  78.5  26   Valenc…     9.40  78.5    26  
    ##  9 Valencia 09ja… 2002  1         9 3        15  8.80  83.6  18.7 Valenc…     8.80  83.6    18.7
    ## 10 Valencia 10ja… 2002  1        10 4        19  9.60  76.2  39.6 Valenc…     9.60  76.2    39.6
    ## # ℹ 3,642 more rows
    ## # ℹ 15 more variables: tmean_l1 <dbl>, rh_l1 <dbl>, pm10_l1 <dbl>, tmean_l2 <dbl>, rh_l2 <dbl>,
    ## #   pm10_l2 <dbl>, tmean_l3 <dbl>, rh_l3 <dbl>, pm10_l3 <dbl>, l03tmean <dbl>, l03rh <dbl>,
    ## #   l01pm10 <dbl>, ns.tmean <dbl[,6]>, ns.rh <dbl[,3]>, ind <int>

# Fit conditional Poisson with space-time-stratified strata

``` r
# space-time-stratified strata #-----------
model.cc.adj.city <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                         data = wdf, 
                         family = quasipoisson, 
                         subset = ind>0, 
                         eliminate = stratum)

summary(model.cc.adj.city)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = wdf, subset = ind > 0)
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
Epi::ci.exp(model.cc.adj.city, subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.007165 1.003522 1.010821

# Stratified analysis by city

``` r
# stratified analysis by city #-----------------------
```

## Valencia

``` r
## Valencia #-----------------------
model.cc.vlc <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data = wdf %>% filter(city == "Valencia"), 
                    family = quasipoisson, 
                    subset = ind>0, 
                    eliminate = stratum)

summary(model.cc.vlc)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = wdf %>% filter(city == "Valencia"),     subset = ind > 0)
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

## London

``` r
## London #-------------
model.cc.ldn <- gnm(all ~ ns.tmean + ns.rh + l01pm10, 
                    data = wdf %>% filter(city == "London"), 
                    family = quasipoisson, 
                    subset = ind>0, 
                    eliminate = stratum)

summary(model.cc.ldn)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = wdf %>% filter(city == "London"),     subset = ind > 0)
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

## Relative Risk

``` r
## relative risk #-----------------------
Epi::ci.exp(model.cc.vlc, subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002945 1.043341

``` r
Epi::ci.exp(model.cc.ldn, subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10   1.00641 1.002558 1.010278

# Interaction analysis by location

``` r
# interaction by location #--------------
model.cc.city.int <- gnm(all ~ ns.tmean + ns.rh + city + l01pm10:city, 
                         data = wdf, 
                         family = quasipoisson, 
                         subset = ind>0, 
                         eliminate = stratum)

summary(model.cc.city.int)
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + city + l01pm10:city, eliminate = stratum, 
    ##     family = quasipoisson, data = wdf, subset = ind > 0)
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
Epi::ci.exp(model.cc.city.int, subset = "l01pm10") 
```

    ##                      exp(Est.)     2.5%    97.5%
    ## cityLondon:l01pm10    1.006701 1.003014 1.010402
    ## cityValencia:l01pm10  1.022566 1.002665 1.042861

``` r
# Likelihood ratio test for effect modification.
anova(model.cc.adj.city, model.cc.city.int, test = "LRT")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: all ~ ns.tmean + ns.rh + l01pm10 - 1
    ## Model 2: all ~ ns.tmean + ns.rh + city + city:l01pm10 - 1
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
    ## 1      2795     2976.7                     
    ## 2      2794     2974.2  1   2.5022   0.1239
