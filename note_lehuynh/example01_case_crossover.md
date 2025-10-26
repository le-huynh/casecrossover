Case-crossover analysis - Example 1
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

ls_wdf
```

    ## $valencia
    ## # A tibble: 1,826 × 31
    ##    city     date      year  month   day dow     all allage1 allage2 tmean    rh  pm10 stratum tmean_l0 rh_l0 pm10_l0
    ##    <chr>    <chr>     <fct> <fct> <int> <fct> <int>   <int>   <int> <dbl> <dbl> <dbl> <fct>      <dbl> <dbl>   <dbl>
    ##  1 Valencia 01jan2002 2002  1         1 2        31       5      26  8.90  85    17.8 2002:1…     8.90  85      17.8
    ##  2 Valencia 02jan2002 2002  1         2 3        29       3      26  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7
    ##  3 Valencia 03jan2002 2002  1         3 4        23       4      19 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2
    ##  4 Valencia 04jan2002 2002  1         4 5        17       2      15  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8
    ##  5 Valencia 05jan2002 2002  1         5 6        19       6      13  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9
    ##  6 Valencia 06jan2002 2002  1         6 0        26       4      22  9.70  87.2  20.5 2002:1…     9.70  87.2    20.5
    ##  7 Valencia 07jan2002 2002  1         7 1        25       6      19  8.80  72.4  48.9 2002:1…     8.80  72.4    48.9
    ##  8 Valencia 08jan2002 2002  1         8 2        21       0      21  9.40  78.5  26   2002:1…     9.40  78.5    26  
    ##  9 Valencia 09jan2002 2002  1         9 3        15       3      12  8.80  83.6  18.7 2002:1…     8.80  83.6    18.7
    ## 10 Valencia 10jan2002 2002  1        10 4        19       1      18  9.60  76.2  39.6 2002:1…     9.60  76.2    39.6
    ## # ℹ 1,816 more rows
    ## # ℹ 15 more variables: tmean_l1 <dbl>, rh_l1 <dbl>, pm10_l1 <dbl>, tmean_l2 <dbl>, rh_l2 <dbl>, pm10_l2 <dbl>,
    ## #   tmean_l3 <dbl>, rh_l3 <dbl>, pm10_l3 <dbl>, l03tmean <dbl>, l03rh <dbl>, l01pm10 <dbl>, ns.tmean <ns[,6]>,
    ## #   ns.rh <ns[,3]>, ind <int>
    ## 
    ## $london
    ## # A tibble: 1,826 × 29
    ##    city   date      year  month   day dow     all   tmean    rh  pm10 stratum  tmean_l0 rh_l0 pm10_l0 tmean_l1 rh_l1
    ##    <chr>  <chr>     <fct> <fct> <int> <fct> <int>   <dbl> <dbl> <dbl> <fct>       <dbl> <dbl>   <dbl>    <dbl> <dbl>
    ##  1 London 01jan2002 2002  1         1 2       199 -0.225   75.7  71.7 2002:1:2  -0.225   75.7    71.7  NA       NA  
    ##  2 London 02jan2002 2002  1         2 3       231  0.0875  77.5  40.2 2002:1:3   0.0875  77.5    40.2  -0.225   75.7
    ##  3 London 03jan2002 2002  1         3 4       210  0.850   81.3  41.8 2002:1:4   0.850   81.3    41.8   0.0875  77.5
    ##  4 London 04jan2002 2002  1         4 5       203  0.538   85.4  50.4 2002:1:5   0.538   85.4    50.4   0.850   81.3
    ##  5 London 05jan2002 2002  1         5 6       224  4.25    93.5  49.4 2002:1:6   4.25    93.5    49.4   0.538   85.4
    ##  6 London 06jan2002 2002  1         6 0       198  7.07    96.4  31.1 2002:1:0   7.07    96.4    31.1   4.25    93.5
    ##  7 London 07jan2002 2002  1         7 1       180  5.19    93.5  48.6 2002:1:1   5.19    93.5    48.6   7.07    96.4
    ##  8 London 08jan2002 2002  1         8 2       188  3.51    81.5  48.6 2002:1:2   3.51    81.5    48.6   5.19    93.5
    ##  9 London 09jan2002 2002  1         9 3       168  3.22    88.3  59.2 2002:1:3   3.22    88.3    59.2   3.51    81.5
    ## 10 London 10jan2002 2002  1        10 4       194  5.32    85.4  48.7 2002:1:4   5.32    85.4    48.7   3.22    88.3
    ## # ℹ 1,816 more rows
    ## # ℹ 13 more variables: pm10_l1 <dbl>, tmean_l2 <dbl>, rh_l2 <dbl>, pm10_l2 <dbl>, tmean_l3 <dbl>, rh_l3 <dbl>,
    ## #   pm10_l3 <dbl>, l03tmean <dbl>, l03rh <dbl>, l01pm10 <dbl>, ns.tmean <ns[,6]>, ns.rh <ns[,3]>, ind <int>

# Example 1

``` r
# example 1 #-----------
```

Exposure-outcome association and environmental time-varying confounders
adjustment

``` r
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
```

### Valencia

``` r
summary(ls_model1[[1]])
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
Epi::ci.exp(ls_model1[[1]], subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002945 1.043341

### London

``` r
summary(ls_model1[[2]])
```

    ## 
    ## Call:
    ## gnm(formula = all ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = data, subset = ind > 0)
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
# Get Relative Risk for PM10.
Epi::ci.exp(ls_model1[[2]], subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10   1.00641 1.002558 1.010278

# Adjustment of sub-population time-invariant covariates

``` r
# sub-population #-----------------------
# reshape data
(df_model2 <- ls_wdf[[1]] %>% 
    pivot_longer(cols = c(allage1, allage2),
                 names_to = "age",
                 values_to = "allage"))
```

    ## # A tibble: 3,652 × 31
    ##    city  date  year  month   day dow     all tmean    rh  pm10 stratum tmean_l0 rh_l0 pm10_l0 tmean_l1 rh_l1 pm10_l1
    ##    <chr> <chr> <fct> <fct> <int> <fct> <int> <dbl> <dbl> <dbl> <fct>      <dbl> <dbl>   <dbl>    <dbl> <dbl>   <dbl>
    ##  1 Vale… 01ja… 2002  1         1 2        31  8.90  85    17.8 2002:1…     8.90  85      17.8    NA     NA      NA  
    ##  2 Vale… 01ja… 2002  1         1 2        31  8.90  85    17.8 2002:1…     8.90  85      17.8    NA     NA      NA  
    ##  3 Vale… 02ja… 2002  1         2 3        29  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7     8.90  85      17.8
    ##  4 Vale… 02ja… 2002  1         2 3        29  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7     8.90  85      17.8
    ##  5 Vale… 03ja… 2002  1         3 4        23 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2     8.5   89.4    21.7
    ##  6 Vale… 03ja… 2002  1         3 4        23 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2     8.5   89.4    21.7
    ##  7 Vale… 04ja… 2002  1         4 5        17  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8    10.8   86.8    43.2
    ##  8 Vale… 04ja… 2002  1         4 5        17  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8    10.8   86.8    43.2
    ##  9 Vale… 05ja… 2002  1         5 6        19  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9     9.60  88.1    19.8
    ## 10 Vale… 05ja… 2002  1         5 6        19  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9     9.60  88.1    19.8
    ## # ℹ 3,642 more rows
    ## # ℹ 14 more variables: tmean_l2 <dbl>, rh_l2 <dbl>, pm10_l2 <dbl>, tmean_l3 <dbl>, rh_l3 <dbl>, pm10_l3 <dbl>,
    ## #   l03tmean <dbl>, l03rh <dbl>, l01pm10 <dbl>, ns.tmean <ns[,6]>, ns.rh <ns[,3]>, ind <int>, age <chr>,
    ## #   allage <int>

``` r
# Fit fixed-effects conditional quasi-Poisson regression adjusted by age
model.cc.adj.age <- gnm(allage ~ factor(age) + ns.tmean + ns.rh + l01pm10, 
                        data = df_model2,
                        family = quasipoisson,
                        eliminate = stratum)

summary(model.cc.adj.age)
```

    ## 
    ## Call:
    ## gnm(formula = allage ~ factor(age) + ns.tmean + ns.rh + l01pm10, 
    ##     eliminate = stratum, family = quasipoisson, data = df_model2)
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
Epi::ci.exp(model.cc.adj.age, subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002817 1.043473

Generate age-time-stratified strata

``` r
df_model2 %>% 
    mutate(age = as.factor(age),
           stratum4 = as.factor(year:month:dow:age)) %>% 
    select(year, month, dow, age, stratum4)
```

    ## # A tibble: 3,652 × 5
    ##    year  month dow   age     stratum4        
    ##    <fct> <fct> <fct> <fct>   <fct>           
    ##  1 2002  1     2     allage1 2002:1:2:allage1
    ##  2 2002  1     2     allage2 2002:1:2:allage2
    ##  3 2002  1     3     allage1 2002:1:3:allage1
    ##  4 2002  1     3     allage2 2002:1:3:allage2
    ##  5 2002  1     4     allage1 2002:1:4:allage1
    ##  6 2002  1     4     allage2 2002:1:4:allage2
    ##  7 2002  1     5     allage1 2002:1:5:allage1
    ##  8 2002  1     5     allage2 2002:1:5:allage2
    ##  9 2002  1     6     allage1 2002:1:6:allage1
    ## 10 2002  1     6     allage2 2002:1:6:allage2
    ## # ℹ 3,642 more rows

``` r
(df_model3 <- df_model2 %>% 
    mutate(age = as.factor(age),
           stratum4 = as.factor(year:month:dow:age)))
```

    ## # A tibble: 3,652 × 32
    ##    city  date  year  month   day dow     all tmean    rh  pm10 stratum tmean_l0 rh_l0 pm10_l0 tmean_l1 rh_l1 pm10_l1
    ##    <chr> <chr> <fct> <fct> <int> <fct> <int> <dbl> <dbl> <dbl> <fct>      <dbl> <dbl>   <dbl>    <dbl> <dbl>   <dbl>
    ##  1 Vale… 01ja… 2002  1         1 2        31  8.90  85    17.8 2002:1…     8.90  85      17.8    NA     NA      NA  
    ##  2 Vale… 01ja… 2002  1         1 2        31  8.90  85    17.8 2002:1…     8.90  85      17.8    NA     NA      NA  
    ##  3 Vale… 02ja… 2002  1         2 3        29  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7     8.90  85      17.8
    ##  4 Vale… 02ja… 2002  1         2 3        29  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7     8.90  85      17.8
    ##  5 Vale… 03ja… 2002  1         3 4        23 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2     8.5   89.4    21.7
    ##  6 Vale… 03ja… 2002  1         3 4        23 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2     8.5   89.4    21.7
    ##  7 Vale… 04ja… 2002  1         4 5        17  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8    10.8   86.8    43.2
    ##  8 Vale… 04ja… 2002  1         4 5        17  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8    10.8   86.8    43.2
    ##  9 Vale… 05ja… 2002  1         5 6        19  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9     9.60  88.1    19.8
    ## 10 Vale… 05ja… 2002  1         5 6        19  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9     9.60  88.1    19.8
    ## # ℹ 3,642 more rows
    ## # ℹ 15 more variables: tmean_l2 <dbl>, rh_l2 <dbl>, pm10_l2 <dbl>, tmean_l3 <dbl>, rh_l3 <dbl>, pm10_l3 <dbl>,
    ## #   l03tmean <dbl>, l03rh <dbl>, l01pm10 <dbl>, ns.tmean <ns[,6]>, ns.rh <ns[,3]>, ind <int>, age <fct>,
    ## #   allage <int>, stratum4 <fct>

``` r
(df_ind4 <- df_model3 %>% 
        group_by(stratum4) %>% 
        summarise(ind4 = sum(all)))
```

    ## # A tibble: 840 × 2
    ##    stratum4          ind4
    ##    <fct>            <int>
    ##  1 2002:1:0:allage1    97
    ##  2 2002:1:0:allage2    97
    ##  3 2002:1:1:allage1    80
    ##  4 2002:1:1:allage2    80
    ##  5 2002:1:2:allage1   116
    ##  6 2002:1:2:allage2   116
    ##  7 2002:1:3:allage1   113
    ##  8 2002:1:3:allage2   113
    ##  9 2002:1:4:allage1   109
    ## 10 2002:1:4:allage2   109
    ## # ℹ 830 more rows

``` r
(df_model3a <- df_model3 %>% 
    left_join(df_ind4,
              by = join_by(stratum4)))
```

    ## # A tibble: 3,652 × 33
    ##    city  date  year  month   day dow     all tmean    rh  pm10 stratum tmean_l0 rh_l0 pm10_l0 tmean_l1 rh_l1 pm10_l1
    ##    <chr> <chr> <fct> <fct> <int> <fct> <int> <dbl> <dbl> <dbl> <fct>      <dbl> <dbl>   <dbl>    <dbl> <dbl>   <dbl>
    ##  1 Vale… 01ja… 2002  1         1 2        31  8.90  85    17.8 2002:1…     8.90  85      17.8    NA     NA      NA  
    ##  2 Vale… 01ja… 2002  1         1 2        31  8.90  85    17.8 2002:1…     8.90  85      17.8    NA     NA      NA  
    ##  3 Vale… 02ja… 2002  1         2 3        29  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7     8.90  85      17.8
    ##  4 Vale… 02ja… 2002  1         2 3        29  8.5   89.4  21.7 2002:1…     8.5   89.4    21.7     8.90  85      17.8
    ##  5 Vale… 03ja… 2002  1         3 4        23 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2     8.5   89.4    21.7
    ##  6 Vale… 03ja… 2002  1         3 4        23 10.8   86.8  43.2 2002:1…    10.8   86.8    43.2     8.5   89.4    21.7
    ##  7 Vale… 04ja… 2002  1         4 5        17  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8    10.8   86.8    43.2
    ##  8 Vale… 04ja… 2002  1         4 5        17  9.60  88.1  19.8 2002:1…     9.60  88.1    19.8    10.8   86.8    43.2
    ##  9 Vale… 05ja… 2002  1         5 6        19  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9     9.60  88.1    19.8
    ## 10 Vale… 05ja… 2002  1         5 6        19  9.20  88.8  15.9 2002:1…     9.20  88.8    15.9     9.60  88.1    19.8
    ## # ℹ 3,642 more rows
    ## # ℹ 16 more variables: tmean_l2 <dbl>, rh_l2 <dbl>, pm10_l2 <dbl>, tmean_l3 <dbl>, rh_l3 <dbl>, pm10_l3 <dbl>,
    ## #   l03tmean <dbl>, l03rh <dbl>, l01pm10 <dbl>, ns.tmean <ns[,6]>, ns.rh <ns[,3]>, ind <int>, age <fct>,
    ## #   allage <int>, stratum4 <fct>, ind4 <int>

``` r
# Fit conditional quasi-Poisson with age-time-stratified strata.
model.cc.str4 <- gnm(allage ~ ns.tmean + ns.rh + l01pm10, 
                     data = df_model3a,
                     family = quasipoisson,
                     subset = ind4>0,
                     eliminate = stratum4)

summary(model.cc.str4)
```

    ## 
    ## Call:
    ## gnm(formula = allage ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum4, 
    ##     family = quasipoisson, data = df_model3a, subset = ind4 >         0)
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
Epi::ci.exp(model.cc.str4, subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.003121 1.043157

``` r
Epi::ci.exp(model.cc.adj.age, subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002817 1.043473

``` r
Epi::ci.exp(ls_model1[[1]], subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.022943 1.002945 1.043341

# Investigation of effect modification by age

``` r
# effect modification by age #---------------------
# age1 = <65 years
model.cc.age1 <- gnm(allage1 ~ ns.tmean + ns.rh + l01pm10, 
                     data = ls_wdf[["valencia"]], 
                     family = quasipoisson, 
                     subset = ind>0, 
                     eliminate = stratum)
summary(model.cc.age1)
```

    ## 
    ## Call:
    ## gnm(formula = allage1 ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = ls_wdf[["valencia"]], subset = ind >         0)
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
Epi::ci.exp(model.cc.age1, subset = "l01pm10") 
```

    ##         exp(Est.)      2.5%    97.5%
    ## l01pm10  1.012693 0.9650332 1.062707

``` r
# age2 = >65 years
model.cc.age2 <- gnm(allage2 ~ ns.tmean + ns.rh + l01pm10, 
                     data = ls_wdf[["valencia"]], 
                     family = quasipoisson, 
                     subset = ind>0, 
                     eliminate = stratum)
summary(model.cc.age2)
```

    ## 
    ## Call:
    ## gnm(formula = allage2 ~ ns.tmean + ns.rh + l01pm10, eliminate = stratum, 
    ##     family = quasipoisson, data = ls_wdf[["valencia"]], subset = ind >         0)
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
Epi::ci.exp(model.cc.age2, subset = "l01pm10") 
```

    ##         exp(Est.)     2.5%    97.5%
    ## l01pm10  1.025051 1.003277 1.047297

## Interaction analysis by age

``` r
model.cc.age.int <- gnm(allage ~ ns.tmean + ns.rh + factor(age) + l01pm10:factor(age), 
                        data = df_model3a, 
                        family = quasipoisson, 
                        subset = ind4>0, 
                        eliminate = stratum4)
summary(model.cc.age.int)
```

    ## 
    ## Call:
    ## gnm(formula = allage ~ ns.tmean + ns.rh + factor(age) + l01pm10:factor(age), 
    ##     eliminate = stratum4, family = quasipoisson, data = df_model3a,     subset = ind4 > 0)
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
Epi::ci.exp(model.cc.age.int, subset = "l01pm10") 
```

    ##                            exp(Est.)      2.5%    97.5%
    ## factor(age)allage1:l01pm10  1.012637 0.9668299 1.060614
    ## factor(age)allage2:l01pm10  1.024988 1.0034458 1.046992

``` r
# Likelihood ratio-test for effect for interaction.
anova(model.cc.str4, model.cc.age.int, test = "LRT")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: allage ~ ns.tmean + ns.rh + l01pm10 - 1
    ## Model 2: allage ~ ns.tmean + ns.rh + factor(age) + factor(age):l01pm10 - 
    ##     1
    ##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
    ## 1      2796     2850.3                     
    ## 2      2795     2850.1  1  0.21419   0.6359
