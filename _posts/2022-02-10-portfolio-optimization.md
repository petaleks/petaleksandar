---
layout: post
title: "Course code in static and dynamic Portfolio Optimization"
categories: misc
---


This code was part of the exercise session I gave in the Course: *International Portfolio Management and Investment Analysis*, autumn 2021.

I did that course together with my supervisor: Prof. Andreas Stephan and my colleage Toni Duras.

In the exercise I replicate parts of the [Eun and Resnik (1988)](https://www.jstor.org/stable/2328331?seq=1#metadata_info_tab_contents) analyses 
of international diversification under flexible exchange rates.
They find that exchange rate risk is a largely non-diversifiable factor which adversely affects the risk-return performance of international portfolios.

I used more recent data and also to take the Swedish perspective instead of being an US investor.

## Read data

``` r
urlfile_ind <-'https://raw.githubusercontent.com/petaleks/petaleksandar/master/data/IFM_2021_indices.csv'
indices <-read.csv(urlfile_ind, sep=";", header = T)
urlfile_erates <-'https://raw.githubusercontent.com/petaleks/petaleksandar/master/data/IFM_2021_Exchange_rates.csv'
erates <-read.csv(urlfile_erates, sep=";", header = T)

head(indices)
```

    ##         Date  DJES50I DKKFXIN FTSE100  HNGKNGI ISTA100 KORCOMP  SPCOMP SWISSMI
    ## 1 2014-12-01 3232.910  757.73 6656.37 23367.45 1306.70 1965.22 2053.44 9146.18
    ## 2 2015-01-01 3146.430  744.44 6566.09 23605.04 1288.01 1915.59 2058.90 8983.37
    ## 3 2015-02-01 3370.113  806.44 6782.55 24484.74 1268.88 1952.68 2020.85 8429.20
    ## 4 2015-03-01 3591.092  879.19 6940.64 24887.44 1360.23 1996.81 2117.39 9055.69
    ## 5 2015-04-01 3714.892  965.33 6809.50 25082.75 1421.12 2028.45 2059.69 9137.26
    ## 6 2015-05-01 3615.586  963.77 6985.95 28133.00 1423.70 2127.17 2108.29 9077.12
    ##   TOKYOSE SDTB30D   MSWRLD SWEDOMX
    ## 1 1421.65   0.066 1729.888 1450.03
    ## 2 1407.51   0.080 1709.680 1464.55
    ## 3 1408.75   0.010 1694.348 1575.95
    ## 4 1524.97  -0.094 1776.646 1686.64
    ## 5 1528.99  -0.250 1738.150 1669.49
    ## 6 1585.61  -0.269 1787.401 1628.04

``` r
colnames(indices)
```

    ##  [1] "Date"    "DJES50I" "DKKFXIN" "FTSE100" "HNGKNGI" "ISTA100" "KORCOMP"
    ##  [8] "SPCOMP"  "SWISSMI" "TOKYOSE" "SDTB30D" "MSWRLD"  "SWEDOMX"

``` r
# choose the portfolio
select_indices = c("DKKFXIN", "FTSE100", "TOKYOSE", "SWISSMI","SWEDOMX" )
prices = indices[ , select_indices]
summary(prices)
```

    ##     DKKFXIN          FTSE100        TOKYOSE        SWISSMI         SWEDOMX    
    ##  Min.   : 744.4   Min.   :5455   Min.   :1254   Min.   : 7688   Min.   :1340  
    ##  1st Qu.: 955.2   1st Qu.:6516   1st Qu.:1511   1st Qu.: 8688   1st Qu.:1526  
    ##  Median : 999.5   Median :7020   Median :1608   Median : 9137   Median :1591  
    ##  Mean   :1080.0   Mean   :6884   Mean   :1612   Mean   : 9381   Mean   :1648  
    ##  3rd Qu.:1118.9   3rd Qu.:7322   3rd Qu.:1721   3rd Qu.: 9952   3rd Qu.:1678  
    ##  Max.   :1810.5   Max.   :7702   Max.   :1986   Max.   :12433   Max.   :2382

``` r
# choose currency
select_currency = c("SDDKKSP", "SDGBPSP", "SDJPYSP", "SDCHFSP")
erates2 = erates[, select_currency]

# Introduce risk-free, take the last one known
rf_serie = indices$SDTB30D/100  # choose the column from indices data frame
# !!!!divide by 100, as it is in percentages
rf_serie = rf_serie[-(1:2)] # erase first two observations to make it equal to returns
l=length(rf_serie)          # for numeric vectors use length() function
rf = rf_serie[l]/12         # convert annual rf to monthly
```

Calculate returns

``` r
# calculate returns on prices
n = dim(prices)[1]
rets = log(prices[2:n,] / prices[1:(n-1),] ) # log- returns
rets = rets[-1,] # remove first row
summary(rets*12)
```

    ##     DKKFXIN           FTSE100            TOKYOSE            SWISSMI        
    ##  Min.   :-1.3519   Min.   :-2.38678   Min.   :-1.47179   Min.   :-1.12545  
    ##  1st Qu.:-0.2448   1st Qu.:-0.23484   1st Qu.:-0.21864   1st Qu.:-0.14524  
    ##  Median : 0.1566   Median : 0.02756   Median : 0.10428   Median : 0.05789  
    ##  Mean   : 0.1228   Mean   : 0.01005   Mean   : 0.05103   Mean   : 0.03756  
    ##  3rd Qu.: 0.4729   3rd Qu.: 0.26795   3rd Qu.: 0.45070   3rd Qu.: 0.28854  
    ##  Max.   : 1.1216   Max.   : 1.45649   Max.   : 1.14125   Max.   : 0.87442  
    ##     SWEDOMX        
    ##  Min.   :-1.96083  
    ##  1st Qu.:-0.21450  
    ##  Median : 0.15686  
    ##  Mean   : 0.06389  
    ##  3rd Qu.: 0.42133  
    ##  Max.   : 1.21481

``` r
# calculate returns on exchange rates 
i = dim(erates2)[1]
fxrets = log(erates2[2:i,] / erates2[1:(i-1),] ) # log-returns

# create data frame where you place share and fx returns
# will need to calculate correlation and covariances out of it
new = data.frame(rets, fxrets)
# add column for swedish currency
fxrets$SEK = 0 

# convert in SEK returns
sekrets = (1+rets) * (1+fxrets) -1
returns = as.matrix(sekrets) # convert to matrix format
summary(returns)
```

    ##     DKKFXIN            FTSE100              TOKYOSE         
    ##  Min.   :-0.10637   Min.   :-0.1901881   Min.   :-0.107690  
    ##  1st Qu.:-0.01885   1st Qu.:-0.0222864   1st Qu.:-0.018019  
    ##  Median : 0.01716   Median : 0.0042134   Median : 0.010043  
    ##  Mean   : 0.01092   Mean   : 0.0004373   Mean   : 0.005824  
    ##  3rd Qu.: 0.04022   3rd Qu.: 0.0236718   3rd Qu.: 0.032702  
    ##  Max.   : 0.09997   Max.   : 0.1146010   Max.   : 0.081703  
    ##     SWISSMI             SWEDOMX         
    ##  Min.   :-0.106142   Min.   :-0.163402  
    ##  1st Qu.:-0.014815   1st Qu.:-0.017875  
    ##  Median : 0.012921   Median : 0.013071  
    ##  Mean   : 0.005044   Mean   : 0.005324  
    ##  3rd Qu.: 0.029103   3rd Qu.: 0.035111  
    ##  Max.   : 0.069936   Max.   : 0.101234

Plot of returns in foreign currency

``` r
# stack data over each other
datum = as.Date(c(indices$Date), "%Y-%m-%d") # change to date format, so it is easier for plotting

ret_M  = data.frame(date = datum[-(1:2)], rets*12)
plotdf <- melt(ret_M[,], id="date")
p = ggplot(plotdf, aes(x=date, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  ggtitle("Returns in foreign currency")+
  theme(axis.text.x = element_text(size=6, angle=90))+
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
```
![](https://github.com/petaleks/petaleksandar/blob/ee34ae0fd86ea66c0a31d36676e5ffc5b554cb45/data/figure-gfm/plot_ret_in_fx-1.png?raw=true)


``` r
namepic = paste0("Returns in foreign currency.jpg")
ggsave(namepic, plot = p)
```

    ## Saving 7 x 5 in image

Plot of returns in local currency

``` r
# stack data over each other
ret_M  = data.frame(date = datum[-(1:2)], returns*12)
plotdf <- melt(ret_M[,], id="date")
p = ggplot(plotdf, aes(x=date, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  ggtitle("Returns in local currency")+
  theme(axis.text.x = element_text(size=6, angle=90))+
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
```
![](https://github.com/petaleks/petaleksandar/blob/ee34ae0fd86ea66c0a31d36676e5ffc5b554cb45/data/figure-gfm/plot_ret_in_SEK-1.png?raw=true)


``` r
namepic = paste0("Returns in local currency.jpg")
ggsave(namepic, plot = p)
```

    ## Saving 7 x 5 in image

## Functions

``` r
# Constraint function
constraint1  <- function(W){
  cons = (sum(W)-1) ^2              # "x = 1" constraint
  return(cons)
}

# Minimum Variance portfolio 
Min_Var <- function(W, SIGMA){
  QFORM =  t(W) %*% SIGMA %*% W  
  goal = sqrt(QFORM * 12)  +  100*constraint1(W) 
  return(goal)
}

# Maximum Sharpe Portfolio
Max_sharp <- function(Wr, SIGMA, MU){
  rf = MU[length(MU)]
  preturn = (t(Wr) %*% MU)*12
  pvar =  t(Wr) %*% SIGMA %*% Wr 
  sharpe = (preturn-rf)/ sqrt(pvar*12)
  goal = -sharpe +  100 * constraint1(Wr)
  return(goal)
}

# Optimize CVar portfolio 
MCVar <- function(W, ret, alpha){
  time_portr = (ret %*% W)*12
  VaR = quantile(time_portr, probs =c(alpha)) # cut-off point for alpha confidence
  # calculate difference R-VaR only for those lower then VaR
  min_vec = sapply(time_portr-VaR, function(x) min(x,0))  
  cvar_vec =min_vec[min_vec <0]
  CVaR = VaR + sum(cvar_vec)/length(cvar_vec)
  goal = -CVaR + 100 * constraint1(W)
  return(goal)
}
```

## Calculate mean vector and covarience matrix

``` r
MU = apply(returns, 2 , mean) # monthly mean return
MU = as.matrix(MU)            # transform into matrix format
MU
```

    ##                 [,1]
    ## DKKFXIN 0.0109232240
    ## FTSE100 0.0004373459
    ## TOKYOSE 0.0058241996
    ## SWISSMI 0.0050437693
    ## SWEDOMX 0.0053239289

``` r
SIGMA   <- cov(as.matrix(returns)) # monthly sample covariance
SIGMA
```

    ##              DKKFXIN     FTSE100      TOKYOSE     SWISSMI     SWEDOMX
    ## DKKFXIN 0.0016432888 0.001018114 0.0009709699 0.001073115 0.001099622
    ## FTSE100 0.0010181141 0.002182548 0.0015171324 0.001214068 0.001826426
    ## TOKYOSE 0.0009709699 0.001517132 0.0017627917 0.001023805 0.001368267
    ## SWISSMI 0.0010731154 0.001214068 0.0010238053 0.001296902 0.001146451
    ## SWEDOMX 0.0010996224 0.001826426 0.0013682666 0.001146451 0.002224341

``` r
#COVAR   <- SIGMA * 12              # annual sample covariance
#sigma_d <- sqrt(diag(COVAR))
#sigma_a <- as.matrix(sigma_d * sqrt(12))

# start with equal weights for k-assets
k = dim(returns)[2] # number of assets
vec_eq = c(rep(1/k,k)) # create vector of equal weights
W <- as.matrix(vec_eq) # transform into matrix format
```

add risk free to return vector and create bordered covariance matrix

``` r
MUr = rbind(MU, rf)                   # monthly mean vector with risk-free
MUr
```

    ##                  [,1]
    ## DKKFXIN  0.0109232240
    ## FTSE100  0.0004373459
    ## TOKYOSE  0.0058241996
    ## SWISSMI  0.0050437693
    ## SWEDOMX  0.0053239289
    ## rf      -0.0006083333

``` r
zeroc <- matrix(c(rep(0,k)),k,1)      # column with k zeros
zeror <- matrix(c(rep(0,k+1)),1, k+1) # row with k+1 zeros
colnames(zeroc)[1] =c("rf")           # add column name for rf
rownames(zeror)[1] =c("rf")           # add row name for rf
SIGMAr = cbind(SIGMA, zeroc)          # bind mu and rf by columns
SIGMAr = rbind(SIGMAr, zeror)         # bind mu and rf by rows
SIGMAr                                # monthly bordered covariance matrix
```

    ##              DKKFXIN     FTSE100      TOKYOSE     SWISSMI     SWEDOMX rf
    ## DKKFXIN 0.0016432888 0.001018114 0.0009709699 0.001073115 0.001099622  0
    ## FTSE100 0.0010181141 0.002182548 0.0015171324 0.001214068 0.001826426  0
    ## TOKYOSE 0.0009709699 0.001517132 0.0017627917 0.001023805 0.001368267  0
    ## SWISSMI 0.0010731154 0.001214068 0.0010238053 0.001296902 0.001146451  0
    ## SWEDOMX 0.0010996224 0.001826426 0.0013682666 0.001146451 0.002224341  0
    ## rf      0.0000000000 0.000000000 0.0000000000 0.000000000 0.000000000  0

Calculate correlation

``` r
rcorr(as.matrix(new)) 
```

    ##         DKKFXIN FTSE100 TOKYOSE SWISSMI SWEDOMX SDDKKSP SDGBPSP SDJPYSP SDCHFSP
    ## DKKFXIN    1.00    0.51    0.51    0.66    0.61   -0.25    0.14   -0.13   -0.02
    ## FTSE100    0.51    1.00    0.61    0.64    0.72   -0.29   -0.06   -0.20   -0.17
    ## TOKYOSE    0.51    0.61    1.00    0.68    0.70   -0.24    0.32   -0.53   -0.29
    ## SWISSMI    0.66    0.64    0.68    1.00    0.69   -0.14    0.28   -0.32   -0.31
    ## SWEDOMX    0.61    0.72    0.70    0.69    1.00   -0.13    0.36   -0.18   -0.02
    ## SDDKKSP   -0.25   -0.29   -0.24   -0.14   -0.13    1.00    0.48    0.60    0.55
    ## SDGBPSP    0.14   -0.06    0.32    0.28    0.36    0.48    1.00    0.20    0.32
    ## SDJPYSP   -0.13   -0.20   -0.53   -0.32   -0.18    0.60    0.20    1.00    0.73
    ## SDCHFSP   -0.02   -0.17   -0.29   -0.31   -0.02    0.55    0.32    0.73    1.00
    ## 
    ## n= 81 
    ## 
    ## 
    ## P
    ##         DKKFXIN FTSE100 TOKYOSE SWISSMI SWEDOMX SDDKKSP SDGBPSP SDJPYSP SDCHFSP
    ## DKKFXIN         0.0000  0.0000  0.0000  0.0000  0.0226  0.1989  0.2374  0.8707 
    ## FTSE100 0.0000          0.0000  0.0000  0.0000  0.0095  0.5916  0.0736  0.1231 
    ## TOKYOSE 0.0000  0.0000          0.0000  0.0000  0.0341  0.0042  0.0000  0.0078 
    ## SWISSMI 0.0000  0.0000  0.0000          0.0000  0.2026  0.0112  0.0041  0.0046 
    ## SWEDOMX 0.0000  0.0000  0.0000  0.0000          0.2395  0.0008  0.1107  0.8384 
    ## SDDKKSP 0.0226  0.0095  0.0341  0.2026  0.2395          0.0000  0.0000  0.0000 
    ## SDGBPSP 0.1989  0.5916  0.0042  0.0112  0.0008  0.0000          0.0740  0.0031 
    ## SDJPYSP 0.2374  0.0736  0.0000  0.0041  0.1107  0.0000  0.0740          0.0000 
    ## SDCHFSP 0.8707  0.1231  0.0078  0.0046  0.8384  0.0000  0.0031  0.0000

``` r
# first matrix is correlation matrix
# second matrix are p-values for the significance of the correlation coeficents
```

Calculate total variances

``` r
cov_df = data.frame(cov(new)) # covariance matrix
write.csv(cov_df, "cov_df.csv")
# calculate total var
total_VAR = sum(cov(new)*12*((1/5)^2))
total_VAR
```

    ## [1] 0.01610799

``` r
# variance on the diagonals
vars = diag(cov(new)*12)*((1/5)^2)
vars
```

    ##      DKKFXIN      FTSE100      TOKYOSE      SWISSMI      SWEDOMX      SDDKKSP 
    ## 0.0008330430 0.0008154415 0.0011428215 0.0006208591 0.0010676835 0.0001036600 
    ##      SDGBPSP      SDJPYSP      SDCHFSP 
    ## 0.0002977352 0.0004309153 0.0002527578

``` r
total_vars =sum(vars)
total_vars/ total_VAR
```

    ## [1] 0.3454755

``` r
# covariances
covars = total_VAR - total_vars
covars
```

    ## [1] 0.01054308

``` r
# 0.01054308
covars/total_VAR
```

    ## [1] 0.6545245

``` r
# 0.6545245
```

Equally weighted portfolio

``` r
p_ret_e =  mean(returns %*% W)*12 
p_var_e = ( t(W) %*% SIGMA %*% W ) * 12  
p_stdv_e = sqrt(p_var_e)      
sharpe_EQ = (p_ret_e-rf)/ p_stdv_e 
p_ret_e  # "Annual return of equally weighted portfolio"
```

    ## [1] 0.06612592

``` r
p_var_e  # "Annual variance of equally weighted portfolio"
```

    ##            [,1]
    ## [1,] 0.01614039

``` r
p_stdv_e # annual st.dev of equally weighted portfolio
```

    ##           [,1]
    ## [1,] 0.1270448

``` r
sharpe_EQ # "Annual sharpe ratio of equally weighted portfolio"
```

    ##           [,1]
    ## [1,] 0.5252811

# Static optimization

## Minimum variance

### Static optimization with SS (short sale allowed)

``` r
W <- as.matrix(vec_eq)
MinVar_SS_results <- optim( par=W, Min_Var, SIGMA=SIGMA)
optw_MV1 = MinVar_SS_results$par
# sum(optw_MV1) # check if they add up to approximately 1
p_ret_opt = (t(optw_MV1) %*% MU)*12               
p_var_opt = (t(optw_MV1) %*% SIGMA %*% optw_MV1)*12  # annual variance of MV portfolio (SS)
p_stdv_opt = sqrt(p_var_opt)   
sharpe_MV1 = (p_ret_opt-rf)/p_stdv_opt 
optw_MV1 # "optimal weights of MV portfolio (SS)"
```

    ##             [,1]
    ## [1,]  0.21493340
    ## [2,] -0.11481445
    ## [3,]  0.27281289
    ## [4,]  0.57353814
    ## [5,]  0.05292898

``` r
as.numeric(p_ret_opt)   # "annual return of MV portfolio (SS)"
```

    ## [1] 0.08473263

``` r
as.numeric(p_stdv_opt)  # "annual st.dev of MV portfolio (SS)"
```

    ## [1] 0.1187085

``` r
as.numeric(sharpe_MV1) # "annual Sharpe ratio of MV portfolio (SS)"
```

    ## [1] 0.7189121

### Static optimization - no SS (no Short sale)

``` r
MinVar_NoSS_results <- optim( par=W, Min_Var, SIGMA=SIGMA, 
                            lower=c(rep(0.001,k)), upper=c(rep(0.999,k)), method=c("L-BFGS-B"))
optw_MV2 = MinVar_NoSS_results$par
p_ret_opt = (t(optw_MV2) %*% MU)*12
p_var_opt = (t(optw_MV2) %*% SIGMA %*% optw_MV2)*12
p_stdv_opt = sqrt(p_var_opt) 
sharpe_MV2 = (p_ret_opt-rf)/p_stdv_opt
optw_MV2               # "optmal weights of MV portfolio (no SS)"
```

    ##           [,1]
    ## [1,] 0.2319831
    ## [2,] 0.0010000
    ## [3,] 0.2293858
    ## [4,] 0.5360369
    ## [5,] 0.0010000

``` r
as.numeric(p_ret_opt)  # "annual return of MV portfolio (no SS)"
```

    ## [1] 0.07895279

``` r
as.numeric(p_stdv_opt) # "annual st.dev of MV portfolio (no  SS)"
```

    ## [1] 0.11903

``` r
as.numeric(sharpe_MV2) # "annual Sharpe ratio of MV portfolio (no SS)"
```

    ## [1] 0.6684125

## Maximum Sharpe

### Short Sales allowed

``` r
Wr <- as.matrix(c(rep(1/(k+1),k+1)))    # create Weight matrix with k+1 assets
MaxSharp_results <- optim( par=Wr, Max_sharp, SIGMA=SIGMAr, MU=MUr)
optw_MS = MaxSharp_results$par
optw_MS # optimal weights of MS portfolio (SS) 
```

    ##            [,1]
    ## [1,]  2.5678705
    ## [2,] -3.7204617
    ## [3,]  1.7339991
    ## [4,] -0.3158146
    ## [5,]  1.9095375
    ## [6,] -1.1753890

``` r
sum(optw_MS)
```

    ## [1] 0.9997418

``` r
p_ret_opt = (t(optw_MS) %*% MUr)*12    
p_var_opt = (t(optw_MS) %*% SIGMAr %*% optw_MS)*12
p_stdv_opt = sqrt(p_var_opt) 
sharpe_MS1 = (p_ret_opt-rf)/p_stdv_opt             # annual Sharpe ratio of MS portfolio (SS)
p_ret_opt # portfolio return of MS portfolio (SS)
```

    ##           [,1]
    ## [1,] 0.5497179

``` r
p_var_opt # portfolio variance of MS portfolio (SS)
```

    ##           [,1]
    ## [1,] 0.2004866

``` r
p_stdv_opt # portfolio s.deviation of MS portfolio (SS)
```

    ##           [,1]
    ## [1,] 0.4477573

``` r
sharpe_MS1 #  annual Sharpe ratio of MS portfolio (SS)
```

    ##          [,1]
    ## [1,] 1.229073

### No Short Sales allowed

``` r
Wr <- as.matrix(c(rep(1/(k+1),k+1)))
MaxSharp_results2 <- optim( par=Wr, Max_sharp, SIGMA=SIGMAr, MU=MUr, 
                           lower=c(rep(0.001,k+1)), upper=c(rep(0.999,k+1)), method=c("L-BFGS-B"))
optw_MS2 = MaxSharp_results2$par
optw_MS2 # optimal portfolio weights of MS portfolio (no SS)
```

    ##           [,1]
    ## [1,] 0.9949833
    ## [2,] 0.0010000
    ## [3,] 0.0010000
    ## [4,] 0.0010000
    ## [5,] 0.0010000
    ## [6,] 0.0010000

``` r
p_ret_opt = (t(optw_MS2) %*% MUr)*12     
p_var_opt = (t(optw_MS2) %*% SIGMAr %*% optw_MS2)*12
p_stdv_opt = sqrt(p_var_opt)  
sharpe_MS2 = (p_ret_opt-rf)/p_stdv_opt    # calculate annual Sharpe ratio of MS portfolio (no SS)
p_ret_opt  # portfolio return of MS portfolio (no SS)
```

    ##           [,1]
    ## [1,] 0.1306134

``` r
p_stdv_opt # portfolio st.deviation of MS portfolio (no SS)
```

    ##           [,1]
    ## [1,] 0.1400777

``` r
sharpe_MS2 # annual Sharpe ratio of MS portfolio (no SS)
```

    ##           [,1]
    ## [1,] 0.9367776

## Minimize Conditional VaR - no Short Sales allowed

``` r
CVar_results <- optim( par=W, MCVar, ret=returns, alpha=0.05,
                            lower=c(rep(0.001,k+1)), upper=c(rep(0.999,k+1)), method=c("L-BFGS-B"))
optw_CV = CVar_results$par
optw_CV
```

    ##             [,1]
    ## [1,] 0.659190383
    ## [2,] 0.001000000
    ## [3,] 0.001000000
    ## [4,] 0.325996743
    ## [5,] 0.008651878

``` r
p_ret_opt = (t(optw_CV) %*% MU)*12     
p_var_opt = (t(optw_CV) %*% SIGMA %*% optw_CV)*12
p_stdv_opt = sqrt(p_var_opt)  
sharpe_CV = (p_ret_opt-rf)/p_stdv_opt # calculate annual Sharpe ratio of CV portfolio (no SS)
p_ret_opt # portfolio return of optimized CV portfolio (no SS)
```

    ##           [,1]
    ## [1,] 0.1067647

``` r
p_var_opt # portfolio variance of optimized CV portfolio (no SS)
```

    ##            [,1]
    ## [1,] 0.01603704

``` r
p_stdv_opt # portfolio st.deviation of optimized CV portfolio (no SS)
```

    ##           [,1]
    ## [1,] 0.1266375

``` r
sharpe_CV  # annual Sharpe ratio of optimized CV portfolio (no SS)
```

    ##           [,1]
    ## [1,] 0.8478776

Plot the returns from static weights, no short sales allowed

``` r
# take all period and apply the weights
ret_rf = as.matrix(cbind(returns, rf))  # add risk-free

# 1.equaly weighted return
EW_r = rowMeans(ret_rf)
# 2. Minimum variance return
MinV_r  = returns %*% optw_MV2
# 3. Maximum Sharpe return
MaxS_r  = ret_rf %*% optw_MS2
# 3. Maximum Sharpe return
CV_r  = returns %*% optw_CV

# bind the results by c=column and store them in the big matrix
OPT_STAT = data.frame( MinV_r, MaxS_r, CV_r)
summary(OPT_STAT)
```

    ##      MinV_r              MaxS_r              CV_r          
    ##  Min.   :-0.106505   Min.   :-0.10628   Min.   :-0.105796  
    ##  1st Qu.:-0.010750   1st Qu.:-0.01881   1st Qu.:-0.016813  
    ##  Median : 0.009403   Median : 0.01713   Median : 0.013578  
    ##  Mean   : 0.006579   Mean   : 0.01088   Mean   : 0.008897  
    ##  3rd Qu.: 0.032149   3rd Qu.: 0.03999   3rd Qu.: 0.033399  
    ##  Max.   : 0.069987   Max.   : 0.09965   Max.   : 0.089014

Statistics on monthly base

``` r
means = colMeans(OPT_STAT)
means # "Means of returns for different strategies"                      
```

    ##     MinV_r     MaxS_r       CV_r 
    ## 0.00657940 0.01088445 0.00889706

``` r
vari = apply(OPT_STAT, 2, var)    # variance of returns for different strategies
stdev = sqrt(vari)   
stdev # "St.deviation of returns for different strategies"
```

    ##     MinV_r     MaxS_r       CV_r 
    ## 0.03436100 0.04043696 0.03655708

``` r
sharpe = (means - rf) / stdev   
sharpe # "Sharpe ratio of returns for different strategies"
```

    ##    MinV_r    MaxS_r      CV_r 
    ## 0.2091829 0.2842147 0.2600151

Plot Cumulative Returns from static optimization

``` r
CUM_STAT = cumsum(OPT_STAT)
ret_M  = data.frame(date = datum[-(1:2)], CUM_STAT)
plotdf <- melt(ret_M[,], id="date")
p = ggplot(plotdf, aes(x=date, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  ggtitle("Cumulative Returns from static optimization")+
  theme(axis.text.x = element_text(size=6, angle=90))+
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
```
![](https://github.com/petaleks/petaleksandar/blob/ee34ae0fd86ea66c0a31d36676e5ffc5b554cb45/data/figure-gfm/plot_cum_ret-1.png?raw=true)


``` r
namepic = paste0("Cumulative Returns from static optimization.jpg")
ggsave(namepic, plot = p)
```

    ## Saving 7 x 5 in image

# Dynamic optimization

Prepare data for loop

``` r
rollw = 24         
rollw # "length of rolling window"
```

    ## [1] 24

``` r
n = dim(returns)[1]    # number of total time-periods
nom_rolls = n - rollw  
nom_rolls # "number of iterations"
```

    ## [1] 57

``` r
# create empty storage data frame for returns calculated using optimized weights
OPT_mat <- data.frame(EW_r=double(),
                     MinV_r=double(),
                     MaxS_r=double(),
                     CV_r=double(),
                     stringsAsFactors=FALSE)

# create empty storage data frame for optimized MV weights
MV_weights <- data.frame(matrix(ncol = k, nrow = 0)) 
colnames(MV_weights) = colnames(returns) # change the names of columns so they show up good on figures

# create empty storage data frame for optimized MS weights
MS_weights <- data.frame(matrix(ncol = k+1, nrow = 0))
colnames(MS_weights) = c(colnames(returns), "rf")

# create empty storage data frame for optimized CV weights
CV_weights <- data.frame(matrix(ncol = k, nrow = 0))
colnames(CV_weights) = colnames(returns)

alpha = 0.05 # define confidence interval for CVaR
```

Run the loop ( it takes a 15 sec)

``` r
# start the loop
for (r in 1:nom_rolls){
  
  # choose the rolling window with the data
  wind = returns[r:(r-1 +rollw),]
  
  # find MU and sigma for each new rolling period
  MU = as.matrix(apply(wind, 2 , mean))   # monthly mean return
  SIGMA   <- cov(as.matrix(wind))         # monthly sample covariance
  
  
  # run the solver for minimum variance (no SS)
  MWopt_results <- optim( par=W, Min_Var, SIGMA=SIGMA,
                          lower=c(rep(0.001,k)), upper=c(rep(0.999,k)), method=c("L-BFGS-B"))
  MWopt = MWopt_results$par  # extract weights from the list object
  MV_weights[r,] = t(MWopt)  # store as separate row in the matrix of MV weights
  
  # Introduce risk-free, take the last one known
  rf = rf_serie[(r-1 +rollw)]/12 # change to monthly risk free rate
  
  # create bordered matrix
  zeroc <- matrix(c(rep(0,k)),k,1) # column with zeros
  zeror <- matrix(c(rep(0,k+1)),1, k+1) # row with zeros
  colnames(zeroc)[1] =c("rf")
  rownames(zeror)[1] =c("rf")
  SIGMAr = cbind(SIGMA, zeroc)
  SIGMAr = rbind(SIGMAr, zeror)
  MUr = rbind(MU, rf)
  
  # run the solver for maximum Sharpe (no SS)
  MSopt_results <- optim( par=Wr, Max_sharp, SIGMA=SIGMAr, MU=MUr, 
                          lower=c(rep(0.001,k)), upper=c(rep(0.999,k)), method=c("L-BFGS-B"))
  MSopt = MSopt_results$par
  MSopt
  MS_weights[r,] = t(MSopt)
  
  # run the solver for CvaR
  CVar_results <- optim( par=W, MCVar, ret=wind, alpha=alpha,
                          lower=c(rep(0.001,k+1)), upper=c(rep(0.999,k+1)), method=c("L-BFGS-B"))
  CV_opt = CVar_results$par
  CV_weights[r,] = t(CV_opt)
  
  
  # take next period and apply the weights
  new_per = data.frame(returns[r +rollw,])
  new_per = as.matrix(new_per)                # make it in matrix form
  new_per_rf = as.matrix(rbind(new_per, rf))  # add risk-free
  
  # 1.equaly weighted return
  EW_r = colMeans(new_per_rf)
  # 2. Minimum variance return
  MinV_r  = t(MWopt) %*% new_per
  # 3. Maximum Sharpe return
  MaxS_r  = t(MSopt) %*% new_per_rf
  # 3. Maximum Sharpe return
  CV_r  = t(CV_opt) %*% new_per
  
  # bind the results by c=column and store them in the big matrix
  OPT_mat[r,] = cbind(EW_r, MinV_r, MaxS_r, CV_r)
  
}
```

``` r
# returns on Swedish index is in last column (k column)
# need to slect rows from row rollw+1=21 till row n=81 and column k
OPTIM_mat = data.frame(OPT_mat, SWEDOMXr = returns[(rollw+1):n, k] ) # add swedish index
summary(OPTIM_mat) 
```

    ##       EW_r               MinV_r              MaxS_r         
    ##  Min.   :-0.081829   Min.   :-0.090575   Min.   :-0.081008  
    ##  1st Qu.:-0.007178   1st Qu.:-0.010105   1st Qu.:-0.012334  
    ##  Median : 0.005787   Median : 0.014431   Median : 0.014112  
    ##  Mean   : 0.004998   Mean   : 0.006903   Mean   : 0.007409  
    ##  3rd Qu.: 0.024643   3rd Qu.: 0.028720   3rd Qu.: 0.026529  
    ##  Max.   : 0.061708   Max.   : 0.071360   Max.   : 0.067518  
    ##       CV_r              SWEDOMXr        
    ##  Min.   :-0.090556   Min.   :-0.163402  
    ##  1st Qu.:-0.007426   1st Qu.:-0.017821  
    ##  Median : 0.014754   Median : 0.013071  
    ##  Mean   : 0.009255   Mean   : 0.006835  
    ##  3rd Qu.: 0.030424   3rd Qu.: 0.038125  
    ##  Max.   : 0.069508   Max.   : 0.101234

``` r
# extract the data, so one can use it in regression
# in order to calculate jensen alpha, beta  and other performance measures
write.csv(OPTIM_mat, "Returns from different strategies.csv")
```

## Plot dynamic weights

Plot the cumulative returns from dynamic strategies

``` r
# stack data over each other
CUMOPT = cumsum(OPTIM_mat) # create cumulative returns
OPT_M  = data.frame(date = datum[(rollw +1):n], CUMOPT) # add column with dates
plotdf <- melt(OPT_M[,], id="date")
p = ggplot(plotdf, aes(x=date, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  ggtitle("Cumulative Returns from dynamic optimization")+
  theme(axis.text.x = element_text(size=6, angle=90))+
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
```
![](https://github.com/petaleks/petaleksandar/blob/ee34ae0fd86ea66c0a31d36676e5ffc5b554cb45/data/figure-gfm/plot_cumret_dyn-1.png?raw=true)

``` r
namepic = paste0("Cumulative returns from dynamic optimization.jpg") # define name of figure
ggsave(namepic, plot = p)                  # save it
```

    ## Saving 7 x 5 in image

Plot the Minimum Variance weights

``` r
MVw  = data.frame(date = datum[(rollw +1):n], MV_weights) # add column with dates
plotdf <- melt(MVw[,], id="date")
p = ggplot(plotdf, aes(x=date, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  ggtitle("Minimum Variance weights")+
  theme(axis.text.x = element_text(size=6, angle=90))+
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
```
![](https://github.com/petaleks/petaleksandar/blob/ee34ae0fd86ea66c0a31d36676e5ffc5b554cb45/data/figure-gfm/dyn_minV_we-1.png?raw=true)


``` r
namepic = paste0("Minimum Variance weights.jpg")
ggsave(namepic, plot = p)
```

    ## Saving 7 x 5 in image

Plot the Maximum Sharpe weights

``` r
MSw  = data.frame(date = datum[(rollw +1):n], MS_weights) # add column with dates
plotdf <- melt(MSw[,], id="date")
p = ggplot(plotdf, aes(x=date, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  ggtitle("Maximum Sharpe weights")+
  theme(axis.text.x = element_text(size=6, angle=90))+
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
```
![](https://github.com/petaleks/petaleksandar/blob/ee34ae0fd86ea66c0a31d36676e5ffc5b554cb45/data/figure-gfm/dyn_maxS_we-1.png?raw=true)


``` r
namepic = paste0("Maximum Sharpe weights.jpg")
ggsave(namepic, plot = p)
```

    ## Saving 7 x 5 in image

Plot the CvaR weights

``` r
CVw  = data.frame(date = datum[(rollw +1):n], CV_weights) # add column with dates
plotdf <- melt(CVw[,], id="date")
p = ggplot(plotdf, aes(x=date, y=value, group=variable)) +
  geom_line(aes(color=variable))+
  ggtitle("CvaR weights")+
  theme(axis.text.x = element_text(size=6, angle=90))+
  theme(plot.title = element_text(hjust = 0.5))
plot(p)
```
![](https://github.com/petaleks/petaleksandar/blob/ee34ae0fd86ea66c0a31d36676e5ffc5b554cb45/data/figure-gfm/plot_dyn_cvar_weights-1.png?raw=true)

``` r
namepic = paste0("CvaR weights.jpg")
ggsave(namepic, plot = p)
```

    ## Saving 7 x 5 in image

## CAPM Regression

``` r
# formula for different strategies
rel_EW = as.formula("EW_r   ~ SWEDOMXr") 
rel_MV = as.formula("MinV_r ~ SWEDOMXr") 
rel_MS = as.formula("MaxS_r ~ SWEDOMXr") 
rel_CV = as.formula("CV_r   ~ SWEDOMXr") 

# regression for different strategies
# store them in lm (linear model) R object
lin_EW = lm(rel_EW, data=OPTIM_mat)
lin_MV = lm(rel_MV, data=OPTIM_mat)
lin_MS = lm(rel_MS, data=OPTIM_mat)
lin_CV = lm(rel_CV, data=OPTIM_mat)
```

Create summary object to get more statistics

``` r
# 
# see the summary
sumlin_EW = summary(lin_EW)
sumlin_EW
```

    ## 
    ## Call:
    ## lm(formula = rel_EW, data = OPTIM_mat)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.039302 -0.007108 -0.002476  0.010014  0.021081 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.00137    0.00168   0.816    0.418    
    ## SWEDOMXr     0.53079    0.03446  15.404   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01256 on 55 degrees of freedom
    ## Multiple R-squared:  0.8118, Adjusted R-squared:  0.8084 
    ## F-statistic: 237.3 on 1 and 55 DF,  p-value: < 2.2e-16

``` r
sumlin_MV = summary(lin_MV)
sumlin_MV
```

    ## 
    ## Call:
    ## lm(formula = rel_MV, data = OPTIM_mat)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.064611 -0.014549 -0.001519  0.015222  0.044943 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.00354    0.00299   1.184    0.242    
    ## SWEDOMXr     0.49197    0.06134   8.020 8.08e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02235 on 55 degrees of freedom
    ## Multiple R-squared:  0.5391, Adjusted R-squared:  0.5307 
    ## F-statistic: 64.33 on 1 and 55 DF,  p-value: 8.076e-11

``` r
sumlin_MS = summary(lin_MS)
sumlin_MS
```

    ## 
    ## Call:
    ## lm(formula = rel_MS, data = OPTIM_mat)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.058542 -0.018860 -0.003708  0.023783  0.051570 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 0.004353   0.003457   1.259    0.213    
    ## SWEDOMXr    0.447191   0.070922   6.305 5.14e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02584 on 55 degrees of freedom
    ## Multiple R-squared:  0.4196, Adjusted R-squared:  0.409 
    ## F-statistic: 39.76 on 1 and 55 DF,  p-value: 5.145e-08

``` r
sumlin_CV = summary(lin_CV)
sumlin_CV
```

    ## 
    ## Call:
    ## lm(formula = rel_CV, data = OPTIM_mat)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.068736 -0.013770 -0.000197  0.016993  0.049688 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 0.006076   0.003401   1.786   0.0796 .  
    ## SWEDOMXr    0.465164   0.069774   6.667 1.32e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02543 on 55 degrees of freedom
    ## Multiple R-squared:  0.4469, Adjusted R-squared:  0.4369 
    ## F-statistic: 44.44 on 1 and 55 DF,  p-value: 1.325e-08

Extract parameters from the summary object, for each strategy

``` r
alpha_EW = sumlin_EW$coefficients[1]
beta_EW  = sumlin_EW$coefficients[2]
sigma_e_EW = sqrt(var(sumlin_EW$residuals)) # standard deviation of the residuals

alpha_MV = sumlin_MV$coefficients[1]
beta_MV  = sumlin_MV$coefficients[2]
sigma_e_MV = sqrt(var(sumlin_MV$residuals)) # standard deviation of the residuals

alpha_MS = sumlin_MS$coefficients[1]
beta_MS  = sumlin_MS$coefficients[2]
sigma_e_MS = sqrt(var(sumlin_MS$residuals)) # standard deviation of the residuals

alpha_CV = sumlin_CV$coefficients[1]
beta_CV  = sumlin_CV$coefficients[2]
sigma_e_CV = sqrt(var(sumlin_CV$residuals)) # standard deviation of the residuals
```

## Statistics from dynamic strategies

Sharpe ratios

``` r
# Sharpe ratio= (mu_return - rf)/sigma_return

means = colMeans(OPTIM_mat)
means # "Means of returns for different strategies"                             
```

    ##        EW_r      MinV_r      MaxS_r        CV_r    SWEDOMXr 
    ## 0.004998066 0.006902521 0.007409207 0.009255042 0.006834958

``` r
vari = apply(OPTIM_mat, 2, var)    # variance of returns for different strategies
stdev = sqrt(vari)    
stdev # "St.deviation of returns for different strategies"
```

    ##       EW_r     MinV_r     MaxS_r       CV_r   SWEDOMXr 
    ## 0.02868763 0.03262955 0.03361941 0.03388360 0.04869706

``` r
sharpe = (means - rf) / stdev   
sharpe # "Sharpe ratio of returns for different strategies"
```

    ##      EW_r    MinV_r    MaxS_r      CV_r  SWEDOMXr 
    ## 0.1954292 0.2301857 0.2384795 0.2910958 0.1528489

Appraisal ratio

``` r
# Appraisal ratio = alpha / sigma_residuals
# Use alphas and sigma_e from regression output to calculate appraisal ratio
app_ratio_EW = alpha_EW / sigma_e_EW
app_ratio_MV = alpha_MV / sigma_e_MV
app_ratio_MS = alpha_MS / sigma_e_MS
app_ratio_CV= alpha_CV / sigma_e_CV
cbind(app_ratio_EW, app_ratio_MV, app_ratio_MS, app_ratio_CV)
```

    ##      app_ratio_EW app_ratio_MV app_ratio_MS app_ratio_CV
    ## [1,]       0.1101    0.1597991    0.1699393    0.2411093

Treynor ratios

``` r
# Treynor ratio = (mu_return - rf)/beta

beta_vec = as.matrix(c(beta_EW, beta_MV, beta_MS, beta_CV, 1)) # create vector of betas
beta_vec = t(beta_vec) # transpose it to be row, not column vector
# note that beta for the index is 1
# Use this vector with betas to calculate Treynor ratio
Treynors = (means - rf) / beta_vec ; 
colnames(Treynors) <- names(means) 
Treynors # see the Treynors
```

    ##            EW_r     MinV_r     MaxS_r       CV_r    SWEDOMXr
    ## [1,] 0.01056237 0.01526694 0.01792868 0.02120406 0.007443291

Information ratios

``` r
# Information ratios = E(Rp-Rb)/ SD(Rp-Rb)

# create difference between portfolio returns for each of the 4 strategies 
# and the benchmark (Swedish stock index) 
OPTIMIZE_mat = data.frame(EW_rp_rb   = OPTIM_mat$EW_r- OPTIM_mat$SWEDOMXr,
                          MinV_rp_rb = OPTIM_mat$MinV_r- OPTIM_mat$SWEDOMXr,
                          MaxS_rp_rb = OPTIM_mat$MaxS_r- OPTIM_mat$SWEDOMXr,
                          CV_rp_rb   = OPTIM_mat$CV_r- OPTIM_mat$SWEDOMXr )

head(OPTIMIZE_mat) # see first 6 rows
```

    ##        EW_rp_rb    MinV_rp_rb   MaxS_rp_rb      CV_rp_rb
    ## 1 -0.0159143088 -0.0064492041 -0.014585235 -0.0092259830
    ## 2  0.0075043101  0.0028337930  0.019068422 -0.0047068358
    ## 3  0.0005565496 -0.0009458905 -0.008214167 -0.0009434538
    ## 4 -0.0106109407  0.0108851757 -0.027111922  0.0215905522
    ## 5  0.0124297767  0.0114488067  0.005795693  0.0245035268
    ## 6 -0.0101766098 -0.0048380293 -0.013203620 -0.0047573387

``` r
summary(OPTIMIZE_mat) # see the summary
```

    ##     EW_rp_rb           MinV_rp_rb           MaxS_rp_rb        
    ##  Min.   :-0.054323   Min.   :-6.966e-02   Min.   :-7.568e-02  
    ##  1st Qu.:-0.021712   1st Qu.:-1.489e-02   1st Qu.:-2.304e-02  
    ##  Median :-0.002521   Median : 2.834e-03   Median :-3.574e-05  
    ##  Mean   :-0.001837   Mean   : 6.756e-05   Mean   : 5.742e-04  
    ##  3rd Qu.: 0.014957   3rd Qu.: 1.254e-02   3rd Qu.: 2.018e-02  
    ##  Max.   : 0.081573   Max.   : 1.217e-01   Max.   : 1.241e-01  
    ##     CV_rp_rb       
    ##  Min.   :-0.07174  
    ##  1st Qu.:-0.02239  
    ##  Median : 0.00246  
    ##  Mean   : 0.00242  
    ##  3rd Qu.: 0.02314  
    ##  Max.   : 0.13488

``` r
diff_means = colMeans(OPTIM_mat)
diff_means # "Difference between return for different strategies and benchmark"                           
```

    ##        EW_r      MinV_r      MaxS_r        CV_r    SWEDOMXr 
    ## 0.004998066 0.006902521 0.007409207 0.009255042 0.006834958

``` r
diif_vari = apply(OPTIM_mat, 2, var)    
diif_vari # "Variance of difference between return for different strategies and benchmark for different strategies"
```

    ##         EW_r       MinV_r       MaxS_r         CV_r     SWEDOMXr 
    ## 0.0008229799 0.0010646874 0.0011302645 0.0011480987 0.0023714035

``` r
te = sqrt(diif_vari)              
te # "Tracking error for different strategies"
```

    ##       EW_r     MinV_r     MaxS_r       CV_r   SWEDOMXr 
    ## 0.02868763 0.03262955 0.03361941 0.03388360 0.04869706

``` r
inform_ratios = diff_means / te   
inform_ratios # "Information ratio of returns for different strategies"
```

    ##      EW_r    MinV_r    MaxS_r      CV_r  SWEDOMXr 
    ## 0.1742237 0.2115420 0.2203848 0.2731422 0.1403567

#### Readings
Eun, Cheol S., and Bruce G. Resnick. *Exchange rate uncertainty, forward contracts, and international portfolio selection*. Journal of Finance (1988): 197-215.

Body, Marcus, Kane *Investments*, chap. 25, McGraw-Hill.

#### R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.
