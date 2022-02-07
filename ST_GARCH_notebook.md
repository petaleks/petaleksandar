---
layout: page
title: Spatial-Temporal GARCH
permalink: /ST_GARCH/
--- 


### Required libraries
```r
library(ggplot2)
library(reshape2) 
#for sp models
library(sp)
library(raster)
library(rgdal)
library(sf)
library(spData)
library(spDataLarge)
library(spatialreg)
library(spdep)
library(plyr)
library(tidyr)
# for numerical optimization and calculating Hessian/Jacobian
library(numDeriv)
# for direct download from github
library(RCurl)
library(foreign)
```

### Get residuals from house prices in Iowa per year/month with their coordinates

``` r
urlfile<-'https://raw.githubusercontent.com/petaleks/first_repo/master/SAR_dnb_resiudals.csv'
RESID <-read.csv(urlfile)
head(RESID)
```

    ##        ym         X        Y          eps
    ## 1  2010_5 -93.61975 42.05403  0.065484080
    ## 2  2006_5 -93.61846 42.05341 -0.054902706
    ## 3  2006_8 -93.61818 42.05333 -0.003273209
    ## 4 2009_10 -93.61890 42.05304 -0.067943595
    ## 5  2006_6 -93.61926 42.05311 -0.070721464
    ## 6 2006_10 -93.61644 42.05252  0.044184872

### Model specifications

Classical ARCH model structure

-   conditional mean
    *μ*<sub>*t*</sub> = *μ*(*θ*, *x*<sub>*t*</sub>) = *E*(*Y*<sub>*t*</sub>\|*x*<sub>*t*</sub>)
-   conditional variance
    *σ*<sub>*t*</sub><sup>2</sup> = *σ*<sub>*t*</sub><sup>2</sup>(*θ*, *x*<sub>*t*</sub>) = *E*((*Y*<sub>*t*</sub> − *μ*<sub>*t*</sub>)<sup>2</sup>\|*x*<sub>*t*</sub>)

**Spatio-temporal ARCH (1,1)** model in matrix format

*σ*<sub>*n*<sub>*t*</sub>, *t*</sub><sup>2</sup> = *ω* + *α W̄*<sub>*t*</sub>*ϵ*<sub>*n*<sub>*t*</sub>, *t* − 1</sub><sup>2</sup>


**Spatio-temporal GARCH (1,1)** model in matrix format

*σ*<sub>*n*<sub>*t*</sub>, *t*</sub><sup>2</sup> = *ω* + *α W̄*<sub>*t*</sub>*ϵ*<sub>*n*<sub>*t*</sub>, *t* − 1</sub><sup>2</sup> + *γ W̄*<sub>*t*</sub>*σ*<sub>*n*<sub>*t*</sub>, *t* − 1</sub><sup>2</sup>


**Assymetric Spatio-Temporal GARCH (1,1)** model in matrix format 

*σ*<sub>*n*<sub>*t*</sub>, *t*</sub><sup>2</sup> = *ω* + *α W̄*<sub>*t*</sub>*ϵ*<sub>*n*<sub>*t*</sub>, *t* − 1</sub><sup>2</sup> + *γ W̄*<sub>*t*</sub>*σ*<sub>*n*<sub>*t*</sub>, *t* − 1</sub><sup>2</sup> + *δ I*<sub>(*ϵ*<sub>*n*<sub>*t*</sub>, *t* − 1</sub> &lt; 0)</sub>*W̄*<sub>*t*</sub>*ϵ*<sub>*n*<sub>*t*</sub>, *t* − 1</sub><sup>2</sup>

## Help functions

### function for trace

``` r
tr <- function(A) sum(diag(A))
```

### function for creating spatial-temporal matrices

``` r
create_matrices <- function(treshold, RESID_new) {
  
  # function for summing null in nbd objects
  sum_null <- function(nbdobj){
    c=0
    for (l in nbdobj) {
      if (is.null(l)) {
        c=c+1
      }
    }
    return(c)
  }
  
  ym = as.character(RESID_new$ym)
  table(ym)
  vecym = unique(ym)
  koord = cbind(RESID_new$X,RESID_new$Y)
  rownames(koord) = 1:dim(koord)[1]
  resids = RESID_new$eps
  mean_resid = mean(resids^2)
  
  
  # EXTRACT matrix
  {  EXTRACT <- data.frame(Characters=character(),
                           Doubles=double(),
                           Doubles=double(),
                           Doubles=double(),
                           Doubles=double(),
                           stringsAsFactors=FALSE)
  }
  colnames(EXTRACT) =c("ym", "X", "Y", "eps" )
  
  
  for (j in 2:length(vecym)) {
    
    indeks = ym == vecym[j-1]
    indeks1 = ym == vecym[j]
    
    coords  <- as.matrix(koord[indeks,])
    coords1 <- as.matrix(koord[indeks1,])
 
    eps = resids[indeks]
    eps1 = resids[indeks1]
    
    # 2) double matrices
    
    # bind 2 datasets
    all_coord = rbind(coords, coords1)
    all_eps = c(eps, eps1) # create common error-vector and weight matrix for 2 periods
    
    IDs_all <- rownames(all_coord)
    # IDs_all <- 1:dim(all_coord)[1]
    dnb_all <- dnearneigh(all_coord, 0, treshold, row.names=IDs_all,longlat=TRUE)
    nbd_all <- nbdists(dnb_all, all_coord,longlat=TRUE);
   
    gl_all     <- lapply(nbd_all, function(x) 1/x);
    lw_nb_all  <- nb2listw(dnb_all, glist=gl_all, zero.policy=F) # class(lw_nb) #"listw" "nb"
    all_W      <- nb2mat(lw_nb_all$neighbours,glist=gl_all, style="W",zero.policy=F)
    
    name_all = paste0("./", treshold,"km_", j,"_all_mat.csv")
    write.csv(all_W, file=name_all)
    
    
    # create first Single matrix
    if(j==2){
      m = length(eps)
      W = diag(m)
      if ( class(coords)[1]== "numeric" ){
        coords = (t(as.matrix(coords)))    
      }
      name = paste0("./", treshold,"km_", j-1,"_mat.csv")
      write.csv(W, file=name)
      
      # store coords and epsilons
      add1 = data.frame( rep(vecym[j-1],length(eps)),  coords, eps  ) 
      colnames(add1) =c("ym", "X", "Y", "eps" )
      EXTRACT = rbind(EXTRACT,add1)
    }
    
    
    # create second Single matrix
    m1 = length(eps1) # size is determined already to sample j
    W1 = diag(m1)
    if ( class(coords1)[1]== "numeric" ){
      coords1 = (t(as.matrix(coords1)))    }
    name = paste0("./", treshold,"km_", j,"_mat.csv")
    write.csv(W1, file=name)
    
    add2= data.frame(  rep(vecym[j],  length(eps1)), coords1, eps1 )
    colnames(add2) =c("ym", "X", "Y", "eps" )
    EXTRACT = rbind(EXTRACT,add2)
    
  }
  colnames(EXTRACT) = c("ym", "X", "Y", "eps" )
  rownames(EXTRACT) = 1: dim(EXTRACT)[1]
  name_all = paste0("./Extract_RES_", treshold,"km.csv")
  write.csv(EXTRACT, file=name_all)
}
```

### function for identifying points with no neighbours

``` r
delete_null_coords <- function(nbdobj, koordinati, epsiloni){
  c = 0
  nulllist = vector()
  for (l in nbdobj) {
    c=c+1
    if (is.null(l)) {
      nulllist = c(nulllist, c)
    }
  }
  koordinati =data.frame(koordinati)
  
  if (length(nulllist) > 0) {
    if ( dim(koordinati)[1] >1 ) {
      redovi = rownames( koordinati[nulllist,] )
      koordinati = koordinati[-nulllist,]
    }
    if ( dim(koordinati)[1] == 1 ) {
      redovi = rownames(koordinati[nulllist])
      koordinati = koordinati[-nulllist]
    }
    epsiloni   = epsiloni[-nulllist]
  } else {redovi = ""}
  
  return(list(redovi, koordinati,epsiloni))
}
```

### function for deleting points with no neighbours

``` r
no_null_neighbours = function(resdf, treshold){
  ym = as.character(resdf$ym)
  table(ym)
  vecym = unique(ym)
  koord = data.frame(cbind(resdf$X,resdf$Y))
  resids = resdf$eps
  clean_red =c()
  for (j in 2:(length(vecym)-1)) {
    
    indeks = ym == vecym[j-1]
    indeks1 = ym == vecym[j]
    indeks2 = ym == vecym[j+1]
    
    coords  <- as.matrix(koord[indeks,])
    coords1 <- as.matrix(koord[indeks1,])
    coords2 <- as.matrix(koord[indeks2,])

    sum(indeks)
    sum(indeks1)
    sum(indeks2)
    eps = resids[indeks]
    eps1 = resids[indeks1]
    eps2 = resids[indeks2]
    
    
    # a) bind 1-2 datasets
    all_coorda = rbind(coords, coords1)
    all_epsa = c(eps, eps1) # create common error-vector and weight matrix for 2 periods
    IDs_alla <- 1:dim(all_coorda)[1]
    dnb_alla <- dnearneigh(all_coorda, 0, treshold, row.names=IDs_alla,longlat=TRUE)
    nbd_alla <- nbdists(dnb_alla, all_coorda,longlat=TRUE);
    
    # b) bind 2-3 datasets
    all_coordb = rbind(coords1, coords2)
    all_epsb = c(eps1, eps2) # create common error-vector and weight matrix for 2 periods
    IDs_allb <- 1:dim(all_coordb)[1]
    dnb_allb <- dnearneigh(all_coordb, 0, treshold, row.names=IDs_allb,longlat=TRUE)
    nbd_allb <- nbdists(dnb_allb, all_coordb,longlat=TRUE);
    
    # coordinates without the neighbours for given treshold
    cleana = delete_null_coords(nbd_alla, all_coorda, all_epsa)
    cleanb = delete_null_coords(nbd_allb, all_coordb, all_epsb)
    redovia = cleana[[1]]
    redovib = cleanb[[1]]
    
    # update the vector with coordinates without the neighbours for given treshold
    clean_red = unique(c(clean_red, redovia, redovib))
    
   
  }
  
  clean_red = as.numeric(clean_red)
  clean_red = clean_red[!is.na(clean_red)]
  if( length(clean_red)==0) {
    new_resid=resdf
    clean_red = NULL
  } else {
    new_resid = resdf[-clean_red,]
  }
  
  return(list(new_resid, clean_red))
}
```

``` r
out = no_null_neighbours(RESID, treshold=2)
RESID_new = data.frame(out[[1]])
no_neighbours = out[[2]]
dim(RESID_new)
```

    ## [1] 2681    4

``` r
length(no_neighbours)
```

    ## [1] 10

``` r
c = 0
while( length(no_neighbours) > 0){
  c=c+1
  out = no_null_neighbours(RESID_new, treshold=2)
  RESID_new = data.frame(out[[1]])
  no_neighbours = out[[2]]
  print(dim(RESID_new)[1])
  print(length(no_neighbours))
  if  ( c == 20 | length(no_neighbours)  == 0) break; 
}
```

    ## [1] 2680
    ## [1] 1
    ## [1] 2680
    ## [1] 0

``` r
ggplot(RESID, aes(x=X, y=Y)) +geom_point()
```

![](ST_GARCH_rmarkdown_notebook_files/figure-gfm/plot1-1.png)<!-- -->

``` r
# Create W matrix
IDs <-  1:dim(RESID_new)[1]
koord = cbind(RESID_new$X, RESID_new$Y)
# k-Nearest neighbors
knb <- knn2nb(knearneigh(koord,k=1,longlat=T))
# Take a distance threshold such that all observations have at least one neighbor
treshold <- max(unlist(nbdists(knb, koord,longlat=T))) # longlat=T means we calculate in kilometres
# maximum distance to 1st neighbour is the threshold distance around centroid
dnb <- dnearneigh(koord, 0, treshold, row.names=IDs,longlat=TRUE)
# from 0km to the threshold
nbd <- nbdists(dnb, koord,longlat=TRUE);# nbd
gl     <- lapply(nbd, function(x) 1/x); # gl
lw_nb  <- nb2listw(dnb, glist=gl,zero.policy=T) # class(lw_nb) #"listw" "nb"
W      <- nb2mat(lw_nb$neighbours, glist=gl, style="W",zero.policy=T)
# W      <- listw2mat(lw_nb)
treshold
```

    ## [1] 0.398399

## Iterate functions

### function for creating SARCH model time series

``` r
create_SARCH_serie <- function(pars, RS, modelot, Wmat) {
  
  if (modelot == "SARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    gamma1   <- 0
    gamma2   <- 0
    delta1   <- 0
    delta2   <- 0
  }
  
  if (modelot == "ASARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    gamma1   <- 0
    gamma2   <- 0
    delta1   <- 0
    delta2   <- pars[3]
  }
  
  # print(c(omega, alpha1, alpha2, gamma1, gamma2, delta1, delta2 ))
  ym = as.character(RS$ym)
  coords   <- cbind(RS$X,RS$Y)
  eps = RS$eps
  N = length(eps)
  
  
  # 0. create dummy for negative shocks
  asim_dummy = rep(0, length(eps))
  asim_dummy[eps <= 0] <- 1
  
  # effect of past squared errors (spatial interpolation)
  weps = Wmat %*% (eps^2)     # spatially interpolated errors
  
  
  h1 <-as.vector(omega              + alpha2*weps + delta2*asim_dummy*weps)  # ASARCH
  xi <- eps / sqrt(h1)
  
  HM = data.frame(cbind( as.character(ym),coords,eps,h1, xi, weps))
  colnames(HM) =c("ym", "X", "Y", "eps", "h", "xi","weps")
  rownames(HM) = 1:dim(HM)[1]
  # calculate LL
  HM$eps  = gsub(",",".", HM$eps )
  HM$eps  = as.numeric(HM$eps)
  HM$h  = gsub(",",".", HM$h )
  HM$h  = as.numeric(HM$h)
  HM$xi  = gsub(",",".", HM$xi )
  HM$xi  = as.numeric(HM$xi)
  HM$weps  = gsub(",",".", HM$weps )
  HM$weps  = as.numeric(HM$weps)
 
  #name = paste0("./static_",modelot,"_HM_mat.csv")
  #write.csv(HM, file=name)
  vec_par = pars
  
  return(list(HM, vec_par))
}
```

### function for creating assymetric ST-GARCH model time series

``` r
create_ASTAGARCH_past_serie <- function(pars,RESID, modelot, treshold) {
  
  if (modelot == "ARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "SARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "ASARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- pars[3]
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "GARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- pars[3]
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "SARCH-GARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- pars[3]
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "ARCH-SGARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- pars[3]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "SARCH-SGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- pars[3]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "AGARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- pars[3]
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- pars[4]
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "STAGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- pars[3]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- pars[4]
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  # HM matrix
  {  HM <- data.frame(Characters=character(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      Doubles=double(),
                      stringsAsFactors=FALSE)
  }
  colnames(HM) =c("ym", "X", "Y", "eps", "h", "xi","leps", "weps","lh", "wh", "aleps", "aweps", "alh", "awh")
  
  ime =  paste0("./Extract_RES_", treshold,"km.csv")
  extract = data.frame(read.csv(ime, sep=",", header = T)); extract$X.1 =NULL
  vecym = unique(extract$ym)
  reziduali = extract$eps
  mean_resid = mean(reziduali^2) # variance
  koordinati = data.frame(cbind(extract$X, extract$Y))
  
  # big matrix for storing small matrices
  NT = length(reziduali)
  WW = data.frame(matrix(rep(0, NT*NT), nrow = NT, ncol= NT))
  colnames(WW)= 1:NT
  WWpp = data.frame(matrix(rep(0, NT*NT), nrow = NT, ncol= NT))
  colnames(WWpp)= 1:NT
  
  for (i in 2:length(vecym)) {
    
    name = paste0("./", treshold,"km_", i-1,"_mat.csv")
    W = read.csv(name, sep=",", header = T); W$X =NULL
    W = as.matrix(W)
    dim(W)[1]; 
    
    name1 = paste0("./", treshold,"km_", i,"_mat.csv")
    W1 = read.csv(name1, sep=",", header = T); W1$X =NULL
    W1 = as.matrix(W1)
    dim(W1)[1]
    
    name_all = paste0("./", treshold,"km_", i,"_all_mat.csv")
    all_W = read.csv(name_all, sep=",", header = T); all_W$X =NULL
    all_W =as.matrix(all_W)
    dim(all_W)[1]
    
    month_iminus = vecym[i-1]
    month_i      = vecym[i]
    
    indeks = extract$ym == month_iminus
    indeks1 = extract$ym == month_i
    
    # create common epsilon-vector 
    
    eps = reziduali[indeks]
    eps1 =reziduali[indeks1]
    all_eps =c(eps, eps1)
    
    coords = koordinati[indeks ,]
    coords1 = koordinati[indeks1 ,]
    
    # create common h-vector and weight matrix for 2 periods
    if (i==2){
      all_h  = c( rep(mean_resid,dim(W)[1]), rep(0, dim(W1)[1]))  # initiate the vector with zeros in new period, as later formula will fulfill this
    }
    if (i > 2){
      all_h  = c(h, rep(mean_resid, dim(W1)[1]))  
      # initiate the vector with zeros in new period, as later formula will fulfill this
    }
    
    
    # create W, T and WT matrices (MAKE ENDOGENIOUS matrices)
    first   <- c(rep(0, dim(W)[1]), rep(1, dim(W1)[1]))
    second  <- c(rep(1, dim(W)[1]), rep(0, dim(W1)[1])) # effect of old on new only
    third   <- c(rep(0, dim(W)[1]), rep(1, dim(W1)[1])) # effect of old (and new) on new only
    T_past    = first %o% second  # outer product to get 1/0 zero matrix
    T_current = first %o% third   # outer product to get 1/0 zero matrix
    T_pastpresent = T_past+T_current
    
    WT = T_past * all_W        # Hadamard product to make the W matrix -> directed W matrix
    rowsum = matrix(apply(WT,1,sum)) # make WT matrix row standardized
    WTs = matrix(rep(0, dim(WT)[1]* dim(WT)[1]), nrow = dim(WT)[1], ncol= dim(WT)[2] )
    WTs[rowsum > 0,] = WT[rowsum > 0,] / rowsum[rowsum > 0]
    
    WTcur = T_current * all_W        # Hadamard product to make the W matrix -> directed W matrix
    rowsum = matrix(apply(WTcur,1,sum)) # make WT matrix row standardized
    WTscur = matrix(rep(0, dim(WTcur)[1]* dim(WTcur)[1]), nrow = dim(WTcur)[1], ncol= dim(WTcur)[2] )
    WTscur[rowsum > 0,] = WTcur[rowsum > 0,] / rowsum[rowsum > 0]
    
    WTpp = T_pastpresent * all_W        # Hadamard product to make the W matrix -> directed W matrix
    rowsum = matrix(apply(WTpp,1,sum)) # make WT matrix row standardized
    WTspp = matrix(rep(0, dim(WTpp)[1]* dim(WTpp)[1]), nrow = dim(WTpp)[1], ncol= dim(WTpp)[2] )
    WTspp[rowsum > 0,] = WTpp[rowsum > 0,] / rowsum[rowsum > 0]
    
    #  fill the big matrix with small matrices
    ii= dim(W)[1];ii     # length of first matrix
    jj = dim(W1)[1]; jj  # length of second matrix
    gy = rownames(koordinati[indeks ,])  # place of coordinates in the matrix up-down
    gx = rownames(koordinati[indeks1 ,]) # place of coordinates in the matrix left-right
    gg = c(gy,gx)
    WW[gx, gy] <- WT[c((ii+1):(ii+jj)),(1:ii)] 
    
    # also for past-present matrices
    WWpp[gx, gg] <- WTpp[c((ii+1):(ii+jj)),(1:(ii+jj))] 
    
    # create equally weighted T matrix
    T_equal = matrix(rep(0, dim(T_past)[1] * dim(T_past)[1]), 
                     nrow = dim(T_past)[1],
                     ncol = dim(T_past)[2] )
    rst = matrix(apply(T_past,1,sum))
    T_equal[rst > 0,] = T_past[rst > 0,] / rst[rst > 0]
    
    T_equal_cur = matrix(rep(0, dim(T_current)[1]* dim(T_current)[1]), nrow = dim(T_current)[1], ncol= dim(T_current)[2] )
    rst = matrix(apply(T_current,1,sum))
    T_equal_cur[rst > 0,] = T_current[rst > 0,] / rst[rst > 0]
    
    # 0. create dummy for negative shocks
    asimdummy = rep(0, length(all_eps))
    asimdummy[all_eps <= 0] <- 1
    asim_dummy = asimdummy[ (dim(W)[1]+1) : dim(WT)[1] ]  # extract it
    
    # 1-a. effect of past squared errors (linear interpolation)
    W_eps_eq = T_equal %*% (all_eps^2)             # get linearly interpolataed error
    leps = W_eps_eq[ (dim(W)[1]+1) : dim(WT)[1] ]  # extract it
    aleps = leps * asim_dummy
    
    # 2-a. effect of past squared errors (spatial interpolation)
    W_eps = WTs %*% (all_eps^2)     # spatially interpolataed errors
    weps = W_eps[ (dim(W)[1]+1) : dim(WT)[1] ] # extract it
    aweps = weps * asim_dummy
    
    # 3-a. effect of past variance (equally weighted)
    l_h = T_equal %*% all_h                # linearly interpolated variance
    lh = l_h[ (dim(W)[1]+1) : dim(WT)[1] ] # extract it
    alh = lh * asim_dummy
    
    # 4-a. effect of past variance (spatialy weighted)
    w_h = WTs %*% all_h                    #  spatially interpolated variance
    wh = w_h[ (dim(W)[1]+1) : dim(WT)[1] ] # extract it
    awh = wh * asim_dummy
    
    # 1-b. effect of present squared errors (linear interpolation)
    Seeq = solve( diag(length(all_eps)) -(alpha3)*T_equal_cur -(beta3)*asimdummy*T_equal_cur)
    See = Seeq[ (dim(W)[1]+1) : dim(WT)[1],(dim(W)[1]+1) : dim(WT)[1] ] # extract it
    
    # 2-b. effect of present squared errors (spatial interpolation)
    Sewe = solve( diag(length(all_eps)) -(alpha4)*WTscur -(beta4)*asimdummy*WTscur)
    Sew = Sewe[ (dim(W)[1]+1) : dim(WT)[1],(dim(W)[1]+1) : dim(WT)[1] ] # extract it
    
    # 3-b. effect of present variance (equally weighted)
    Sveq = solve( diag(length(all_h)) -(gamma3)*T_equal_cur -(delta3)*asimdummy*T_equal_cur)
    Sve = Sveq[ (dim(W)[1]+1) : dim(WT)[1],(dim(W)[1]+1) : dim(WT)[1] ] # extract it
    
    # 4-b. effect of present variance (spatialy weighted)
    Svwe = solve( diag(length(all_h)) -(gamma4)*WTscur -(delta4)*asimdummy*WTscur)
    Svw = Svwe[ (dim(W)[1]+1) : dim(WT)[1],(dim(W)[1]+1) : dim(WT)[1] ] # extract it
    
    if (modelot %in% c("ARCH","SARCH", "ASARCH")) {
      h1 <- as.vector( See %*% Sew %*% (omega + alpha1*leps +alpha2*weps +gamma1*lh +gamma2*wh 
                                        + beta1*aleps + beta2 *aweps + delta1*alh + delta2*awh))  # ST-AGARCH
    }
    
    if (!modelot %in% c("ARCH","SARCH","ASARCH")) {
      h1 <- as.vector( Sve %*% Svw %*% (omega + alpha1*leps +alpha2*weps +gamma1*lh +gamma2*wh 
                                        + beta1*aleps + beta2 *aweps + delta1*alh + delta2*awh ))  # ST-AGARCH
    }
    
    xi <- eps1 / sqrt(h1)
    TS = length(xi)
    all = cbind( rep( as.character(vecym[i]), length(eps1)), coords1 ,eps1, h1, xi, leps, weps, lh, wh, aleps, aweps, alh, awh)
    colnames(all) <- colnames(HM)
    HM = rbind(HM, all)
    
    # store for the next step
    h      = h1
    
  }
  
  # standardize matrix of matrices
  # first cut the first n rows and last n columns
  cutoff = length(reziduali[extract$ym == vecym[1]])
  cutoff2 = length(reziduali[extract$ym == vecym[length(vecym)]])
  WW = WW[((cutoff+1):dim(WW)[1]), 1:(dim(WW)[2]-cutoff)   ]
  WWs = WW / rowSums(WW)
  WWs[WWs == "NaN"] <- 0
  
  WWpp = WWpp[((cutoff+1):dim(WWpp)[1]), 1:(dim(WWpp)[2]-cutoff)   ]
  WWpps = WWpp / rowSums(WWpp)
  WWpps[WWpps == "NaN"] <- 0
  
  HM = data.frame(HM)
  rownames(HM) = 1:dim(HM)[1]
  
  HM[, 2:dim(HM)[2]] = apply( HM[, 2:dim(HM)[2]], 2, function(y) as.numeric(as.character(gsub("%", "", y))))
  
  #name = paste0("./",modelot,"_",treshold,"_HM_mat.csv")
  #write.csv(HM, file=name)
  return(list(HM, pars, WWs, WWpps))
}
```

## QML functions

### Spatial ARCH

``` r
SARCH <- function(sarchpars, RES, eig, modelot, Wmat) {

  NT = dim(RES)[1]
  dum = rep(0,NT)
  dum[RES$eps <= 0] = 1 # dummy for negative shocks
  
  Y = RES$eps^2
  Z = as.matrix(data.frame( ones = rep(1,NT)))
  
  if (modelot == "SARCH") {
    omega     <- sarchpars[1]
    alpha     <- sarchpars[2]
    delta     <- 0
  }
  
  if (modelot == "ASARCH") {
    omega     <- sarchpars[1]
    alpha     <- sarchpars[2]
    delta     <- sarchpars[3]
  }

  rezid = Y - alpha* (Wmat %*% Y) - delta*dum*(Wmat %*% Y) - Z %*% omega
  Sigma_e = sqrt(mean(rezid^2))
  
  # logdet1 = log(det(diag(N)-alpha*Wmat -delta*dum*Wmat ) ) # slow log-determinant
  logdet1 = log(prod(1 - alpha*eig -delta*dum*eig))    # fast log-determinant
  
  A = diag(NT) -alpha*Wmat -delta*dum*Wmat
  AY = A %*% Y
  
  LL =  -1 *( -NT/2 *log(2*pi)
              -NT/2 *log( Sigma_e ^2 )
              +logdet1
              -1/2*( t(AY - Z %*% omega) %*% (AY - Z %*% omega) ) / (Sigma_e ^2)  )
  
  exit_pars = c(omega, alpha, delta)
  
  #opt_res = t(matrix(c(exit_pars, -LL)))
  #print(opt_res)
 
  return(LL)
  
}
```

### Spatial - Temporal GARCH: allowing past -> present effects

``` r
past_ASTAGARCH <- function(pars, RES, modelot, treshold) {
  
  if (modelot == "ARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "SARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "ASARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- pars[3]
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "GARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- pars[3]
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "SARCH-GARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- pars[3]
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "ARCH-SGARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- pars[3]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "SARCH-SGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- pars[3]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "AGARCH") {
    omega    <- pars[1]
    alpha1   <- pars[2]
    alpha2   <- 0
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- pars[3]
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- pars[4]
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "STAGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- 0
    gamma1   <- 0
    gamma2   <- pars[3]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- pars[4]
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  pars = c(omega, alpha1, alpha2, alpha3, alpha4,
           gamma1, gamma2 ,gamma3 ,gamma4,
           beta1, beta2, beta3, beta4,
           delta1, delta2, delta3, delta4)
  
  result = create_ASTAGARCH_past_serie(pars,RESID, modelot, treshold)
  
  HM      = result[[1]]
  vec_par = result[[2]]
  WM = as.matrix(result[[3]])
  
  NT = dim(HM)[1]
  
  dum = rep(0,NT)
  dum[HM$eps <= 0] = 1 # dummy for negative shocks
  
  Y = HM$eps^2
  X = data.frame(leps = HM$leps, weps = HM$weps, lh = HM$lh, wh= HM$wh, aleps = HM$aleps, aweps = HM$aweps, alh =HM$alh, awh =HM$awh )
  Z = as.matrix(data.frame( ones = rep(1,NT), X))
  
  betas = c(omega,alpha1,alpha2,gamma1,gamma2, beta1, beta2, delta1,delta2)
  
  rezid   = Y         -(alpha3 + alpha4)*(WM %*% Y) -(beta3+beta4)*dum*(WM %*% Y) -Z %*% betas
  # logdet1 = log(prod(1-(alpha3 + alpha4)*eig        -(beta3+beta4)*dum*eig))  # fast log-determinant
  A       = diag(NT)  -(alpha3 + alpha4)*WM         -(beta3+beta4)*dum*WM
  logdet1 = log(det(A)) # slow-log determinant
  
  
  Sigma_e = sqrt(mean(rezid^2))
  AY = A %*% Y
  
  LL =  -1 *( -NT/2 *log(2*pi)
              -NT/2 *log( Sigma_e ^2 )
              +logdet1
              -1/2*( t(AY - Z %*% betas) %*% (AY - Z %*% betas) ) / (Sigma_e ^2)  )
  
  opt_res = t(matrix(c(pars, -LL)))
  print(opt_res)

  return(LL)
}
```

### Spatial-Temporal GARCH: allowing past + present ->  present effects

``` r
pastpresent_ASTAGARCH <- function(pars, RES, modelot, treshold) {
  
  if (modelot == "SARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- pars[3]
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "ASARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- pars[3]
    gamma1   <- 0
    gamma2   <- 0
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- pars[4]
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "SARCH-SGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- pars[3]
    gamma1   <- 0
    gamma2   <- pars[4]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- 0
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  if (modelot == "STAGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    alpha3   <- 0
    alpha4   <- pars[3]
    gamma1   <- 0
    gamma2   <- pars[4]
    gamma3   <- 0
    gamma4   <- 0
    beta1    <- 0
    beta2    <- pars[5]
    beta3    <- 0
    beta4    <- 0
    delta1   <- 0
    delta2   <- 0
    delta3   <- 0
    delta4   <- 0
  }
  
  pars = c(omega, alpha1, alpha2, alpha3, alpha4,
           gamma1, gamma2 ,gamma3 ,gamma4,
           beta1, beta2, beta3, beta4,
           delta1, delta2, delta3, delta4)
  
  result = create_ASTAGARCH_pastpresent_serie(pars,RESID, modelot, treshold)
  
  HM      = result[[1]]
  vec_par = result[[2]]
  WM = as.matrix(result[[4]])
  
  NT = dim(HM)[1]
  
  dum = rep(0,NT)
  dum[HM$eps <= 0] = 1 # dummy for negative shocks
  
  Y = HM$eps^2
  X = data.frame(leps = HM$leps, weps = HM$weps, lh = HM$lh, wh= HM$wh, aleps = HM$aleps, aweps = HM$aweps, alh =HM$alh, awh =HM$awh )
  Z = as.matrix(data.frame( ones = rep(1,NT), X))
  
  betas = c(omega,alpha1,alpha2,gamma1,gamma2,beta1,beta2,delta1,delta2)
  
  
  #if (modelot %in% c("ARCH","SARCH","ASARCH")) {
  rezid   = Y         -(alpha3 + alpha4)*(WM %*% Y) -(beta3+beta4)*dum*(WM %*% Y) -Z %*% betas
  # logdet1 = log(prod(1-(alpha3 + alpha4)*eig        -(beta3+beta4)*dum*eig))  # fast log-determinant
  A       = diag(NT)  -(alpha3 + alpha4)*WM         -(beta3+beta4)*dum*WM
  logdet1 = log(det(A)) # slow-log determinant
  #}
  

  Sigma_e = sqrt(mean(rezid^2))
  AY = A %*% Y
  
  LL =  -1 *( -NT/2 *log(2*pi)
              -NT/2 *log( Sigma_e ^2 )
              +logdet1
              -1/2*( t(AY - Z %*% betas) %*% (AY - Z %*% betas) ) / (Sigma_e ^2)  )
  
  
  opt_res = t(matrix(c(pars, -LL)))
  print(opt_res)
 
  return(LL)
  
}
```

## GMM functions

### Function for creating moments for SARCH series

``` r
Fmoment_sarch  <- function(pars, RESID, modelot, Wmat){
  
  if (modelot == "SARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[3]
    gamma1   <- 0
    gamma2   <- 0
    delta1   <- 0
    delta2   <- 0
    LT_vol = omega /(1-alpha1-alpha2-gamma1-gamma2-delta1-delta2)
  }
  
  if (modelot == "ASARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[3]
    gamma1   <- 0
    gamma2   <- 0
    delta1   <- 0
    delta2   <- pars[7]
    LT_vol = omega /(1-alpha1-alpha2-gamma1-gamma2-delta1-delta2)
  }
  
  result = create_SARCH_serie(pars, RESID, modelot, Wmat)
  HM      = result[[1]]
  vec_par = result[[2]]
  e  = HM$eps
  h  = HM$h
  xi = HM$xi
  N = length(e)
  
  # new algo
  ex_eps = HM$weps
  
  # 1. Long term volatility = omega / (1-alpha-gamma)
  # F1 = e^2 - LT_vol
  F1 = (e^2-h)*ex_eps
  
  # 2. Mean of standardized residuals is 0
  # F2 = xi
  F2 = (e^2-h) # mean of residuals (realized shocks squared - predicted) is 0
  
  if (modelot %in% c("ASARCH")){
    # 3. Var of standardized residuals is 1
    F3 = (xi-mean(xi))^2 - 1
  }
  
  if (modelot %in% c("SARCH")){FM = cbind(F1, F2)}
  if (modelot %in% c("ASARCH")){FM = cbind(F1, F2, F3)}
  
  return(FM)
  
}
```

``` r
summary(Fmoment_sarch(pars=c(0.01,0.70), RESID=RESID_new, modelot="SARCH", Wmat=W))
```

    ##        F1                   F2           
    ##  Min.   :-0.3250316   Min.   :-0.479592  
    ##  1st Qu.:-0.0002793   1st Qu.:-0.017675  
    ##  Median :-0.0001069   Median :-0.012569  
    ##  Mean   :-0.0005442   Mean   :-0.003938  
    ##  3rd Qu.:-0.0000238   3rd Qu.:-0.002995  
    ##  Max.   : 0.0536273   Max.   : 3.672749

### Function for creating moments for ST-GARCH series

``` r
Fmoment  <- function(pars, RESID, modelot, treshold){
  
  if (modelot == "SARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    gamma1   <- 0
    gamma2   <- 0
    beta1    <- 0
    beta2    <- 0
    delta1   <- 0
    delta2   <- 0
    LT_vol = omega /(1-alpha1-alpha2-gamma1-gamma2-beta1-beta2-delta1-delta2)
  }
  
  if (modelot == "ASARCH") {
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    gamma1   <- 0
    gamma2   <- 0
    beta1    <- 0
    beta2    <- pars[3]
    delta1   <- 0
    delta2   <- 0
    LT_vol = omega /(1-alpha1-alpha2-gamma1-gamma2-beta1-beta2-delta1-delta2)
  }
  
  if (modelot == "SARCH-SGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    gamma1   <- 0
    gamma2   <- pars[3]
    beta1    <- 0
    beta2    <- 0
    delta1   <- 0
    delta2   <- 0
    LT_vol = omega /(1-alpha1-alpha2-gamma1-gamma2-beta1-beta2-delta1-delta2)
  }
  
  if (modelot == "STAGARCH"){
    omega    <- pars[1]
    alpha1   <- 0
    alpha2   <- pars[2]
    gamma1   <- 0
    gamma2   <- pars[3]
    beta1    <- 0
    beta2    <- pars[4]
    delta1   <- 0
    delta2   <- 0
    LT_vol = omega /(1-alpha1-alpha2-gamma1-gamma2-beta1-beta2-delta1-delta2)
  }
  
  result = create_ASTAGARCH_past_serie(pars, RESID, modelot, treshold)
  HM      = result[[1]]
  vec_par = result[[2]]
  e  = HM$eps
  h  = HM$h
  xi = HM$xi
  N = length(e)
  
  ex_eps = HM$weps
  
  # 2 params optimization
  if (length(pars)==2){
    F1 = (e^2-h)*ex_eps
    F2 = (e^2-h) 
    FM = cbind(F1, F2)
  }
  
  if (length(pars)==3){
    # 1.residuals (realized shocks squared - predicted) exogenious to X (previous shocks)
    F1 = (e^2-h)*ex_eps
    # 2. mean of residuals (realized shocks squared - predicted) is 0
    F2 = (e^2-h) #
    # 3. Var of standardized residuals is 1
    F3 = (xi-mean(xi))^2 - 1
    FM = cbind(F1, F2, F3)
  }
  
  if (length(pars)==4){
    # 1.residuals (realized shocks squared - predicted) exogenious to X (previous shocks)
    F1 = (e^2-h)*ex_eps
    # 2. mean of residuals (realized shocks squared - predicted) is 0
    F2 = (e^2-h) #
    # 3. Var of standardized residuals is 1
    F3 = (xi-mean(xi))^2 - 1
    # 4. Long term volatility = LTvol (omega / (1-alpha-gamma-...))
    F4 = e^2 - LT_vol
    FM = cbind(F1, F2, F3, F4)
  }
  return(FM)
}
```

``` r
create_matrices(treshold=2, RESID_new)
summary(Fmoment(pars=c(0.01,0.40, 0.30), RESID=RESID_new, modelot="SARCH-SGARCH", treshold=2))
```

    ##        F1                   F2                  F3          
    ##  Min.   :-0.4968026   Min.   :-0.448900   Min.   : -1.0000  
    ##  1st Qu.:-0.0003298   1st Qu.:-0.020809   1st Qu.: -0.9448  
    ##  Median :-0.0001434   Median :-0.016456   Median : -0.7581  
    ##  Mean   :-0.0006985   Mean   :-0.005958   Mean   : -0.1782  
    ##  3rd Qu.:-0.0000207   3rd Qu.:-0.004968   3rd Qu.: -0.2245  
    ##  Max.   : 0.0873368   Max.   : 3.665603   Max.   :136.0631

### Create average of moments for SARCH series

``` r
Gmoment_sarch <- function(pars, RESID, modelot,Wmat){
  Fs = Fmoment_sarch(pars, RESID, modelot,Wmat)
  Gs = t(matrix(apply(Fs, 2, mean)))
  return(Gs)
}
```

test it

``` r
Gmoment_sarch(pars=c(0.01,0.70), RESID=RESID_new, modelot="SARCH", Wmat=W)
```

    ##               [,1]         [,2]
    ## [1,] -0.0005441925 -0.003938163

### Function for average of moments for ST-GARCH series

``` r
Gmoment <- function(pars, RESID, modelot, treshold){
  Fs = Fmoment(pars, RESID, modelot, treshold)
  Gs = t(matrix(apply(Fs, 2, mean)))
  return(Gs)
}
```

### Function for quadratic form of average moment function for SARCH series

``` r
MinQ_sarch <- function(pars, RESID, modelot, W, A, Wmat){
  Gs = Gmoment_sarch(pars, RESID, modelot, Wmat)
  QFORM = (Gs %*% t(A)) %*% W %*% t(Gs %*% t(A)) 
  #print(pars)
  #print(QFORM)
  opt_res = t(matrix(c(pars, QFORM)))
  name = paste0("./",modelot, "_GMM_opt_par_vec.csv")
  write.csv(opt_res, file=name)
  return(QFORM)
}
```

test it

``` r
MinQ_sarch(pars=c(0.01,0.70), RESID=RESID_new, modelot="SARCH", W = as.matrix(diag(rep(1,2))), A = as.matrix(diag(rep(1,2))), Wmat=W)
```

    ##              [,1]
    ## [1,] 1.580527e-05

### Function for quadratic form of average moment function for ST-GARCH series

``` r
MinQ <- function(pars, RESID, modelot, W, A, treshold){
  Gs = Gmoment(pars, RESID, modelot, treshold)
  QFORM = (Gs %*% t(A)) %*% W %*% t(Gs %*% t(A)) 
  #print(pars)
  #print(QFORM)
  opt_res = t(matrix(c(pars, QFORM)))
  name = paste0("./",modelot, "_GMM_opt_par_vec.csv")
  write.csv(opt_res, file=name)
  return(QFORM)
}
```

## Estimation functions

``` r
SARCH_package <- function(RESID, treshold){
  
  out = no_null_neighbours(RESID, treshold)
  RESID_new = data.frame(out[[1]])
  no_neighbours = out[[2]]
  dim(RESID_new)
  length(no_neighbours)

  c = 0
  while( length(no_neighbours) > 0){
    c=c+1
    out = no_null_neighbours(RESID_new, treshold)
    RESID_new = data.frame(out[[1]])
    no_neighbours = out[[2]]
    print(dim(RESID_new)[1])
    print(length(no_neighbours))
    if  ( c == 20 | length(no_neighbours)  == 0) break;
  }
  
  RESID_new$eps2 = RESID_new$eps^2
  koordinats <- cbind(RESID_new$X, RESID_new$Y)
  Ns = dim(koordinats)[1]
  rownames(koordinats) = 1:Ns
  IDi <-  1:Ns
  dnbs     <- dnearneigh(koordinats, 0, treshold, row.names=IDi,longlat=TRUE)
  nbds     <- nbdists(dnbs, koordinats,longlat=TRUE);
  gls     <- lapply(nbds, function(x) 1/x); # gl
  stat_nb <- nb2listw(dnbs, glist=gls,zero.policy=T)
  #stat_W  <- nb2mat(stat_nb$neighbours,glist=gls, style="W", zero.policy=T)
  m_test = lagsarlm(eps2 ~ 1, data=RESID_new, stat_nb, tol.solve=1.0e-30)
  sum_test = summary(m_test, Nagelkerke=T)

  koef    = t(t(m_test$coefficients))
  se_koef = t(t(m_test$rest.se))
  t_koef   = koef/se_koef
  p_val = 2*pnorm(-abs(t_koef))
  adjR2 = round(sum_test$NK, digits=4)
  sar_AIC = AIC(sum_test)
  LL = as.numeric(sum_test$LL)
  rhoto = sum_test$rho
  p_sar = round(p_val,digits = 4)
  KOEF_SAR = cbind(treshold, round(koef,digits = 4), rhoto, adjR2, LL, p_sar, sum_test$Wald1$p.value)
  return(KOEF_SAR)
  
}
```

``` r
SARCH_qml <- function(RESID, treshold, pars, mo){
  
  out = no_null_neighbours(RESID, treshold)
  RESID_new = data.frame(out[[1]])
  no_neighbours = out[[2]]
  dim(RESID_new)
  length(no_neighbours)

  c = 0
  while( length(no_neighbours) > 0){
    c=c+1
    out = no_null_neighbours(RESID_new, treshold)
    RESID_new = data.frame(out[[1]])
    no_neighbours = out[[2]]
    print(dim(RESID_new)[1])
    print(length(no_neighbours))
    if  ( c == 20 | length(no_neighbours)  == 0) break;
  }
 
  koordinats <- cbind(RESID_new$X, RESID_new$Y)
  Ns = dim(koordinats)[1]
  rownames(koordinats) = 1:Ns
  IDi <-  1:Ns
  dnbs     <- dnearneigh(koordinats, 0, treshold, row.names=IDi,longlat=TRUE)
  nbds     <- nbdists(dnbs, koordinats,longlat=TRUE);# nbd
  gls     <- lapply(nbds, function(x) 1/x); # gl
  stat_nb <- nb2listw(dnbs, glist=gls,zero.policy=T)
  stat_W  <- nb2mat(stat_nb$neighbours,glist=gls, style="W", zero.policy=T)
  stat_eig <- eigenw(stat_nb)

  try( optimiser0 <- optim( pars, SARCH, RES=RESID_new, eig=stat_eig,
                            modelot=mo, Wmat= stat_W, hessian=TRUE), silent =T)
  opt_values0 =optimiser0$par
  ind = opt_values0 != 0
  opt_values0 = opt_values0[1:length(pars)]
  hes0 = optimiser0$hessian
  I0 = solve(hes0)
  se0 = try(sqrt(diag(I0)), silent=T)
  t_koef0   = opt_values0/se0;    t_koef0
  pvaln0 = round(2*pnorm(-abs(t_koef0)), digits = 4); pvaln0
  STATDF  = rbind( c(treshold,opt_values0,-optimiser0$value), c(treshold, pvaln0, ""))
  
  return(STATDF)
}
```

``` r
ASTAGARCH_past_QML = function(RESID, treshold, pars, mo) {
  
  out = no_null_neighbours(RESID, treshold)
  RESID_new = data.frame(out[[1]])
  no_neighbours = out[[2]]
  dim(RESID_new)
  length(no_neighbours)
  
  c = 0
  while( length(no_neighbours) > 0){
    c=c+1
    out = no_null_neighbours(RESID_new, treshold)
    RESID_new = data.frame(out[[1]])
    no_neighbours = out[[2]]
    print(dim(RESID_new)[1])
    print(length(no_neighbours))
    if  ( c == 20 | length(no_neighbours)  == 0) break; 
  }
  
  create_matrices(treshold, RESID_new)
  
  try( optimiser2 <- optim(pars, past_ASTAGARCH, RES=RESID_new, modelot=mo,
                             treshold = treshold, hessian=TRUE), silent =T)

  indb = optimiser2$par
  hes_b = optimiser2$hessian
  Ib = solve(hes_b)
  se_b = try(sqrt(diag(Ib)), silent=T)
  t_koefb = indb/se_b; t_koefb
  pvaln_b = round(2*pnorm(-abs(t_koefb)), digits = 4); pvaln_b

  rezultat = rbind( c(treshold, mo, indb, -optimiser2$value), c(treshold, mo, pvaln_b, ""))

  return(rezultat)
}
```

``` r
ASTAGARCH_pastpresent_QML = function(RESID, treshold, pars, mo) {
  
  out = no_null_neighbours(RESID, treshold)
  RESID_new = data.frame(out[[1]])
  no_neighbours = out[[2]]
  dim(RESID_new)
  length(no_neighbours)
  
  c = 0
  while( length(no_neighbours) > 0){
    c=c+1
    out = no_null_neighbours(RESID_new, treshold)
    RESID_new = data.frame(out[[1]])
    no_neighbours = out[[2]]
    print(dim(RESID_new)[1])
    print(length(no_neighbours))
    if  ( c == 20 | length(no_neighbours)  == 0) break; 
  }
  
  create_matrices(treshold, RESID_new)
  
  try( optimiser2 <- optim(pars, pastpresent_ASTAGARCH, RES=RESID_new, modelot=mo,
                             treshold = treshold, hessian=TRUE), silent =T)

  indb = optimiser2$par
  hes_b = optimiser2$hessian
  Ib = solve(hes_b)
  se_b = try(sqrt(diag(Ib)), silent=T)
  t_koefb = indb/se_b; t_koefb
  pvaln_b = round(2*pnorm(-abs(t_koefb)), digits = 4); pvaln_b

  rezultat = rbind( c(treshold, mo, indb, -optimiser2$value), c(treshold, mo, pvaln_b, ""))

  return(rezultat)
}
```

``` r
ASTAGARCH_past_GMM = function(RESID, treshold, pars, mo) {
  
  out = no_null_neighbours(RESID, treshold)
  RESID_new = data.frame(out[[1]])
  no_neighbours = out[[2]]
  dim(RESID_new)
  length(no_neighbours)
  
  c = 0
  while( length(no_neighbours) > 0){
    c=c+1
    out = no_null_neighbours(RESID_new, treshold)
    RESID_new = data.frame(out[[1]])
    no_neighbours = out[[2]]
    print(dim(RESID_new)[1])
    print(length(no_neighbours))
    if  ( c == 20 | length(no_neighbours)  == 0) break; 
  }
  
  create_matrices(treshold, RESID_new)
  
  {
    # optimize 2 params
    if ( mo %in% c("SARCH")) {
      W = as.matrix(diag(rep(1,2)))
      A = as.matrix(diag(rep(1,2)))
      pars=pars
    }
    
    # optimize 3 params
    if ( mo %in% c("ASARCH","SARCH-SGARCH") ) {
      W = as.matrix(diag(rep(1,3)))
      A = as.matrix(diag(rep(1,3)))
    }
    
    # optimize 4 params
    if (mo %in% c("STAGARCH")) {
      W = as.matrix(diag(rep(1,4)))
      A = as.matrix(diag(rep(1,4)))
      pars=pars
    }
    
    
  }
  
  try(GMMoptimiser <- optim( pars, MinQ, RESID=RESID_new, modelot=mo, W=W, A=A, treshold=treshold,
                             lower=rep(0.00001,length(pars)), upper=rep(0.99999,length(pars)), method=c("L-BFGS-B")), silent=T)
  
  if(! class(GMMoptimiser)=="try-error") {
    opt_values = GMMoptimiser$par
    NT = dim(RESID_new)[1]
    Fm = Fmoment(opt_values, RESID_new, mo, treshold)
    S2 = (t(Fm) %*% Fm) / NT
    D = as.matrix(diag(rep(-1,dim(Fm)[2])))
    W2 = solve(S2)
    CovB    = solve(t(D) %*% solve(S2) %*% D)/NT
    Se = diag(sqrt(CovB)); Se
    Tstat = opt_values[opt_values >0]/Se ; Tstat
    Pval = 2*pnorm(-abs(Tstat)); Pval
    
    GMMoptimiser2 <- optim( opt_values, MinQ, RESID=RESID_new, modelot=mo, W=W2, A=A, treshold=treshold,
                            lower=rep(0.00001,length(pars)), upper=rep(0.99999,length(pars)),
                            method=c("L-BFGS-B"))
    
    if(! class(GMMoptimiser2)=="try-error") {
      opt_values = GMMoptimiser2$par
      NT = dim(RESID_new)[1]
      Fm = Fmoment(opt_values, RESID_new, mo, treshold)
      S2 = (t(Fm) %*% Fm) / NT
      D = as.matrix(diag(rep(-1,dim(Fm)[2])))
      W2 = solve(S2)
      CovB    = solve(t(D) %*% solve(S2) %*% D)/NT
      Se = diag(sqrt(CovB)); Se
      Tstat = opt_values[opt_values >0]/Se ; Tstat
      Pval = 2*pnorm(-abs(Tstat)); Pval
      
      rez = rbind( c( opt_values, -GMMoptimiser2$value, Pval))
      
      rm(pars,opt_values)
      
      return(rez)
    }
  } 
  
}
```

# Estimation

Estimate using package

``` r
out= SARCH_package(RESID, treshold=2)
```
``` r
out
```

    ##             treshold            rhoto  adjR2       LL             
    ## (Intercept)        2 0.0162 0.1602987 0.0015 2798.009 0 0.03791738

Choose between two models: “SARCH” and “ASARCH”. Initialize with these
parameters, respectively: c(0.01, 0.15) or c(0.01, 0.10, 0.10))

``` r
out= SARCH_qml(RESID, treshold=2, pars=c(0.01, 0.15), mo="SARCH")
```
``` r
out
```

    ##      [,1] [,2]                 [,3]                [,4]              
    ## [1,] "2"  "0.0161600361997262" "0.160152593201492" "2798.00892236926"
    ## [2,] "2"  "0"                  "0.0417"            ""

Estimate using QML: past -> present

``` r
out= ASTAGARCH_past_QML(RESID, treshold=2, pars=c(0.01, 0.5), mo="SARCH")
```


``` r
out
```

    ##      [,1] [,2]    [,3]                 [,4]                 [,5]              
    ## [1,] "2"  "SARCH" "0.0186381273335633" "0.0292923924625195" "2690.02483815309"
    ## [2,] "2"  "SARCH" "0"                  "0.4957"             ""






