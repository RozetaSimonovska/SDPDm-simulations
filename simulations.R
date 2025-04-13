rm(list=ls())
##############################################################################
## Functions
###################################
## Data Generating process
dgp <- function(seedst = 12345,          ## seed start
                t = 8,                   ## time
                W = rownor(usa46),       ## spatial weights matrix
                Z = list("rn","rn","ru"),## type of distribution (rm=normal, ru=uniform)
                Zpar = list(c(mu=0.1,sd=0.9),
                            c(-0.05,1),
                            c(-0.2,0.81)),## distribution parameters
                sderr = 1,               ## sd for error term
                betas = c(0,-1,1.4,0.8,0.3,2,-0.6), ## beta coefficients 
                rho = 0.65,              ## spatial coefficient
                SDM = TRUE,              ## spatial Durbin model if TRUE, otherwise SAR
                dynamic = TRUE,          ## dynamic model if TRUE, otherwise static model
                dyncoef = c(0.3,-0.255)  ## dynamic model coefficients
) {
  
  set.seed(seedst)
  
  n <- dim(W)[1]
  k <- length(Z)
  
  ind <- rep(1:n, t)
  tind <- rep(1:t, each=n)
  
  zm <- matrix(NA, nrow = n*t, ncol = k)
  for(i in 1:k){
    if(Z[[i]]=="rn"){
      zm[,i] <- rnorm(n*t, mean = Zpar[[i]][1], sd = Zpar[[i]][2])
    }else if(Z[[i]]=="ru"){
      zm[,i] <- runif(t*n, min = Zpar[[i]][1], max = Zpar[[i]][2])
    }
  }
  
  colnames(zm) <- paste0("X",seq(1,k))
  
  Wnt <- kronecker(diag(t), W)
  
  e <- rnorm(n*t, mean = 0, sd = sderr)
  
  y <- rep(NA, n*t)
  if(SDM){
    y0 <- betas[1] + rowSums(t(t(zm)*betas[2:(k+1)])) + 
      rowSums(Wnt%*%t(t(zm)*betas[(k+2):(2*k+1)])) + e
  }else{
    y0 <- betas[1] + rowSums(t(t(zm)*betas[2:(k+1)])) + e
  }
  for(i in 1:t) {
    yt <- y0[1:n + (i-1)*n]
    if(dynamic){
      if(i>1){
        ylag <- y0[1:n + (i-2)*n]
      }else{
        ylag <- rnorm(n)
      }
      y[1:n + (i-1)*n] <- as.vector(
        solve(diag(n) - rho*W,
              (diag(dyncoef[1],n)+dyncoef[2]*W)%*%ylag+yt) )
    }else{
      y[1:n + (i-1)*n] <- as.vector(solve(diag(n) - rho*W, yt))
    }
  }
  
  datain <- data.frame(ind, tind, zm, y)
  
  datain <- datain[order(datain$ind, datain$tind),]
  
  return(datain)
}
###############################################
## Simulation
sim <- function(seedstart = 12345,            ## seed start 
                nsim = 200,                   ## number of simulations
                ncores = 4,
                Time = 50,                    ## time
                Wmat = rownor(usa46),         ## spatial weight matrix
                modeltype = "sdm",            ## model type ("sdm","sar")
                covDistPar ,                  ## distributional parameters
                rhoval ,                      ## spatial coefficient
                betaval ,                     ## beta coefficients
                dynamicmodel = TRUE,          ## dynamic model if TRUE
                dyncval                       ## dynamic coefficients
){
  seeds <- seq(seedstart, seedstart+nsim*11, 11)
  ncov <- length(covDistPar)
  if(ncov==2){
    fm <- y~X1+X2
    ZDist <- list("rn","ru")
  } else if(ncov==3){
    fm <- y~X1+X2+X3
    ZDist <- list("rn","rn","ru")
  } else{
    stop("Only 2 and 3 covariates are supported!")
  }
  if(modeltype == "sar"){ 
    nbpar <- 1 
    sdmmod <- FALSE
    if((ncov+1)!=length(betaval)) 
      stop("The model or the length of the distribution parameters and the betas differs!")
  } else if(modeltype == "sdm"){ 
    nbpar <- 2 
    sdmmod <- TRUE
    if((2*ncov+1)!=length(betaval)) 
      stop("The model or the length of the distribution parameters and the betas differs!")
  } else { stop("Only 'sar' and 'sdm' models are supported!")  }
  if(dynamicmodel){ 
    dynpar <- 2 
    if(is.null(dyncval)) 
      stop("No values entered for time and space-time lag coefficients!")
  } else { dynpar <- 0 }
  
  datas <- mod <- vector("list", nsim)
  rhos <- rep(NA, nsim)
  coefm <- matrix(NA, nrow = nsim, ncol = nbpar*ncov+dynpar)
  for(i in 1:nsim){
    datas[[i]] <- dgp(seedst = seeds[i],
                      t = Time,
                      W = Wmat, 
                      Z = ZDist,
                      Zpar = covDistPar,
                      betas = betaval, 
                      rho = rhoval,
                      SDM = sdmmod,
                      dynamic = dynamicmodel,
                      dyncoef = dyncval
    )
  }
  if(ncores==1){
    for(i in 1:nsim){
      mod[[i]] <- SDPDmod::SDPDm(
        formula = fm,
        data = datas[[i]],
        W = Wmat,
        index = c("ind","tind"),
        model = modeltype,
        effect = "twoways",
        LYtrans = TRUE,
        dynamic = dynamicmodel,
        tlaginfo = list(ind=NULL, tl=TRUE, stl=TRUE)
      )
    }
  }else if(ncores>1){
    modi <- vector("list",(nsim/ncores))
    for(i in 1:(nsim/ncores)){
      modi[[i]]<- foreach(j = 1:ncores) %dopar%  {
        SDPDmod::SDPDm(
          formula = fm,
          data = datas[[(i-1)*ncores+j]],
          W = Wmat,
          index = c("ind","tind"),
          model = modeltype,
          effect = "twoways",
          LYtrans = T,
          dynamic = dynamicmodel,
          tlaginfo = list(tl=TRUE, stl=TRUE)
        )
      }
    }
    indm<-0
    for(i in 1:(nsim/ncores)){
      for(j in 1:ncores){
        indm<-indm+1
        mod[[indm]]<-modi[[i]][[j]]
      }}
  }
  for(i in 1:nsim){
    ss <- summary(mod[[i]])
    if(dynamicmodel){
      rhos[i] <- ss$rho1
      coefm[i,]<-as.vector(ss$coefficients1)
      
      actualv <- c(rhoval, dyncval, betaval[2:length(betaval)])
    } else{
      rhos[i] <- ss$rho
      coefm[i,]<-as.vector(ss$coefficients)
      
      actualv <- c(rhoval, betaval[2:length(betaval)])
    }
  }
  rm<-mean(rhos)
  cm<-colMeans(coefm)
  
  predictedv <- cbind(rhos, coefm)
  
  bias<-rmse<-rep(NA,length(actualv))
  for(i in 1:length(actualv)){
    bias[i]<-mean(predictedv[,i])-actualv[i]
    rmse[i]<-sqrt(mean((predictedv[,i]-actualv[i])^2))
  }
  
  return(list(seeds= seeds,
              # data=datas,
              mods=mod,
              rho=rhos, 
              beta = coefm,
              meanrho=rm,
              meanbetas=cm,
              actualval=actualv,
              predictedval=predictedv,
              bias=bias,
              rmse=rmse))
}
##############################################################################
library("spdep")
library("SDPDmod")
library("sf")
library("dplyr")
# library(foreach)
# library(doParallel)

ger   <- st_read(system.file("shape/GermanyNUTS3.shp",
                             package = "SDPDmod"),
                 quiet = TRUE)
nc    <- st_read(system.file("shape/nc.shp", 
                             package = "sf"),
                 quiet = TRUE)
gern1 <- ger %>% 
  mutate(NUTS1=substr(NUTS_CODE,1,3)) %>%
  group_by(NUTS1) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Wgn1     <- mOrdNbr(gern1)
Wnc     <- mOrdNbr(nc)
W <- vector("list",3)
W[[1]] <-rownor(Wgn1)
W[[2]] <- rownor(usa46)
W[[3]] <- rownor(Wnc)
times <- c(10,50,100)
r <- 0.4; b <- c(1.4,-0.8,1.1,-0.5); d <- c(0.3,-0.2)
###
s <- vector("list",length(W)*length(times))
nsimulations <- 1000
nvec <- tvec <- rep(NA, length(W)*length(times))
res.b <- res.r <- matrix(NA, 
                         nrow = length(W)*length(times), 
                         ncol = (1+length(b)+2))
###########################
ncrs <- 1
# registerDoParallel(ncrs) 
ind <- 0
start_time <- Sys.time()
for(i in 1:length(W)){
  for(j in 1:length(times)){
    ind <- ind + 1
    print(ind)
    s[[ind]] <- sim(seedstart = 12345+(ind)*(nsimulations+1)*11,
                    ncores = ncrs,
                    nsim = nsimulations,
                    Time = times[j],
                    Wmat = W[[i]],
                    modeltype = "sdm",
                    covDistPar = list(
                      c(-0.05,1),
                      c(-0.2,0.81)),
                    rhoval = r,
                    betaval = c(0,b),
                    dynamicmodel = TRUE,
                    dyncval = d)
    nvec[ind] <- nrow(W[[i]])
    tvec[ind] <- times[j]
    res.b[ind,] <- s[[ind]]$bias
    res.r[ind,] <- s[[ind]]$rmse
  }
}
end_time <- Sys.time()
durtime <- end_time - start_time
durtime

###########
rows.combined <- nrow(res.b) 
cols.combined <- ncol(res.b) + ncol(res.r)
res.br <- matrix(NA, nrow=rows.combined, ncol=cols.combined)
res.br[, seq(1, cols.combined, 2)] <- res.b
res.br[, seq(2, cols.combined, 2)] <- res.r
res <- cbind(nvec,tvec,res.br)
df <- as.data.frame(res) %>% 
  mutate_at(vars(nvec,tvec), list(~round(., 0)))%>% 
  mutate_at(vars(-nvec,-tvec), list(~round(., 4)))
