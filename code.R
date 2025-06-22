### SDPDmod R Package :Spatial Dynamic Panel Data modeling
rm(list=ls())

iopt<-options()


## Package and options
library("SDPDmod")
options(prompt = "R> ", continue = "+  ", width = 70,
  useFancyQuotes = FALSE)


## Data
data("Cigar", package = "plm")
data("usa46", package = "SDPDmod") ## binary contiguity matrix of 46 USA states

W <- rownor(usa46) ## row-normalization

fm <- logc ~ logp + logy

Cigar$logc <- log(Cigar$sales)
Cigar$logp <- log(Cigar$price/Cigar$cpi)
Cigar$logy <- log(Cigar$ndi/Cigar$cpi)
Cigar$lpm <- log(Cigar$pimin/Cigar$cpi)



#############################################################
### Bayesian posterior probabilities                      ###
res1<-blmpSDPD(formula = fm, data = Cigar, W = W,
               index = c("state","year"),
               model = list("sar","sdm","sem","sdem"),
               effect = "individual")
res1

res2<-blmpSDPD(formula = fm, data = Cigar, W = W,
               index = c("state","year"),
               model = list("sar","sdm","sem","sdem"),
               effect = "individual",
               dynamic = TRUE,
               prior = "beta")

#############################################################
### Spatial panel data modeling                          ###
mod1 <- SDPDm(formula = fm,
              data = Cigar,
              index = c("state","year"),
              effect = "individual",
              model = "sar",
              LYtrans = FALSE,
              W = W)
summary(mod1)

mod2<-SDPDm(formula = fm,
            data = Cigar,
            index = c("state","year"),
            effect = "twoways",
            model = "sdm",
            W = W,
            LYtrans = TRUE,
            dynamic = TRUE,
            tlaginfo = list(ind = NULL, tl = TRUE, stl = TRUE))

imp<-impactsSDPDm(mod1)
summary(imp)

#############################################################
### Spatial weights matrices                              ###
library("sf")
ger <- st_read(system.file(dsn = "shape/GermanyNUTS3.shp",
                           package = "SDPDmod"),
               quiet = TRUE)	
data(gN3dist, package = "SDPDmod")

sn <- ger[which(substr(ger[["NUTS_CODE"]],1,3)=="DED"),]
sn_dist <- gN3dist[which(substr(rownames(gN3dist),1,3)=="DED"),
                   which(substr(colnames(gN3dist),1,3)=="DED")]

W_1on<-mOrdNbr(sf_pol=sn, m=1)
W_2on<-mOrdNbr(sf_pol=sn, m=2)
W_len_sh<-SharedBMat(sf_pol=sn)
W_nn<-mNearestN(distMat=sn_dist, m=5)
W_inv<-InvDistMat(distMat=sn_dist, distCutOff=100000, powr=2)
W_exp<-ExpDistMat(distMat=sn_dist, distCutOff=100000, expn=0.001)
W_dd<-DDistMat(distMat=sn_dist, distCutOff=100000, powr=3) 

################################################################
####
##Reset options
options(iopt)


#############################################################
#### Numeric comparison to other software                 ###
##### Bayesian posterior probabilities                 ######
r1<-blmpSDPD(formula = fm, data = Cigar, W = W,
             index = c("state","year"),
             model = list("sar","sdm","sem","sdem"),
             effect = "twoways")
r1
##### Spatial static panel data modeling               ######
m1 <- SDPDm(formula = fm,
            data = Cigar,
            index = c("state","year"),
            effect = "twoways",
            model = "sar",
            LYtrans = FALSE,
            W = W)
summary(m1)
i1 <- impactsSDPDm(m1, sd = 12345)
summary(i1)
m1$sige
m1$likl

library("splm")
library("spdep")
m2 <- spml(formula = fm,
           data = Cigar,
           index = c("state","year"),
           listw = mat2listw(W, style = "W"),
           model = "within",
           effect = "twoways",
           lag = TRUE,
           spatial.error = "none")
summary(m2)
set.seed(12345)
i2 <- impacts(m2, 
              listw = mat2listw(W, style = "W"), 
              time = length(unique(Cigar$year)))
s2 <- summary(i2, zstats=TRUE, short=T)
s2
m2$sigma2
m2$logLik

##### Spatial dynamic panel data modeling              ######

m3<- SDPDm(formula = fm,
           data = Cigar,
           index = c("state","year"),
           effect = "individual",
           model = "sdm",
           dynamic = TRUE,
           tlaginfo = list(ind = NULL, tl = TRUE, stl = TRUE),
           W = W,
           LYtrans = TRUE)
summary(m3)
m3$sige1
m3$likl1
########################################################################
## tables
## tab1
t1<-rbind(sprintf(r1$lmarginal, fmt = '%#.1f'),  
          sprintf(r1$probs, fmt = '%#.4f'))
t1<-as.data.frame(t1)
colnames(t1)<-names(r1$lmarginal)

## tab2
cf1<-sprintf(c(m1$rho, m1$coefficients), fmt = '%#.5f')
se1<-sprintf(c(m1$rho.se,m1$std), fmt = '%#.5f')

cf2<-sprintf(c(m2$coefficients[1], m2$coefficients[2:length(m2$coefficients)]), fmt = '%#.5f')
se2<-sprintf(c(sqrt(diag(m2$vcov))[1],sqrt(diag(m2$vcov))[2:length(m2$coefficients)]), fmt = '%#.5f')

c1<- c(paste0(cf1," (",se1,")"), 
       paste0(sprintf(i1$DIRECT.tab[,1], fmt = '%#.5f')," (",sprintf(i1$DIRECT.tab[,2], fmt = '%#.5f'),")"), 
       paste0(sprintf(i1$INDIRECT.tab[,1], fmt = '%#.5f')," (",sprintf(i1$INDIRECT.tab[,2], fmt = '%#.5f'),")"),
       paste0(sprintf(i1$TOTAL.tab[,1], fmt = '%#.5f')," (",sprintf(i1$TOTAL.tab[,2], fmt = '%#.5f'),")"),
       sprintf(m1$sige, fmt = '%#.4f'),
       sprintf(m1$likl, fmt = '%#.1f')
)
c2<- c(paste0(cf2," (",se2,")"), 
       paste0(sprintf(i2$res$direct, fmt = '%#.5f')," (",sprintf(s2$semat[,1], fmt = '%#.5f'),")"),
       paste0(sprintf(i2$res$indirect, fmt = '%#.5f')," (",sprintf(s2$semat[,2], fmt = '%#.5f'),")"),
       paste0(sprintf(i2$res$total, fmt = '%#.5f')," (",sprintf(s2$semat[,3], fmt = '%#.5f'),")"),
       sprintf(m2$sigma2, fmt = '%#.4f'),
       sprintf(m2$logLik, fmt = '%#.1f')
)
t2<-cbind(c1,c2)
t2<-as.data.frame(t2)

## tab3
cf3<-sprintf(c(m3$rho1, m3$coefficients1), fmt = '%#.5f')
se3<-sprintf(c(m3$rho.se1, m3$std1), fmt = '%#.5f')

c3<- c(paste0(cf3," (",se3,")"), 
       sprintf(m3$sige1, fmt = '%#.4f'),
       sprintf(m3$likl1, fmt = '%#.1f')
)
t3<-as.data.frame(c3)


#############################################################
### Session Info                                          ###
sessionInfo()
