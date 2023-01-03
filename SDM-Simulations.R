#########################################################################################.
# Inferring single- and multi-species distributional aggregation using quadrat sampling #
# This code is used to do numerical and empirical tests                                 # 
#########################################################################################.
# (1) Help functions                                                                 ----

# CECI()       - Clark and Evans Competition Index
# DC()         - Deviation Coefficient or Diffusion Coefficient 
# species.distribution() - Simulation of species distribution using Poisson cluster process
# r2.test      - calculate r2, RMSE, NRMSD based on actual values and predicted values
source("./functions/Spatial distribution Index.R")

# rdirichlet  - random variate generation
# rdirichlet1 - random variate generation (avoid NAs)
# dDir        - probability density (not mass!) function of Dirichlet distribution
# dMDir       - probability mass (or density) function of the multinomial-Dirichlet distribution
# dNBD        - probability mass function for independent negative binomial model 
# dNMD        - probability mass function for negative multinomial model
# rMDir       - simulation of multinomial-Dirichlet distribution
# rMDir1      - simulation of multinomial-Dirichlet distribution by avoiding NAs
# rand        - random matrices generation
# likelihood  - calculation of the negative log-likelihood function for the Dirichlet-Multinomial model
# likelihood0 - calculation of the negative log-likelihood function for the null model: multinomial model
# fit         - fitting of the SDM model
# fitNBD      - fitting of the independent NBD model
# fitNMD      - fitting of the NMD model
# fitMD       - fitting of the ordinary multinomial model

source("./functions/Dirichlet-Multinomial distribution-S1.R")

# ---------------------------------------------------------------------------------- ----  
# (2) Numerical test I                                                               ----

#simulat SDM data, fit SDM and NBD models:
#the first set of simulations
#20 species in 10 quadrats
a=c(0.01,0.05,0.1,0.5,1,5,10,50)
m=s=r=lik=vector()
par(mfrow=c(3,3))
index_res_m <- data.frame()
for(j in 1:length(a))
{
  v=v1=v2=vc=va=v1a=v2a=vector()
  for(i in 1001:2000)
  {
    set.seed(i + j^2)
    dat=rMDir1(20,rep(a[j],1,10),sample(100:10000,20))
    res=fitNBD(dat)
    res1=fit(dat) #fit the standard species-site matrix
    v1=c(v1,res$par[1])
    v=c(v,res1$par[1])
    #square error:
    v1a=c(v1a,(res$par[1]-a[j])^2)
    va=c(va,(res1$par[1]-a[j])^2)
    #AIC values:
    vc=rbind(vc,c(res1$objective*2+2*1,res$objective*2+2*2))
    # DC
    dat2 = apply(dat,2,sum) 
    dc = DC(dat2)
    dc_mean = dc$dc
    # index_res
    index_res_m <- rbind(index_res_m, data.frame("initial_k" = a[j] , "SDM_alpha" = res1$par, "NBD_k" = res$par[1], 
                                                 "NBD_mu" = res$par[2], "DC" = dc_mean))
    
  }#i
  #estimation bias:
  m=rbind(m,c(mean(v),mean(v1))-a[j])
  #sample error:
  s=rbind(s,c(sd(v),sd(v1)))
  #root mean square error:
  r=rbind(r,c(sqrt(mean(va)),sqrt(mean(v1a))))
  #
  lik=rbind(lik,vc)
  #
  #windows()
  boxplot(cbind("SDM (¦Á)" = v,"NBD (k)" = v1),
          main = paste("Preset ¦Á = ",a[j],sep=""),
          xlab = "Fitted value",
          ylab = "Model",
          col = c("orange","grey"),
          border = "black",
          horizontal = TRUE,
          notch = TRUE,
          outline = F )
  abline(v = a[j], col = "red")
}#j

# BIC_SDM < BIC_NBD - 2
scales::percent(length(which(lik[,1]<lik[,2]-2))/(length(a)*1000), 1)

#
colnames(m)=colnames(s)=colnames(r)=c("SDM","NBD")
rownames(m)=rownames(s)=rownames(r)=a
#
m = round(m,4)
s = round(s,4)
r = round(r,4)
#
msr=cbind(m,s,r)
colnames(msr)=c("bias.SDM","bias.NBD",
                "se.SDM","se.NBD",
                "rmse.SDM","rmse.NBD")
rownames(msr)=a
msr
#
plot(lik[,1:2],xlim=c(0,4000),ylim=c(0,4000),
     xlab="BICs for SDM model",ylab="BICs for NBD model")
lines(1:10000,1:10000,col=2)

#very exciting, NBD is not good!!! #as variance is larger!!!


#simulat NBD data, fit SDM and NBD models:
#the second set of simulations
#1 species in 200 quadrats
k=c(0.01,0.05,0.1,0.5,1,5,10,50)
par(mfrow=c(3,3))
m=s=r=lik=vector()
index_res <- data.frame()
for(j in 1:length(k))
{
  v=v1=v2=vc=va=v1a=v2a=vector()
  for(i in 1001:2000)
  {
    set.seed(i + j^2)
    dat=rnbinom(200,size=k[j],mu=sample(100:10000,1))
    res=fitNBD(dat)
    res1=fit(dat)
    v1=c(v1,res$par[1])
    v=c(v,res1$par[1])
    #square error:
    v1a=c(v1a,(res$par[1]-k[j])^2)
    va=c(va,(res1$par[1]-k[j])^2)
    # AIC values:
    # AIC = - 2*logLik + k * edf
    # k = 2 corresponds to the traditional AIC, using k = log(n) provides the BIC (Bayesian IC) instead.
    vc=rbind(vc,c(res1$objective*2+log(length(dat))*1,res$objective*2+log(length(dat))*2))
    # DC
    dc = DC(dat)
    
    # index_res
    index_res <- rbind(index_res, data.frame("initial_k" = k[j] , "SDM_alpha" = res1$par, "NBD_k" = res$par[1], 
                                             "NBD_mu" = res$par[2], "DC" = dc$dc))
    
  }#i
  
  m=rbind(m,c(mean(v),mean(v1))-k[j])
  s=rbind(s,c(sd(v),sd(v1)))
  #root mean square error:
  r=rbind(r,c(sqrt(mean(va)),sqrt(mean(v1a))))
  lik=rbind(lik,vc)
  #
  #windows()
  boxplot(cbind("SDM (¦Á)" = v,"NBD (k)" = v1),
          main = paste("Preset k = ",k[j],sep=""),
          xlab = "Fitted value",
          ylab = "Model",
          col = c("orange","grey"),
          border = "black",
          horizontal = TRUE,
          notch = TRUE,
          outline = F )
  abline(v = k[j], col = "red")
}#j

# BIC_SDM < BIC_NBD 
scales::percent(length(which(lik[,1]<lik[,2]-2))/(length(k)*1000), 1)

#
colnames(m)=colnames(s)=colnames(r)=c("SDM","NBD")
rownames(m)=rownames(s)=rownames(r)=k
#
m = round(m,4)
s = round(s,4)
r = round(r,4)
#
msr=cbind(m,s,r)
colnames(msr)=c("bias.SDM","bias.NBD",
                "se.SDM","se.NBD",
                "rmse.SDM","rmse.NBD")
rownames(msr)=k
msr
#
plot(lik[,c(1,2)],xlim=c(0,4500),ylim=c(0,4500),
     xlab="BICs for SDM model",ylab="BICs for NBD model")
lines(1:10000,1:10000,col=2)
#
#very exciting, NBD is not good!!! as variance is larger!!!



par(mfrow=c(2,2))
index_res_1 <- by(index_res[,"SDM_alpha"],index_res[,"initial_k"], mean)
index_res_2 <- by(index_res[,"NBD_k"],index_res[,"initial_k"], mean)
index_res_3 <- by(index_res[,"DC"], index_res[,"initial_k"], mean)
index_res_33 <- by(index_res[,"DC"],index_res[,"initial_k"], sd)
index_res_333 <- table(index_res$initial_k)
index_res_3333 <- index_res_33/sqrt(index_res_333)



x1 <- summary(lm(log(index_res_3) ~ log(index_res_1)))
x2 <- summary(lm(log(index_res_3) ~ log(index_res_2)))
mod=lm(log(index_res_3) ~ log(index_res_1))
mod2=lm(log(index_res_3) ~ log(index_res_2))
# plot(log(index_res_1),log(index_res_3), xlab = "log(SDM_alpha)", ylab = "log(DC)", main = paste0("Random abundance data (1 species in 200 quadrats)"))
# lines(log(index_res_1),predict(mod,newdata=list(log(index_res_1))),col= "blue")
# text(1.9,12,paste0("adj R-squred = ", round(x1$adj.r.squared,3)),pch=1,cex=1)


# plot(log(index_res_2),log(index_res_3), xlab = "log(NBD_k)", ylab = "log(DC)", main =  paste0("Random abundance data (1 species in 200 quadrats)"))
# lines(log(index_res_2),predict(mod2,newdata=list(log(index_res_2))),col="red")
# text(1.9,12,paste0("adj R-squred = ", round(x2$adj.r.squared,3)),pch=1,cex=1)

index_res_m1 <- by(index_res_m[,"SDM_alpha"],index_res_m[,"initial_k"], mean)
index_res_m2 <- by(index_res_m[,"NBD_k"],index_res_m[,"initial_k"], mean)
index_res_m3 <- by(index_res_m[,"DC"], index_res_m[,"initial_k"], mean)
index_res_m33 <- by(index_res_m[,"DC"],index_res_m[,"initial_k"], sd)
index_res_m333 <- table(index_res_m$initial_k)
index_res_m3333 <- index_res_m333/sqrt(index_res_3333)


x3 <- summary(lm(log(index_res_m3) ~ log(index_res_m1)))
x4 <- summary(lm(log(index_res_m3) ~ log(index_res_m2)))
mod=lm(log(index_res_m3) ~ log(index_res_m1))
mod2=lm(log(index_res_m3) ~ log(index_res_m2))
# plot(log(index_res_m1),log(index_res_m3), xlab = "log(SDM_alpha)", ylab = "log(DC)", main = paste0("Random abundance data (20 species in 10 quadrats)"))
# lines(log(index_res_m1),predict(mod,newdata=list(log(index_res_m1))),col= "blue")
# text(1.9,8,paste0("adj R-squred = ", round(x3$adj.r.squared,3)),pch=1,cex=1)

# plot(log(index_res_m2),log(index_res_m3), xlab = "log(NBD_k)", ylab = "log(DC)", main =  paste0("Random abundance data (20 species in 10 quadrats)"))
# lines(log(index_res_m2),predict(mod2,newdata=list(log(index_res_m2))),col="red")
# text(0,8,paste0("adj R-squred = ", format(round(x4$adj.r.squared, digits = 3), nsmall = 3)),pch=1,cex=1)

df1 <- data.frame("DC" = unlist(list(index_res_3)), "SE" = unlist(list(index_res_3333)),"SDM_alpha" = unlist(list(index_res_1)), "type" = rep("Single", length(index_res_1)) )  
df2 <- data.frame("DC" = unlist(list(index_res_m3)), "SE" = unlist(list(index_res_m3333)),"SDM_alpha" = unlist(list(index_res_m1)), "type" = rep("Multiple", length(index_res_m1)) )
df_sdm <- rbind(df1,df2)
  
df3 <- data.frame("DC" = unlist(list(index_res_3)), "SE" = unlist(list(index_res_3333)),"NBD_k" = unlist(list(index_res_2)), "type" = rep("Single", length(index_res_2)) )  
df4 <- data.frame("DC" = unlist(list(index_res_m3)),"SE" = unlist(list(index_res_m3333)), "NBD_k" = unlist(list(index_res_m2)), "type" = rep("Multiple", length(index_res_m2)) )
df_nbd <- rbind(df3,df4)


###Plot data
library(ggplot2)
p1 <- ggplot(df_sdm, aes(x = log(SDM_alpha), y = log(DC), colour = type)) +
    xlab("log(SDM_alpha)") +
    ylab("log(DC)") +
    stat_smooth(method = 'lm', formula = 'y~x', se=T) +
   geom_point(size=4, pch=21,color = "black", stroke=1.5, aes(fill=type)) +
   theme_test()   
   # geom_errorbar(aes(ymin=log(DC)-log(SE), ymax=log(DC)+log(SE)),colour= "blue", width=.03,size=0.75)
p1
 
p2 <- ggplot(df_nbd, aes(x = log(NBD_k), y = log(DC), colour = type)) +
   xlab("log(NBD_k)") +
   ylab("log(DC)") +
   stat_smooth(method = 'lm', formula = 'y~x') +
   geom_point(size=4, pch=21,color = "black", stroke=1.5, aes(fill=type)) +
   theme_test() 
   # geom_errorbar(aes(ymin=log(DC)-log(SE), ymax=log(DC)+log(SE)),colour= "blue", width=.03,size=0.75)
p2
 



#plot the bubble graph
#each row represents a species
#each column represents a site
bubble<-function(ssmat,factor=1)
{
  spsn=dim(ssmat)[1]
  siten=dim(ssmat)[2]
  mat=ssmat/sum(ssmat)
  mat=mat/max(mat)
  ##########################
  xy=expand.grid(1:spsn,1:siten)
  # windows()
  par(cex.axis=1.5,cex.lab=1.5)
  plot(xy,xlim=c(1,(spsn+1)),ylim=c(1,(siten+1)),axes=FALSE,xlab="Species",ylab="Sites",type="n")
  axis(1,c(1:spsn)+.5,1:spsn,tck=FALSE,col=NA)
  axis(2,c(1:siten)+.5,1:siten,tck=FALSE,col=NA)
  ##########################
  for(i in 1:siten)
  {
    id=which(xy[,2]==i)
    z=mat[,i]
    points(xy[id,]+.5,cex=z*factor,col=adjustcolor("red",alpha.f=.5),pch=16)
  }#i
  ##########################
  for(i in 1:(siten+1))
  {
    lines(1:(spsn+1),rep(i,1,spsn+1),col="gray")
  }#
  #
  for(i in 1:(spsn+1))
  {
    lines(rep(i,1,siten+1),1:(siten+1),col="gray")
  }#
  ##########################
  #
}#end





#making Fig. 1

set.seed(3)
dat=rMDir1(8,rep(.1,1,10),sample(500:1000,8))
bubble(dat,factor=5)
mtext(expression(paste(alpha,"=0.1",sep="")),cex=1.5)

#
set.seed(4)
dat=rMDir1(8,rep(1,1,10),sample(500:1000,8))
bubble(dat,factor=5)
mtext(expression(paste(alpha,"=1",sep="")),cex=1.5)

#
set.seed(1123)
dat=rMDir1(8,rep(10,1,10),sample(500:1000,8))
bubble(dat,factor=5)
mtext(expression(paste(alpha,"=10",sep="")),cex=1.5)
#


# ---------------------------------------------------------------------------------- ----  
# (3) Numerical test II                                                              ----

### example
par(mfrow=c(1,1))
test1 <- species.distribution(delta = 5, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
plot(test1, xlim = c(0, 200), ylim = c(0, 200), col=2, pch = 1, cex = 0.5, xlab = "Coordinate X",ylab = "Coordinate Y")
# abline(v= seq( 0, 200, by=20),h=seq(0, 200,by=20),col = grey(3/8),lty = "dotted")
CECI(test1, 200 * 200)


test2 <- species.distribution(delta = 50, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
CECI(test2, 200 * 200)

test3 <- species.distribution(delta = 500, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
CECI(test3, 200 * 200)

test4 <- species.distribution(delta = 5000, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
CECI(test4, 200 * 200)

points(test2, col =1, pch = 22,cex = 0.5)
points(test3, col =3, pch =4,cex = 0.5)
points(test4, col =5, pch =3,cex = 0.5)


par(mfrow=c(2,2))
library(spatstat) 

xy <- as.ppp(test1, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
DC(q1)
plot(q1, main = "")
plot(xy, col=2, pch = 1, cex = 0.5, alpha = 0.05, add = TRUE)

xy <- as.ppp(test2, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
DC(q1)
plot(q1, main = "")
plot(xy, col=1, pch = 1, cex = 0.5, add = TRUE)

xy <- as.ppp(test3, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
DC(q1)
plot(q1, main = "")
plot(xy, col=3, pch = 1, cex = 0.5, add = TRUE)

xy <- as.ppp(test4, W = c(0,200,0,200))
q1 <- quadratcount(xy, nx=10, ny=10)
DC(q1)
plot(q1, main = "")
plot(xy, col=5, pch = 1, cex = 0.5, add = TRUE)


# Dispersion Test for Spatial Point Pattern Based on Quadrat Counts
# Performs a test of Complete Spatial Randomness for a given point pattern, based on quadrat counts. 
# Alternatively performs a goodness-of-fit test of a fitted inhomogeneous Poisson model. By default performs chi-squared tests; can also perform Monte Carlo based tests.
quadrat.test(q1, alternative="two.sided",method="MonteCarlo",conditional=TRUE,nsim=999)
quadrat.test(q1, method="Chisq")


# formal test start

# Simulation I - sinle species
xy = rbind(c(4,4),c(5,5),c(8,8), c(10,10),c(20,20),c(25,25),c(40,40),c(50,50),c(100,100))
par(mfrow=c(3,3))
index_res_single <- data.frame()
library(spatstat) 
for(j in 1:dim(xy)[1])
{
  sdm_bic = nbd_bic = sdm_nrmsd = nbd_nrmsd = vector()
  for(i in 1:4000)
  {
    if(i  <= 1000){
      delta = 5
    }else if(i <= 2000){
      delta = 50
    }else if(i <= 3000){
      delta = 500
    }else if(i <= 4000){
      delta = 5000
    }
    
    # because sometimes dat could not be simulated correctly.
    for (ii in 1:3) {
      cat("It is the", ii, "runs", "\n") 
      dat=species.distribution(delta = delta, lambda = 3, prob = 0.01, cnum = 200, rnum = 200)
      AI = CECI(dat, 200 * 200)
      xy_ppp <- as.ppp(dat, W = c(0,200,0,200))
      q1 <- quadratcount(xy_ppp, nx=xy[j,1], ny=xy[j,2])
      # print(q1)
      if (length(which(q1>0)) == 1) {
        next 
      } else {
        break 
      } 
    }
    
    q2 <- matrix(q1, nrow = xy[j,1], ncol = xy[j,2])
    DC_mean = DC(q2)
    

    ## NBD start
    # convert into 1 row and n quadrats
    dat=as.vector(q2)
    res=fitNBD(dat)
    
    # model performance 
    nbd_num  = rnbinom(length(dat),size=res$par[1],mu=res$par[2])
    NBD_res <- r2.test(y_actual = sort(dat),y_predicted = sort(nbd_num))
    NBD_res
    
  
    ## SDM start
    res1 = fit(dat) #fit the standard species-site matrix
    p = rdirichlet1(1,alpha = rep(res1$par,length(dat)))
    sdm_num = rmultinom(1,size=sum(dat),prob=p)

    SDM_res <- r2.test(y_actual = sort(dat),y_predicted = sort(sdm_num))
    SDM_res

    # 
    # KS test
    library("dgof")
    nbd_ks <- dgof::ks.test(nbd_num, dat, simulate.p.value=TRUE, B=10000)
    sdm_ks <- dgof::ks.test(sdm_num, dat, simulate.p.value=TRUE, B=10000)

    
    # # index for plot
    sdm_bic = c(sdm_bic, as.numeric(format(round(res1$objective*2+log(length(dat))*1, digits = 4), nsmall = 4)))
    nbd_bic = c(nbd_bic, as.numeric(format(round(res$objective*2+log(length(dat))*2, digits = 4), nsmall = 4)))
    
    # index_res
    index_res_single <- rbind(index_res_single, data.frame("scale" = paste0(xy[j,1],"X",xy[j,2]),
                                                 "delta" = delta,
                                                 "DC"   = as.numeric(format(round(DC_mean$dc, digits = 4), nsmall = 4)),
                                                 "AI"   = as.numeric(format(round(as.numeric(AI), digits = 4), nsmall = 4)),
                                                 "SDM_alpha" = as.numeric(format(round(res1$par, digits = 4), nsmall = 4)), 
                                                 "SDM_Rsqure" = as.numeric(format(round(SDM_res$rsqure, digits = 4), nsmall = 4)),
                                                 "SDM_RMSE" = as.numeric(format(round(SDM_res$RMSE, digits = 4), nsmall = 4)),
                                                 "SDM_NRMSD" = as.numeric(format(round(SDM_res$NRMSD, digits = 4), nsmall = 4)),
                                                 "SDM_BIC" = as.numeric(format(round(res1$objective*2+log(length(dat))*1, digits = 4), nsmall = 4)),
                                                 "SDM_ks_D" = as.numeric(format(round(sdm_ks[["statistic"]][["D"]], digits = 4), nsmall = 4)),
                                                 "SDM_ks_p" = as.numeric(format(round(sdm_ks[["p.value"]], digits = 4), nsmall = 4)),
                                                 "NBD_k" = as.numeric(format(round(res$par[1], digits = 4), nsmall = 4)), 
                                                 "NBD_mu" = as.numeric(format(round(res$par[2], digits = 4), nsmall = 4)), 
                                                 "NBD_Rsqure" = as.numeric(format(round(NBD_res$rsqure, digits = 4), nsmall = 4)),
                                                 "NBD_RMSE" = as.numeric(format(round(NBD_res$RMSE, digits = 4), nsmall = 4)),
                                                 "NBD_NRMSD" = as.numeric(format(round(NBD_res$NRMSD, digits = 4), nsmall = 4)),
                                                 "NBD_BIC" = as.numeric(format(round(res$objective*2+log(length(dat))*2, digits = 4), nsmall = 4)),
                                                 "NBD_ks_D" = as.numeric(format(round(nbd_ks[["statistic"]][["D"]], digits = 4), nsmall = 4)),
                                                 "NBD_ks_p" = as.numeric(format(round(nbd_ks[["p.value"]], digits = 4), nsmall = 4))
                                                 ))
    
  }#i
  
  boxplot(cbind("SDM" = array(sdm_bic),"NBD" = array(nbd_bic)),
          main = paste("Spatial scale = ",paste0(xy[j,1]," X ",xy[j,2]),sep=""),
          xlab = "BIC",
          ylab = "Model",
          col = c("orange","grey"),
          border = "black",
          horizontal = TRUE,
          notch = TRUE,
          outline = F )

}#j

# the simulation will take ca. 30 min, then save results
saveRDS(index_res_single,"./output/randomized_single_species_distribution.rds")


# read in 
res_single <- readRDS("./output/randomized_single_species_distribution.rds") 


library(ggplot2)
library(colorspace)  
library(ggpmisc)


res_single$scale <- factor(res_single$scale, levels = c("4X4","5X5","8X8","10X10","20X20",
                                                        "25X25","40X40","50X50","100X100"))


# Simulation II - multiple species
# Spnum means how many species you want simulate.
# A means the areal size of study area
library(dgof)

Mul_SDM_comp <- function(x = 4, y = 4, delta = 5, lambda = 3, prob = 0.01, cnum = 200, rnum = 200, Spnum = 10)
{
  options(warn = -1)
  dat = q1 = q2 = list()
  for (u in 1:Spnum) {
    dat[[u]]=species.distribution(delta = delta, lambda = lambda, prob = prob, cnum = cnum, rnum = rnum)
    xy_ppp <- spatstat::as.ppp(dat[[u]], W = c(0,cnum,0,rnum))
    q1[[u]] <- spatstat::quadratcount(xy_ppp, nx=x, ny=y)
    q2[[u]] <- matrix(q1[[u]], nrow = 1, ncol = x*y )
  }
  
  # AI is the Clark and Evans competition index
  # AI
  dat = Reduce(rbind,dat)
  AI_mean = CECI(dat, cnum * rnum)
  
  # convert into n row and n quadrats
  dat = Reduce(rbind,q2)
  
  # DC is the Deviation Coefficient
  dat2 = apply(dat,2,sum) 
  DC_mean = DC(dat2)
  
  # NBD start
  res=fitNBD(dat)
  
  # model performance 
  nbd_num  = rnbinom(length(dat),size=res$par[1],mu=res$par[2])
  NBD_res <- r2.test(y_actual = sort(dat),y_predicted = sort(nbd_num))
  NBD_res
  
  ## SDM start
  res1 = fit(dat) #fit the standard species-site matrix
  
  # model performance 
  sdm_num = rMDir1(dim(dat)[[1]],rep(res1$par,1,dim(dat)[[2]]),apply(dat,1,sum))
  
  
  all.equal(apply(sdm_num, 1, sum), apply(dat, 1, sum))
  
  SDM_res <- r2.test(y_actual = sort(dat),y_predicted = sort(sdm_num))
  SDM_res
  
  if(SDM_res$NRMSD > NBD_res$NRMSD){print("NBD better")}else{print("SDM better")}
  
  # Kolmogorov?CSmirnov test
  nbd_ks <- dgof::ks.test(nbd_num, c(dat), simulate.p.value=TRUE, B=10000)
  sdm_ks <- dgof::ks.test(sdm_num, dat, simulate.p.value=TRUE, B=10000)
  
  # index_res
  index_res_multi  = data.frame("scale" = paste0(x,"X",y),
                                "delta" = delta,
                                "DC"   = as.numeric(format(round(DC_mean$dc, digits = 4), nsmall = 4)),
                                "AI"   = as.numeric(format(round(as.numeric(AI_mean), digits = 4), nsmall = 4)),
                                "SDM_alpha" = as.numeric(format(round(res1$par, digits = 4), nsmall = 4)), 
                                "SDM_Rsqure" = as.numeric(format(round(SDM_res$rsqure, digits = 4), nsmall = 4)),
                                "SDM_RMSE" = as.numeric(format(round(SDM_res$RMSE, digits = 4), nsmall = 4)),
                                "SDM_NRMSD" = as.numeric(format(round(SDM_res$NRMSD, digits = 4), nsmall = 4)),
                                "SDM_BIC" = as.numeric(format(round(res1$objective*2+log(length(dat))*1, digits = 4), nsmall = 4)),
                                "SDM_ks_D" = as.numeric(format(round(sdm_ks[["statistic"]][["D"]], digits = 4), nsmall = 4)),
                                "SDM_ks_p" = as.numeric(format(round(sdm_ks[["p.value"]], digits = 4), nsmall = 4)),
                                "NBD_k" = as.numeric(format(round(res$par[1], digits = 4), nsmall = 4)), 
                                "NBD_mu" = as.numeric(format(round(res$par[2], digits = 4), nsmall = 4)), 
                                "NBD_Rsqure" = as.numeric(format(round(NBD_res$rsqure, digits = 4), nsmall = 4)),
                                "NBD_RMSE" = as.numeric(format(round(NBD_res$RMSE, digits = 4), nsmall = 4)),
                                "NBD_NRMSD" = as.numeric(format(round(NBD_res$NRMSD, digits = 4), nsmall = 4)),
                                "NBD_BIC" = as.numeric(format(round(res$objective*2+log(length(dat))*2, digits = 4), nsmall = 4)),
                                "NBD_ks_D" = as.numeric(format(round(nbd_ks[["statistic"]][["D"]], digits = 4), nsmall = 4)),
                                "NBD_ks_p" = as.numeric(format(round(nbd_ks[["p.value"]], digits = 4), nsmall = 4))
  )
  
  return(index_res_multi)
  
}

# simple test
Mul_SDM_comp(x = 4, y = 4, delta = 5, lambda = 3, prob = 0.01, cnum = 200, rnum = 200, Spnum = 10)

# combinations
xy = rbind(c(4,4),c(5,5),c(8,8), c(10,10),c(20,20),c(25,25),c(40,40),c(50,50),c(100,100))
delta = c(5,50,500,5000)
input_dat <- data.frame()
for(j in 1:dim(xy)[1]){
  for (i in 1:length(delta)) {
    for (rr in 1:1000) {
      input_dat = rbind(input_dat,data.frame("x" = xy[j,1], "y" = xy[j,1], "delta" = delta[i]))
    }
  }
}

# loop
# Parallel loop start...
# library
t1 <- Sys.time()
library(doParallel)
detectCores()
cl <- makeCluster(96) # server 46 cores 
registerDoParallel(cl)
Res <- NULL
Res <- foreach(x = 1:dim(input_dat)[1], .packages=c("dgof","MASS","MCMCpack","dirmult"),.combine = rbind) %dopar% {
  if(!dir.exists("output")){dir.create("output")}
  r <- Mul_SDM_comp(x = input_dat[x,"x"], 
                    y = input_dat[x,"y"], 
                    delta = input_dat[x,"delta"], 
                    lambda = 3, 
                    prob = 0.01, 
                    cnum = 200, 
                    rnum = 200, 
                    Spnum = 10)
  saveRDS(r, paste0("./output/",x,".rds"))
  r
}

t2 <- Sys.time()-t1
t2

stopCluster(cl)

# the simulation will take ca. 5 hours, then save results
saveRDS(index_res_single,"./output/randomized_10_species_distribution.rds")


# read in 
res_multi <- readRDS("./output/randomized_10_species_distribution.rds")

delta_bic_single <- (res_single$SDM_BIC - res_single$NBD_BIC)/res_single$SDM_BIC * 100
delta_bic_multiple <- (res_multi$SDM_BIC - res_multi$NBD_BIC)/res_multi$SDM_BIC * 100

library(vioplot)

par(mfrow = c(2,2))

# data
dat1 <- delta_bic_single
# [-which(delta_bic_single < -50)] 

# Histogram
hist(dat1, probability =T, xlab ="¦¤BIC%", ylab = "Density", col = "grey",
     axes = FALSE, main = "Randomized spatial-point distribution data - Single species")

# Axis
axis(1)
# axis(2)

# Density
# lines(density(dat1), col = "red", lwd = 1)

# Add boxplot
par(new = TRUE)
# boxplot(dat1, horizontal = TRUE, axes = FALSE,
#        lwd = 1, col = rgb(0, 1, 1, alpha = 0.15))

vioplot(dat1, horizontal = TRUE, yaxt = "n", axes = FALSE,
        col = rgb(0, 1, 1, alpha = 0.15))

#
dat2 <- delta_bic_multiple

# Histogram
hist(dat2, probability = TRUE, xlab ="¦¤BIC%", ylab = "Density", col = "grey",
     axes = FALSE, main = "Randomized spatial-point distribution data - Multiple species")

# Axis
axis(1)
# axis(2)

# Add boxplot
par(new = TRUE)
vioplot(dat2, horizontal = TRUE, yaxt = "n", axes = FALSE,
        col = rgb(1, 1, 0, alpha = 0.15))

delta_D_single <- res_single$SDM_ks_D - res_single$NBD_ks_D
delta_D_multiple <- res_multi$SDM_ks_D - res_multi$NBD_ks_D

# Histogram
hist(delta_D_single, probability =T, xlab ="¦¤D-statistic", ylab = "Density", col = "grey",
     axes = FALSE, main = "Randomized spatial-point distribution data - Single species")

# Axis
axis(1)
# axis(2)

# Add boxplot
par(new = TRUE)

vioplot(delta_D_single, horizontal = TRUE, yaxt = "n", axes = FALSE,
        col = rgb(0, 1, 1, alpha = 0.15))


# Histogram
hist(delta_D_multiple, probability = TRUE, xlab ="¦¤D-statistic", ylab = "Density", col = "grey",
     axes = FALSE, main = "Randomized spatial-point distribution data - Multiple species")

# Axis
axis(1)
# axis(2)


# Add boxplot
par(new = TRUE)
vioplot(delta_D_multiple, horizontal = TRUE, yaxt = "n", axes = FALSE,
        col = rgb(1, 1, 0, alpha = 0.15))

# ---------------------------------------------------------------------------------- ---- 
# (4) Empirical test II                                                              ----

# BCI data - spatial aggregation index dynamic with time squences 
# Spatial sacle effects

for (i in 1: 8) {
  load(paste0("./input/bci.tree",i,".rdata")) 
}

BCI_List <- list()
BCI_List[[1]] <- bci.tree1
BCI_List[[2]] <- bci.tree2
BCI_List[[3]] <- bci.tree3
BCI_List[[4]] <- bci.tree4
BCI_List[[5]] <- bci.tree5
BCI_List[[6]] <- bci.tree6
BCI_List[[7]] <- bci.tree7
BCI_List[[8]] <- bci.tree8

rm(bci.tree1, bci.tree2, bci.tree3, bci.tree4, bci.tree5, bci.tree6, bci.tree7, bci.tree8)
BCI_data = BCI_List
load("./input/bci.spptable.rdata")
BCI_spp = bci.spptable


# stand density index calculator
BCI_SDI <- function(scale = c(20,20), census = 1, recruitment = T, BCI_data)
{
  if(recruitment == F){
    ids1 = which(BCI_data[[1]]$status=="A" & BCI_data[[1]]$dbh >= 10)
    ids2 = which(is.na(BCI_data[[1]]$gx) | is.na(BCI_data[[1]]$gy))
    bci = BCI_data[[census]]
    bci = bci[ids1,]
    bci = bci[-ids2,]
    bci = bci[which(bci$status=="A" & bci$dbh >= 10),]
  
  } else {
    #BCI dataset:
    bci = BCI_data[[census]]
    ids = which(bci$status=="A" & bci$dbh >= 10)
  if(length(ids)>0)
  {
    bci=bci[ids,]
  }
  ids = which(is.na(bci$gx) | is.na(bci$gy))
  if(length(ids)>0)
  {
    bci=bci[-ids,]
  }
 }
  #bci has x range 0-1000, y range 0-500:
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 500
  #
  
  x.size= scale[1] #different scale size
  y.size= scale[2] #different scale size
  
  xint=seq(xmin,xmax,by=x.size)
  yint=seq(ymin,ymax,by=y.size)
  xint=xint[-length(xint)]
  yint=yint[-length(yint)]
  #
  xyint=expand.grid(xint,yint)
  xyint=cbind(xyint,xyint[,1]+x.size,xyint[,2]+y.size)
  colnames(xyint)=c("x1","y1","x2","y2")
  
  sps = as.character(unique(bci$sp))
  tx=bci$gx
  ty=bci$gy
  #
  ids_records <- list()
  for(j in 1:dim(xyint)[1])
  {
    ids=which(tx>xyint[j,1] & tx<=xyint[j,3] & ty>xyint[j,2] & ty<=xyint[j,4])
    if(length(ids)>0)
    {
      ids_records[[j]] = ids
    }#if
  }#j
  
  SDI <- SR <- abundance <-  NULL
  
  dc = seq(100, 3500, 100)
  for(xx in 1:dim(xyint)[1])
  {
    # SDI calculation
    # Diameter class (100 mm)
    bci_sel = bci[ids_records[[xx]],]
    Di = Ni = NULL
    for (qq in 1: length(dc)) {
      D_all = bci_sel$dbh[bci_sel$dbh > (dc[qq]-100) & bci_sel$dbh <= dc[qq]]
      D_mean = sqrt(sum((D_all/10)^2)/length(D_all))
      Di = c(Di,D_mean)
      Ni = c(Ni,length(D_all)*(10000/(scale[1]*scale[2])))
    }
    SR = c(SR, length(unique(bci_sel$sp)))
    abundance = c(abundance, length(bci_sel$sp))
    SDI = c(SDI, sum(Ni*(Di/25)^1.6, na.rm = T))
    
   }
  return(list("SDI" = SDI, "SR" = SR, "abundance" = abundance))
}

# test
BCI_SDI(scale = c(20,20), census = 1, recruitment = T, BCI_data = BCI_data)


# family
BCI_SDM <- function(scale = c(20,20), recruitment = T, family = "Asteraceae", census = 1, BCI_data, BCI_spp)
{
  if(recruitment == F){
    ids1 = which(BCI_data[[1]]$status=="A" & BCI_data[[1]]$dbh >= 10)
    ids2 = which(is.na(BCI_data[[1]]$gx) | is.na(BCI_data[[1]]$gy))
    bci = BCI_data[[census]]
    bci = bci[ids1,]
    bci = bci[-ids2,]
    bci = bci[which(bci$status=="A" & bci$dbh >= 10),]
    
  } else {
    #BCI dataset:
    bci = BCI_data[[census]]
    ids = which(bci$status=="A" & bci$dbh >= 10)
    if(length(ids)>0)
    {
      bci=bci[ids,]
    }
    ids = which(is.na(bci$gx) | is.na(bci$gy))
    if(length(ids)>0)
    {
      bci=bci[-ids,]
    }
  }
  
  #bci has x range 0-1000, y range 0-500:
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 500
  #
  out=vector()
  
  x.size= scale[1] #different scale size
  y.size= scale[2] #different scale size
  
  xint=seq(xmin,xmax,by=x.size)
  yint=seq(ymin,ymax,by=y.size)
  xint=xint[-length(xint)]
  yint=yint[-length(yint)]
  #
  xyint=expand.grid(xint,yint)
  xyint=cbind(xyint,xyint[,1]+x.size,xyint[,2]+y.size)
  colnames(xyint)=c("x1","y1","x2","y2")
  
  # bci species
  bci_spp <- as.character(unique(bci$sp))
  
  #create species-site matrix
  if(family == "ALL"){
    BCI_spp_sel = BCI_spp
    sps = as.character(unique(BCI_spp_sel$sp))
    sps = sps[sps %in% bci_spp]
  }else{
    BCI_spp_sel = subset(BCI_spp, BCI_spp$Family == family)
    sps = as.character(unique(BCI_spp_sel$sp))
    sps = sps[sps %in% bci_spp]
  }
  
 
  # AI
  AI = CECI(bci[bci$sp %in% sps,c("gx","gy")], 1000 * 500)
  
  
  mat=matrix(0,nrow=length(sps),ncol=dim(xyint)[1])
  #
  for(i in 1:length(sps))
  {
    id=which(bci$sp==sps[i])
    tx=bci$gx[id]
    ty=bci$gy[id]
    #
    for(j in 1:dim(xyint)[1])
    {
      ids=which(tx>xyint[j,1] & tx<=xyint[j,3] & ty>xyint[j,2] & ty<=xyint[j,4])
      if(length(ids)>0)
      {
        mat[i,j]=length(ids)
      }#if
    }#j
  }#i
  
  row.names(mat) = sps
  
  quardat_SR <- quardat_abundance <-  NULL
  for(xx in 1:dim(mat)[2])
  {
    mat_sel <- mat[which(mat[, xx] > 0), xx]
    quardat_SR = c(quardat_SR,  length(mat_sel))
    quardat_abundance = c(quardat_abundance,  sum(mat_sel))
  }
  
  res=fit(mat)
  # DC
  dat2 = apply(mat,2,sum) 
  dc = DC(dat2)
  dc_mean = dc$dc
  
  return(list("family" = family, "total species" = sps, "total SR" = length(sps), 
              "quardat SR" = quardat_SR,
              "quardat Abundance" = quardat_abundance,
              "SDM_alpha" = round(res$par,3),
              "DC" = dc_mean, 
              "AI" =  AI,
              "Abundance_table" = mat))
}

BCI_SDM_simple <- function(scale = c(20,20), recruitment = T, family = "Asteraceae", census = 1, BCI_data, BCI_spp)
{
  if(recruitment == F){
    ids1 = which(BCI_data[[1]]$status=="A" & BCI_data[[1]]$dbh >= 10)
    ids2 = which(is.na(BCI_data[[1]]$gx) | is.na(BCI_data[[1]]$gy))
    bci = BCI_data[[census]]
    bci = bci[ids1,]
    bci = bci[-ids2,]
    bci = bci[which(bci$status=="A" & bci$dbh >= 10),]
    
  } else {
    #BCI dataset:
    bci = BCI_data[[census]]
    ids = which(bci$status=="A" & bci$dbh >= 10)
    if(length(ids)>0)
    {
      bci=bci[ids,]
    }
    ids = which(is.na(bci$gx) | is.na(bci$gy))
    if(length(ids)>0)
    {
      bci=bci[-ids,]
    }
  }
  
  #bci has x range 0-1000, y range 0-500:
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 500
  #
  out=vector()
  
  x.size= scale[1] #different scale size
  y.size= scale[2] #different scale size
  
  xint=seq(xmin,xmax,by=x.size)
  yint=seq(ymin,ymax,by=y.size)
  xint=xint[-length(xint)]
  yint=yint[-length(yint)]
  #
  xyint=expand.grid(xint,yint)
  xyint=cbind(xyint,xyint[,1]+x.size,xyint[,2]+y.size)
  colnames(xyint)=c("x1","y1","x2","y2")
  
  # bci species
  bci_spp <- as.character(unique(bci$sp))
  
  #create species-site matrix
  if(family == "ALL"){
    BCI_spp_sel = BCI_spp
    sps = as.character(unique(BCI_spp_sel$sp))
    sps = sps[sps %in% bci_spp]
  }else{
    BCI_spp_sel = subset(BCI_spp, BCI_spp$Family == family)
    sps = as.character(unique(BCI_spp_sel$sp))
    sps = sps[sps %in% bci_spp]
  }
  
  mat=matrix(0,nrow=length(sps),ncol=dim(xyint)[1])
  #
  for(i in 1:length(sps))
  {
    id=which(bci$sp==sps[i])
    tx=bci$gx[id]
    ty=bci$gy[id]
    #
    for(j in 1:dim(xyint)[1])
    {
      ids=which(tx>xyint[j,1] & tx<=xyint[j,3] & ty>xyint[j,2] & ty<=xyint[j,4])
      if(length(ids)>0)
      {
        mat[i,j]=length(ids)
      }#if
    }#j
  }#i
  
  row.names(mat) = sps
  
  quardat_SR <- quardat_abundance <-  NULL
  for(xx in 1:dim(mat)[2])
  {
    mat_sel <- mat[which(mat[, xx] > 0), xx]
    quardat_SR = c(quardat_SR,  length(mat_sel))
    quardat_abundance = c(quardat_abundance,  sum(mat_sel))
  }
  
  res=fit(mat)
  # DC
  dat2 = apply(mat,2,sum) 
  dc = DC(dat2)
  dc_mean = dc$dc
  
  return(list("family" = family, "total species" = sps, "total SR" = length(sps), 
              "quardat SR" = quardat_SR,
              "quardat Abundance" = quardat_abundance,
              "SDM_alpha" = round(res$par,3),
              "DC" = dc_mean, 
              "Abundance_table" = mat))
}

# test
BCI_SDM(scale = c(20,20), recruitment = T, family = "Asteraceae", census = 1, BCI_data,  BCI_spp)
BCI_SDM_simple(scale = c(20,20), recruitment = T, family = "Asteraceae", census = 1, BCI_data,  BCI_spp)

# please note that calculating AI index for such a large number of trees is time-consuming, 
# therefore, We have embedded parallel operations in the CECI function
# these simulation will take ca. 12 hours
ALP <- data.frame()
Abundance_matrix = NULL
for (cen in 1:8) {
    alp <- BCI_SDM(scale = c(20,20), recruitment = T, family = "ALL", census = cen, BCI_data, BCI_spp)
    ALP <- rbind(ALP, data.frame("census" = cen, "SDM_alpha" = alp$SDM_alpha, "SR" = alp$`total SR`, "DC" =  alp$DC, "AI" = alp$AI))
    Abundance_matrix <- cbind(Abundance_matrix, alp$`quardat Abundance`)
    cat("Census",cen, "finished...","\n")
}
 
colnames(Abundance_matrix) <- 1:8
saveRDS(Abundance_matrix,"./output/Abundance_matrix.rds")
write.csv(ALP,"./output/ALP.csv")
write.csv(Abundance_matrix,"./output/Abundance_matrix.csv")


for (i in 1:8) {
  nam <- paste0("SDI",i)
  sdi = BCI_SDI(scale = c(20,20), census = i, BCI_data, recruitment = T)
  assign(nam, sdi$SDI) 
}
SDI = cbind(SDI1, SDI2, SDI3,SDI4,SDI5,SDI6,SDI7,SDI8)
rm(SDI1);rm(SDI2);rm(SDI3);rm(SDI4);rm(SDI5);rm(SDI6);rm(SDI7);rm(SDI8)
final_SDI = cbind(xyint,SDI)
final_SDI$ID = 1:1250
write.csv(final_SDI,"./output/BCI_SDI.csv")



BCI_species <- merge(data.frame("species" = BCI_data[[1]]$sp),BCI_spp[,c(1,5)], by.x = "species",by.y = "sp")
BCI_species <- unique(BCI_species$Family)

ALP <- data.frame()
for(fam in 60:length(BCI_species)){
  cat("Family",fam, "started...","\n")
  for (cen in 1:8) {
    alp <- BCI_SDM_simple(scale = c(20,20), recruitment = T, family = BCI_species[fam], census = cen, BCI_data, BCI_spp)
    saveRDS(alp, paste0("./outcome/family_res/",BCI_species[fam],"_",cen,".rds"))
    ALP <- rbind(ALP, data.frame("census" = cen,"family" = BCI_species[fam], "SDM_alpha" = alp$SDM_alpha, "SR" = alp[["total SR"]],"DC" = alp[["DC"]]))
  }
  cat("Family",BCI_species[fam], "finished...","\n")
}

write.csv(ALP,"./output/all_family_alpha.csv")


AM <- readRDS("./output/Abundance_matrix.rds")
AM <- data.frame("Abundance" = c(AM[,1],AM[,2],AM[,3],AM[,4],AM[,5],AM[,6],AM[,7],AM[,8]),
                 "Census" = rep(c(1982,1985,1990,1995,2000,2005,2010,2015),each=1250))
AM$Census <- as.factor(AM$Census)


library(ggridges)
library(RColorBrewer)
ggplot(AM, aes(x = `Abundance`, y = `Census`, fill = ..density..)) + 
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.00, size = 0.3) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8,'Spectral')))(32))+
  theme_bw()


# ---------------------------------------------------------------------------------- ----
# EOF
