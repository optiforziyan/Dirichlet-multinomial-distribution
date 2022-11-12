################################################
library(MASS)
#simulation of species distribution using Poisson cluster process
#generate a matrix with cnum and rnum, species presence in each cell was recorded as 1
#using bivariate probability function to generate offsprings
#Poisson distribution to generate number of parental points
#geometric distribution to generate number of offsprings for each parental point
#uniform distribution to obtain locations of parental points
#cnum determines number of columns, while rnum determines number of rows
#lambda is the parameter for Poisson distribution, prob is the parameter for geometric distribution
species.distribution<-function(delta,lambda,prob,cnum,rnum)
{
  pnum<-rpois(1,lambda)+1 # avoiding zeros, Poission distribution,lambda =  parameter intensity ¦Ñ
  offnum<-rgeom(pnum,prob)+1 # avoiding zeros,The Geometric Distribution, prob = 0.01
  # generate the locations of parental points
  x<-sample(cnum,pnum,replace=FALSE)
  y<-sample(rnum,pnum,replace=FALSE)
  #
  #now generating offspring locations
  sig<-matrix(c(delta,0,0,delta),2,2)
  pp<-vector()
  for(i in 1:pnum)
  {
    this<-mvrnorm(offnum[i],c(x[i],y[i]),sig)#bivariate normal distribution,delta link to ¦Ò2 
    this<-rbind(this,c(x[i],y[i])) #including parental point is tricky, otherwise sometimes this is NULL
    #exluding out-of-bound points
    idd<-which(this[,1]>cnum | this[,2]>rnum | this[,1]<1 | this[,2]<1)
    if(length(idd)!=0)
    {
      this<-this[-idd,]
    }
    pp<-rbind(pp,this)
  }
  #
  return(pp)
}#end




### The Cluster-Random-Regular Continuum
#  Clark, P. J., & Evans, F. C. (1954). Distance to nearest neighbor as a measure of spatial relationships in populations. Ecology, 35, 445-453.
#  R is the Clark and Evans competition index
#  the input mat contains x and y coordinates 

CECI <- function(mat, area)
{
  # if(length(x) != length(y)){
  #   stop("x and y should have the same number")	
  # }
  n = dim(mat)[1] # numbers of points
  area = area   # study area (square meter)
  Rho = n/area  # density of points per square meter
  # minimal distance between a target point to its neighbours
  if(n < 1000){
  dd <- as.matrix(dist(mat, p=2))
  diag(dd) <- NA
  dmin = apply(dd, 1, min, na.rm = TRUE)
  } else{
  
  t1 <- Sys.time()
  library(doParallel)
  detectCores()
  cl <- makeCluster(10)
  registerDoParallel(cl)
  Res <- foreach(x = 1:n,.combine = c) %dopar% {
    d <- sqrt((mat[x,1] - mat[,1])^2 + (mat[x,2] - mat[,2])^2)
    dmin = sort(d)[2]
    dmin
  }
  t2 <- Sys.time()-t1
  t2
  stopCluster(cl)
  dmin = Res
  }
  # R is the Clark and Evans competition index
  R = (1/n) * sum(dmin) * 2 * sqrt(Rho)
  names(R) = "Clark & Evans Competition Index"
  if(R < 1){
    # The average nearest neighbor distance decreases as points become more tightly clustered; hence, the closer to 0 the value for R becomes.
    print("Clustered pattern detected") 
  } else if(R == 1){
    # The closer the points are to being randomly dispersed, the more similar are the values of and the closer to 1 the value for R becomes.
    print("Randomly dispersed pattern detected") 
  } else if(R > 1){
    # The average nearest neighbor distance decreases as points become more uniformity; hence, the R value greater than 1.
    print("Regularly dispersed pattern detected") 
  }
  return(R)
}



####
# DC means deviation coefficient or diffusion coefficient 
# G. E. BLACKMAN, Statistical and Ecological Studies in the Distribution of Species in Plant Communities: I. Dispersion as a Factor in the Study of Changes in Plant Populations, Annals of Botany, Volume 6, Issue 2, April 1942, Pages 351¨C370, https://doi.org/10.1093/oxfordjournals.aob.a088411
DC <- function(dat)
{ 
  N  = length(dat)
  dat <- array(dat)
  dat <- data.frame(table(dat))
  dat$dat <- as.numeric(as.character(dat$dat))
  #variance
  V = (sum(dat$dat^2 * dat$Freq) - (((sum(dat$dat * dat$Freq))^2)/N))/(N-1)
  #mean
  X = sum(dat$dat * dat$Freq)/N
  # deviation coefficient
  dc = V/X
  # names(dc) = "Diffusion coefficient"
  # S is the standard error
  # t value
  t = (dc - 1) / sqrt(2/(N - 1))
  # p value 
  p_value = round(2 * pt(-abs(t),df = N-1, lower.tail = T),3)
  # p_value = 2 * pt(abs(t),df = N-1, lower.tail = F)
  if(dc < 1){
    # The average nearest neighbor distance decreases as points become more tightly clustered; hence, the closer to 0 the value for R becomes.
    print("Regularly dispersed pattern detected") 
  } else if(dc == 1){
    # The closer the points are to being randomly dispersed, the more similar are the values of and the closer to 1 the value for R becomes.
    print("Randomly dispersed pattern detected") 
  } else if(dc > 1){
    # The average nearest neighbor distance decreases as points become more uniformity; hence, the R value greater than 1.
    print("Clustered pattern detected") 
  }
  return(list("dc" = dc,"t" = t, "p" = p_value))
}

### calculate r2, RMSE, NRMSD based on actual values and predicted values
r2.test<-function(y_actual,y_predicted)
{
  avr_y_actual <- mean(y_actual)
  ss_total <- sum((y_actual - avr_y_actual)^2)
  ss_residuals <- sum((y_actual - y_predicted)^2)
  rsquare <- 1 - ss_residuals / ss_total
  #return(rsquare) 
  n1<-length(y_actual)
  n2<-length(y_predicted)#
  meansquare<-ss_residuals/(n1-2)
  #return(meansquare)#MS
  RMSE<-sqrt(ss_residuals/n1) # 
  NRMSD<-RMSE/(max(y_actual)-min(y_actual))# NRMSD (normalized root mean square error deviation) 
  return(list("rsqure" = rsquare, "RMSE" = RMSE, "NRMSD"= NRMSD))
}
