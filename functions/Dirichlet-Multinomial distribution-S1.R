
# pacakga necessary
library(MCMCpack)
library(dirmult)
#

##############
#random variate generation
rdirichlet<-function (n = 1, alpha) 
{
    Gam <- matrix(0, n, length(alpha))
    for (i in 1:length(alpha)) Gam[, i] <- rgamma(n, shape = alpha[i])
    Gam/rowSums(Gam)
}


#################
#avoid NAs:
rdirichlet1<-function (n = 1, alpha) 
{
rv=rdirichlet(n,alpha)
if(n==1)
{
rv=t(as.matrix(rv))
}
#
ids=which(is.na(rv[,1]))
if(length(ids)>0)
{
for(i in 1:length(ids))
{
rv[ids[i],]=rep(0,1,length(alpha))
rv[ids[i],sample(1:length(alpha),1)]=1
}#i
}#
#
return(rv)
}



##################################################
#probability density (not mass!) function of Dirichlet distribution
dDir<-function(x,alpha)
{
x=as.vector(x)
if(length(x)==length(alpha))
{
num=length(x)
#
part1=lgamma(sum(alpha))
part2=sum(lgamma(alpha))
part3=sum((alpha-1)*log(x))
#
f=part1-part2+part3
#
return(exp(f))
}else
{
return(NA)
}
}#
##################################################


##############################################################################
#probability mass (or density) function of the multinomial-Dirichlet distribution
dMDir<-function(x,alpha,log=FALSE)
{
x=as.vector(x)
if(length(x)==length(alpha))
{
num=length(x)
n=sum(x)
#
p1=lgamma(n+1)
p2=lgamma(sum(alpha))
p3=lgamma(n+sum(alpha))
#
p4=sum(lgamma(x+alpha))
p5=sum(lgamma(x+1))
p6=sum(lgamma(alpha))
#
f=p1+p2-p3+p4-p5-p6
#
if(log==FALSE)
{
return(exp(f))
}else
{
return(f)
}
#
}else
{
return(NA)
}
}#
##################################################



##################################################################
#probability mass function for independent negative binomial model
#two parameters mu and k (or size)
#dNBD is equal to dnbinom(10,size=.1,mu=10)
dNBD<-function(x,mu,k,log=FALSE)
{
p1=lgamma(k+x)
p2=lgamma(x+1)
p3=lgamma(k)
p4=k*log(k/(k+mu))
p5=x*log(mu/(k+mu))
v=p1-p2-p3+p4+p5
######
if(log==TRUE)
{
return(v)
}else
{
return(exp(v))
}
}#
#


##################################################################
#probability mass function for negative multinomial model
#two parameters mu and k (or size)
dNMD<-function(x,mu,k,log=FALSE)
{
n=sum(x)
q=length(x)
#
p1=lgamma(k+n)
p2=sum(lgamma(x+1))
p3=lgamma(k)
p4=k*log(k/(k+mu))
p5=sum(x*log(mu/q/(k+mu)))
v=p1-p2-p3+p4+p5
######
if(log==TRUE)
{
return(v)
}else
{
return(exp(v))
}
}#
#





##################################################
#simulation of multinomial-Dirichlet distribution
#N is the total number of organisms, can be a single value or a vector
rMDir<-function(n,alpha,N)
{
p=MCMCpack::rdirichlet(n,alpha)
#
if(length(N)==1)
{
N=rep(N,1,n)
}#
#
mat=vector()
for(i in 1:n)
{
v=rmultinom(1,size=N[i],prob=p[i,])
mat=rbind(mat,as.vector(v))
}#
#
return(mat)
}#
#
#rMDir(10,rep(.1,1,8),rpois(10,200))
#rMDir(10,rep(1,1,8),rpois(10,200))



##################################################
#simulation of multinomial-Dirichlet distribution by avoiding NAs
#N is the total number of organisms, can be a single value or a vector
rMDir1<-function(n,alpha,N)
{
p=rdirichlet1(n,alpha)
#
if(length(N)==1)
{
N=rep(N,1,n)
}#
#
mat=vector()
for(i in 1:n)
{
v=rmultinom(1,size=N[i],prob=p[i,])
mat=rbind(mat,as.vector(v))
}#
#
return(mat)
}#
#

##################################################
#random matrices generation
#rowSum is kept without chnage!
#method="proportion" means the relative occurrence rate is kept without change
#method="constant" means the relative occurrence rate is equal
#method="random" means the relative occurrence rate is randomly drawn from a uniform distribution
rand<-function(mat,reptime=1000,method="proportion")
{
if(is.vector(mat))
{
mat=t(as.matrix(mat))
}#
##################
out=array(0,dim=c(dim(mat)[1],dim(mat)[2],reptime))
one=function(v)
{
n=sum(v)
#return(rmultinom(1,size=n,prob=rep(1,1,length(v)))) #even distribution
if(method=="random")
{
p=runif(length(v),0,1)
}
if(method=="constant")
{
p=rep(1,1,length(v))
}
if(method=="proportion")
{
p=v
}
#
p=p/sum(p)
return(rmultinom(1,size=n,prob=p))
}#
##################
for(i in 1:reptime)
{
out[,,i]=t(apply(mat,1,one))
}#i
#
return(out)
}#
#########



######################################################
#calculation of the negative log-likelihood function
#for the Dirichlet-Multinomial model
#each row represents a species
likelihood<-function(mat,alpha)
{
one=function(v)
{
return(dMDir(v,alpha,log=TRUE))
}#
#########
all=apply(mat,1,one)
#########
return(-sum(all))
}#
#


######################################################
#calculation of the negative log-likelihood function
#for the null model: multinomial model
#each row represents a species
likelihood0<-function(mat,method="proportion")
{
one=function(v)
{
if(method=="random")
{
p=runif(length(v),0,1)
}
if(method=="constant")
{
p=rep(1,1,length(v))
}
if(method=="proportion")
{
p=v
}
#
p=p/sum(p)
########################
return(dmultinom(v,size=sum(v),prob=p,log=TRUE))
}#
#########
all=apply(mat,1,one)
#########
return(-sum(all))
}#
#


######################################################
#fitting of the SDM model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fit<-function(mat)
{
if(is.vector(mat))
{
mat=t(as.matrix(mat))
}#
##################
n=dim(mat)[2] #number of sites
likelihood<-function(pars)
{
alpha=pars[1]
alpha=rep(alpha,1,n)
#
one=function(v)
{
return(dMDir(v,alpha,log=TRUE))
}#
#
all=apply(mat,1,one)
return(-sum(all))
}#
##################
res=nlminb(runif(1),likelihood,lower=1e-20,upper=1e+20)
#
return(res)
}#
#







######################################################
#fitting of the independent NBD model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fitNBD<-function(mat)
{
x=as.vector(mat) #use vector data for fitting
likelihood<-function(pars)
{
k=pars[1]
mu=pars[2]
#
all=dnbinom(x,size=k,mu=mu,log=TRUE)
return(-sum(all,na.rm=TRUE))
}#
##################
res=nlminb(runif(2),likelihood,lower=rep(1e-20,1,2),upper=rep(1e+20,1,2))
#
return(res)
}#
#




######################################################
#fitting of the NMD model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fitNMD<-function(mat)
{
if(is.vector(mat))
{
mat=t(as.matrix(mat))
}#
##################
n=dim(mat)[2] #number of sites
likelihood<-function(pars)
{
k=pars[1]
mu=pars[2]
#
one=function(v)
{
return(dNMD(v,k=k,mu=mu,log=TRUE))
}#
#
all=apply(mat,1,one)
return(-sum(all))
}#
##################
res=nlminb(runif(2),likelihood,lower=rep(0.000001,1,2),upper=rep(100000,1,2))
#
return(res)
}#
#









######################################################
#fitting of the ordinary multinomial model
#mat is a species-site matrix
#if it is a vector, convert it into a matrix
fitMD<-function(mat)
{
if(is.vector(mat))
{
mat=t(as.matrix(mat))
}#
##################
n=dim(mat)[2] #number of sites
likelihood<-function()
{
#
one=function(v)
{
return(dmultinom(v,size=sum(v),prob=rep(1/length(v),1,length(v)),log=TRUE))
}#
#
all=apply(mat,1,one)
return(-sum(all))
}#
##################
res=likelihood()
#
return(res)
}#
#







