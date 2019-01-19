#### Empirical Dynamic Quantile for Visualization of High-Dimensional Time Series
#### This program has three parts.
#### They are
#### (1) EDQts: Compute empirical dynamic quantile for a given probability "p" based on the weighted algoritm proposed in the 
####            article by Pena, Tsay and Zamar.
#### (2) EDQsim: Simple programs to generate high-dimensional Vector AR(1) or Vector MA(1) or Vector GARCH(1,1) time series.
####
#### (3) EDQplot: plot the observed time series with selected EDQs.
####
#### Created on January 15, 2019.
#### Special Notes:
#### (1) EDQts: This function uses the R code: WDQTS written originally by Daniel Pena, R. Tsay and R. Zamar. 
#### (2) EDQsim was written originally by R. Tsay and stores individual time series in each column.
####  

"EDQts" <- function(x,p=0.5,h=30){
### Input variables:
###  x: an T by m matrix: store m time series in each column and there are T data points in each series.
###  p: a probability, the quantile series of which is to be computed. Defaulty is p=0.5
###  h: the number of time series used in the algorithm used. The larger h is the longer to compute. Default = 30.
###
### Output:
### The column number which stores the EDQ of interest.
if(!is.matrix(x))x <- as.matrix(x)
if(p <= 0) p <- 0.5
if(h < 1)h <- 30

### The following are internal functions used in "EDQts".

### Below is the WDQTS program by R. Zamar.
#### The individual time seies are stored in each row. 
###############################################
"WDQTS" <- function(x,p,h) {
  #This function computes the empirical dynamic p quantile of a set of time series  using an iterative algorithm based on weights 
  ## First the program computes the nearest h series to the series of pointwise quantiles with the L1 metric.
  # Second, it uses these h series as initial values for an iterative algorithm based on weights. 
  # Third, the solution that minimizes the objective function  is reported. 
  # Dependency 
  #  This program uses the subrutines RO and WEPT
  # Input:
  # x : The series must be in the rows of a matrix x of data, in columns the periods 1:T
  # p: the quantile to be computed 
  # h:  the number of series used as starting points in the algorithm
  #Output
  ## ifinal: The index of the series for the q quantile
  
  dd=dim(x)
  N=dd[1]
  T=dd[2]
  q=seq(1,T,1)
  for (tt in 1:T) {
    q[tt]=quantile(x[,tt],p)     }
  
  #q is the pointwise pth quantile of the series  
  #We compute the index of the  closest series to q, i1, and also the indeces of the h closest series to q, that are in orden.
  MA=x-x
  Mq=matrix(rep(q,N),N,T,byrow=TRUE)
  MA=abs(x-Mq)
  vdq=apply(MA,1,sum)
  i1=which.min(vdq)
  orden=sort(vdq,index.return=TRUE)$ix[1:h]
  # Now finds the best series among the h iterative solutions and gives the index of the best found series  in ifinal
  Msol=matrix(0,h,2)
  
  for (k in 1:h)  {
  
    ic=orden[k]
    yy=x[ic,]
    pesos=RO(x,yy,q, p)
    fin= WOPT(x,ic,pesos,q,p)
    
    Msol[k,]=c(fin[[1]],fin[[2]])
  }
  ifinal2=which.min(Msol[,2])[1]
  ifinal=Msol[ifinal2,1]
  
  out=ifinal
  return(out)
}
################################################################

RO=function(x,yy,q,p){
  #This function computes the corresponding T weights for the L1 distance between the series  yy and a given timewise quantile
  # Input:
  # x : The series must be in the rows of a matrix x of data, in columns the periods 1:T
  # yy: a given time series
  # q:  the timewise quantile of order p for the given data 
  # p: the quantile to be computed 
  
  #Output
  ## pesos: The weights for the given time series
  
  dd=dim(x)
  N=dd[1]
  T=dd[2]
  Mb=matrix(0,N,T)
  for (i in 1:N){
    ss1=sign(x[i,]-yy) 
    ss2=sign(x[i,]-q)
    a1=abs(x[i,]-yy)
    a2=abs(x[i,]-q)
    vd1=ifelse (ss1>0,p*a1,(1-p)*a1)
    vd2=ifelse (ss2>0,p*a2,(1-p)*a2)
    vd=vd1-vd2
    
    Mb[i,]=ifelse(abs(yy-q)>0, vd/abs(yy-q),0)   
  }
  pesos=apply(Mb,2,mean)
  out=pesos
  return(out)
}
################################################################
WOPT=function(x,i1,pesos,q,p) {
    # This function iterates from an initial time series and weights and try to find another series that improves the objective function
  #This function computes the corresponding T weights for the L1 distance between the series  yy and a given timewise quantile
  # Input:
  # x : The series must be in the rows of a matrix x of data, in columns the periods 1:T
  # i1: the index of the time series to be used as starting value
  # pesos: The weights for the given time series
  # q:  the timewise quantile of order p for the given data 
  # p: the quantile to be computed 
  
  #Output
  # A list containing:
  # i1: the index of the final time series
  # Fob: The value of the objective funcion for this series
   
  dd=dim(x)
  N=dd[1]
  T=dd[2]
  pesos1=pesos
  dpes1=sum(abs(x[i1,]-q)*pesos1)
  #print(c(i1,dpes1))
  dpes=rep(0,N)
  for (i in 1:N) {
    dpes[i]=sum(abs(x[i,]-q)*pesos1)
  }
  i2=which.min(dpes)
  pesos2=RO(x,x[i2,],q, p)
  dpes2=sum(abs(x[i2,]-q)*pesos2)
  
  while (dpes2<dpes1) {
    # print(c(i2,dpes2))
    i1=i2
    dpes1=dpes2
    pesos=RO(x,x[i1,],q, p)
    for (i in 1:N) {
      dpes[i]=sum(abs(x[i,]-q)*pesos)
    }
    i2=which.min(dpes)
    pesos2=RO(x,x[i2,],q, p)
    dpes2=sum(abs(x[i2,]-q)*pesos2)
  }
  VV=rep(0,N)
  for (j in 1:N) {
    ss1=sign(x[j,]-q) 
    a1=abs(x[j,]-q)
    vd1=ifelse (ss1>0,p*a1,(1-p)*a1)
    
    VV[j]=sum(vd1)  
    
  }
  VVT=sum(VV)
  FOb=dpes1*N+VVT
  #print(c(i1,dpes1))
  #the output is the index and the value of the objective funcion 
  out=list(i1=i1,FOb=FOb)
  return(out)
}

### 
loc <- WDQTS(t(x),p=p,h=h)
return(loc)

#### ENd of "EDQts"
}
##############################       

"HDgen" <- function(m=100,T=100,rho=0.5,AR=TRUE,band=1){
### Generate High-dimensional VAR(1) or VMA(1) model
### inout variables:
### m: dimension of the time series: Defauty is 100. m must be greater than 1.
### T: sample size. Default is 100.
### AR: indicator for VAR(1) model. 
###           AR=TRUE, VAR(1) model is used; AR=FALSE, VMA(1) model is used.
### rho: starting correlation coefficient used to generate contemporaneous cross-dependence between time series: Default is 0.5.
### band: Coefficient matrix used is a band matrix. band denotes the width of the banded matrix. Default us 1.
###
### Output variables:
###  et: The white noise series used in the simulation.
### Xt: The generated time series
### Sigma: The covariance matrix of the white noise series
### Coef: The coefficient matrix used in the simulation.
###
if(m < 2) m <- 2
Sig <- diag(rep(1,m))
for (i in 2:m){
 for (j in 1:(i-1)){
  Sig[i,j] <- rho^(i-j)
  Sig[j,i] <- Sig[i,j]
 }
}
### Compute the square-root matrix of Sig
m1 <- eigen(Sig)
Pmtx <- m1$vectors
L <- diag(sqrt(m1$values))
A1 <- Pmtx%*%L%*%t(Pmtx)
et <- matrix(rnorm(m*T),T,m)
et <- et%*%A1
###
### Generating coefficient matrix
mm1 <- m*(m-1)/2
 phi <- runif(m,min=-0.7,max=0.7)
 Phi <- diag(phi)
 wrk <- runif(mm1,min=-0.5,max=0.5)
 jdx <- 0
if(band > 0){
 for (i in 1:(m-1)){
   iend <- min(m,i+band)
   for (j in (i+1):iend){
    jdx <- jdx+1
    Phi[i,j] <- wrk[jdx]
    }
  }
}
X <- et
phi0 <- matrix(runif(m,min=-1,max=1),m,1)
if(AR){
  for (i in 2:T){
   wk <- phi0+Phi%*%matrix(c(X[i-1,]),m,1)
   X[i,] <- X[i,]+c(wk)
   }
  }else{
 X1 <- rbind(rep(0,m),X[-T,])
 X <- X - X1%*%t(Phi)+matrix(rep(phi0,m),T,m,byrow=T)
 }

HDgen <- list(et=et,Xt=X,Sigma=Sig,Coef=Phi)
}

###
"sqrtMtx" <- function(S){
### Obtain the square-root matrix of a positive-definite matrix
S <- (S+t(S))/2
m1 <- eigen(S)
P <- m1$vectors
L <- diag(sqrt(m1$values))
Sroot <- P%*%L%*%t(P)
Sroot
}


###
"HDgarch" <- function(m=100,T=100,rho=0.5,band=1){
### Generate univariate GARCH(1,1) model with Gaussian innovations
### Multivariate series are obtained by non-singular linear transformation
###
### Input variables are the same as "HDgen".
### Output variables:
### et: The underlying Gaussian white noises used
### Xt: The uncorrelated GARCH(1,1) series
### Yt: The generated vector GARCH(1,1) series. Yt has cross-dependence
### par: GARCH(1,1) parameters used in the simulation.
###
#### Use additional 50 data points to mitigate the effect of initial values.
nobe <- T+50
et <- matrix(rnorm(m*nobe),nobe,m)
X <- NULL; par <- NULL
coef <- runif(m,min=0.75,max=0.98)  ## The persistent parameter
cnst <- runif(m,min=0.02,max=0.15)  ## Constant term
for (i in 1:m){
  ab <- coef[i]
  c0 <- cnst[i]
  beta <- runif(1,min=(ab/1.5),max=ab)
  alpha <- ab-beta
  sig <- c0/(1-ab)
  x <- et[1,i]*sqrt(sig)
  for (t in 2:nobe){
   sig <- c0+alpha*x[t-1]^2+beta*sig
   x <- c(x,et[t,i]*sqrt(sig))
  }
  par <- rbind(par,c(c0,alpha,beta,ab))
  X <- cbind(X,x)
 }
Sig <- diag(rep(1,m))

for (i in 1:(m-1)){
  iend <- min(m,i+band)
  for (j in (i+1):iend){
   Sig[i,j] <- rho^(j-i)  ## make sure j > i
   Sig[j,i] <- Sig[i,j]
 }
}
wk <- sqrtMtx(Sig)
Yt <- X%*%wk
Yt <- Yt[51:nobe,]
cnst <- runif(m,min=-1,max=1)   #### Add means to the observed GARCH series Yt.
Yt <- Yt + matrix(rep(cnst,m),T,m,byrow=TRUE)

HDgarch <- list(et=et,Yt=Yt, Xt=X,par=par)

}

###########
"EDQplot" <- function(x,prob=c(0.05,0.5,0.95),h=30,loc=NULL,color=c("yellow","red","green","blue")){
### Plot the observed time series and selected EDQs
### x: a T-by-m data matrix of m series with T observations. 
### prob: the probability vector contains probabilities to compute EDQ
### loc: locations of the EDQ. If loc is not null, then prob is not used.
### h: number of series used in the algorithm to compute EDS. Default is 30.
### color: colors for plotting the EDQ
if(!is.matrix(x))x <- as.matrix(x)
np <- length(prob)
if(is.null(loc)){
 idx <- NULL
 for (i in 1:np){
  ii <- EDQts(x,p=prob[i],h=h)
  idx <- c(idx,ii)
 }
}else{idx <- loc}
ts.plot(x)
nc <- length(color)
for (i in 1:np){
 ii <- idx[i]
 if(i <= nc) {
        j <- i}else{ j <- i%%nc
         if(j==0)j <- nc
      }
 lines(x[,ii],col=color[j],lwd=2)
 }
###
}