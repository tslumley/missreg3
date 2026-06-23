##################################################################################
datagen <- function(N, sigma, ypars, zpars, nwave=8, index=FALSE)
##################################################################################
{

##### Models to be considered:
### Y-model: rinter + wave(nwave) + rcon(nwave) + rbin(nwave)
### Z-model: inter + Y1
### with unstratified case-control sampling;

a <- rnorm(N, mean=0, sd=sigma)  #an individual intercept pp

# Generate time-variant y and Xs (wave, rbinary)
y <- wave <- rcon <- rbin <- matrix(NA, N, nwave)

for (j in 1:nwave) wave[,j] <- rep(j, N)
for (i in 1:N) {
   rcon[i,] <- rnorm(nwave, 0, 1)  #random variable of Normal(0, 1)
   rbin[i,] <- rbinom(nwave, 1, 0.5)  #random values of 1/0 pp
}

if (index) cat("Distribution of y:\n")
for (j in 1:nwave) {
  Xj <- cbind(rep(1,N), wave[,j], rcon[,j], rbin[,j])
  if (dim(Xj)[2] != length(ypars)) stop("ypars length differs from the design matrix!")

  etaj <- a + Xj %*% ypars
  pyj <- plogis(etaj)
  y[,j] <- 1*(runif(N)<pyj)

  if (index) print(sum(y[,j][y[,j]==1])/N)
}

X2 <- cbind(rep(1,N), y[,1])
if (dim(X2)[2] != length(zpars)) stop("zpars length differs from the design matrix!")

eta2 <- X2 %*% zpars
pz <- plogis(eta2)
z <- 1*(runif(N)<pz)

if (index) {
 cat("\n Distribution of z:\n")
 print(sum(z[z==1]))
 print(sum(z[z==1])/N)
 cat("\n")
}

#### Give population counts
NMat <- matrix(c(table(z)),2,1)
dimnames(NMat) <- list(sort(unique(z)),"N")

ID <- rep(0, N)
ID[z==1] <- 1
ctsamp <- sample((1:N)[z==0], sum(z[z==1]))
ID[ctsamp] <- 1

wave2 <- c(t(wave[ID==1,]))
rcon2 <- c(t(rcon[ID==1,]))
rbin2 <- c(t(rbin[ID==1,]))
y2 <- c(t(y[ID==1,]))
z2 <- rep(z[ID==1], rep(nwave,length(ID[ID==1])))
ClusID <- rep((1:N)[ID==1], rep(nwave,length(ID[ID==1])))

simdata <- data.frame(id=ClusID, wave=wave2, proband=z2, adhd=y2, rcon=rcon2, rbin=rbin2)
simdata$obstype <- rep("retro", dim(simdata)[1])
simdata$probandS <- 2 - simdata$proband #used for ystrata as 1/2 for case/control

list(simdata=simdata, NMat=NMat)
}


##################################################################################
ADHDdatagen <- function(N, sigma, ypars, zpars, nwave=8, index=FALSE)
##################################################################################
{

##### based on ADHD data structure (with nwave as input; sex as stratification)
wave <- 1:nwave
age <- rnorm(N, mean=5.2, sd=0.77)  
sex <- rbinom(N, 1, 0.2)  # 1=female and 0=male
a <- rnorm(N, mean=0, sd=sigma)  #random intercept

y <- time <- matrix(NA, N, nwave)
if (index) cat("Distribution of y:\n")

for (j in 1:nwave) {
  time[,j] <- rep(wave[j], N)
  inter <- rep(1,N)
  Xj <- cbind(inter, time[,j], age, sex)
  if (dim(Xj)[2] != length(ypars)) stop("ypars length differs from X-matrix!")

  etaj <- a + Xj %*% ypars
  pyj <- plogis(etaj)
  y[,j] <- 1*(runif(N)<pyj)

  if (index) print(sum(y[,j][y[,j]==1])/N)
}

inter2 <- rep(1, N)
X2 <- cbind(inter2, y[,1])
if (dim(X2)[2] != length(zpars)) stop("zpars length differs from X2-matrix!")

eta2 <- X2 %*% zpars
pz <- plogis(eta2)
z <- 1*(runif(N)<pz)

if (index) {
 cat("\n Distribution of z:\n")
 print(sum(z[z==1]))
 print(sum(z[z==1])/N)
 cat("\n")
}

#### Give population counts
NMat <- matrix(c(ftable(sex~z)),2,2)
dimnames(NMat) <- list(sort(unique(z)), sort(unique(sex)))

ID <- rep(0, N)
ID[z==1] <- 1
sampM <- sample((1:N)[z==0 & sex==0], sum(z[z==1 & sex==0]))
sampF <- sample((1:N)[z==0 & sex==1], sum(z[z==1 & sex==1]))
ID[sampM] <- 1
ID[sampF] <- 1

age2 <- rep(age[ID==1], rep(length(wave),length(ID[ID==1])))
sex2 <- rep(sex[ID==1], rep(length(wave),length(ID[ID==1])))
time2 <- c(t(time[ID==1,]))
y2 <- c(t(y[ID==1,]))
z2 <- rep(z[ID==1], rep(length(wave),length(ID[ID==1])))
ClusID <- rep((1:N)[ID==1], rep(length(wave),length(ID[ID==1])))

simdata <- data.frame(id=ClusID, wave=time2, proband=z2, adhd=y2, age=age2, sexF=sex2)
simdata$obstype <- rep("retro", dim(simdata)[1])
simdata$probandS <- 2 - simdata$proband #used for ystrata as 1/2 for case/control

list(simdata=simdata, NMat=NMat)
}



####################################################################################
clmodelfn <- function(theta, ndev=2, z, y, X, clusid, intraclusid, modelfn, Nz, ...) 
####################################################################################
{

### WITHIN EACH X-STRATUM:
### calculate loglk, score, and inf of the conditional likelihood function;

### theta == (alpha, beta);
### z is the categorical sampling variable of length n, with K categories (K=2 for binary);
### y is the binary response of length n, which values are varied within each cluster;
### X is the design matrix of dim[n,p+1] with intercept and p covariates of interest;
### clusid is the ID of clusters in the data;
### intraclusid is the ID of individual observations within a cluster;
### modelfn relates to the link function of h(z_k | y_1; alpha);
### Nz is population counts of z-categories within x-stratum, sorted by nature order of z-values;

### Notes: 
# only y_1 within each cluster is associated with h(z_k | y_1; alpha);
 

theta <- as.vector(theta)
ntheta <- length(theta)

alpha <- theta[1:2] # a model with y_1 only
beta <- theta[3:ntheta] #the model of interest
if (length(beta) != dim(X)[2]) stop("Dimention of X differs from theta for beta coefficients")

mat <- cbind(clusid, intraclusid, z, y, X)
miss <- (1:dim(mat)[1])[apply(mat,1,function(x) any(is.na(x)))]

mydat <- data.frame(clusid, intraclusid, z, y, X)[-miss,,drop=FALSE]
mydat <- mydat[with(mydat, order(clusid, intraclusid)), ]

mycluster <- mydat[!duplicated(mydat$clusid), ]
ncluster <- dim(mycluster)[1]

K <- length(unique(mycluster$z))
if (K != length(Nz)) stop("Nz has different length from the number of z-categories")
nz <- as.numeric(table(mycluster$z))
pz <- nz/Nz
 
eta1 <- alpha[1]+alpha[2]*mycluster$y
lhzy <- modelfn(y=mycluster$z,eta=eta1)
hzy <- modelfn(y=mycluster$z,eta=eta1,report="dens")


# All 2^(Jmax) permutations of (0,1) with the cluster size of Jmax 
#library(gtools)

Js <- unique(as.numeric(table(mydat$clusid)))
if (length(Js) != 1) warning("Variable cluster sizes are found!")
Jmax <- max(Js)

T <- A <- rep(0, ncluster)
dA1 <- matrix(0, ncluster, length(alpha))
dA2 <- matrix(0, ncluster, length(beta))
d2A11 <- array(0, c(ncluster, length(alpha), length(alpha)))
d2A12 <- array(0, c(ncluster, length(alpha), length(beta)))
d2A21 <- array(0, c(ncluster, length(beta), length(alpha)))
d2A22 <- array(0, c(ncluster, length(beta), length(beta)))

Xcols <- (dim(mydat)[2]-dim(X)[2]+1):dim(mydat)[2]

for (i in 1:ncluster)  {
   mydati <- mydat[mydat$clusid==mycluster$clusid[i],,drop=FALSE]
   Ji <- dim(mydati)[1] 

   ystar0 <- permutations(n=2,r=Ji,v=c(0,1),repeats.allowed=TRUE)
   ytotal0 <- apply(ystar0, 1, sum)

   T[i] <- sum(mydat$y[mydat$clusid==mycluster$clusid[i]]) 
   ystari <- ystar0[ytotal0==T[i],,drop=FALSE]
   
   Xi <- mydat[mydat$clusid==mycluster$clusid[i],Xcols,drop=FALSE]

   for (j in 1:dim(ystari)[1]) {
      fyxj <- apply(ystari[j,]*Xi,2,sum) %*% beta
      etaj <- rep(alpha[1]+alpha[2]*ystari[j,1],K)
      zj <- sort(unique(mycluster$z))
      hzyj <- modelfn(y=zj, eta=etaj, report="dens")

      A[i] <- A[i] + exp(fyxj)*sum(pz*hzyj$fy)
      dA1[i,1] <- dA1[i,1] + exp(fyxj)*sum(pz*hzyj$dfy)
      dA1[i,2] <- dA1[i,2] + exp(fyxj)*sum(pz*hzyj$dfy*ystari[j,1])

      dbeta <- apply(ystari[j,]*Xi,2,sum)
      dA2[i,] <- dA2[i,] + dbeta*exp(fyxj)*sum(pz*hzyj$fy)

      detaj <- matrix(c(1, ystari[j,1], ystari[j,1], ystari[j,1]^2),2,2)
      for (k in 1:K) 
       d2A11[i,,] <- d2A11[i,,] + as.vector(exp(fyxj)*pz[k]*hzyj$d2fy[k,1,1])*detaj
      
      d2A22[i,,] <- d2A22[i,,] + ((dbeta*exp(fyxj)) %*% t(dbeta))*sum(pz*hzyj$fy)

      d2A12[i,1,] <- d2A12[i,1,] + dbeta*exp(fyxj)*sum(pz*hzyj$dfy)
      d2A12[i,2,] <- d2A12[i,2,] + dbeta*exp(fyxj)*sum(pz*hzyj$dfy*ystari[j,1])
   }
}

loglk <- 0
score <- if(ndev >= 1) rep(0, ntheta) else NULL
inf <- if(ndev >= 2) matrix(0, ntheta, ntheta) else NULL

mydatX <- mydat[,Xcols,drop=FALSE]
loglk <- sum(nz*log(pz)) + sum(lhzy$logfy) + apply(mydat$y*mydatX,2,sum) %*% beta - sum(log(A))   

if (ndev >= 1) {
  score[1] <- sum(hzy$dfy/hzy$fy) - sum(dA1[,1]/A)
  score[2] <- sum(hzy$dfy*mycluster$y/hzy$fy) - sum(dA1[,2]/A)

  score[3:ntheta] <- apply(mydat$y*mydatX,2,sum) - apply(dA2/A,2,sum)
} 

if (ndev == 2) {
  deta <- array(c(rep(1,ncluster), mycluster$y, mycluster$y, mycluster$y^2), 
	    c(ncluster, length(alpha), length(alpha)))

  tmp1 <- array(NA, c(ncluster, length(alpha), length(alpha)))
  for (l in 1:length(alpha)) tmp1[,,l] <- dA1*dA1[,l]
  
  part11 <- (as.vector(hzy$d2fy)/hzy$fy-(as.vector(hzy$dfy)/hzy$fy)^2)*deta - (d2A11/A-tmp1/A^2)
  inf[1:length(alpha), 1:length(alpha)] <- -apply(part11,c(2,3),sum)

  tmp2 <- array(NA, c(ncluster, length(beta), length(beta)))
  for (l in 1:length(beta)) tmp2[,,l] <- dA2*dA2[,l]

  part22 <- -(d2A22/A-tmp2/A^2)
  inf[(length(alpha)+1):ntheta, (length(alpha)+1):ntheta] <- -apply(part22,c(2,3),sum)

  tmp12 <- array(NA, c(ncluster, length(alpha), length(beta)))
  for (l in 1:length(beta)) tmp12[,,l] <- dA1*dA2[,l]

  part12 <- -(d2A12/A-tmp12/A^2)
  inf[1:length(alpha), (length(alpha)+1):ntheta] <- -apply(part12,c(2,3),sum)
  inf[(length(alpha)+1):ntheta, 1:length(alpha)] <- -t(apply(part12,c(2,3),sum)) 

} 

list(loglk=loglk, score=as.vector(score), inf=inf)
} 


###################################################################################
clmodelfn2 <- function(theta, ndev=2, z, y, X, clusid, intraclusid, modelfn, NzMat,  
		       xStrat=NULL, ...) 
###################################################################################
{

### calculate loglk, score, and inf of the conditional likelihood function;

theta <- as.vector(theta)
ntheta <- length(theta)

loglk <- 0
score <- if(ndev >= 1) rep(0, ntheta) else NULL
inf <- if(ndev >= 2) matrix(0, ntheta, ntheta) else NULL

if(is.null(xStrat)) xStrat <- rep(1,dim(X)[1])
maxxStrat <- max(xStrat)

for (i in 1:maxxStrat){
     zi <- z[xStrat==i]
     yi <- y[xStrat==i]
     Xi <- X[xStrat==i,,drop=FALSE]
     clusidi <- clusid[xStrat==i]
     intraclusidi <- intraclusid[xStrat==i]
     Nzi <- NzMat[,i]

     w <- clmodelfn(theta=theta, ndev=ndev, z=zi, y=yi, X=Xi, clusid=clusidi,
		    intraclusid=intraclusidi, modelfn=modelfn, Nz=Nzi)

     loglk <- loglk + w$loglk
     if (ndev >= 1) score <- score + as.vector(w$score)
     if (ndev >= 2) inf <- inf + w$inf
} 

list(loglk=loglk, score=as.vector(score), inf=inf)
} 


########################################################################################
condclusbin <- function(formula, zvar, clusid, intraclusid, data, NzMat, start, modelfn, 
			xstrata=NULL, devcheck=FALSE, control=mlefn.control(...), 
			control.inner=mlefn.control.inner(...), ...)
########################################################################################
{
# NzMat: categories of z and xstrata must be in default nature order!


ymodel <- model.frame(formula, data, na.action=na.pass)
y <- ymodel[,1]

#X <- model.matrix(formula, data, na.action=na.pass)[,-1,drop=FALSE]
X<- model.apply(formula, data)$X[,-1,drop=FALSE]

z <- model.frame(zvar, data, na.action=na.pass)[,1]
clusid <- model.frame(clusid, data, na.action=na.pass)[,1]
intraclusid <- model.frame(intraclusid, data, na.action=na.pass)[,1]

xStrat <- if (is.null(xstrata)) rep(1,dim(X)[1]) else {
    xstratform <- as.formula(paste("~","-1 + ",paste(xstrata,sep="",collapse=":")))
    data1 <- as.data.frame(lapply(data, factor))
    xstratmat <- model.matrix(xstratform, model.frame(xstratform,data1,
                              na.action=function(x)x))
    for(i in 1:ncol(xstratmat)) xstratmat[,i] <- i * xstratmat[,i]
    apply(xstratmat,1,sum)
}
names(xStrat) <- NULL
nxStrat <- max(xStrat)

if (nxStrat != dim(NzMat)[2])
   stop("Number of xstrata is different from dim(NzMat)[2]!")

xsvar <- data[,match(xstrata,names(data))]
xkey <- as.character(sort(unique(xsvar)))
if (!all(xkey==dimnames(NzMat)[[2]])) warning("dimnames(NzMat)[[2]] are not in default categories!")

if (devcheck) {
  divchk(theta=start, loglkfn=clmodelfn2, nderivs=2, z=z, y=y, X=X, clusid=clusid, 
	 intraclusid=intraclusid, modelfn=modelfn, NzMat=NzMat, xStrat=as.numeric(xStrat))
  stop("Derivatives check on clmodelfn")
}

res <-  mlefn(theta=start, loglkfn=clmodelfn2, z=z, y=y, X=X, clusid=clusid, 
	      intraclusid=intraclusid, modelfn=modelfn, NzMat=NzMat, 
	      xStrat=as.numeric(xStrat), control=control, control.inner=control.inner)  

names(res$theta) <- names(res$score) <- c("zinter", names(ymodel)[1], dimnames(X)[[2]])
dimnames(res$inf) <- list(names(res$theta),names(res$theta))
covmat <- solve(res$inf); dimnames(covmat) <- dimnames(res$inf)
dd <- sqrt(diag(covmat))
correlation <- covmat/outer(dd, dd)

# ------------------------------------------------------------------------------------
ans<-list(coefficients=res$theta, loglk=res$loglk, score=res$score, inf=res$inf, cov=covmat, 
	  dd=dd, cor=correlation)
class(ans)<-"condclusbin"
ans
}


























