#################################
###  SECTION 2: KEY FUNCTIONS ###
#################################

#############################################################################
RClusProspModInf3 <- function(theta,nderivs=2,y,x,wts,modelfn,Jis=NULL,
                             inxStrat=1, gamma, yStrat, nzval0=20, ...){
#############################################################################
# New dataframe: 
#1. theta - include betas, logsig1, logsig2, and rho;
#2. x - include fixed effects Xs and one random effect variable V;

v <- x[,dim(x)[2],]

## ignore first two cols (inter2/3) and last col (v) in x-array;
x <- x[,-c(1:2,dim(x)[2]),,drop=FALSE]  

nn <- dim(x)[1]
Jmax <-  dim(x)[3]
if (is.null(Jis)) {
    if (any(is.na(x))) stop("RClusProspModInf: Missing values in x when none expected")
    Jis <- rep(Jmax,nn)
}
else {
    Jis <- Jis[Jis[,2]==inxStrat,1]
    if (min(Jis <- as.vector(Jis)) < 1 | max(Jis) > Jmax) stop(paste(
               "RClusProspModInf: Jis- values out of range",range(Jis)))
}

if (nrow(y)!=nn) stop(paste(
               "RClusProspModInf: dimension problems, nrows(y,x)",nrow(y),nn))

if (!is.null(yStrat)) yStrat <- as.vector(yStrat[yStrat[,2]==inxStrat,1])
      else { yStrat <- rep(2, nn); yStrat[y[,1]==1] <- 1 }

if (!is.null(gamma)) {
   gamma <- gamma[,inxStrat]
   if (length(gamma)!=2 || any(gamma<0) || any(gamma>1))
     stop(paste("RClusProspModInf: Wrong gamma values in xStrat", inxStrat))

   pgamma <- rep(1, nn)
   pgamma[apply(2-y,1,sum,na.rm=TRUE)==1 & yStrat==1] <- gamma[1]
   pgamma[apply(2-y,1,sum,na.rm=TRUE)> 1 & yStrat==1] <- gamma[2]
   pgamma[apply(2-y,1,sum,na.rm=TRUE)==1 & yStrat!=1] <- 1-gamma[1]
   pgamma[apply(2-y,1,sum,na.rm=TRUE)> 1 & yStrat!=1] <- 1-gamma[2]
}
 
ntheta <- length(theta)
betas  <- theta[1:(ntheta-3)]
logsig1 <- theta[ntheta-2]
sigma1  <- exp(logsig1)
logsig2 <- theta[ntheta-1]
sigma2  <- exp(logsig2)
rtan <- theta[ntheta]  ##!!
rho <- tanh(rtan) ####!!

eta <- matrix(0,nrow=nn,ncol=Jmax)
for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      eta[ind,j] <- x[ind, ,j]%*%betas
}

loglk <- 0
score <- if (nderivs >= 1) rep(0,ntheta)  else NULL
inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta) else  NULL
dlogldtheta <- NULL

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)
zval1 <- zval2 <- zval
intwts1 <- intwts2 <- intwts

Py <- rep(0,nn)
if (nderivs >= 1) {
    t <- 0
    storePjzs <- as.list(1:nzval^2) 
    storePygz <- array(0,c(nn,nzval,nzval))

    temscore <- matrix(0,nrow=nn,ncol=ntheta)
    temscore2 <- rep(0, length=nn)
}

for (k in 1:length(zval2)){
  for (i in 1:length(zval1)){

    if (nderivs >=1) {
      t <- t + 1
      storePjzs[[t]] <- as.list(1:Jmax)
    }

    Pygz <- rep(1,nn)
    for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      Cjz <- sigma1*zval1[i]+sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*v[ind,j]+eta[ind,j]
      Pjz <- modelfn(y[ind,j], Cjz, nderivs=nderivs, report="vals")
      if (nderivs >= 1) storePjzs[[t]][[j]] <- Pjz
      Pygz[ind] <- Pygz[ind]*Pjz$fy
    }
    Py <- Py + Pygz*intwts1[i]*intwts2[k]

    if (nderivs >=1) {
      storePygz[,i,k] <- Pygz

      for (j in 1:Jmax){
        ind <- (1:nn)[Jis >= j]
	dev1_beta <- matrix(x[ind,,j], nrow=length(ind))
	dev1_ls1 <- sigma1*zval1[i]
	dev1_ls2 <- sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*v[ind,j]
	dev1_rtan <- sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(1-rho^2)*v[ind,j]    ##!!
	dev2_rtan <- sigma2*zval2[k]/((1-rho^2)*sqrt(1-rho^2))*(1-rho^2)^2*v[ind,j] +   
		     sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(2*rho)*(1-rho^2)*v[ind,j]      ##!!
	Xj <- cbind(dev1_beta, dev1_ls1, dev1_ls2, dev1_rtan)   ##!!
        Pjz <- storePjzs[[t]][[j]]

        temscore[ind, ]<-temscore[ind, ]+intwts1[i]*intwts2[k]*Pygz[ind]*(as.vector(Pjz$dfy)/Pjz$fy)*Xj
	temscore2[ind]<-temscore2[ind]+intwts1[i]*intwts2[k]*Pygz[ind]*(as.vector(Pjz$dfy)/Pjz$fy)*dev2_rtan  ##!!
      }
    }

  }
}

if (is.null(gamma)) loglk <- sum(wts*log(Py))
  else loglk <- sum(wts*log(pgamma*Py))

if (nderivs >= 1){
     dlogldtheta <- temscore/Py
     score <- apply(wts*dlogldtheta,2,sum)

     dlogldrtan <- temscore2/Py           ##!!
     rtanscore <- sum(wts*dlogldrtan)     ##!!

}

if (nderivs >=2) {
  inf <- t(dlogldtheta) %*% (wts*dlogldtheta)
  correctmat <- matrix(0,nrow=ntheta,ncol=ntheta)
  t <- 0

  for (k in 1:length(zval2)){
    for (i in 1:length(zval1)){
      Pygz <- storePygz[,i,k]
      multvec <- wts*(Pygz/Py)*intwts1[i]*intwts2[k]
      t <- t + 1

      for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= j]
	 dev1_beta <- matrix(x[ind,,j], nrow=length(ind))
	 dev1_ls1 <- sigma1*zval1[i]
	 dev1_ls2 <- sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*v[ind,j]
	 dev1_rtan <- sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(1-rho^2)*v[ind,j]            ##!!
	 Xj <- cbind(dev1_beta, dev1_ls1, dev1_ls2, dev1_rtan)                         ##!!
         Pjz <- storePjzs[[t]][[j]]

         correctmat <- correctmat + t(Xj) %*% (multvec[ind]*as.vector(Pjz$d2fy)/Pjz$fy*Xj)

         if (j < Jmax) for (l in (j+1):Jmax){ ##***THINKING HERE FOR RAGGED ****
            indl <- (1:nn)[Jis >= l]
	    dev1_beta <- matrix(x[indl,,l], nrow=length(indl))
	    dev1_ls1 <- sigma1*zval1[i]
	    dev1_ls2 <- sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*v[indl,l]
	    dev1_rtan <- sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(1-rho^2)*v[indl,l]   ##!!
	    Xl <- cbind(dev1_beta, dev1_ls1, dev1_ls2, dev1_rtan)                           ##!!
            Plz <- storePjzs[[t]][[l]] 

            jrowuse <- (1:length(ind))[table(factor(indl,levels=ind))==1]
            Xjjrowuse <- if (length(jrowuse)>1) Xj[jrowuse,] 
                         else matrix(Xj[jrowuse,], nrow=1)
            temp <- t(Xjjrowuse) %*% ((multvec[indl])*(as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*
                           (as.vector(Plz$dfy)/Plz$fy)*Xl)
            correctmat <- correctmat + temp + t(temp)
          }
      }
    }
  }
    
  inf <- inf - correctmat

  inf[(ntheta-2),(ntheta-2)] <- inf[(ntheta-2),(ntheta-2)] - score[ntheta-2]
  inf[(ntheta-1),(ntheta-1)] <- inf[(ntheta-1),(ntheta-1)] - score[ntheta-1]
  inf[ntheta,ntheta] <- inf[ntheta,ntheta] + rtanscore 				##!!

  inf[(ntheta-1),ntheta] <- inf[(ntheta-1),ntheta] - score[ntheta]
  inf[ntheta,(ntheta-1)] <- inf[ntheta,(ntheta-1)] - score[ntheta]
}

indcaseclus <- (1:nn)[yStrat==1]
extra <- list(dlogldtheta=dlogldtheta,indcaseclus=indcaseclus) 
list(loglk=loglk, score=as.vector(score), inf=inf, extra=extra)
}


###################################################################################	
MEtaProspModInf3 <- function(theta, nderivs=2, y, x, wts=1, modelfn, inxStrat,
		   	     Jis=NULL, npar, nzval0=20, ...){
###################################################################################

# npar is a vector of length two given the numbers of parameters in each of the models


theta <- as.vector(theta)
ntheta <- length(theta)

theta1 <- theta[1:npar[1]]
theta2 <- theta[(npar[1]+1):sum(npar)]
ntheta1 <- length(theta1)
ntheta2 <- length(theta2)

y1 <- y[,-dim(y)[2],drop=FALSE]
y2 <- as.vector(y[,dim(y)[2]])

x1 <- x[,1:npar[1],,drop=FALSE]  
x2 <- x[,(npar[1]+1):sum(npar),1,drop=FALSE]  
off.set2 <- matrix(0, dim(x2)[1], dim(x2)[3])

loglk <- 0
score <- if (nderivs >= 1) rep(0,ntheta)  else NULL
inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta) else  NULL


part1 <- RClusProspModInf3(theta=theta1, nderivs=nderivs, y=y1, x=x1, wts=wts, modelfn=modelfn, 
			  Jis=Jis, inxStrat=inxStrat, gamma=NULL, nzval0=nzval0, ...)

part2 <- MEtaProspModInf(theta=theta2, nderivs=nderivs, y=y2, x=x2, wts=wts, modelfn=modelfn, 
			 off.set=off.set2, nzval0=nzval0)


loglk <- part1$loglk+part2$loglk
if (nderivs >= 1) score <- c(part1$score, part2$score)
if (nderivs >= 2) {
	inf[1:npar[1],1:npar[1]] <- part1$inf
 	inf[(npar[1]+1):sum(npar), (npar[1]+1):sum(npar)] <- part2$inf
}

list(loglk=loglk, score=as.vector(score), inf=inf) 
}



###############################################################################################
Clusbinlogis3 <- function(y, eta, V, X1, X2, nderivs=2, modelfn, sigma1, sigma2, rtan, Jis=NULL,
			  nzval0=20, ...){
###############################################################################################

### Calculate joint density and its derivatives with respect to theta;
### with dimensions of nn, nn*ntheta and nn*ntheta*ntheta for fy, dfy and d2fy;
### served for MEtaStratModInf2 function only;

### "y" is the response of interest with dimension of (nn, Jmax+1) with the last
### column corresponding to the case-control response as a binary variable

### "eta" has a dimention of (nn, Jmax+1) with the last column corresponding to
###  the linear predictor for binary model part (z|y,x; theta2);
### sigma1/sigma2/rtan are not part of eta1 matrix so input separately
### Jis indicates the size of each cluster, maximum equal to Jmax

fy <- dfy <- d2fy <- NULL

nn <- dim(eta)[1]
Jmax <-  dim(eta)[2]-1

if (is.null(Jis)) Jis <- rep(Jmax,nn)
else {
  if (any(Jis>Jmax)) stop("Clusbinlogis2: Maximum of Jis should not beyond Jmax in eta")
  if (length(Jis) != nn) stop("Clusbinlogis2: Dimension problem in Jis - different from nrow(eta)")
}

if (nrow(y)!=nn) stop("Clusbinlogis2: Dimension problem in y - different from nrow(eta)")


y1 <- y[,-dim(y)[2],drop=FALSE]  # response of interest
y2 <- as.vector(y[,dim(y)[2]])   # case-control response (sampling)

eta1 <- eta[,-dim(eta)[2],drop=FALSE]
eta2 <- eta[,dim(eta)[2],drop=FALSE]
rho <- tanh(rtan)   ##!!

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)
zval1 <- zval2 <- zval
intwts1 <- intwts2 <- intwts

ntheta1 <- dim(X1)[2]+3
ntheta2 <- dim(X2)[2]
ntheta <- ntheta1+ntheta2

Py <- Pz <- rep(0,nn)
if (nderivs >=1) {
    dPy <- matrix(0,nrow=nn,ncol=ntheta1)
    dPyr <- rep(0, length=nn)
    dPz <- matrix(0,nrow=nn,ncol=ntheta2)
}
if (nderivs >= 1) {
    tt <- 0
    storePjzs <- as.list(1:nzval^2) 
    storePygz <- array(0,c(nn,nzval,nzval))
}

for (k in 1:length(zval2)){
  for (i in 1:length(zval1)){

    if (nderivs >=1) {
      tt <- tt + 1
      storePjzs[[tt]] <- as.list(1:Jmax)
    }

    Pygz <- rep(1,nn)
    for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      Cjz <- sigma1*zval1[i]+sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*V[ind,j]+eta1[ind,j]
      Pjz <- modelfn(y1[ind,j],Cjz,nderivs=nderivs,report="vals")
      if (nderivs >= 1) storePjzs[[tt]][[j]] <- Pjz
      Pygz[ind] <- Pygz[ind]*Pjz$fy
    }
    Py <- Py + Pygz*intwts1[i]*intwts2[k]

    if (nderivs >= 1) {
      storePygz[,i,k] <- Pygz

      for (j in 1:Jmax){
        ind <- (1:nn)[Jis >= j]
	dev1_beta <- matrix(X1[ind,,j], nrow=length(ind))
	dev1_ls1 <- sigma1*zval1[i]
	dev1_ls2 <- sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*V[ind,j]
	dev1_rtan <- sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(1-rho^2)*V[ind,j]   ##!!
	dev2_rtan <- sigma2*zval2[k]/((1-rho^2)*sqrt(1-rho^2))*(1-rho^2)^2*V[ind,j] +
		     sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(2*rho)*(1-rho^2)*V[ind,j]     ##!!

	Xj <- cbind(dev1_beta, dev1_ls1, dev1_ls2, dev1_rtan)   ##!!
        Pjz <- storePjzs[[tt]][[j]]

        dPy[ind, ] <- dPy[ind, ] + intwts1[i]*intwts2[k]*Pygz[ind]*(as.vector(Pjz$dfy)/Pjz$fy)*Xj 
        dPyr[ind] <- dPyr[ind] + intwts1[i]*intwts2[k]*Pygz[ind]*(as.vector(Pjz$dfy)/Pjz$fy)*dev2_rtan  ##!! 
      }
    }

  }
}


if (nderivs >= 2) {
  d2Py_a <- d2Py_b <- array(0,c(nn, ntheta1, ntheta1))
  tt <- 0   

  for (k in 1:length(zval2)){
    for (i in 1:length(zval1)){
      tt <- tt + 1
      Pygz <- storePygz[,i,k]

      for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= j]
	 dev1_beta <- matrix(X1[ind,,j], nrow=length(ind))
	 dev1_ls1 <- sigma1*zval1[i]
	 dev1_ls2 <- sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*V[ind,j]
	 dev1_rtan <- sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(1-rho^2)*V[ind,j]   ##!!
	 Xj <- cbind(dev1_beta, dev1_ls1, dev1_ls2, dev1_rtan)   ##!!
         Pjz <- storePjzs[[tt]][[j]]

         Xjnew <- array(t(Xj),c(ntheta1,1,length(ind)))
         tmp <- apply(Xjnew,3,function(x) x %*% t(x))     
         Xj2 <- array(tmp, c(ntheta1, ntheta1, length(ind)))
         Xj2 <- aperm(Xj2,c(3,1,2))
         d2Py_a[ind,,] <- d2Py_a[ind,,] + intwts1[i]*intwts2[k]*Pygz[ind]*as.vector(Pjz$d2fy)/Pjz$fy*Xj2 

         if (j < Jmax) for (l in (j+1):Jmax){ 
            indl <- (1:nn)[Jis >= l]
	    dev1_beta <- matrix(X1[indl,,l], nrow=length(indl))
	    dev1_ls1 <- sigma1*zval1[i]
	    dev1_ls2 <- sigma2*(rho*zval1[i]+sqrt(1-rho^2)*zval2[k])*V[indl,l]
	    dev1_rtan <- sigma2*(zval1[i]-rho/sqrt(1-rho^2)*zval2[k])*(1-rho^2)*V[indl,l]    ##!!
	    Xl <- cbind(dev1_beta, dev1_ls1, dev1_ls2, dev1_rtan)   ##!!
            Plz <- storePjzs[[tt]][[l]]

            jrowuse <- (1:length(ind))[table(factor(indl,levels=ind))==1]
            Xjjrowuse <- if (length(jrowuse)>1) Xj[jrowuse,] 
                         else matrix(Xj[jrowuse,], nrow=1)

	    tmp1 <- tmp2 <- NULL
 	    for (m in 1:ntheta1) {
		tmp1 <- cbind(tmp1,Xjjrowuse[,m]*Xl)
		tmp2 <- cbind(tmp2,Xl[,m]*Xjjrowuse)
	    }
	    Xjj2 <- array(tmp1, c(length(indl),ntheta1,ntheta1))
	    Xl2 <- array(tmp2, c(length(indl),ntheta1,ntheta1))
	 
	    part1 <- intwts1[i]*intwts2[k]*Pygz[indl]*(as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*(as.vector(Plz$dfy)/Plz$fy)*Xjj2
	    part2 <- intwts1[i]*intwts2[k]*Pygz[indl]*(as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*(as.vector(Plz$dfy)/Plz$fy)*Xl2
            d2Py_b[indl,,] <- d2Py_b[indl,,] + part1 + part2  
          }
      }

    }
  }

  d2Py <- d2Py_a + d2Py_b


  d2Py[,ntheta1-2,ntheta1-2] <- d2Py[,ntheta1-2,ntheta1-2] + dPy[,ntheta1-2]
  d2Py[,ntheta1-1,ntheta1-1] <- d2Py[,ntheta1-1,ntheta1-1] + dPy[,ntheta1-1]
  d2Py[,ntheta1,ntheta1] <- d2Py[,ntheta1,ntheta1] - dPyr

  d2Py[,ntheta1-1,ntheta1] <- d2Py[,ntheta1-1,ntheta1] + dPy[,ntheta1]
  d2Py[,ntheta1,ntheta1-1] <- d2Py[,ntheta1,ntheta1-1] + dPy[,ntheta1]
}


zpart <- modelfn(y2,eta2,nderivs=nderivs,report="vals")
Pz <- zpart$fy
if (nderivs >= 1) dPz <- as.vector(zpart$dfy)*X2
if (nderivs >= 2) {
        Xnew <- array(t(X2),c(ntheta2,1,nn))
        tmp <- apply(Xnew,3,function(x) x %*% t(x))
        Xnew2 <- array(tmp, c(ntheta2, ntheta2, nn))
        Xnew2 <- aperm(Xnew2,c(3,1,2))

	d2Pz <- as.vector(zpart$d2fy)*Xnew2
}


fy <- Py*Pz
if (nderivs >= 1) dfy <- cbind(dPy*Pz, Py*dPz)
if (nderivs >= 2) {
	d2fy <- array(0, c(nn, ntheta, ntheta))
	d2fy[,1:ntheta1,1:ntheta1] <- d2Py*Pz
	d2fy[,(ntheta1+1):ntheta,(ntheta1+1):ntheta] <- Py*d2Pz

	tmp1 <- tmp2 <- NULL
	for (t in 1:ntheta1) tmp1 <- cbind(tmp1,dPy[,t]*dPz)
 	for (t in 1:ntheta2) tmp2 <- cbind(tmp2,dPz[,t]*dPy)
	part1 <- array(tmp1, c(nn,ntheta2,ntheta1))
	part2 <- array(tmp2, c(nn,ntheta1,ntheta2))

	d2fy[,(ntheta1+1):ntheta,1:ntheta1] <- part1
	d2fy[,1:ntheta1,(ntheta1+1):ntheta] <- part2

}

list(fy=fy, dfy=dfy, d2fy=d2fy)
}


###########################################################################################
MEtaStratModInf3 <- function(theta,nderivs=2,ymat,xlist,Bcounts,ptildes,modelfn,inxStrat,
			     Jis=NULL, npar, nzval0=20, ...)
###########################################################################################
{ 
# New dataframe: 
#1. theta1 - include betas, logsig1, logsig2, and rho;
#2. x - x1 is for Y-model including the fixed effects Xs and one random effect variable V;
#       x2 is for Z-model;
#3. npar - a vector of length two that gives the numbers of parameters in Y- and Z- models;


theta <- as.vector(theta)
ntheta <- length(theta)

theta1 <- theta[1:npar[1]]
theta2 <- theta[(npar[1]+1):sum(npar)]
ntheta1 <- length(theta1)
ntheta2 <- length(theta2)

betas  <- theta1[1:(ntheta1-3)]
logsig1 <- theta1[ntheta1-2]
sigma1  <- exp(logsig1)
logsig2 <- theta1[ntheta1-1]
sigma2  <- exp(logsig2)
rtan <- theta1[ntheta1]   ##!!


if (length(xlist) != 4) stop("xlist has wrong dimension!")
else if (is.null(names(xlist))) names(xlist) <- c("orig","1","2","xstrata")

strata <- xlist[[match("xstrata",names(xlist))]]
X1  <- (xlist[[match("orig",names(xlist))]])[strata==inxStrat,3:npar[1],1, drop=FALSE]  # reduce to cluster size of 1
X2 <- (xlist[[match("orig",names(xlist))]])[strata==inxStrat,(npar[1]+1):sum(npar),1]
X2_y1 <- (xlist[[match("1",names(xlist))]])[strata==inxStrat,(npar[1]+1):sum(npar),1]
X2_y2 <- (xlist[[match("2",names(xlist))]])[strata==inxStrat,(npar[1]+1):sum(npar),1] 

V <- matrix(X1[,dim(X1)[2],], ncol=dim(X1)[3])  # a matrix of dimention [nn,Jmax];
X1 <- X1[,-dim(X1)[2],,drop=FALSE]  

nn <- dim(X1)[1]
Jmax <-  dim(X1)[3]  #Jmax=1
Jis <- rep(Jmax, nn)

eta1 <- matrix(0,nrow=nn,ncol=Jmax)
for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      eta1[ind,j] <- X1[ind, ,j]%*%betas
}

etarray <- array(0, c(nn,Jmax+1,3))
etarray[,,1] <- cbind(eta1, X2 %*% theta2)
etarray[,,2] <- cbind(eta1, X2_y1 %*% theta2)
etarray[,,3] <- cbind(eta1, X2_y2 %*% theta2)

y11 <- cbind(rep(1,nn), rep(1,nn))
y21 <- cbind(rep(2,nn), rep(1,nn))
y12 <- cbind(rep(1,nn), rep(2,nn))
y22 <- cbind(rep(2,nn), rep(2,nn))

ptildes <- as.vector(ptildes)
if (length(ptildes) != 2) stop("length(ptildes) != 2")
nyStrat <- length(ptildes)


score <- inf <- Qstar <- dQstar <- SptQstar <- dSptQstar <- NULL

loglk <- 0
if (nderivs >= 1) score <- rep(0,ntheta)  
if (nderivs >= 2) inf <- matrix(0,ntheta,ntheta) 

#nzval <- nzval0
#ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
#zval <- ghq20$zeros*sqrt(2)
#intwts <- ghq20$weights/sqrt(pi)

w11 <- Clusbinlogis3(y11,eta=etarray[,,2], V=V, X1=X1, X2=X2_y1, nderivs=nderivs, modelfn=modelfn, 
		     sigma1=sigma1, sigma2=sigma2, rtan=rtan, Jis=Jis, nzval0=nzval0, ...)
w21 <- Clusbinlogis3(y21,eta=etarray[,,3], V=V, X1=X1, X2=X2_y2, nderivs=nderivs, modelfn=modelfn, 
		     sigma1=sigma1, sigma2=sigma2, rtan=rtan, Jis=Jis, nzval0=nzval0, ...)
w12 <- Clusbinlogis3(y12,eta=etarray[,,2], V=V, X1=X1, X2=X2_y1, nderivs=nderivs, modelfn=modelfn, 
		     sigma1=sigma1, sigma2=sigma2, rtan=rtan, Jis=Jis, nzval0=nzval0, ...)
w22 <- Clusbinlogis3(y22,eta=etarray[,,3], V=V, X1=X1, X2=X2_y2, nderivs=nderivs, modelfn=modelfn, 
		     sigma1=sigma1, sigma2=sigma2, rtan=rtan, Jis=Jis, nzval0=nzval0, ...)

pz1 <- w11$fy+w21$fy
pz2 <- w12$fy+w22$fy
Qstar <- cbind(pz1, pz2)
SptQstar <- apply(ptildes*t(Qstar), 2, sum)

loglk <- sum(Bcounts*log(SptQstar))

if (nderivs >= 1) {
    dQstar <- array(0,c(nn,nyStrat,ntheta))
    dQstar[,1,] <- w11$dfy+w21$dfy
    dQstar[,2,] <- w12$dfy+w22$dfy

    dSptQstar <- ptildes[1]*(w11$dfy+w21$dfy) + ptildes[2]*(w12$dfy+w22$dfy)
    dlogdtheta <- dSptQstar/SptQstar
    score <- apply(Bcounts*dlogdtheta,2,sum)
}

if (nderivs >= 2) {
    inf <- t(dlogdtheta) %*% (Bcounts*dlogdtheta)
    
    d2SptQstar <- ptildes[1]*(w11$d2fy+w21$d2fy) + ptildes[2]*(w12$d2fy+w22$d2fy)
    inf2 <- apply(Bcounts*(d2SptQstar/SptQstar),c(2,3),sum)
    inf <- inf - inf2
}    


list(loglk=loglk, score=as.vector(score), inf=inf, Qstar=Qstar, dQstar=dQstar,
     SptQstar=SptQstar, dSptQstar=dSptQstar)
}


