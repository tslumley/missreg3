#############################################################################
RClusProspModInf <- function(theta,nderivs=2,y,x,wts,modelfn,Jis=NULL,
                             inxStrat=1, gamma, yStrat, nzval0=20, ...){
#############################################################################
x <- x[,-1,,drop=FALSE]  ## ignore the first variable (inter2) in XArray

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
betas  <- theta[-ntheta]
logsig <- theta[ntheta]
sigma  <- exp(logsig)

eta <- matrix(0,nrow=nn,ncol=Jmax)
for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      eta[ind,j] <- x[ind, ,j]%*%betas
}

loglk <- 0
score <- if (nderivs >= 1) rep(0,ntheta)  else NULL
inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta) else  NULL
dlogldtheta <- NULL

#zval <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
#nzval <- length(zval)
#intwts <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
#intwts[3] <- 8/15

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)

Py <- rep(0,nn)
if (nderivs >=1) temscore <- matrix(0,nrow=nn,ncol=ntheta)
if (nderivs >= 1) {
    storePjzsi <- as.list(1:nzval)
    storePygz <- matrix(0,nrow=nn,ncol=nzval)
}

for (i in 1:length(zval)){
    Pygz <- rep(1,nn)
    if (nderivs >= 1) storePjzsi[[i]] <- as.list(1:Jmax)
    for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      Pjz <- modelfn(y[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
      if (nderivs >= 1) storePjzsi[[i]][[j]] <- Pjz
      Pygz[ind] <- Pygz[ind]*Pjz$fy
    }
    if (nderivs >= 1) storePygz[,i] <- Pygz
    Py <- Py + Pygz*intwts[i]
    if (nderivs >=1) {
      for (j in 1:Jmax){
        ind <- (1:nn)[Jis >= j]
        Xj <- if (length(ind)>1) cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
              else matrix(c(x[ind, ,j],sigma*zval[i]),nrow=1)
      # Pjz <- modelfn(y[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
        Pjz <- storePjzsi[[i]][[j]]
        temscore[ind, ] <- temscore[ind, ] + Pygz[ind]*intwts[i]*
                           (as.vector(Pjz$dfy)/Pjz$fy)*Xj
      }
    }
}

if (is.null(gamma)) loglk <- sum(wts*log(Py))
  else loglk <- sum(wts*log(pgamma*Py))

if (nderivs >= 1){
     dlogldtheta <- temscore/Py
     score <- apply(wts*dlogldtheta,2,sum)
}

if (nderivs >=2) {
    inf <- t(dlogldtheta) %*% (wts*dlogldtheta)
    correctmat <- matrix(0,nrow=ntheta,ncol=ntheta)
    for (i in 1:length(zval)){
      Pygz <- storePygz[,i]
      multvec <- wts*(Pygz/Py)*intwts[i]
      for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= j]
         Xj <- if (length(ind)>1) cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
               else  matrix(c(x[ind, ,j],sigma*zval[i]),nrow=1)
       # Pjz <- modelfn(y[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
         Pjz <- storePjzsi[[i]][[j]]
         correctmat <- correctmat + t(Xj) %*% (multvec[ind]*as.vector(Pjz$d2fy)/Pjz$fy*Xj)
         if (j < Jmax) for (k in (j+1):Jmax){ ##***THINKING HERE FOR RAGGED ****
            indk <- (1:nn)[Jis >= k]
            Xk <- if (length(indk)>1)  cbind(x[indk, ,k],sigma*rep(zval[i],length(indk)))
                  else matrix(c(x[indk, ,k],sigma*zval[i]),nrow=1)
          # Pkz <- modelfn(y[indk,k],sigma*zval[i]+eta[indk,k],nderivs=nderivs,report="vals")
            Pkz <- storePjzsi[[i]][[k]]
            jrowuse <- (1:length(ind))[table(factor(indk,levels=ind))==1]
            Xjjrowuse <- if (length(jrowuse)>1) Xj[jrowuse,] 
                         else matrix(Xj[jrowuse,], nrow=1)
            temp <- t(Xjjrowuse) %*% ((multvec[indk])*(as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*
                           (as.vector(Pkz$dfy)/Pkz$fy)*Xk)
            correctmat <- correctmat + temp + t(temp)
          }
      }
    }
    inf <- inf - correctmat
    inf[ntheta,ntheta] <- inf[ntheta,ntheta] - score[ntheta]
}

indcaseclus <- (1:nn)[yStrat==1]
extra <- list(dlogldtheta=dlogldtheta,indcaseclus=indcaseclus) 
list(loglk=loglk, score=as.vector(score), inf=inf, extra=extra)
}


#################################################################################
ProspInf <- function(theta,nderivs=2,ProspModInf,x,y,Aposn,Acounts, 
                     xStrat=rep(1,dim(x)[1]),off.set=matrix(0,dim(x)[1],dim(x)[3]), 
                     control.inner=mlefn.control.inner(...), ...){
#################################################################################
## this function accounts for strata when ProspModInf function is used 
## it gives identical results with or without stratification (prospective sample)

n <- length(y[,1])
maxxStrat <- max(xStrat)
if ((dim(x)[1] != n) || dim(x)[2] != length(theta) || (length(xStrat) != n) ||
       max(Aposn)>n || min(Aposn)<1 )
    stop("ProspInf: data dimension problems")

loglk <- 0
ntheta <- length(theta)
score <- if (nderivs >= 1) rep(0,ntheta)  else NULL
inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta) else  NULL

for (i in 1:maxxStrat){
    Aposni <- Aposn[xStrat[Aposn]==i]
    if (length(Aposni)<1) next()
    Acountsi <- Acounts[xStrat[Aposn]==i]

    w <- ProspModInf(theta,nderivs=2,y=y[Aposni,,drop=FALSE],x=x[Aposni,,, drop=FALSE],
                     wts=Acountsi,inxStrat=i,off.set=off.set[Aposni,,drop=FALSE], ...)
    error <- 0
    if (!is.null(w$error)) if (w$error != 0) {
	error <- 1
        loglk <- score <- inf <- NULL
     	return(list(loglk=loglk,score=score,inf=inf,error=error))  }
 
    loglk <- loglk + w$loglk
    if (nderivs >= 1) score <- score + as.vector(w$score)
    if (nderivs >= 2) inf <- inf + w$inf
}

list(loglk=loglk,score=as.vector(score),inf=inf,error=error)
}


###################################################################################
RProbandStratModInf <- function(theta,nderivs=2,ptildes,x, Bcounts, modelfn, off.set=0,
                       nzval0=20, ...){
###################################################################################
# this function is used when Y-stratum is defined by case-control status of the proband
# and this stratification is also applied in the loglikelihood;

x <- x[,-1,,drop=FALSE]  ## ignore the first variable (inter2) in XArray

ntheta <- length(theta)
betas  <- theta[-ntheta]
logsig <- theta[ntheta]
sigma  <- exp(logsig)
eta <- matrix(x[, ,1]%*%betas,ncol=1)

loglk <- 0
score <- inf <- NULL
Qstar <-  dQstar <- SptQstar <- dSptQstar <- NULL

#zval <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
#nzval <- length(zval)
#intwts <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
#intwts[3] <- 8/15

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)

P1y <- rep(0,nrow(x[, ,1]))
dP1y <- if (nderivs >=1) matrix(0,nrow=nrow(x[, ,1]),ncol=ntheta) else NULL
for (i in 1:length(zval)){
    P1z <- modelfn(rep(1,nrow(x[, ,1])),sigma*zval[i]+eta[,1],nderivs=nderivs,report="vals")
    P1y <- P1y + P1z$fy*intwts[i]
    if (nderivs >=1) {
         X1 <- cbind(x[, ,1],sigma*rep(zval[i],nrow(x[, ,1])))
         dP1y <- dP1y    +   intwts[i]*(as.vector(P1z$dfy)*X1)
    }
}
SptQstar  <- (ptildes[1]-ptildes[2])*P1y + ptildes[2]  # = ptildes[1]*P1y + ptildes[2]*(1-P1y)
if (nderivs >=1) dSptQstar <- (ptildes[1]-ptildes[2])*dP1y

loglk <- sum(Bcounts*log(SptQstar))
if (nderivs >= 1){
    score <- apply(Bcounts*dSptQstar/SptQstar,2,sum)
    dQstar <- array(0,c(nrow(x[, ,1]),2,ntheta)); dQstar[,1,] <- dP1y ;  dQstar[,2,] <- - dP1y
}
if (nderivs >=2) {
    tempmat <- dSptQstar/SptQstar
    inf <- t(tempmat) %*% (Bcounts*tempmat)
    correctmat <- matrix(0,nrow=ntheta,ncol=ntheta)
    correctcorner <- 0
    for (i in 1:length(zval)){
         X1 <- cbind(x[, ,1],sigma*rep(zval[i],nrow(x[, ,1])))
         P1z <- modelfn(rep(1,nrow(x[, ,1])),sigma*zval[i]+eta[,1],nderivs=nderivs,report="vals")
         multvec <- Bcounts*intwts[i]
         correctmat <- correctmat + t(X1) %*% (multvec*(ptildes[1]-ptildes[2])*
                       (as.vector(P1z$d2fy)/SptQstar)*X1)
         correctcorner <- correctcorner + sum(multvec*(ptildes[1]-ptildes[2])*
                          zval[i]*P1z$dfy*sigma/SptQstar)
    }
    inf <- inf - correctmat
    inf[ntheta,ntheta] <- inf[ntheta,ntheta] - correctcorner
}
list(loglk=loglk,score=as.vector(score),inf=inf,Qstar=cbind(P1y,1-P1y),dQstar=dQstar,
     SptQstar=SptQstar,dSptQstar=dSptQstar)
}


#############################################################################
R1StratModInf <- function(theta,nderivs=2,ptildes,x, Bcounts, modelfn, off.set=0, 
                       Jis=NULL, inxStrat=1, gamma,  nzval0=20, ...){
#############################################################################
# this function is used when Y-stratum is defined by case-control status of the proband
# but with a new stratification function (with gamma) applied in the loglikelihood;

x <- x[,-1,,drop=FALSE]  ## ignore the first variable (inter2) in XArray

gamma <- gamma[,inxStrat]
if (length(gamma)!=2 || any(gamma<0) || any(gamma>1))
  stop(paste("R1StratModInf: Wrong gamma values in xStrat", inxStrat))

nn <- dim(x)[1]
Jmax <-  dim(x)[3]
if (is.null(Jis)) {
    if (any(is.na(x))) stop("R1StratModInf: Missing values in x when none expected")
    Jis <- rep(Jmax,nn)
}
else {
    Jis <- Jis[Jis[,2]==inxStrat,1]
    if (min(Jis <- as.vector(Jis)) < 1 | max(Jis) > Jmax) stop(paste(
               "R1StratModInf: Jis- values out of range",range(Jis)))
}

ntheta <- length(theta)
betas  <- theta[-ntheta]
logsig <- theta[ntheta]
sigma  <- exp(logsig)

eta <- matrix(0,nrow=nn,ncol=Jmax)
for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      eta[ind,j] <- x[ind, ,j]%*%betas
}

loglk <- 0
score <- inf <- NULL
Qstar <-  dQstar <- SptQstar <- dSptQstar <- NULL

#zval <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
#nzval <- length(zval)
#intwts <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
#intwts[3] <- 8/15

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)

if (nderivs >=1) 
    temscore0 <- temscore1 <- matrix(0,nrow=nn,ncol=ntheta) 
else temscore0 <- temscore1 <- NULL
if (nderivs >= 2)
    teminf0 <- teminf1 <- array(0,c(nn,ntheta,ntheta))
else teminf0 <- teminf1 <- NULL

if (nderivs >= 1) {
    storePjzsi1 <- as.list(1:Jmax)
    storePygz1 <- array(0,c(nn,Jmax,nzval))
    storePjzsi0 <- as.list(1:nzval)
    storePygz0 <- matrix(0,nrow=nn,ncol=nzval)
}

y0 <- matrix(0,nrow=nn,ncol=Jmax)
Py0 <- rep(0,nn)
for (i in 1:length(zval)){  
  if (nderivs >= 1) storePjzsi0[[i]] <- as.list(1:Jmax)
  Pygz <- rep(1,nn)
  
  for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      Pjz <- modelfn(y0[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
      if (nderivs >= 1) storePjzsi0[[i]][[j]] <- Pjz
      Pygz[ind] <- Pygz[ind]*Pjz$fy
  }
  Py0 <- Py0 + Pygz*intwts[i]

  if (nderivs >= 1) {
     storePygz0[,i] <- Pygz
     for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= j]
         Xj <- cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
         Pjz <- storePjzsi0[[i]][[j]]
         temscore0[ind, ] <- temscore0[ind, ] + Pygz[ind]*intwts[i]*
                            (as.vector(Pjz$dfy)/Pjz$fy)*Xj
     }
  }

  if (nderivs >= 2) {
     for (j in 1:Jmax) {
       ind <- (1:nn)[Jis >= j]
       Xj <- cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
       Pjz <- storePjzsi0[[i]][[j]]
       part0 <- diag(Pygz[ind]*intwts[i]*as.vector(Pjz$d2fy)/Pjz$fy)
       id <- 1
       for (t in ind) {
          teminf0[t,,] <- teminf0[t,,] + t(Xj) %*% (part0[,id]*Xj)
          id <- id + 1 }
       if (j < Jmax) for (k in (j+1):Jmax) {
          indk <- (1:nn)[Jis >= k]
          Xk <- cbind(x[indk, ,k],sigma*rep(zval[i],length(indk)))
          Pkz <- storePjzsi0[[i]][[k]]
          jrowuse <- (1:length(ind))[table(factor(indk,levels=ind))==1]
          partk0 <- diag(Pygz[indk]*intwts[i]*(as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*
                    (as.vector(Pkz$dfy)/Pkz$fy))
          id <- 1
          for (tt in indk) {
             temp <- t(Xj[jrowuse,]) %*% (partk0[,id]*Xk)
             teminf0[tt,,] <- teminf0[tt,,] + temp + t(temp) 
             id <- id + 1 }
       }
    }
  }
}

Py1 <- rep(0, nn)
for (l in 1:Jmax) {
  y1 <- y0;  y1[,l] <- 1
  if (nderivs >= 1) storePjzsi1[[l]] <- as.list(1:nzval)
  for (i in 1:length(zval)){
    Pygz <- rep(1,nn)
    id0 <- (1:nn)[Jis < l]; Pygz[id0] <- 0
    if (nderivs >= 1) storePjzsi1[[l]][[i]] <- as.list(1:Jmax)
    for (j in 1:Jmax){
       ind <- (1:nn)[Jis >= max(j,l)]
       Pjz <- modelfn(y1[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
       if (nderivs >= 1) storePjzsi1[[l]][[i]][[j]] <- Pjz  
       Pygz[ind] <- Pygz[ind]*Pjz$fy
    }
    Py1 <- Py1 + Pygz*intwts[i]

    if (nderivs >= 1) {
       storePygz1[,l,i] <- Pygz
       for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= max(j,l)]
         Xj <- cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
         Pjz <- storePjzsi1[[l]][[i]][[j]]
         temscore1[ind, ] <- temscore1[ind, ] + Pygz[ind]*intwts[i]*
                            (as.vector(Pjz$dfy)/Pjz$fy)*Xj
       }
    }

    if (nderivs >= 2) {
      for (j in 1:Jmax) {
        ind <- (1:nn)[Jis >= max(j,l)]
        Xj <- cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
        Pjz <- storePjzsi1[[l]][[i]][[j]]
        part1 <- diag(Pygz[ind]*intwts[i]*as.vector(Pjz$d2fy)/Pjz$fy)
        id <- 1
        for (t in ind) {
          teminf1[t,,] <- teminf1[t,,] + t(Xj) %*% (part1[,id]*Xj)
          id <- id + 1 }
        if (j < Jmax) for (k in (j+1):Jmax) {
          indk <- (1:nn)[Jis >= max(k,l)]
          Xk <- cbind(x[indk, ,k],sigma*rep(zval[i],length(indk)))
          Pkz <- storePjzsi1[[l]][[i]][[k]]
          jrowuse <- (1:length(ind))[table(factor(indk,levels=ind))==1]
          partk1 <- diag(Pygz[indk]*intwts[i]*(as.vector(Pjz$dfy)/
                  Pjz$fy)[jrowuse]*(as.vector(Pkz$dfy)/Pkz$fy))
          id <- 1
          for (tt in indk) {
             temp <- t(Xj[jrowuse,]) %*% (partk1[,id]*Xk)
             teminf1[tt,,] <- teminf1[tt,,] + temp + t(temp)
             id <- id + 1 }
        }
      }
    }  
  }
}


Py <- gamma[1]*Py1 + gamma[2]*(1-Py1-Py0)  
temscore <- gamma[1]*temscore1 - gamma[2]*(temscore1+temscore0)
teminf <- gamma[1]*teminf1 - gamma[2]*(teminf1+teminf0)

Qstar <- cbind(Py,1-Py)
SptQstar  <- ptildes[1]*Qstar[,1] + ptildes[2]*Qstar[,2]

error <- 0
if (min(SptQstar) <= 0) {
   error <- 1
   print("R1StratModInf: WARNING -- -ve values to log, loglk evaluation failed")
}

loglk <- sum(Bcounts*log(SptQstar))
if (nderivs >= 1){
    dQstar <- array(0,c(nrow(x[, ,1]),2,ntheta))
    dQstar[,1,] <- temscore; dQstar[,2,] <- -temscore
    dSptQstar <- ptildes[1]*dQstar[,1,] + ptildes[2]*dQstar[,2,]
    tempmat <- dSptQstar/SptQstar
    score <- apply(Bcounts*tempmat,2,sum)
}

if (nderivs >=2) {
  inf <-  t(tempmat) %*% (Bcounts*tempmat)

  correctmat <- (ptildes[1]-ptildes[2])*(teminf/SptQstar)
  inf1 <- apply(Bcounts*correctmat, c(2,3), sum)
  inf <- inf - inf1

  inf[ntheta,ntheta] <- inf[ntheta,ntheta] - score[ntheta]
}

list(loglk=loglk,score=as.vector(score),inf=inf,Qstar=Qstar,dQstar=dQstar,
     SptQstar=SptQstar,dSptQstar=dSptQstar, error=error)
}


#############################################################################
R0StratModInf <- function(theta,nderivs=2,ptildes,x, Bcounts, modelfn, off.set=0, 
                          Jis=NULL, inxStrat=1, nzval0=20, ...){
#############################################################################
# this function is used when Y-stratum is defined by control-status of all individuals
# in a family, along with this stratification applied in the loglikelihood;

x <- x[,-1,,drop=FALSE]  ## ignore the first variable (inter2) in XArray

nn <- dim(x)[1]
Jmax <-  dim(x)[3]
if (is.null(Jis)) {
    if (any(is.na(x))) stop("R0StratModInf: Missing values in x when none expected")
    Jis <- rep(Jmax,nn)
}
else {
    Jis <- Jis[Jis[,2]==inxStrat,1]
    if (min(Jis <- as.vector(Jis)) < 1 | max(Jis) > Jmax) stop(paste(
               "R0StratModInf: Jis- values out of range",range(Jis)))
}

y <- matrix(0,nrow=nn,ncol=Jmax)

ntheta <- length(theta)
betas  <- theta[-ntheta]
logsig <- theta[ntheta]
sigma  <- exp(logsig)

eta <- matrix(0,nrow=nn,ncol=Jmax)
for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      eta[ind,j] <- x[ind, ,j]%*%betas
}

loglk <- 0
score <- inf <- NULL
Qstar <-  dQstar <- SptQstar <- dSptQstar <- NULL

#zval <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
#nzval <- length(zval)
#intwts <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
#intwts[3] <- 8/15

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)

Py <- rep(0,nn)
if (nderivs >=1) temscore <- matrix(0,nrow=nn,ncol=ntheta) else temscore <- NULL
if (nderivs >= 1) {
    storePjzsi <- as.list(1:nzval)
    storePygz <- matrix(0,nrow=nn,ncol=nzval)
}
for (i in 1:length(zval)){
    Pygz <- rep(1,nn)
    if (nderivs >= 1) storePjzsi[[i]] <- as.list(1:Jmax)
    for (j in 1:Jmax){
          ind <- (1:nn)[Jis >= j]
          Pjz <- modelfn(y[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
          if (nderivs >= 1) storePjzsi[[i]][[j]] <- Pjz
          Pygz[ind] <- Pygz[ind]*Pjz$fy
    }
    if (nderivs >= 1) storePygz[,i] <- Pygz
    Py <- Py + Pygz*intwts[i]
    if (nderivs >=1) {
      for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= j]
         Xj <- cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
       # Pjz <- modelfn(y[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
         Pjz <- storePjzsi[[i]][[j]]
         temscore[ind, ] <- temscore[ind, ] + Pygz[ind]*intwts[i]*
                            (as.vector(Pjz$dfy)/Pjz$fy)*Xj
      }
    }
}
Qstar <- cbind(1-Py,Py)
SptQstar  <- ptildes[1]*Qstar[,1] + ptildes[2]*Qstar[,2]

loglk <- sum(Bcounts*log(SptQstar))
if (nderivs >= 1){
    dQstar <- array(0,c(nrow(x[, ,1]),2,ntheta))
    dQstar[,1,] <- - temscore; dQstar[,2,] <- temscore
    dSptQstar <- ptildes[1]*dQstar[,1,] + ptildes[2]*dQstar[,2,]
    tempmat <- dSptQstar/SptQstar
    score <- apply(Bcounts*tempmat,2,sum)
}

if (nderivs >=2) {
    inf <-  t(tempmat) %*% (Bcounts*tempmat)
    correctmat <- matrix(0,nrow=ntheta,ncol=ntheta)
    for (i in 1:length(zval)){
         Pygz <- storePygz[,i]
         multvec <- (ptildes[2]-ptildes[1])*Bcounts*(Pygz/SptQstar)*intwts[i]
#        multvec <- Bcounts*(Pygz/Py)*intwts[i]
         for (j in 1:Jmax){
            ind <- (1:nn)[Jis >= j]
            Xj <- cbind(x[ind, ,j],sigma*rep(zval[i],length(ind)))
          # Pjz <- modelfn(y[ind,j],sigma*zval[i]+eta[ind,j],nderivs=nderivs,report="vals")
            Pjz <- storePjzsi[[i]][[j]]
            correctmat <- correctmat + t(Xj) %*% (multvec[ind]*as.vector(Pjz$d2fy)/Pjz$fy*Xj)
            if (j < Jmax) for (k in (j+1):Jmax){ ##***THINKING HERE FOR RAGGED ****
               indk <- (1:nn)[Jis >= k]
               Xk <- cbind(x[indk, ,k],sigma*rep(zval[i],length(indk)))
             # Pkz <- modelfn(y[indk,k],sigma*zval[i]+eta[indk,k],nderivs=nderivs,report="vals")
               Pkz <- storePjzsi[[i]][[k]]
               jrowuse <- (1:length(ind))[table(factor(indk,levels=ind))==1]
               temp <- t(Xj[jrowuse,]) %*% ((multvec[indk])*(as.vector(Pjz$dfy)/
                       Pjz$fy)[jrowuse]*(as.vector(Pkz$dfy)/Pkz$fy)*Xk)
               correctmat <- correctmat + temp + t(temp)
            }
         }
    }
    inf <- inf - correctmat
    inf[ntheta,ntheta] <- inf[ntheta,ntheta] - score[ntheta]
}

list(loglk=loglk,score=as.vector(score),inf=inf,Qstar=Qstar,dQstar=dQstar,SptQstar=SptQstar,
     dSptQstar=dSptQstar)
}
