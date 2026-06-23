##############################################################################	
MLdirectInf <- function(theta,nderivs=2,deltamat=NULL, modelfn,hmodelfn,x,y,xStrat,
                    Aposn,Acounts,Bposn,Bcounts,hvalue,Cmult,hxStrat,
                    off.set=matrix(0,dim(x)[1],dim(x)[3]),extra=NULL,
                    control.inner=mlefn.control.inner(...), ...){
##############################################################################
    nhvalue <- nrow(as.matrix(hvalue))     # allow for both matrices and vectors
    Cmult <- as.vector(Cmult)
    theta <- as.vector(theta)
    ntheta <- length(theta)
    lA <- length(Aposn)
    lB <- length(Bposn)
    n <- nrow(as.matrix(y))     # allow for both matrices and vectors
#   DIMENSION CHECKING
    if ( dim(x)[1] != n | dim(x)[2] != ntheta | max(Aposn,Bposn)>n | min(Aposn,Bposn)<1 |
       lA != length(Acounts) | lB != length(Bcounts) |
       length(Cmult) != nhvalue )
        stop("MLdirectInf: Data-dimension problems")

nxStrat <- max(xStrat)
nBposns <- rep(0,nxStrat)
for (i in 1:nxStrat) nBposns[i] <- length(Bposn[xStrat[Bposn]==i])

# The incoming deltamat (if non-null) is a starting value used in mlefn's 1st call to MLdirect only
# For subsequent calls we use the updated version from the last call to MLdirect that
# mlefn grabs from MLdirect on return and passes back

if (!is.null(extra$deltamat))   deltamat <- extra$deltamat
else if (is.null(deltamat)){ # calculate uniform starting values for deltas, put in deltamat
      deltamat <- matrix(NA,nrow=max(nBposns),ncol=nxStrat)
      for (i in 1:nxStrat){
         Bposni <- Bposn[xStrat[Bposn]==i]
         Cmulti <- Cmult[hxStrat==i]
         if (length(Bposni)<1 || length(Cmulti)<1 || sum(abs(Cmulti))==0 ) next()
         if (length(Bposni) < 3) stop(paste("MLdirectInf: Less than 3 distinct x-obs in xStrat=",i))
         delta <- Bcounts[xStrat[Bposn]==i]
         deltamat[,i][1:nBposns[i]] <- delta/sum(delta)
      }
}
if (any(dim(deltamat) != c(max(nBposns),nxStrat)))
             stop(paste("Dimension problems with deltamat in MLdirect, dim=",max(nBposns),nxStrat))

loglk <- 0
ntheta <- length(theta)
score <- if (nderivs >= 1) rep(0,ntheta)  else NULL
inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta) else  NULL

dimx <- dim(x)


for (i in 1:nxStrat){
     Aposni <- Aposn[xStrat[Aposn]==i]; Acountsi <- Acounts[xStrat[Aposn]==i]
     Bposni <- Bposn[xStrat[Bposn]==i]; Bcountsi <- Bcounts[xStrat[Bposn]==i]
     if (length(Aposni)<1 & length(Bposni)<1) next()
     hvaluei <- if (is.matrix(hvalue)) hvalue[hxStrat==i,] else hvalue[hxStrat==i]
     Cmulti <- Cmult[hxStrat==i]
     if(length(Cmulti)>0  && sum(abs(Cmulti))>0  && length(Bcountsi)>0  && sum(abs(Bcountsi))>0){
        if (length(hvaluei) < 2) stop(
               paste("MLdirectInf: Less than 2 distinct hvalue-obs in xStrat=",i))
        deltai <- deltamat[,i][1:nBposns[i]]
     }
     else  #  Do not wish to estimate deltai = pr(x | x-stratum=i)
          deltai <- Cmulti <- hvaluei <- Bposni <- Bcountsi <- numeric(0)

     w <- ML2directInf(theta,nderivs=nderivs, modelfn=modelfn,hmodelfn=hmodelfn, x=x,y=y,
                     Aposn=Aposni,Acounts=Acountsi,Bposn=Bposni,Bcounts=Bcountsi,
                     hvalue=hvaluei,Cmult=Cmulti,delta=deltai,
                     off.set=off.set, inxStrat=i, control.inner=control.inner, ...)

     if (!is.null(w$error))
         return(list(extra=list(deltamat=deltamat),error=1))

     if(length(Cmulti)>0) 
         deltamat[,i][1:nBposns[i]] <- w$delta  # update deltas for next call to MLdirectInf

     loglk <- loglk + w$loglk
     if (nderivs >= 1) score <- score + as.vector(w$score)
     if (nderivs >= 2) inf <- inf + w$inf
}
list(loglk=loglk,score=as.vector(score),inf=inf,extra=list(deltamat=deltamat))
}


##############################################################################	
ML2directInf <- function(theta,nderivs=2,modelfn,hmodelfn,x,y,
                        Aposn,Acounts,Bposn,Bcounts,hvalue,Cmult,delta,
                        off.set=matrix(0,dim(x)[1],dim(x)[3]), 
                        inxStrat, control.inner=mlefn.control.inner(...), ...){
##############################################################################

# inxStrat may be missing. It is only used for printing a diagnostic

    nhvalue <- nrow(as.matrix(hvalue))     # allow for both matrices and vectors
    theta <- as.vector(theta)
    ntheta <- length(theta)
    lA <- length(Aposn)
    lB <- length(Bposn)
    neta <- dim(x)[3]
    n <- nrow(as.matrix(y))     # allow for both matrices and vectors

#   DIMENSION CHECKING
    if ( dim(x)[1] != n || dim(x)[2] != ntheta || max(Aposn)>n || (lB>0 && max(Bposn)>n) ||
       min(Aposn)<1 || (lB>0 && min(Bposn)<1) ||
       length(delta) != lB || lA != length(Acounts) || lB != length(Bcounts) ||
       length(Cmult) != nhvalue )
        stop("ML2directInf: Data-dimension problems")

    eta <- matrix(0,dim(x)[1],neta)
    use <- unique(c(Aposn,Bposn))
    for (j in 1:neta) eta[use,j] <- x[use, ,j]%*%theta
    eta[use,] <- eta[use,] + off.set[use,]  # only calc eta values that will be used

#   Update the deltas's

    if (lB>0) {
        rho <- log(delta[-lB]/delta[lB])
        if (control.inner$guide=="auto") control.inner$guide <- "uphill"
        w <- mlefn(rho,rhodirectInf,eta=eta[Bposn, ,drop=FALSE],Bcounts=Bcounts,hvalue=hvalue,Cmult=Cmult,
                  hmodelfn=hmodelfn,control=control.inner, ...)
        if (control.inner$print.progress>= 1 && w$error != 0 && !missing(inxStrat))
                   cat("mlefn error",w$error, "from x-Stratum",inxStrat,"\n")
        rho <- w$theta
        delta <- c(exp(rho),1) / (1 + sum(exp(rho)))
    }
    yuse <- if (is.matrix(y)) y[Aposn, ,drop=FALSE] else y[Aposn]
    wmod <- modelfn(yuse,eta[Aposn, ,drop=FALSE],nderivs=nderivs, ...)
    loglk <- sum(Acounts*as.vector(wmod$logfy))
    if (lB>0) loglk <- loglk + w$loglk

    score <- if (nderivs >= 1) rep(0,ntheta)  else  NULL
    inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta)   else NULL

    if (nderivs >= 1) {
         if (lA > 0) for (j in 1:neta)  score <- score +
                           as.vector(t(x[Aposn, ,j]) %*% (Acounts*wmod$dlogfy[,j]))
         if (lB > 0) {
             indhvalue <- rep(1:nhvalue,lB)
             indeta <- rep(Bposn,nhvalue*rep(1,lB))
             inddelta <- rep(1:lB,nhvalue*rep(1,lB))
             hvalsend <- if (is.matrix(hvalue)) hvalue[indhvalue,]  else hvalue[indhvalue]
             wonly <- hmodelfn(hvalsend,eta[indeta, ,drop=FALSE],nderivs=nderivs, ...)
             Dvec <- apply(matrix(wonly$Qstar*delta[inddelta],nhvalue,lB),1,sum)
             wonly$dQstar <- wonly$dQstar*delta[inddelta]/Dvec[indhvalue]
             dftheta <- matrix(0,length(indhvalue),ntheta)
             for (k in 1:neta)  dftheta <- dftheta + wonly$dQstar[,k]*x[indeta, ,k]
             wonly$dQstar <- NULL   # RELEASE THE STORAGE
             score <- score + apply(Cmult[indhvalue]*dftheta,2,sum)
             }
         }
    if (nderivs >= 2){
        if (lA > 0){
           # 1st calculate the I(theta,theta) matrix
           for (j in 1:neta)
                 inf <- inf -
                      t(x[Aposn, ,j])%*%((Acounts*wmod$d2logfy[,j,j])*x[Aposn, ,j])
            if (neta > 1)
                 for (j in 1:(neta-1))  for (k in (j+1):neta){
                       temp1 <- t(x[Aposn, ,j])%*%((Acounts*wmod$d2logfy[,j,k])*x[Aposn, ,k])
                       inf <- inf - (temp1 + t(temp1))
                       }
            }
        if (lB > 0){
           Rmatlong <- wonly$Qstar*delta[inddelta]/Dvec[indhvalue]
           dfthetasm <- apply(array(dftheta,c(nhvalue,lB,ntheta)),c(1,3),sum)
#                Combines  "dftheta3 <- array(dftheta,c(nhvalue,lB,ntheta))"
#                and "dfthetasm <- apply(dftheta3,c(1,3),sum)" # Collapse dftheta over B dimension
           Irhotheta <-  - apply(array(Cmult[indhvalue]*dftheta,c(nhvalue,lB,ntheta)),c(2,3),sum) +
                          apply(array(
                           Cmult[indhvalue]*Rmatlong*dfthetasm[indhvalue,],c(nhvalue,lB,ntheta)),c(2,3),sum)
           Irhotheta <- Irhotheta[-lB, ,drop=FALSE]
           rm(dftheta, Rmatlong)    # RELEASE THE STORAGE

           # Find part of inf[theta,theta] from hvalue part of likelihood
           wonly$d2Qstar <- wonly$d2Qstar*delta[inddelta]/Dvec[indhvalue]
           infadd <- t(Cmult*dfthetasm) %*% dfthetasm

           for (j in 1:neta)
                 infadd <- infadd - t(Cmult[indhvalue]*wonly$d2Qstar[,j,j]*x[indeta, ,j]) %*% x[indeta, ,j]
           if (neta > 1)
                 for (j in 1:(neta-1))  for (k in (j+1):neta){
                       temp1 <- t(Cmult[indhvalue]*wonly$d2Qstar[,j,k]*x[indeta, ,j]) %*% x[indeta, ,k]
                       infadd <- infadd - (temp1 + t(temp1))
                       }

           soln <- qr.solve.cw(w$inf,Irhotheta,messg="ML2directInf:")  # w$inf is Irhorho
           if (soln$err != 0)
              return(list(loglk=loglk,score=as.vector(score),inf=NULL,delta=delta,error=1))
           inf <- inf + infadd - t(Irhotheta) %*% soln$soln
           }
        }
    list(loglk=loglk,score=as.vector(score),inf=inf,delta=delta)
    }



##############################################################################	
rhodirectInf <- function(rho,nderivs=2,eta,Bcounts,hvalue,Cmult,hmodelfn, ...){
##############################################################################	

    Cmult <- as.vector(Cmult)
    nB <- sum(Bcounts)
    N <- nB + sum(Cmult)
    nhvalue <- nrow(as.matrix(hvalue))
    nrho <- length(rho)
    lB <- nrow(eta)

#   DIMENSION CHECK
    if (lB != (nrho+1) | length(Bcounts) != lB | length(Cmult) != nhvalue)
                            stop("rhodirectInf: Input dimension problems")

    denom <- 1 + sum(exp(rho))
    delta <- c(exp(rho)/denom,1/denom)
    indhvalue <- rep(1:nhvalue,lB)
    indeta <- rep(1:lB,nhvalue*rep(1,lB))
    hvalsend <- if (is.matrix(hvalue)) hvalue[indhvalue,]  else hvalue[indhvalue]
    fdelta <- matrix(
              hmodelfn(hvalue=hvalsend,eta=eta[indeta, ,drop=FALSE],
                       nderivs=0, ...)$Qstar*delta[indeta],  nrow=nhvalue,ncol=lB)
    Dvec <- apply(fdelta,1,sum)

    loglk <- sum(Bcounts*log(delta)) + sum(Cmult*log(Dvec))

    NULL -> score -> inf -> Rmat
    if (nderivs >= 1){
            deltachk <- delta[-lB]
            Rmat <- fdelta[,-lB]/Dvec
            temp <- apply(Cmult*Rmat,2,sum)
            score <- Bcounts[-lB] + temp - N*deltachk 
             }
    if (nderivs >= 2)
          inf <- N*(diag(deltachk) - deltachk%*%t(deltachk)) - diag(temp) +
                      t(Cmult*Rmat) %*% Rmat
    list(loglk=loglk,score=as.vector(score),inf=inf)
}
