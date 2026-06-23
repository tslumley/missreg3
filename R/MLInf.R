##############################################################################	
MLInf <- function(theta,nderivs=2,ProspModInf,StratModInf,x,y,
                   Aposn,Acounts,Bposn,Bcounts,rmat,Qmat,xStrat=rep(1,dim(x)[1]),
                   extra=NULL,off.set=matrix(0,dim(x)[1],dim(x)[3]),
                    control.inner=mlefn.control.inner(...), ...){
##############################################################################
n <- length(y[,1])
maxxStrat <- max(xStrat)
if ((dim(x)[1] != n) || dim(x)[2] != length(theta) || (length(xStrat) != n) ||
       max(Aposn,Bposn)>n || min(Aposn,Bposn)<1 ||
       (length(Bposn)>0 && min(ncol(rmat),ncol(Qmat)) < maxxStrat))
    stop("MLInf: data dimension problems")

if (!is.null(extra$Qmat)) Qmat <- extra$Qmat
          # replace by updated version from last call to MLInf (passed back from mlefn)
loglk <- 0
ntheta <- length(theta)
score <- if (nderivs >= 1) rep(0,ntheta)  else NULL
inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta) else  NULL

for (i in 1:maxxStrat){
     Aposni <- Aposn[xStrat[Aposn]==i]
     Bposni <- Bposn[xStrat[Bposn]==i]
     if (length(Aposni)<1 & length(Bposni)<1) next()

    if (length(Bposni)>0){
        rvec <- rmat[,i][!is.na(rmat[,i])]
        Qs <- Qmat[,i][!is.na(Qmat[,i])]
        Qs[rvec != 0 & (Qs == 0 | Qs == 1)] <- 0.5
    }
    else rvec <- Qs <- numeric(0)

    if (length(Bposni)>0 && length(rvec) != length(Qs)) stop(paste("MLInf: In x-stratum",i,
                            "rmat[,i] and Qmat[,i] have diff non-NA lengths",
                             "(", length(rvec),"vs",length(Qs),")"))
    w <- ML2Inf(theta,nderivs=nderivs,ProspModInf=ProspModInf,StratModInf=StratModInf,x=x,y=as.matrix(y),
           Aposn=Aposni,Acounts=Acounts[xStrat[Aposn]==i],Bposn=Bposni, Bcounts=Bcounts[xStrat[Bposn]==i], 
           rvec=rvec,Qs=Qs,havexStrat=TRUE,inxStrat=i, off.set=off.set, control.inner=control.inner, ...)
    if (!is.null(w$error)) return(list(extra=list(Qmat=Qmat),error=1))  
 

    if (length(Bposni)>0) Qmat[,i][1:length(w$Qs)] <- w$Qs # capture updated Qs for next iteration
    loglk <- loglk + w$loglk
    if (nderivs >= 1) score <- score + as.vector(w$score)
    if (nderivs >= 2) inf <- inf + w$inf
}

list(loglk=loglk,score=as.vector(score),inf=inf,extra=list(Qmat=Qmat))
}

###########################################################################################	
ML2Inf <- function(theta,nderivs=2,ProspModInf,StratModInf,x,y,Aposn,Acounts,Bposn,Bcounts,
                   rvec,Qs,usage="thetaonly",thetaparts=0,paruse="auto",inxStrat,
                   off.set=matrix(0,dim(x)[1],dim(x)[3]), control.inner=mlefn.control.inner(...), ...){
###########################################################################################

# inxStrat may be missing. It is only used for printing a diagnostic
# Qs may be missing, if usage="combined" (Values of usage "thetaonly", "combined", "Qfixed")

# Set up and dimension checks
theta <- as.vector(theta)
if (usage == "combined"){  # using for both theta params and either xis or rhos
        # Check thetaparts
        if (length(thetaparts)!=2 | sum(thetaparts) != length(theta)) stop("ML2Inf: thetaparts error")
        thetain <- theta
        theta <- theta[1:thetaparts[1]]
        thetarest <- thetain[thetaparts[1]+(1:thetaparts[2])]
}
n <- length(y[,1])
nn <- dim(x)[1]
nyStrat <- length(rvec)
ntheta <- length(theta)
nQs <- if (usage !="combined") length(Qs) else nyStrat
lA <- length(Aposn)
lB <- length(Bposn)


# Data Checking
if (nn != n | max(Aposn,Bposn)>n | min(Aposn,Bposn)<1 | nQs != length(rvec) |
             lA != length(Acounts) | lB != length(Bcounts) )
       stop("ML2Inf: data dimension problems")

rplus <- sum(rvec)
nB <- maxabsA <- maxabsB <- 0
if (lA >= 1) maxabsA <- max(abs(Acounts))
if (lB >= 1){
         nB <- sum(Bcounts)
         maxabsB <- max(abs(Bcounts))
}

# Update the Q's if they are not being treated as fixed

if (maxabsB>0 & usage != "Qfixed"){
    nzeros <- length(rvec[rvec==0])
    if (paruse != "xis" & paruse != "rhos")  paruse <- if (nzeros > 1) "xis" else "rhos"
    if (control.inner$guide=="auto"){ #Set up hill-climbing direction for inner mlefn call
          if (all(rvec >= 0)) control.inner$guide <- "downhill"
          else if (all(rvec <= 0)) control.inner$guide <- "uphill"
          else control.inner$guide <- "no"  # Straight Newton Raphson
    }
    if (any(rvec != 0) & maxabsB>0){
        if (paruse=="xis"){   # xi parameters
              if (usage == "thetaonly"){
                  xi <- qlogis(Qs[rvec != 0])
                  w <- mlefn(xi,xiInf,theta4xiInf=theta,StratModInf=StratModInf,
                       x=x[Bposn, , ,drop=FALSE],Bcounts=Bcounts,off.set=off.set[Bposn,,drop=FALSE],
                        rvec=rvec, inxStrat=inxStrat, control=control.inner, ...)
                  if (control.inner$print.progress>= 1 && w$error != 0 && !missing(inxStrat))
                               cat("Error from x-Stratum",inxStrat,"\n")
                  xi <- w$theta
              }
              else {
                  xi <- thetarest
                  w <- xiInf(xi,nderivs,theta4xiInf=theta,StratModInf=StratModInf,x=x[Bposn, , ,drop=FALSE],
                        Bcounts=Bcounts,rvec=rvec, inxStrat=inxStrat, off.set=off.set[Bposn,,drop=FALSE], ...)
                  if (!is.null(w$error)) if (w$error != 0) return(list(Qs=Qs, error=1))
              }
              Qs[rvec != 0] <- plogis(xi)
              nxi <- length(xi)
         }
         else {   # 0 or 1 zero, use rho parameters
              if (usage == "thetaonly"){
                  rho <- log(Qs[-nQs]/Qs[nQs])
                  w <- mlefn(rho,rhoInf,theta4rhoInf=theta,StratModInf=StratModInf,
                         x=x[Bposn, , ,drop=FALSE],Bcounts=Bcounts,off.set=off.set[Bposn,,drop=FALSE],
                         rvec=rvec, inxStrat=inxStrat, control=control.inner, ...)
                  if (control.inner$print.progress>= 1 && w$error != 0 && !missing(inxStrat)) 
                         cat("Error from x-Stratum",inxStrat,"\n")
                  denom <- 1 + sum(exp(w$theta))
                  Qs <- c(exp(w$theta)/denom,1/denom)
              }
              else {
                  rho <- thetarest
                  w <- rhoInf(rho,nderivs,theta4rhoInf=theta,StratModInf=StratModInf,
                           x=x[Bposn, , ,drop=FALSE],Bcounts=Bcounts,rvec=rvec,
                           inxStrat=inxStrat, off.set=off.set[Bposn,,drop=FALSE], ...)
                  if (!is.null(w$error)) if (w$error != 0) return(list(Qs=Qs, error=1))
                  denom <- 1 + sum(exp(rho))
                  Qs <- c(exp(rho)/denom,1/denom)
              }
         }
    }
    if (any(rvec != 0)) lenparams <- if (paruse=="xis") nxi  else (nQs-1)
}

# Calculate the ptildes's
if (maxabsB>0){
    ptildes <- (nB + rplus)*rep(1,nyStrat)    #  Qs[r] may be 0 or NA where rvec[r]=0
    if (any(rvec != 0)) ptildes[rvec != 0] <- ptildes[rvec != 0] - rvec[rvec != 0]/Qs[rvec != 0]
}

# Calculate the loglikelihood contributions
    loglk <- score <- inf <- NULL

    if (maxabsA>0) prospmod <- ProspModInf(
               theta=theta,nderivs=nderivs, y=y[Aposn,,drop=FALSE], x=x[Aposn, , ,drop=FALSE],
               wts=Acounts, inxStrat=inxStrat, off.set=off.set[Aposn,,drop=FALSE], ...)
    if (!is.null(prospmod$error)) 
           if (prospmod$error != 0)  return(list(Qs=Qs, error=2))

    if (maxabsB>0 && any(rvec != 0))  {
          stratmod <- StratModInf(
                    theta=theta,nderivs=nderivs,ptildes=ptildes,x=x[Bposn, , ,drop=FALSE],
                           Bcounts=Bcounts, inxStrat=inxStrat, off.set=off.set[Bposn,,drop=FALSE], ...)
          if (!is.null(stratmod$error)) 
                    if (stratmod$error != 0)  return(list(Qs=Qs, error=3))  } 

    loglk <- 0
    if (maxabsA>0) loglk <- loglk + prospmod$loglk
    if (maxabsB>0 && any(rvec != 0)) loglk <- loglk - stratmod$loglk +
                                       sum(rvec[rvec != 0]*log(Qs[rvec != 0]))

    if (nderivs >= 1){
         score <- rep(0,ntheta)
         if (maxabsA>0) score <- prospmod$score
         if (maxabsB>0 && any(rvec != 0)) score <- score - stratmod$score
         if (maxabsB>0 && usage == "combined" && any(rvec != 0)) score<- c(score,w$score)
    }
    if (nderivs >= 2){ # Initialization
         inf <-  matrix(0,ntheta,ntheta)
         if (maxabsA>0) inf <- inf + prospmod$inf
         if (maxabsB>0 && any(rvec != 0)) inf <- inf - stratmod$inf
         if (maxabsB>0 && any(rvec != 0) && usage != "Qfixed"){
              nonzero <- (1:nyStrat)[rvec != 0]
              Ithetaxi <- matrix(0,ntheta,lenparams)
              if (usage != "Qfixed") {
                   if (paruse=="xis"){   # xi parameters
                       avec <- rvec[nonzero]*exp(-xi)
                       for (r in 1:nxi) Ithetaxi[,r] <- avec[r]*apply(Bcounts*(stratmod$dQstar[,nonzero[r],] -
                            stratmod$Qstar[,nonzero[r]]*stratmod$dSptQstar/stratmod$SptQstar)/stratmod$SptQstar,2,sum)
                   }
                   else {
                       for (r in 1:(nQs-1)){
                           ar <- (rvec[r]/Qs[r])*stratmod$Qstar[,r] - (nB+rplus)*Qs[r]
                           Ithetaxi[,r] <- apply(Bcounts*((rvec[r]/Qs[r])*stratmod$dQstar[,r,] -
                                   ar*stratmod$dSptQstar/stratmod$SptQstar)/stratmod$SptQstar,2,sum)
                       }
                   }
              }
          }
     }
 if (maxabsB>0 && nderivs >= 2 && any(rvec != 0)){
       if (usage == "thetaonly") inf <- inf - Ithetaxi%*%solve(w$inf,t(Ithetaxi)) # w$inf is Ixixi
       else if (usage == "combined") inf <- rbind(cbind(inf,Ithetaxi),cbind(t(Ithetaxi),w$inf))
      #else if (usage == "Qfixed") leave inf as it is
 }
 
 list(loglk=loglk,score=as.vector(score),inf=inf,Qs=Qs)
 }

##############################################################################	
MEtaProspModInf <- function(theta,nderivs=2,y,x,wts=1,modelfn,off.set=0, ...){
##############################################################################

# inf function for models containing several linear predictors

    theta <- as.vector(theta)
    ntheta <- length(theta)
    neta <- dim(x)[3]
    nn <- dim(x)[1]
    eta <- NULL
    for (j in 1:neta){
         xj <- x[, ,j]
         # if (nn==1) xj <- matrix(xj,nrow=1,ncol=ntheta) # dim-collapse problem otherwise
         if (nn==1 | dim(x)[2]==1) xj <- matrix(xj,nrow=nn,ncol=ntheta) # dim-collapse problem otherwise
         eta <- cbind(eta, xj%*%theta)
    }
    eta <- eta + off.set

    error <- 0
    loglk <- score <- inf <- NULL

    w <- modelfn(y,eta,nderivs=nderivs, ...)
    if (!is.null(w$error)) if (w$error != 0) {
     	error <- 1
        return(list(loglk=loglk, score=score, inf=inf, error=error)) }

    loglk <- sum(wts*as.vector(w$logfy))
    score <- if (nderivs >= 1) rep(0,ntheta)  else  NULL
    inf <- if (nderivs >= 2) matrix(0,ntheta,ntheta)   else NULL
    if (nderivs >= 1)  for (j in 1:neta){ # Loop over linear Predictors
         xj <- x[, ,j]
         if (nn==1) xj <- matrix(xj,nrow=1,ncol=ntheta) # dim-collapse problem otherwise
         score <- score + as.vector(t(xj) %*% (wts*w$dlogfy[,j]))
         if (nderivs >= 2){
              inf <-  inf -  t(x[, ,j])%*%((wts*w$d2logfy[,j,j])*x[, ,j])
              if (j < neta)  for (k in (j+1):neta){
                   xk <- x[, ,k]
                   if (nn==1) xk <- matrix(xk,nrow=1,ncol=ntheta) # dim-collapse problem otherwise
                   temp <- t(xj)%*%((wts*w$d2logfy[,j,k])*xk)
                   inf <- inf - (temp + t(temp))
              }
         }
    }
    list(loglk=loglk, score=as.vector(score), inf=inf, error=error)
}

##############################################################################
MEtaStratModInf <- function(theta,nderivs=2,x,Bcounts,ptildes,stratfn,off.set=0, ...){
##############################################################################

# inf function for models containing several linear predictors

score <- inf <- Qstar <- dQstar <- SptQstar <- dSptQstar <- NULL

theta <- as.vector(theta)
ntheta <- length(theta)
neta <- dim(x)[3]
nn <- dim(x)[1]

ptildes <- as.vector(ptildes)
nyStrat <- length(ptildes)

eta <- NULL
for (j in 1:neta){
     xj <- x[, ,j]
     if (nn==1 | dim(x)[2]==1) xj <- matrix(xj,nrow=nn,ncol=ntheta) # dim-collapse problem otherwise
     eta <- cbind(eta, xj%*%theta)
}
eta <- eta + off.set

w <- stratfn(eta,nderivs=nderivs,ptildes=ptildes, ...)
SptQstar <- apply(ptildes*t(w$Qstar) ,2,sum)
loglk <- 0
loglk <- loglk + sum(Bcounts*as.vector(log(SptQstar)))

inf <-  if (nderivs >= 2)  matrix(0,ntheta,ntheta)        else  NULL

if (nderivs >= 1) {
    score <- rep(0,ntheta)
    dQstar <- array(0,c(nn,nyStrat,ntheta))
    dSptQstar <- matrix(0,nrow=nn,ncol=ntheta)
    dSptQstarXX <- apply(ptildes*aperm(w$dQstar,c(2,1,3)),c(2,3),sum)
    for (j in 1:neta){# Loop over linear predictors
         xj <- x[, ,j]
        if (nn==1) xj <- matrix(xj,1,ntheta)  # dim-collapse problem otherwise
        for (g in 1:nyStrat) dQstar[,g,] <- dQstar[,g,] + w$dQstar[,g,j]*xj

if (nderivs >= 2){# Calc contrib of jth linear predictor to information matrix
         inf <- inf   - t(xj)%*% ((Bcounts*((w$SptQstar2[,j,j] - dSptQstarXX[,j]^2/SptQstar)/SptQstar)) * xj)
         if (j < neta)  for (k in (j+1):neta){
              xk <- x[, ,k]
              if (nn==1) xk <- matrix(xk,nrow=1,ncol=ntheta) # dim-collapse problem otherwise
              temp2 <- t(xj)%*%((Bcounts*
                          ((w$SptQstar2[,j,k]-dSptQstarXX[,j]*dSptQstarXX[,k]/SptQstar)/SptQstar))*xk)
         inf <- inf  - (temp2 + t(temp2))
         }
}

    }
#    for (g in 1:nyStrat) dSptQstar <- dSptQstar + ptildes[g]*dQstar[,g,]
    dSptQstar <- apply(ptildes*aperm(dQstar,c(2,1,3)),c(2,3),sum)
    score <- score +  apply((Bcounts*dSptQstar)/SptQstar,2,sum)
}


list(loglk=loglk,score=as.vector(score),inf=inf,Qstar=w$Qstar,dQstar=dQstar,SptQstar=SptQstar,dSptQstar=dSptQstar)
}

###########################################################################################
MEtaStratModInf.spml2<-function(theta,nderivs=2,xlist,Bcounts,ptildes,off.set=0,inxStrat, ...)
###########################################################################################
{
# a special StratModInf function for Spml2 method (bivbin2stg data)
# inf function for models containing several linear predictors

score <- inf <- Qstar <- dQstar <- SptQstar <- dSptQstar <- NULL

if (length(xlist) != 4) stop("xlist has wrong dimension!")
else if (is.null(names(xlist))) names(xlist) <- c("orig","1","0","xstrata")

strata <- xlist[[match("xstrata",names(xlist))]]
X  <- (xlist[[match("orig",names(xlist))]])[strata==inxStrat,,,drop=FALSE]
X1 <- (xlist[[match("1",names(xlist))]])[strata==inxStrat,,,drop=FALSE]
X0 <- (xlist[[match("0",names(xlist))]])[strata==inxStrat,,,drop=FALSE] 

theta <- as.vector(theta)
ntheta <- length(theta)
neta <- dim(X)[3]
nn <- dim(X)[1]

ptildes <- as.vector(ptildes)
if (length(ptildes) != 2) stop("length(ptildes) != 2")
nyStrat <- length(ptildes)

etarray <- array(0, c(nn,neta,3))
for (i in 1:3)
  for (j in 1:neta){
     xarray <- xlist[[i]][strata==inxStrat,,,drop=FALSE]
     xj <- xarray[, ,j]
     if (is.vector(xj)) xj <- matrix(xj,nrow=nn,ncol=ntheta) # dim-collapse problem otherwise
     etarray[,j,i] <- xj %*% theta
  }

offset <- array(off.set, c(dim(off.set)[1],dim(off.set)[2],3))
etarray <- etarray + offset

y11 <- cbind(rep(1,nn),rep(1,nn))
y10 <- cbind(rep(1,nn),rep(0,nn))
y01 <- cbind(rep(0,nn),rep(1,nn))
y00 <- cbind(rep(0,nn),rep(0,nn))
w11 <- spml2(y11,eta=etarray[,,2], nderivs=2)
w10 <- spml2(y10,eta=etarray[,,3], nderivs=2)
w01 <- spml2(y01,eta=etarray[,,2], nderivs=2)
w00 <- spml2(y00,eta=etarray[,,3], nderivs=2)

py1 <- w11$fy+w10$fy
qy1 <- w01$fy+w00$fy
Qstar <- cbind(py1, qy1)
SptQstar <- apply(ptildes*t(Qstar), 2, sum)

if (nderivs >= 2) {
  SptQstar1 <- array(0,c(dim(w11$dfy), 2))
  SptQstar2 <- array(0,c(dim(w11$d2fy), 2))
  SptQstar1[,,1] <- ptildes[1]*w11$dfy + ptildes[2]*w01$dfy
  SptQstar1[,,2] <- ptildes[1]*w10$dfy + ptildes[2]*w00$dfy
  SptQstar2[,,,1] <- ptildes[1]*w11$d2fy + ptildes[2]*w01$d2fy
  SptQstar2[,,,2] <- ptildes[1]*w10$d2fy + ptildes[2]*w00$d2fy
  
  ## Find out cols including Y2 ##
  y2col <- NULL
  for (j in 2:ntheta) {
    if (all(X1[,j,1]==0) & all(X0[,j,1]==0)) next()
      else  if (any((X1[,j,1]-X0[,j,1])!=0)) y2col <- c(y2col,j)
  }
}

loglk <- sum(Bcounts*log(SptQstar))

score <- if (nderivs >= 1) rep(0,ntheta)   else NULL
inf <-  if (nderivs >= 2)  matrix(0,ntheta,ntheta)  else  NULL

if (nderivs >= 1) {
   dQstar <- array(0,c(nn,nyStrat,ntheta))
   dSptQstar <- matrix(0,nrow=nn,ncol=ntheta)
      
   for (j in 1:neta) {# Loop over linear predictors
     x1j <- x1j2 <- matrix(X1[,,j],nn,ntheta)
     x0j <- x0j2 <- matrix(X0[,,j],nn,ntheta)
     
     dQstar[,1,] <- dQstar[,1,] + w11$dfy[,j]*x1j + w10$dfy[,j]*x0j 
     dQstar[,2,] <- dQstar[,2,] + w01$dfy[,j]*x1j + w00$dfy[,j]*x0j 
     
     if (nderivs >= 2){# Calc contrib of jth linear predictor to information matrix
       ## for y2==1 part (x1) ##
       inf1 <-  t(x1j) %*%((Bcounts*((SptQstar2[,j,j,1] - 
               (SptQstar1[,j,1]^2+SptQstar1[,j,1]*SptQstar1[,j,2])/
               SptQstar)/SptQstar)) * x1j)    
       if (!is.null(y2col)) {
          cross <-  t(x1j) %*%((Bcounts*((SptQstar2[,j,j,1] - 
                 SptQstar1[,j,1]^2/SptQstar)/SptQstar)) * x1j) 
          inf1[y2col,y2col] <- cross[y2col,y2col]  
       } 
       inf <- inf - inf1
                   
       if (j < neta)  for (k in (j+1):neta){
          x1k <- matrix(X1[,,k],nn,ntheta)
          temp2 <- t(x1j)%*%((Bcounts*((SptQstar2[,j,k,1]-
                   (SptQstar1[,j,1]*SptQstar1[,k,1]+
                    SptQstar1[,j,1]*SptQstar1[,k,2])/SptQstar)/SptQstar))*x1k)
          inf <- inf - (temp2 + t(temp2))
       }

       ## for y2==0 part (x0) ##
       inf2 <-  t(x0j) %*%((Bcounts*((SptQstar2[,j,j,2] - 
              (SptQstar1[,j,2]^2+SptQstar1[,j,2]*SptQstar1[,j,1])/
              SptQstar)/SptQstar)) * x0j) 
       if (!is.null(y2col))  {      
           cross <-  t(x0j) %*%((Bcounts*((SptQstar2[,j,j,2] - 
                  SptQstar1[,j,2]^2/SptQstar)/SptQstar)) * x0j) 
           inf2[y2col,y2col] <- cross[y2col,y2col]
       }
       inf <- inf - inf2
                       
       if (j < neta)  for (k in (j+1):neta){
          x0k <- matrix(X0[,,k],nn,ntheta)
          temp2 <- t(x0j)%*%((Bcounts*((SptQstar2[,j,k,2]-
                   (SptQstar1[,j,2]*SptQstar1[,k,2]+
                    SptQstar1[,j,2]*SptQstar1[,k,1])/SptQstar)/SptQstar))*x0k)
          inf <- inf - (temp2 + t(temp2))
       }
     }
   }
   dSptQstar <- apply(ptildes*aperm(dQstar,c(2,1,3)),c(2,3),sum)
   score <- score +  apply((Bcounts*dSptQstar)/SptQstar,2,sum)
}

list(loglk=loglk,score=as.vector(score),inf=inf,Qstar=Qstar,dQstar=dQstar,
     SptQstar=SptQstar,dSptQstar=dSptQstar)
}

########################################################################################
MEtaStratModInf.spml2locsc <- function(theta,nderivs=2, x, Bcounts, ptildes, off.set=0, 
                                  inxStrat=1, X4, y2posiny1, nthetaOrig, strata, ...){
########################################################################################

score <- inf <- Qstar <-  dQstar <- SptQstar <- dSptQstar <- NULL

theta <- as.vector(theta)
ntheta <- length(theta)
neta <- dim(x)[3]
nn <- dim(x)[1]
nt <- nthetaOrig

ptildes <- as.vector(ptildes)
nyStrat <- length(ptildes)
if (nyStrat != 2)  stop("length(ptildes) != 2")

eta <- NULL
for (j in 1:neta){
     xj <- x[, ,j]
     if (nn==1 | dim(x)[2]==1) xj <- matrix(xj,nrow=nn,ncol=ntheta) # dim-collapse problem otherwise
     eta <- cbind(eta, xj%*%theta)
}
eta <- eta + off.set

## this variable is added-in and cannot be subset by xstrata in 'MLInf'
X4 <- X4[strata==inxStrat,,drop=FALSE]

## Currently use Normal approximation for eps and f(eps), same as for 'rclusbin' 
eps <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
neps <- length(eps)
feps <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
feps[3] <- 8/15

P1y <- rep(0,nn)
dP1y <- if (nderivs >=1) matrix(0,nrow=nn,ncol=ntheta) else NULL
X1 <- x[,,1]
X2 <- x[,,2]; X3 <- x[,,3]
y2term <- as.vector(X4 %*% theta[y2posiny1])

for (i in 1:neps){
    y2i <- eta[,2]+exp(eta[,3])*eps[i]
    X1[,y2posiny1] <- y2i*X4
    eta1 <- as.vector(X1 %*% theta + off.set[,1])
 
    P1z <- binlogistic(rep(1,nn), eta1, nderivs=nderivs, report="vals")
    P1y <- P1y + P1z$fy*feps[i]

    if (nderivs >=1) {     
       #dP1y <- dP1y + feps[i]*(as.vector(P1z$dfy)*X1) + feps[i]*y2term*(as.vector(P1z$dfy)*X2) +
       #        feps[i]*eps[i]*exp(eta[,3])*y2term*(as.vector(P1z$dfy)*X3)
       dP1y <- dP1y + feps[i]*as.vector(P1z$dfy)*(X1 + y2term*X2 +
               eps[i]*exp(eta[,3])*y2term*X3)
    }
}
Qstar <- cbind(P1y, 1-P1y)
SptQstar  <- (ptildes[1]-ptildes[2])*P1y + ptildes[2]  # = ptildes[1]*P1y + ptildes[2]*(1-P1y)
if (nderivs >=1) dSptQstar <- (ptildes[1]-ptildes[2])*dP1y

loglk <- sum(Bcounts*log(SptQstar))
if (nderivs >= 1){
    score <- apply(Bcounts*dSptQstar/SptQstar,2,sum)
    dQstar <- array(0,c(nn,nyStrat,ntheta)); dQstar[,1,] <- dP1y ;  dQstar[,2,] <- - dP1y
}
if (nderivs >=2) {
    tempmat <- dSptQstar/SptQstar
    inf <- t(tempmat) %*% (Bcounts*tempmat)

    correctmat <- matrix(0,nrow=ntheta,ncol=ntheta)
    X1 <- x[,,1]
    for (i in 1:neps){
         y2i <- eta[,2]+exp(eta[,3])*eps[i]
         X1[,y2posiny1] <- X4i <- y2i*X4
         eta1 <- as.vector(X1 %*% theta + off.set[,1])
         P1z <- binlogistic(rep(1,nn), eta1, nderivs=nderivs, report="vals")
         allpart <- Bcounts*feps[i]*(ptildes[1]-ptildes[2])/SptQstar

         correctmat <- correctmat + t(X1) %*% (allpart*as.vector(P1z$d2fy)*X1) +
                t(X1) %*% (allpart*as.vector(P1z$d2fy)*y2term*X2) +
                t(X1) %*% (allpart*as.vector(P1z$d2fy)*eps[i]*exp(eta[,3])*y2term*X3) +
                t(X2) %*% (allpart*as.vector(P1z$d2fy)*y2term^2*X2) +
                t(X2) %*% (allpart*as.vector(P1z$d2fy)*eps[i]*exp(eta[,3])*y2term^2*X3)+
                t(X3) %*% (allpart*as.vector(P1z$d2fy)*(eps[i]*exp(eta[,3])*y2term)^2*X3)+
                t(X3) %*% (allpart*as.vector(P1z$dfy)*eps[i]*exp(eta[,3])*y2term*X3)

         correctmat[y2posiny1,] <- correctmat[y2posiny1,] + 
                t(X4) %*% (allpart*as.vector(P1z$dfy)*X2) + 
                t(X4) %*% (allpart*as.vector(P1z$dfy)*eps[i]*exp(eta[,3])*X3) 
    }
    
    ## add those symetric values for different blocks of X
    tempmat <- t(correctmat)
    tempmat[1:nt[1],1:nt[1]] <- 0
    tempmat[(nt[1]+1):(nt[1]+nt[2]),(nt[1]+1):(nt[1]+nt[2])] <- 0
    tempmat[(nt[1]+nt[2]+1):ntheta, (nt[1]+nt[2]+1):ntheta] <- 0
    correctmat <- correctmat + tempmat

    inf <- inf - correctmat
}
list(loglk=loglk,score=as.vector(score),inf=inf,Qstar=Qstar,dQstar=dQstar,
     SptQstar=SptQstar,dSptQstar=dSptQstar)
}

##############################################################################	
rhoInf <- function(rho,nderivs=2,theta4rhoInf,StratModInf,x,Bcounts,rvec, ...){
##############################################################################	

    rplus <- sum(rvec)
    nB <- sum(Bcounts)
    nrho <- length(rho)
    nn <- dim(x)[1]
#   DIMENSION CHECK
    if (length(Bcounts) != nn | nrho != (length(rvec)-1)) stop("rhoInf: Input dimension problems")

    denom <- 1 + sum(exp(rho))
    Qs <- c(exp(rho)/denom,1/denom)
    nptildes <- nrho+1
    ptildes <- as.vector(  (nB + rplus)*rep(1,nptildes) - rvec/Qs  )

    w <- StratModInf(theta=theta4rhoInf,nderivs=0,x=x,Bcounts=Bcounts,ptildes=ptildes, ...)
    if (!is.null(w$error)) if (w$error != 0) return(list(error=1))

    loglk <- sum(rvec[rvec != 0]*log(Qs[rvec != 0])) - w$loglk
    NULL -> score -> inf
    if (nderivs >= 1){
          aMat <- (w$Qstar %*% diag(nB+rplus-ptildes) -
                  (nB+rplus)*matrix(rep(1,nn),ncol=1)%*%matrix(Qs,nrow=1))/w$SptQstar
          score <- (rvec-(nB+rplus)*Qs - apply(Bcounts*aMat,2,sum))[1:nrho]
          }
    if (nderivs >= 2){
                inf <- (nB+rplus)*(1-2*sum(Bcounts/w$SptQstar))*(diag(Qs) - Qs%*%t(Qs)) -
                         diag(apply(Bcounts*aMat,2,sum)) - t(Bcounts*aMat) %*% aMat
                inf <- inf[1:nrho,1:nrho]
                }
    list(loglk=loglk,score=as.vector(score),inf=inf)
    }

##############################################################################	
xiInf <- function(xi,nderivs=2,theta4xiInf,StratModInf,x,Bcounts,rvec, ...){
##############################################################################	
# Inf function for the xi parameters -- to be called by ML2Inf

    xi <- as.vector(xi)
    nyStrat <- length(rvec)
    nonzero <- (1:nyStrat)[rvec != 0]
    nxi <- length(xi)
    nn <- dim(x)[1]

#   DIMENSION CHECK
    if (length(Bcounts) != nn | length(nonzero) != nxi) stop("xiInf: Input dimension problems")

    rplus <- sum(rvec)
    nB <- sum(Bcounts)
    Qs <- rep(1,nyStrat)
    Qs[nonzero] <- plogis(xi)  
    ptildes <- (nB + rplus) - rvec/Qs

    w <- StratModInf(theta=theta4xiInf,nderivs=0,x=x,Bcounts=Bcounts,ptildes=ptildes, ...)
    if (!is.null(w$error)) if (w$error != 0) return(list(error=1))

    loglk <- sum(rvec[nonzero]*log(Qs[nonzero])) -
                         sum(Bcounts*as.vector(log(w$SptQstar)))
    NULL -> score -> inf
    if (nderivs >= 1){
            avec <- rvec[nonzero]*exp(-xi)
            tempmat  <-  if (nxi > 1) (w$Qstar[,nonzero,drop=FALSE]/w$SptQstar) %*% diag(avec)
                           else (w$Qstar[,nonzero,drop=FALSE]/w$SptQstar) * (avec)
            tempvec  <- apply(Bcounts*tempmat,2,sum)
            score <- avec*Qs[nonzero] - tempvec
            }
    if (nderivs >= 2){
           tempdiag <- if (nxi >1) diag(rvec[nonzero]*Qs[nonzero]*(1-Qs[nonzero]) - tempvec)
                        else (rvec[nonzero]*Qs[nonzero]*(1-Qs[nonzero]) - tempvec)
           inf <- tempdiag - t(Bcounts*tempmat) %*% tempmat
           }
    list(loglk=loglk,score=as.vector(score),inf=inf)
    }
