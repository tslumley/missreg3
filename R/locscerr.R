#*****************************************
#************  locscerr.R  ***************
#*****************************************

########################################################################################	
locscale <- function(y,eta,nderivs=2, errdist, report="logs", ...){
########################################################################################	
    NULL -> logfy -> dlogfy -> d2logfy -> fy -> dfy -> d2fy
    y <- as.vector(y)

    sigma <- exp(eta[,2])
    eps <- (y - eta[,1])/sigma
    w <- errdist(eps,nderivs,report=report, ...)

    if (report=="logs"){
       logfy <- w$logfeps - eta[,2]
       if (nderivs >= 1) dlogfy <-  - cbind(w$dlogfeps/sigma, 1 + eps*w$dlogfeps)
       if (nderivs >= 2){
           d2logfy <- array(0,c(nrow(eta),2,2))
           d2logfy[,1,1] <- w$d2logfeps/sigma^2
           d2logfy[,1,2] <- d2logfy[,2,1] <- (w$dlogfeps + eps*w$d2logfeps)/sigma
           d2logfy[,2,2] <- eps*w$dlogfeps + eps^2*w$d2logfeps
           }
    }
    if (report != "logs"){
       fy <- w$feps/sigma
       if (nderivs >= 1) dfy <-  - cbind(w$dfeps/sigma, w$feps + eps*w$dfeps)/sigma
       if (nderivs >= 2){
           d2fy <- array(0,c(nrow(eta),2,2))
           d2fy[,1,1] <- w$d2feps/sigma^3
           d2fy[,1,2] <- d2fy[,2,1] <- (2*w$dfeps + eps*w$d2feps)/sigma^2
           d2fy[,2,2] <- (w$feps + 3*eps*w$dfeps + eps^2*w$d2feps)/sigma
           }
    }
    list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy,fy=fy,dfy=dfy,d2fy=d2fy,error=NULL)
}


#################################################################################################
locscale.log <- function(theta, nderivs=2, modelfn=locscale, errdist=errdist, y=y, x=x, wts=w, ...) 
#################################################################################################
{ # for the function locsc2stgwtd 

mx <- rbind(x[,,1],x[,,2])
eta <- matrix(mx %*% as.vector(theta),ncol=2)	
z <- modelfn(y,eta,nderivs=2,errdist=errdist, ...)

loglk <- sum(wts*as.vector(z$logfy))
score <- inf <- NULL

if (nderivs >=1) score <- t(mx) %*% (wts*as.vector(z$dlogfy)) 
if (nderivs >=2) {
	mx1 <- rbind(x[,,1],x[,,2],x[,,1],x[,,2]) #match z$d2logfy sequence of 2nd dimention
	mx2 <- rbind(x[,,1],x[,,1],x[,,2],x[,,2]) #match z$d2logfy sequence of 3nd dimention
	d2logfy <- c(-z$d2logfy[,1,1],-z$d2logfy[,2,1],-z$d2logfy[,1,2],-z$d2logfy[,2,2])
   	inf <- t(mx1) %*% (wts*d2logfy*mx2)  
}

list(loglk=loglk,score=as.vector(score),inf=inf)
}


########################################################################################	
locscstrat1 <- function(eta,nderivs=2,errdistcdf,ptildes,yCuts,havexStrat=FALSE,inxStrat=1, ...)
########################################################################################	
{
    if (havexStrat) {
        yCuts <- yCuts[,inxStrat]
        yCuts <- yCuts[!is.na(yCuts)]
        if (length(yCuts) < 3) stop(
            cat("Illegal yCuts column for x-stratum",inxStrat,"Column is","\n",
                 yCuts,"\n"))
    }

    NULL -> dQstar -> SptQstar2
    n <- nrow(eta)
    nyStrat <- length(ptildes)
    Qstar <- matrix(0,n,nyStrat)

    if (nderivs >= 1) dQstar <- array(0,c(n,nyStrat,2))
    if (nderivs >= 2) SptQstar2 <- array(0,c(n,2,2))

    sigma <- exp(eta[,2])

    epsleft <- (yCuts[1] - eta[,1])/sigma
    left <- errdistcdf(epsleft,nderivs=nderivs)
    for (g in 1:(length(yCuts)-1)){
       epsright <- (yCuts[g+1] - eta[,1])/sigma
       right <- errdistcdf(epsright,nderivs=nderivs)
       Qstar[,g] <- right$F -left$F
       if (nderivs >= 1){
             dQstar[,g,1] <- -(right$dF -left$dF)/sigma
             dQstar[,g,2] <- - (epsright*right$dF - epsleft*left$dF)
             }
       if (nderivs >= 2){
             SptQstar2[,1,1] <- SptQstar2[,1,1] + ptildes[g]*(right$d2F -left$d2F)/sigma^2
             SptQstar2[,1,2] <- SptQstar2[,1,2] +  ptildes[g]*(right$dF - left$dF +
                             epsright*right$d2F - epsleft*left$d2F)/sigma
             SptQstar2[,2,2] <- SptQstar2[,2,2] +  ptildes[g]*(epsright*right$dF - epsleft*left$dF +
                             epsright^2*right$d2F - epsleft^2*left$d2F)
             }
       left <- right
       epsleft <- epsright
       }
       if (nderivs >= 2)   SptQstar2[,2,1] <- SptQstar2[,1,2]

    list(Qstar=Qstar,dQstar=dQstar,SptQstar2=SptQstar2,error=NULL)
}

########################################################################################	
hylocscale <- function(hvalue,eta, nderivs=2, errdist,  ...){
########################################################################################

# hmodel function for direct method with location and scale models where h(y,x) = y
	
    NULL -> Qstar -> dQstar -> d2Qstar
    hvalue <- as.vector(hvalue)

    sigma <- exp(eta[,2])
    eps <- (hvalue - eta[,1])/sigma
    w <- errdist(eps,nderivs,report="dist", ...)

       Qstar <- w$feps/sigma
       if (nderivs >= 1) dQstar <-  - cbind(w$dfeps/sigma, w$feps + eps*w$dfeps)/sigma
       if (nderivs >= 2){
           d2Qstar <- array(0,c(nrow(eta),2,2))
           d2Qstar[,1,1] <- w$d2feps/sigma^3
           d2Qstar[,1,2] <- d2Qstar[,2,1] <- (2*w$dfeps + eps*w$d2feps)/sigma^2
           d2Qstar[,2,2] <- (w$feps + 3*eps*w$dfeps + eps^2*w$d2feps)/sigma
           }
    list(Qstar=Qstar,dQstar=dQstar,d2Qstar=d2Qstar,error=NULL)
}


########################################################################################	
logisterr <- function(eps,nderivs=2, report="logs", ...) {
########################################################################################	
    NULL -> logfeps -> dlogfeps -> d2logfeps -> feps -> dfeps -> d2feps

    if (report == "logs"){
        exp.eps <- exp(eps)
        logfeps <- eps - 2*log(1+exp.eps)
        if (nderivs >= 1)  dlogfeps <- 1-2*plogis(eps)
        if (nderivs >= 2)  d2logfeps <- -2*dlogis(eps)
        }
    if (report != "logs"){
        feps <- dlogis(eps)
        if (nderivs >= 1){
              p <- plogis(eps)
              q <- plogis(-eps)
              dfeps <- feps*(q-p)
              }
        if (nderivs >= 2)  d2feps <- feps*(1-6*feps)
        }

   
    list(logfeps=logfeps,dlogfeps=dlogfeps,d2logfeps=d2logfeps,
                   feps=feps,dfeps=dfeps,d2feps=d2feps,error=NULL)
}

########################################################################################	
logistcdf <- function(eps,nderivs=2, ...) {
########################################################################################	
    F <- plogis(eps)
    NULL -> dF -> d2F
    if (nderivs >= 1)  dF <- dlogis(eps)
    if (nderivs >= 2)  d2F <- dF*(1-2*F)

    list(F=F,dF=dF,d2F=d2F,error=NULL)
}

########################################################################################	
normerr <- function(eps,nderivs=2, report="logs", ...) {
########################################################################################	
    NULL -> logfeps -> dlogfeps -> d2logfeps -> feps -> dfeps -> d2feps

    if (report == "logs"){
        logfeps <- log(1/sqrt(2*pi)) - eps^2/2
        if (nderivs >= 1)  dlogfeps <- -eps
        if (nderivs >= 2)  d2logfeps <- -1
        }
    if (report != "logs"){
        feps <- dnorm(eps)
        if (nderivs >= 1)
              dfeps <- -eps*feps
        if (nderivs >= 2)  d2feps <- -feps - eps*dfeps
        }
   
    list(logfeps=logfeps,dlogfeps=dlogfeps,d2logfeps=d2logfeps,
                   feps=feps,dfeps=dfeps,d2feps=d2feps,error=NULL)
}
########################################################################################	
normcdf <- function(eps,nderivs=2, ...) {
########################################################################################	
    F <- pnorm(eps)
    NULL -> dF -> d2F
    if (nderivs >= 1)  dF <- dnorm(eps)
    if (nderivs >= 2)  d2F <- - eps*dF

    list(F=F,dF=dF,d2F=d2F,error=NULL)
}

########################################################################################	
terr <- function(eps,nderivs=2, errmodpars=4, report="logs", ...) {
########################################################################################
    NULL -> logfeps -> dlogfeps -> d2logfeps -> feps -> dfeps -> d2feps
    dfs <- errmodpars[1]
    if (report == "logs"){
        logfeps <- lgamma((dfs+1)/2)-log(sqrt(dfs*pi))-lgamma(dfs/2)-((dfs+1)/2)*log(1+(eps^2)/dfs)
        if (nderivs >= 1)  dlogfeps <- -(eps*(dfs+1))/(dfs+eps^2)
        if (nderivs >= 2)  d2logfeps <- (-(dfs+1)/((dfs+eps^2)^2))*(dfs-eps^2)
        }
    if (report != "logs"){
        feps <- dt(eps,dfs)
        if (nderivs >= 1)
              dfeps <- feps*((-eps*(dfs+1))/(dfs+eps^2))
        if (nderivs >= 2)  d2feps <- dfeps*((eps^2*(-2-dfs)+dfs)/((dfs+eps^2)*eps))
        }
   
    list(logfeps=logfeps,dlogfeps=dlogfeps,d2logfeps=d2logfeps,
                   feps=feps,dfeps=dfeps,d2feps=d2feps,error=NULL)
}

########################################################################################	
tcdf <- function(eps,nderivs=2, errmodpars=4, ...) {
########################################################################################	
    dfs <- errmodpars[1]
    F <- pt(eps,dfs)
    NULL -> dF -> d2F
    if (nderivs >= 1)  dF <- dt(eps,dfs)
    if (nderivs >= 2)  d2F <- dF*((-eps*(dfs+1))/(dfs+eps^2))

    list(F=F,dF=dF,d2F=d2F,error=NULL)
}
