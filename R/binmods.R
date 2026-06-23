########################################################################################	
binlogistic <- function(y,eta,nderivs=2, report="logs", ...){
########################################################################################	
    NULL -> logfy -> dlogfy -> d2logfy -> fy -> dfy -> d2fy
    y <- as.vector(y)

    eta <- as.vector(eta)
    n <- length(y)
    if (length(eta) != n) stop(
                    "binlogistic: length(y) != length(eta)")
    mult <- ifelse(y==1,1,-1)

    p <- plogis(eta)
    q <- plogis(-eta)  # doing this way to avoid truncation error with 1-p
    if (report=="logs" | report=="both"){
       logfy <- ifelse(y==1,log(p),log(q))
       if (nderivs >= 1) dlogfy <- matrix(ifelse(y==1,q,-p),nrow=n,ncol=1)
       if (nderivs >= 2) d2logfy <- array(-p*q,c(n,1,1))
       }
    else { # (report != "logs")
       fy <- ifelse(y==1,p,q)
       if (nderivs >= 1) {
            mpq <-  mult*p*q
            dfy <-  matrix(mpq,nrow=n,ncol=1)
            }
       if (nderivs >= 2) d2fy <- array(mpq*(q-p),c(n,1,1))
       }   
list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy,fy=fy,dfy=dfy,d2fy=d2fy,error=NULL)
}

########################################################################################	
binlogisticstrat <- function(eta, nderivs=2, ptildes, ...)
########################################################################################	
{
    NULL -> dQstar -> SptQstar2
    eta <- as.vector(eta)
    n <- length(eta)
    if (length(ptildes) != 2) stop("binlogisticstrat: length(ptildes) != 2")
    p <- plogis(eta)
    q <- plogis(-eta)

    Qstar <- cbind(p,q)
    if (nderivs >= 1) {
              dQstar <- array(0,c(n,2,1))
              dQstar[, ,1] <- p*q %o% c(1,-1)
              }
    if (nderivs >= 2) SptQstar2 <- array(p*q*(q-p)*(ptildes[1]-ptildes[2]),c(n,1,1))

    list(Qstar=Qstar,dQstar=dQstar,SptQstar2=SptQstar2,error=NULL)
}


########################################################################################	
hbinlogistic <- function(hvalue,eta, nderivs=2, ...){
########################################################################################

# hmodel function for direct method with binary regression models where h(y,x) = y
# It is interpreting hvalue=2 as Y1 = 1 (i.e. is a case)
	
    NULL -> Qstar -> dQstar -> d2Qstar

    eta <- eta[,1]
    n <- length(hvalue)
    if (length(eta) != n) stop("hpalmgrn1: length(hvalue[,1]) != length(eta)")
    mult <- ifelse(hvalue==2,1,-1)
    Qstar <- plogis(mult*eta)
    if (nderivs >= 1){
        p <- plogis(eta)
        q <- plogis(-eta)
        mpq <-  mult*p*q
        dQstar <-matrix(mpq,n,1)
        }
    if (nderivs >= 2)  d2Qstar <- array(mpq*(q-p),c(n,1,1))

    list(Qstar=Qstar,dQstar=dQstar,d2Qstar=d2Qstar,error=NULL)
}

########################################################################################	
binprobit <- function(y,eta,nderivs=2, report="logs", ...){
########################################################################################	
    NULL -> logfy -> dlogfy -> d2logfy -> fy -> dfy -> d2fy
    y <- as.vector(y)

    eta <- as.vector(eta)
    n <- length(y)
    if (length(eta) != n) stop(
                    "binprobit: length(y) != length(eta)")
    mult <- ifelse(y==1,1,-1)

    p <- pnorm(eta)
    q <- pnorm(-eta)  # doing this way to avoid truncation error with 1-p
    if (nderivs >= 1) dne <- dnorm(eta)

    if (report=="logs" | report=="both"){
       logfy <- ifelse(y==1,log(p),log(q))
       if (nderivs >= 1) {
           r1 <- dne/p
           r0 <- dne/q
           dlogfy <- matrix(ifelse(y==1,r1,-r0),nrow=n,ncol=1)
       }
       if (nderivs >= 2)   d2logfy <- array(ifelse(y==1,-r1*(eta+r1),r0*(eta-r0) ) , c(n,1,1))
    }
    else { # (report != "logs")
       fy <- ifelse(y==1,p,q)
       if (nderivs >= 1)  dfy <-  matrix(mult*dne,nrow=n,ncol=1)
       if (nderivs >= 2)  d2fy <- array(-mult*eta*dne, c(n,1,1))
   }   
list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy,fy=fy,dfy=dfy,d2fy=d2fy,error=NULL)
}

########################################################################################	
binprobitstrat <- function(eta, nderivs=2, ptildes, ...)
########################################################################################	
{
    NULL -> dQstar -> SptQstar2
    eta <- as.vector(eta)
    n <- length(eta)
    if (length(ptildes) != 2) stop("binprobitstrat: length(ptildes) != 2")
    p <- pnorm(eta)
    q <- pnorm(-eta)
    if (nderivs >= 1) dne <- dnorm(eta)

    Qstar <- cbind(p,q)
    if (nderivs >= 1) {
              dQstar <- array(0,c(n,2,1))
              dQstar[, ,1] <- dne %o% c(1,-1)
              }
    if (nderivs >= 2) SptQstar2 <- array(-eta*dne*(ptildes[1]-ptildes[2]),c(n,1,1))

    list(Qstar=Qstar,dQstar=dQstar,SptQstar2=SptQstar2,error=NULL)
}

########################################################################################	
bincloglog <- function(y,eta,nderivs=2, report="logs", ...){
########################################################################################	
    NULL -> logfy -> dlogfy -> d2logfy -> fy -> dfy -> d2fy
    f1memx <- function(x) {x*(1-x/2*(1-x/3*(1-x/4)))}
    y <- as.vector(y)

    eta <- as.vector(eta)
    n <- length(y)
    if (length(eta) != n) stop(
                    "bincloglog: length(y) != length(eta)")
    mult <- ifelse(y==1,1,-1)

    eeta <- exp(eta)
    q <- exp(-eeta)
    p <- 1-q
    f1memx <- function(x) {x*(1-x/2*(1-x/3*(1-x/4)))}
    if (any(eta< -7)) p[eta< -7] <- f1memx(eeta[eta < -7]) # to avoid cancellation error with 1-q

    if (report=="logs" | report=="both"){
       logfy <- ifelse(y==1,log(p),-eeta)
       if (nderivs >= 1) {
           rr <- eeta*q/p
           dlogfy <- matrix(ifelse(y==1,rr,-eeta),nrow=n,ncol=1)
       }
       if (nderivs >= 2)   d2logfy <- array(ifelse(y==1,rr*(1-eeta-rr),-eeta) , c(n,1,1))
    }
    else { # (report != "logs")
       fy <- ifelse(y==1,p,q)
       if (nderivs >= 1)  dfy <-  matrix(mult*eeta*q,nrow=n,ncol=1)
       if (nderivs >= 2)  d2fy <- array(mult*eeta*q*(1-eeta), c(n,1,1))
   }   
list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy,fy=fy,dfy=dfy,d2fy=d2fy,error=NULL)
}

########################################################################################	
bincloglogstrat <- function(eta, nderivs=2, ptildes, ...)
########################################################################################	
{
   NULL -> dQstar -> SptQstar2
    eta <- as.vector(eta)
    n <- length(eta)
    if (length(ptildes) != 2) stop("bincloglogstrat: length(ptildes) != 2")
    eeta <- exp(eta)
    q <- exp(-eeta)
    p <- 1-q
    f1memx <- function(x) {x*(1-x/2*(1-x/3*(1-x/4)))}
    p[eta< -4] <- f1memx(eeta[eta < -4]) # to avoid cancellation error with 1-q
    Qstar <- cbind(p,q)
    if (nderivs >= 1) {
              dQstar <- array(0,c(n,2,1))
              dQstar[, ,1] <- eeta*q %o% c(1,-1)
              }
    if (nderivs >= 2) SptQstar2 <- array(eeta*q*(1-eeta)*(ptildes[1]-ptildes[2]),c(n,1,1))
    list(Qstar=Qstar,dQstar=dQstar,SptQstar2=SptQstar2,error=NULL)
}
