####################################################################	
lpalmgrn <- function(y,eta,nderivs=2, ...){
# 
#  calculates the derivatives of log pr(Y1=y1,Y2=y2)
#  Here we take account of whether each Yi=1 or != 1
####################################################################	

    w <- palmgrn2(y,eta,nderivs)
    logfy <- log(w$fy)
    NULL -> dlogfy -> d2logfy
    if (nderivs >= 1) dlogfy <- w$dfy/w$fy
    if (nderivs >= 2){
        d2logfy <- array(0,c(nrow(eta),3,3))
        for (j in 1:3) for (k in 1:3) d2logfy[,j,k] <-
                              w$d2fy[,j,k]/w$fy - dlogfy[,j]*dlogfy[,k]
        }
    list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy)
}


####################################################################	
palmgrn2 <- function(y,eta,nderivs=2){
# 
#  calculates the derivatives of pr(Y1=y1,Y2=y2)
#  Here we take account of whether each Yi=1 or != 1
####################################################################	
   w <- palmgrn1(eta,nderivs)
   zeta <- cbind(
          1*(y[,1] !=1)*(y[,2] !=1),
          1*(y[,1] == 1)*(y[,2] != 1) - 1*(y[,1] != 1)*(y[,2] != 1),
          1*(y[,1] != 1)*(y[,2] == 1) - 1*(y[,1] != 1)*(y[,2] != 1),
          ifelse(y[,1] == y[,2],1,-1))
   fy <- zeta[,1] + zeta[,2]*w$ps[,1] + zeta[,3]*w$ps[,2] + zeta[,4]*w$p11
   NULL -> dfy -> d2fy
   if (nderivs>=1) {
        p1d1 <- w$ps[,1]*w$qs[,1]    # Deriv dp1/d eta1.  Here p1=pr(Y1=1)
        p2d2 <- w$ps[,2]*w$qs[,2]    # Deriv dp2/d eta2.  Here p2=pr(Y2=1)
        dfy <- cbind(
                (zeta[,2] + zeta[,4]*w$d1)*p1d1  ,
                (zeta[,3] + zeta[,4]*w$d2)*p2d2  ,
                 zeta[,4]*w$d3*w$OR  )
        }
   if (nderivs>=2) {
        d2fy <- array(0,c(dim(eta),3))
        p1d11 <- w$ps[,1]*w$qs[,1]*(w$qs[,1]-w$ps[,1])  # Deriv d^2 p1/d eta1^2
        p2d22 <- w$ps[,2]*w$qs[,2]*(w$qs[,2]-w$ps[,2])  # Deriv d^2 p2/d eta2^2
        d2fy[,1,1] <-  zeta[,4]*w$d11*p1d1^2  +  (zeta[,2] + zeta[,4]*w$d1)*p1d11 
        d2fy[,1,2] <-  d2fy[,2,1] <-  zeta[,4]*w$d12*p1d1*p2d2
        d2fy[,1,3] <-  d2fy[,3,1] <-  zeta[,4]*w$d13*p1d1*w$OR
        d2fy[,2,2] <-  zeta[,4]*w$d22*p2d2^2  +  (zeta[,3] + zeta[,4]*w$d2)*p2d22
        d2fy[,2,3] <-  d2fy[,3,2] <-  zeta[,4]*w$d23*p2d2*w$OR
#       d2fy[,3,3] <-  zeta[,4]*w$d33*w$OR^2
        d2fy[,3,3] <-  zeta[,4]*w$OR*(w$d3 + w$d33*w$OR)
        }
   list(fy=fy,dfy=dfy,d2fy=d2fy)
}


####################################################################	
palmgrn1 <- function(eta,nderivs=2){
# 
#  calculates the derivatives of p11 = pr(Y1=1,Y2=1).
#
####################################################################	
#Also puts out cbind(pr(Y1=1),pr(y2=1)), cbind(pr(Y1!=1),pr(y2!=1)) & Odds Ratio 

    exp.eta <- exp(eta)
    ps <- exp.eta[,c(1,2)]/(1+exp.eta[,c(1,2)])  # cbind(pr(Y1=1),pr(y2=1))
    qs <- 1/(1+exp.eta[,c(1,2)])                 # cbind(pr(Y1!=1),pr(y2!=1)) 
    OR <- exp.eta[,3]                            # odds ratio

    notindep <- function(ps,qs,OR){
          delta <-   OR - 1
          b <- 1 + (ps[,1]+ps[,2])*delta

          #dinside <- b^2 - 4*delta*(delta+1)*ps[,1]*ps[,2]
          # replace by more stable

          pminussq <- (ps[,1]-ps[,2])^2
          dinside <- 1 + (2*( pminussq + ps[,1]*qs[,1] + ps[,2]*qs[,2] )
                                       + pminussq*delta)*delta
          rinside <- sqrt(dinside)
          h1 <- (b-rinside)/(2*delta)
          h2 <- (b+rinside)/(2*delta)
          pmin <- ifelse(ps[,1]<ps[,2],ps[,1],ps[,2])
          mcond <- ifelse(ps[,1]+ps[,2]-1>0, ps[,1]+ps[,2]-1, 0)
          ifelse( h1<pmin & h1>mcond, h1,h2)
          }

     p11 <- ifelse(OR==1,  ps[,1]*ps[,2],  notindep(ps,qs,OR))

# Calc for 1st derivs

  if (nderivs >= 1){
     p12 <- ps[,1] - p11
     p21 <- ps[,2] - p11
     p22 <- 1 - ps[,1] - ps[,2] + p11

     D <- 1/p11 + 1/p12 + 1/p21 + 1/p22

     d1 <- (1/p22 + 1/p12) / D
     d2 <- (1/p22 + 1/p21) / D
     d3 <- 1 / (D*OR)
     }
  else  NULL -> d1 -> d2 -> d3

# Calc for 2nd derivs

  if (nderivs >= 2){
     delta <-   OR - 1
     denom <- p11*p22*D
     d11 <- 2*d1*delta*(d1-1) / denom
     d12 <- ( OR - delta*(d1+d2-2*d1*d2) ) / denom
     d22 <- 2*d2*delta*(d2-1) / denom
     d13 <- ( delta*d3*(2*d1-1) + p21 - d1*(p21+p12) ) / denom
     d23 <- ( delta*d3*(2*d2-1) + p12 - d2*(p21+p12) ) / denom
     d33 <- 2*d3*(delta*d3 - (p12+p21)) / denom
     }
  else   NULL -> d11 -> d12 -> d13 -> d22 -> d23 -> d33

    list(ps=ps,qs=qs,OR=OR,p11=p11,d1=d1,d2=d2,d3=d3,d11=d11,d12=d12,
              d13=d13,d22=d22, d23=d23,d33=d33)
}


####################################################################	
palmgrn1.strat <- function(eta,nderivs=2, ptildes, ...){
# 
#  calculates the derivatives of pr(Y1=1)
####################################################################	
    if (ncol(eta) != 3) stop("palmgrn1.strat: ncol(eta) != 3")

    NULL -> dQstar -> SptQstar2
    n <- nrow(eta)
    if (length(ptildes) != 2) stop("palmgrn1.strat: length(ptildes) != 2")
    p <- plogis(eta[,1])
    q <- plogis(-eta[,1])

    Qstar <- cbind(p,q)
    if (nderivs >= 1) {
              dQstar <- array(0,c(n,2,3))
              dQstar[, ,1] <- p*q %o% c(1,-1)
              }
    if (nderivs >= 2){
             SptQstar2 <- array(0,c(n,3,3))
             SptQstar2[,1,1] <- p*q*(q-p)*(ptildes[1]-ptildes[2])
             }
    list(Qstar=Qstar,dQstar=dQstar,SptQstar2=SptQstar2,error=NULL)

}


####################################################################	
palm11fn <- function(n, eta){
####################################################################	
   if (n != nrow(eta)) stop("palm11fn:  n != nrow(eta)")

   w <- palmgrn1(eta,nderivs=0)
   list(p11=w$p11,ps=w$ps,qs=w$qs)
}


########################################################################################	
hpalmgrn1 <- function(hvalue,eta, nderivs=2, ...){
########################################################################################

# hmodel function for direct method with Palmgren models where h(y,x) = y1
# It is interpreting hvalue=2 as Y1 = 1 (i.e. is a case)
	
    NULL -> Qstar -> dQstar -> d2Qstar

    eta <- eta[,1]
#    n <- length(hvalue[,1])
    n <- length(hvalue)
    if (length(eta) != n) stop("hpalmgrn1: length(hvalue[,1]) != length(eta)")
#    mult <- ifelse(hvalue[,1]==1,1,-1)
    mult <- ifelse(hvalue==2,1,-1)
    Qstar <- plogis(mult*eta)
    if (nderivs >= 1){
        p <- plogis(eta)
        q <- plogis(-eta)
        mpq <-  mult*p*q
        dQstar <-matrix(0,n,3)
        dQstar[,1] <-  mpq
        }
    if (nderivs >= 2){
        d2Qstar <- array(0,c(n,3,3))
        d2Qstar[,1,1] <- mpq*(q-p)
        }

    list(Qstar=Qstar,dQstar=dQstar,d2Qstar=d2Qstar,error=NULL)
}


####################################################################	
bahadur <- function(y,eta,nderivs=2){
# 
#  calculates the derivatives of pr(Y1=y1,Y2=y2) usin the Bahadur model
#  Here we take account of whether each Yi=1 or != 1
####################################################################	
    if (ncol(eta) != 3) stop("bahadur: ncol(eta) != 3")
 
    exp.eta <- exp(eta)
    ps <- exp.eta[,c(1,2)]/(1+exp.eta[,c(1,2)])  # cbind(pr(Y1=1),pr(y2=1))
    qs <- 1/(1+exp.eta[,c(1,2)])                 # cbind(pr(Y1!=1),pr(y2!=1)) 
    rho <- (exp.eta[,3]-1)/(exp.eta[,3] +1)      # correlation


   pry1 <- ifelse(y[,1]==1,ps[,1],qs[,1])
   pry2 <- ifelse(y[,2]==1,ps[,2],qs[,2])

   m1y1 <- ifelse(y[,1]==1,-1,1)
   m1y2 <- ifelse(y[,2]==1,-1,1)
   m1y12 <- m1y1*m1y2
   sqrtterm <- sqrt(ps[,1]*qs[,1]*ps[,2]*qs[,2])

   fy <- pry1*pry2 + m1y12*rho*sqrtterm
   #if (any(fy<=0)) fy[fy<=0] <- 1e-08  #this is dangerous !

   NULL -> dfy -> d2fy
   if (nderivs >= 1){
      dfy <- matrix(0,nrow(eta),3)
      temp <- m1y12*rho*sqrtterm
      dfy[,1] <- -pry2*m1y1*ps[,1]*qs[,1] +
                       0.5*temp*(qs[,1]-ps[,1])
      dfy[,2] <- -pry1*m1y2*ps[,2]*qs[,2] +
                       0.5*temp*(qs[,2]-ps[,2])
      p3 <- exp.eta[,3]/(1+exp.eta[,3])
      q3 <- 1/(1+exp.eta[,3])
      drhodeta <- 2*p3*q3    # deriv wrt eta
      dfy[,3] <- m1y12*sqrtterm*drhodeta
      }
   if (nderivs >= 2){
      d2fy <- array(0,c(nrow(eta),3,3))
      temp <- 0.25*m1y12*rho*sqrtterm
      d2fy[,1,1] <- (-pry2)*m1y1*ps[,1]*qs[,1]*(qs[,1]-ps[,1]) +
                        temp*((qs[,1]-ps[,1])^2-4*ps[,1]*qs[,1])
      d2fy[,2,2] <- (-pry1)*m1y2*ps[,2]*qs[,2]*(qs[,2]-ps[,2]) +
                        temp*((qs[,2]-ps[,2])^2-4*ps[,2]*qs[,2])
      d2fy[,1,2] <- d2fy[,2,1] <- (m1y12)*ps[,1]*qs[,1]*ps[,2]*qs[,2] +
                        temp*(qs[,1]-ps[,1])*(qs[,2]-ps[,2])
      d2fy[,1,3] <- d2fy[,3,1] <- 0.5*m1y12*sqrtterm*(qs[,1]-ps[,1])*drhodeta
      d2fy[,2,3] <- d2fy[,3,2] <- 0.5*m1y12*sqrtterm*(qs[,2]-ps[,2])*drhodeta
      d2fy[,3,3] <- dfy[,3]*(q3-p3)
      }
   list(fy=fy,dfy=dfy,d2fy=d2fy,ps=ps,qs=qs)
}


####################################################################	
lbahad <- function(y,eta,nderivs=2, ...){
# 
#  calculates the derivatives of log pr(Y1=y1,Y2=y2) for Bahadur model
####################################################################	

   error <- 0
   w <- bahadur(y,eta,nderivs)
   if (min(w$fy, na.rm=TRUE) <= 0){
         error <- 1
         print("lbahad: WARNING -- -ve values to log, evaluation failed")
         }
    NULL -> logfy -> dlogfy -> d2logfy
    if (error != 0) return(list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy,error=error))

    logfy <- log(w$fy)
    if (nderivs >= 1) dlogfy <- w$dfy/w$fy
    if (nderivs >= 2){
        d2logfy <- array(0,c(nrow(eta),3,3))
        for (j in 1:3) for (k in 1:3) d2logfy[,j,k] <-
                              w$d2fy[,j,k]/w$fy - dlogfy[,j]*dlogfy[,k]
        }
    list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy,error=error)
}


####################################################################	
bahad11fn <- function(n, eta){
####################################################################	
   if (n != nrow(eta)) stop("palm11fn:  n != nrow(eta)")

   w <- bahadur(matrix(1,n,3),eta,nderivs=0)
   list(p11=w$fy,ps=w$ps,qs=w$qs)
}


####################################################################	
copula1 <- function(eta,nderivs=2){
# 
#  calculates the derivatives of p00 = pr(Y1=0,Y2=0).
#  Frank's bivariate family
####################################################################
alpha <- eta[,3]
hfn  <- function(alpha,t){
if(any(alpha!=0) && any(t!=0))
  res <- log((exp(-alpha)-1)/(exp(-alpha*t)-1))
else
  res <- -log(t)
res
}
hinv <- function(alpha,t){-(1/alpha)*log(exp(-t)*(exp(-alpha)-1)+1)}
dh   <- function(alpha,t){alpha*exp(-alpha*t)*(exp(-alpha*t)-1)^{-1}}
d2h  <- function(alpha,t){dh(alpha,t)*(dh(alpha,t)-alpha)}
dhdalpha    <- function(alpha,t){dh(alpha,t)*t/alpha-exp(-alpha)/(exp(-alpha)-1)}
d2hdalphadt <- function(alpha,t){(1/alpha)*dh(alpha,t)-t*dh(alpha,t)+
               t*(1/alpha)*(dh(alpha,t)^2)}
d2hdalpha2  <- function(alpha,t){d2hdalphadt(alpha,t)*(t/alpha)-(t/alpha^2)*dh(alpha,t)+
               (exp(-alpha)/(exp(-alpha)-1))*(1-(exp(-alpha))/(exp(-alpha)-1))}
Cfn <- function(alpha,us){
 m <- ncol(us)
 hs <- matrix(NA,nrow(us),m)
 for(i in 1:m)
	hs[,i] <- hfn(alpha,us[,i])
 hinv(alpha,apply(hs,1,sum))
}
dudtheta   <- function(u){exp(u)/(1+exp(u))^2}            # for logistic
d2udtheta2 <- function(u){exp(u)*(1-exp(u))/(1+exp(u))^3} # for logistic
#---------------------------------------------------------------------------------
ps <- plogis(eta[,1:2]) 
us <- qs <- plogis(-eta[,1:2]) 

Cs <- Cfn(alpha,us)
dh1 <- dh(alpha,us[,1])
dh2 <- dh(alpha,us[,2])
d2h1 <- d2h(alpha,us[,1])
d2h2 <- d2h(alpha,us[,2])
dhdalpha1 <- dhdalpha(alpha,us[,1])
dhdalpha2 <- dhdalpha(alpha,us[,2])
frac <- (1/dh(alpha,Cs))
du1dtheta <- dudtheta(eta[,1])
du2dtheta <- dudtheta(eta[,2])
d2u1dtheta2 <- d2udtheta2(eta[,1])
d2u2dtheta2 <- d2udtheta2(eta[,2])

# Calc for 1st derivs

if (nderivs >= 1){
dCudeta <- function(eta,alpha,us,whicheta){
  # us has two columns, whicheta is either 1 or 2
  (1/dh(alpha,Cfn(alpha,us)))*dh(alpha,us[,whicheta])*-dudtheta(eta)
}

dCudalpha <- function(alpha,us){
  # us has two columns
  (1/dh(alpha,Cfn(alpha,us)))*(-dhdalpha(alpha,Cfn(alpha,us)) + 
              dhdalpha(alpha,us[,1]) + dhdalpha(alpha,us[,2]))
}

d1 <- dCudeta(eta[,1],alpha,qs,1)
d2 <- dCudeta(eta[,2],alpha,qs,2)
d3 <- dCudalpha(alpha,qs)
d1s <- cbind(d1,d2,d3)
}
else  d1 <- d2 <- d3 <- d1s <- NULL

# Calc for 2nd derivs

if (nderivs >= 2){
d11 <- frac*(d2h1 * (du1dtheta^2) + dh1 * (-d2u1dtheta2) - d2h(alpha,Cs)*d1^2)
d12 <- frac*(-d2h(alpha,Cs)*d1*d2)
d13 <- frac*(d2hdalphadt(alpha,us[,1]) * (-dudtheta(eta[,1])) - d1*(d2h(alpha,Cs)*d3 + 
       d2hdalphadt(alpha,Cs)))
d22 <- frac*(d2h2 * (du2dtheta^2) + dh2 * (-d2u2dtheta2) - d2h(alpha,Cs)*d2^2)
d23 <- frac*(d2hdalphadt(alpha,us[,2]) * (-dudtheta(eta[,2])) - d2*(d2h(alpha,Cs)*d3 + 
       d2hdalphadt(alpha,Cs)))
# Only one to fix now...  
d33 <- frac*(d2hdalpha2(alpha,us[,1]) + d2hdalpha2(alpha,us[,2]) - d3*(d2h(alpha,Cs)*d3 + 
       2*d2hdalphadt(alpha,Cs)) - d2hdalpha2(alpha,Cs))
d2s <- cbind(d11,d22,d33,d12,d23,d13)
}
else  d11 <- d12 <- d13 <- d22 <- d23 <- d33 <- d2s <- NULL

p00 <- Cfn(alpha,1/(1+exp(eta[,1:2])))

list(ps=ps,qs=qs,d1s=d1s,d1=d1,d2=d2,d3=d3,d2s=d2s,d11=d11,d12=d12,d13=d13,d22=d22,
     d23=d23,d33=d33,p00=p00)
}


####################################################################	
copula2 <- function(y,eta,nderivs=2){
# 
#  calculates the derivatives of pr(Y1=y1,Y2=y2)
#  Here we take account of whether each Yi=1 or != 1
####################################################################	
   w <- copula1(eta,nderivs)
   zeta <- cbind(
          1*(y[,1] == 1 & y[,2] == 1),
          1*(y[,1] == 0 & y[,2] == 1) - 1*(y[,1] == 1 & y[,2] == 1),
          1*(y[,1] == 1 & y[,2] == 0) - 1*(y[,1] == 1 & y[,2] == 1),
          ifelse(y[,1] == y[,2],1,-1))
   fy <- zeta[,1] + zeta[,2]*w$qs[,1] + zeta[,3]*w$qs[,2] + zeta[,4]*w$p00
   NULL -> dfy -> d2fy
   if (nderivs>=1) {
        q1d1 <- -w$ps[,1]*w$qs[,1]    # Deriv dq1/d eta1.  Here q1=pr(Y1!=1)
        q2d2 <- -w$ps[,2]*w$qs[,2]    # Deriv dq2/d eta2.  Here q2=pr(Y2!=1)
        dfy <- cbind(
                zeta[,2]*q1d1 + zeta[,4]*w$d1  ,  
                zeta[,3]*q2d2 + zeta[,4]*w$d2  ,   
                zeta[,4]*w$d3)                    
        }
   if (nderivs>=2) {
        d2fy <- array(0,c(dim(eta),3))
        q1d11 <- w$ps[,1]*w$qs[,1]*(w$ps[,1]-w$qs[,1])  # Deriv d^2 q1/d eta1^2
        q2d22 <- w$ps[,2]*w$qs[,2]*(w$ps[,2]-w$qs[,2])  # Deriv d^2 q2/d eta2^2
        d2fy[,1,1] <-  zeta[,4]*w$d11 + zeta[,2]*q1d11 
        d2fy[,1,2] <-  d2fy[,2,1] <-  zeta[,4]*w$d12
        d2fy[,1,3] <-  d2fy[,3,1] <-  zeta[,4]*w$d13
        d2fy[,2,2] <-  zeta[,4]*w$d22 + zeta[,3]*q2d22

        d2fy[,2,3] <-  d2fy[,3,2] <-  zeta[,4]*w$d23
        d2fy[,3,3] <-  zeta[,4]*w$d33
        }
   list(fy=fy,dfy=dfy,d2fy=d2fy)
}


####################################################################	
lcopula <- function (y, eta, nderivs = 2, ...) 
####################################################################	
{
    w <- copula2(y, eta, nderivs)
    logfy <- log(w$fy)
    d2logfy <- dlogfy <- NULL
    if (nderivs >= 1) 
        dlogfy <- w$dfy/w$fy
    if (nderivs >= 2) {
        d2logfy <- array(0, c(nrow(eta), 3, 3))
        for (j in 1:3) for (k in 1:3) d2logfy[, j, k] <- w$d2fy[, 
            j, k]/w$fy - dlogfy[, j] * dlogfy[, k]
    }
    list(logfy = logfy, dlogfy = dlogfy, d2logfy = d2logfy)
}


####################################################################	
biv1pbc.strat <- function(eta,nderivs=2, ptildes, ...){
# 
# two strata, (Y1=1) and (Y1!=1)
# calculates derivatives of pr(Y1=1) for Palmgren, Bahadur and Copula 
# same as "palmgrn1.strat" in contents !
####################################################################	
    if (ncol(eta) != 3) stop("biv1pbc.strat: ncol(eta) != 3")

    NULL -> dQstar -> SptQstar2
    n <- nrow(eta)
    if (length(ptildes) != 2) stop("biv1pbc.strat: length(ptildes) != 2")
    p <- plogis(eta[,1])
    q <- plogis(-eta[,1])

    Qstar <- cbind(p,q)
    if (nderivs >= 1) {
              dQstar <- array(0,c(n,2,3))
              dQstar[, ,1] <- p*q %o% c(1,-1)
              }
    if (nderivs >= 2){
             SptQstar2 <- array(0,c(n,3,3))
             SptQstar2[,1,1] <- p*q*(q-p)*(ptildes[1]-ptildes[2])
             }

    list(Qstar=Qstar,dQstar=dQstar,SptQstar2=SptQstar2,error=NULL)

}


####################################################################	
biv00pbc.strat <- function(eta, nderivs=2, ptildes, method, ...){
#
#  2 strata, (Y1!=1 and Y2!=1) and the rest
#  calculates derivatives of pr(Y1!=1 and Y2!=1) for Palmgren, Bahadur and Copula
####################################################################
   if (ncol(eta) != 3) stop("biv00pbc.strat: ncol(eta) != 3")
   
   NULL -> dQstar -> SptQstar2
   n <- length(eta[,1])
   y00 <- matrix(0,n,2)
   
   if (method == "palmgren")  w00 <- palmgrn2(y00,eta,nderivs)
   else if (method == "bahadur")  w00 <- bahadur(y00,eta,nderivs)
        else if (method == "copula") w00 <- copula2(y00,eta,nderivs)
             else stop("biv00pbc.strat: this method is not implemented!")
      
   Qstar <- cbind(w00$fy, 1-w00$fy)

   if (nderivs >= 1) {
     dQstar <- array(NA,c(n,2,3))
     for (i in 1:3) dQstar[,,i] <- w00$dfy[,i] %o% c(1,-1)
   }

   if (nderivs >=2) {
     SptQstar2 <- array(NA,c(n,3,3))
     SptQstar2 <- w00$d2fy*(ptildes[1]-ptildes[2])
   }

   list(Qstar=Qstar,dQstar=dQstar,SptQstar2=SptQstar2,error=NULL)
}


####################################################################	
lspml2 <- function(y,eta,nderivs=2, ...){
# 
#  calculates the derivatives of log pr(Y1=y1,Y2=y2)
#  Here we take account of whether each Yi=1 or != 1
####################################################################	

    w <- spml2(y,eta,nderivs)    
    logfy <- log(w$fy)
    
    NULL -> dlogfy -> d2logfy
    if (nderivs >= 1) dlogfy <- w$dfy/w$fy
    if (nderivs >= 2){
        d2logfy <- array(0,c(dim(eta),dim(eta)[2]))
        for (j in 1:dim(eta)[2]) 
          for (k in 1:dim(eta)[2]) 
             d2logfy[,j,k] <- w$d2fy[,j,k]/w$fy - dlogfy[,j]*dlogfy[,k]
    }

    list(logfy=logfy,dlogfy=dlogfy,d2logfy=d2logfy)
}


####################################################################	
spml2 <- function(y,eta,nderivs=2){
# 
#  calculates the derivatives of pr(Y1=y1,Y2=y2)
#  Here we take account of whether each Yi=1 or != 1
####################################################################	
    NULL -> fy -> dfy -> d2fy
   
    exp.eta <- exp(eta) 
    ey1 <- exp.eta[,1]; ey2 <- exp.eta[,2]
    py1 <- ifelse(y[,1]==1, ey1/(1+ey1), 1/(1+ey1))
    py2 <- ifelse(y[,2]==1, ey2/(1+ey2), 1/(1+ey2))
    ps <- cbind(py1, py2); rps <- cbind(py2, py1)
    qs <- 1-ps; rqs <- cbind(qs[,2], ps[,1])
    mult <- ifelse(y==1, 1, -1); rmult <- cbind(mult[,2], mult[,1])
    
    fy <- py1*py2
                
    if (nderivs >= 1) {
       part1 <- mult*ps*qs
       dfy <- matrix(part1*rps,dim(eta)[1],dim(eta)[2])
    }
    
    if (nderivs >= 2) {
       d2fy <- array(0,c(dim(eta),dim(eta)[2]))
       part2 <- mult*(qs-ps)
       part3 <- rmult*rps*rqs
       d1 <- part1*part2*rps
       d2 <- part1*part3
       d2fy[,1,1] <- d1[,1]; d2fy[,2,2] <- d1[,2]
       d2fy[,1,2] <- d2[,1]; d2fy[,2,1] <- d2[,2] #should be same
    }

    list(fy=fy, dfy=dfy, d2fy=d2fy)
}
