####################################################################	
lspml2locsc <- function(y,eta,nderivs=2, errdist,  ...){
# 
#  calculates the derivatives of log pr(Y1=y1,Y2=y2)
#  Here we consider logistic model for Y1|Y2 but linear model for Y2
####################################################################	

    w <- spml2locsc(y,eta,nderivs,errdist)   
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
spml2locsc <- function(y,eta,nderivs=2,errdist, ...){
# 
#  calculates the derivatives of pr(Y1=y1,Y2=y2)
#  Here we consider logistic model for Y1|Y2 but linear model for Y2
#  eta1 for Y1|Y2-model, eta2 and eta3 for Y2-model
####################################################################	
    NULL -> fy -> dfy -> d2fy

    ww <- locscale(y[,2],eta[,-1],nderivs,errdist,report="dens") 
    
    ey1 <- exp(eta[,1])
    py1 <- ifelse(y[,1]==1, ey1/(1+ey1), 1/(1+ey1))
    py2 <- ww$fy
    multy1 <- ifelse(y[,1]==1, 1, -1)
    
    fy <- py1*py2
                
    if (nderivs >= 1) {
       dfy1 <- multy1*py1*(1-py1)
       dfy2 <- ww$dfy
       dfy <- cbind(dfy1*py2, py1*dfy2)
    }
    
    if (nderivs >= 2) {
       d2fy <- array(0,c(dim(eta),dim(eta)[2]))
       part1 <- multy1*(1-2*py1)
       d2fy[,1,1] <- dfy1*part1*py2
       d2fy[,1,-1] <- dfy1*dfy2 
       d2fy[,2,1] <- d2fy[,1,2] #==dfy1*dfy2[,1]
       d2fy[,3,1] <- d2fy[,1,3] #==dfy1*dfy2],2]
     
       d2fy2 <- ww$d2fy
       d2fy[,-1,-1] <- py1*d2fy2
    }

    list(fy=fy, dfy=dfy, d2fy=d2fy)
}


