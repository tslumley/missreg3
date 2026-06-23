#################################
### SECTION 1: FRONT FUNCTION ###
#################################


#na.keep = function(x) x


################################################
model.apply <- function(formula,data,weights=NULL){
################################################
  mf <- call <- match.call()
  mf[[1]] <- as.name("model.frame")
  names(mf)[1] <- "model"
  mf$na.action <- as.name("na.keep")

  #mf <- eval(mf[c("model","formula","data","weights","na.action")], sys.frame(sys.parent()))
  mf <- eval(mf, sys.frame(sys.parent()))
  Terms <- attr(mf,"terms")
  X <- model.matrix(Terms,mf)
  if(!is.matrix(X)) X <- matrix(X)
  weights <- model.extract(mf, "weights")
  y <- model.extract(mf,"response")

  terms1 <- attr(terms(formula),"term.labels")
  order1 <- attr(terms(formula),"order")
  assign1 <- attr(X,"assign"); fnames1 <- names(attr(X,"contrasts"))

  list(y=y,X=X,weights=weights,Terms=Terms,terms1=terms1,order1=order1,
	assign1=assign1,fnames1=fnames1)
}
# model.apply(y~x1d+x2+x1d*x2,datcsv)


###########################################################################################
array3form2 <- function(Xdat, ClusInd.name=NULL, IntraClus.name=NULL, MaxInClus=NULL, 
			rmsingletons=FALSE, freqformat=TRUE)
###########################################################################################
{

# This is a new function to convert a data frame into right shape of Array;
# IntraClus ID (if provided) is kept in the original order regardless of missing rows;
# If MaxInClus is supplied, we take only the first MaxInClus individuals in any cluster

# Takes cluster data in which rows are individuals and clustermembership is given by the
# column vector ClusInd and makes a 3-D array in which the 1st dim indexes clusters, the 2nd
# dim is the cols of XMat and the 3rd indexes individuals within the cluster.
# The 3rd index runs from 1 to the maximum number of individuals in a cluster
# (or MaxInClus if supplied).  Array positions for non-existent individuals are
# packed out with missing values (NA)

if (!is.data.frame(Xdat)) stop("array3form2: Xdat should be a data frame")
n <- dim(Xdat)[1]
Xdat0 <- Xdat

pos1 <- pos2 <- NULL
if (is.null(ClusInd.name)) ClusInd <- 1:n
  else {
    	pos1 <- match(ClusInd.name,names(Xdat))
	ClusInd <- Xdat[,pos1]  }

if (is.null(IntraClus.name)) {
  	if (length(unique(ClusInd))==n) IntraClus <- rep(1,n)
  	else stop("IntraClus.name should be provided when Cluster size >1") }
  else {
	pos2 <- match(IntraClus.name,names(Xdat))
	IntraClus <- Xdat[,pos2]  }

pos <- c(pos1,pos2)
if (length(pos)>0) Xdat <- Xdat[,-pos]

id1 <- sort(unique(ClusInd))
id2 <- sort(unique(IntraClus))
id <- expand.grid(x2=id2,x1=id1)
newid <- paste(id$x1,id$x2)
origid <- paste(ClusInd, IntraClus)

newdat0 <- matrix(NA,length(newid),dim(Xdat)[2])
dimnames(newdat0) <- list(1:dim(newdat0)[1], dimnames(Xdat)[[2]])
newdat0 <- t(newdat0)
Xdat1 <- t(Xdat)

for(i in 1:length(origid)) {
   tmpid <- which(origid[i]==newid)
   if (length(tmpid)>0) newdat0[,tmpid] <- Xdat1[,i]
}
newdat0 <- t(newdat0)

newdat <- data.frame(id$x1, id$x2, newdat0)
dimnames(newdat) <- list(1:dim(newdat)[1], c(ClusInd.name, IntraClus.name, dimnames(Xdat)[[2]]))

for (i in 1:dim(Xdat)[2]) {
   vclass <- class(Xdat[,i])
   if ((vclass %in% c("character","factor"))==FALSE) newdat[,2+i] <- as.numeric(newdat0[,i])
}

# Now turn it into the desired array with (nClus, nXvar, Jmax);
XArray <- array(t(newdat),c(dim(newdat)[2],length(id2),length(id1)))
XArray <- aperm(XArray,c(3,1,2))

Jis <- table(newdat[,1])
nClus <- length(Jis)
Jmax <- max(Jis)
ncols <- dim(newdat)[2]

ClusLab <- unique(newdat[,1])  # Start to store Cluster labels from ClusInd
if (rmsingletons) {
    XArray <- XArray[Jis>1, , ,drop=FALSE]
    ClusLab <- ClusLab[Jis>1]
    Jis <- Jis[Jis>1]
    nClus <- length(Jis)
}
if (!is.null(MaxInClus)){
    Jmax <- min(length(id2),MaxInClus)
    Jis[Jis>Jmax] <- Jmax
    XArray <- XArray[, ,1:Jmax,drop=FALSE]
}
counts <- rep(1,nClus)
if (freqformat) {
    # Turn XArray into nClus rows so we can condense it into a frequency format
    XMatWide <- matrix(aperm(XArray,c(2,3,1)),nrow=nClus,ncol=Jmax*ncols,byrow=TRUE)
    z <- weightform(XMatWide)
    Jis <- Jis[z$ind]; counts <- z$counts; nClus <- length(z$counts); ClusLab <- 1:nClus
    XArray <- array(XMatWide[z$ind,],c(nClus,ncols,Jmax))
}
dimnames(XArray) <- list(ClusLab,dimnames(newdat)[[2]],id2)

list(XArray=XArray,Jis=Jis,counts=counts,newdat=newdat)

} # Endfn array3form2



#################################
###  SECTION 2: KEY FUNCTIONS ###
#################################


######################################################################################
Clusbinlogis <- function(y, eta, nderivs=2, modelfn, sigma, Jis=NULL, nzval0=20, ...){
######################################################################################

### "y" is the response of interest with dimension of (nn, Jmax+1) with the last
### column corresponding to the case-control response as a binary variable

### "eta" has a dimention of (nn, Jmax+1) with the last column corresponding to
###  the linear predictor for binary model part (z|y,x; theta2);

### Innerfn is for Y1-model only with three possible links
### sigma is not part of eta matrix so input separately
### Jis indicates the size of each cluster, maximum equal to Jmax


nn <- dim(eta)[1]
Jmax <-  dim(eta)[2]-1

if (is.null(Jis)) Jis <- rep(Jmax,nn)
else {
  if (any(Jis>Jmax)) stop("Clusbinlogis: Maximum of Jis should not beyond Jmax in eta")
  if (length(Jis) != nn) stop("Clusbinlogis: Dimension problem in Jis - different from nrow(eta)")
}

if (nrow(y)!=nn) stop("Clusbinlogis: Dimension problem in y - different from nrow(eta)")


y1 <- y[,-dim(y)[2],drop=FALSE]
y2 <- as.vector(y[,dim(y)[2]])

eta1 <- eta[,-dim(eta)[2],drop=FALSE]
eta2 <- eta[,dim(eta)[2],drop=FALSE]


# ------------------  Random Intercept Model part  -----------------------
# (logit, probit, cloglog links can be incorparated in this model)

dlogldtheta <- NULL

#zval <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
#nzval <- length(zval)
#intwts <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
#intwts[3] <- 8/15

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)

Py <- dPy <- d2Py <- rep(0,nn) 

if (nderivs >= 1) {
    storePjzsi <- as.list(1:nzval)
    storePygz <- matrix(0,nrow=nn,ncol=nzval)
}

for (i in 1:length(zval)){
    Pygz <- rep(1,nn)
    if (nderivs >= 1) storePjzsi[[i]] <- as.list(1:Jmax)

    for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      Pjz <- modelfn(y1[ind,j],sigma*zval[i]+eta1[ind,j],nderivs=nderivs,report="vals")
      if (nderivs >= 1) storePjzsi[[i]][[j]] <- Pjz
      Pygz[ind] <- Pygz[ind]*Pjz$fy
    }
    Py <- Py + Pygz*intwts[i]  ###
    lPy <- log(Py)  ###

    if (nderivs >=1) {
      storePygz[,i] <- Pygz
      dPygz <- rep(0,nn)

      for (j in 1:Jmax){
        ind <- (1:nn)[Jis >= j]
        Pjz <- storePjzsi[[i]][[j]]
 	dPygz[ind] <- dPygz[ind] + as.vector(Pjz$dfy)/Pjz$fy
      }
      dPy <- dPy + Pygz*intwts[i]*dPygz
      dlPy <- dPy/Py
    }
}

if (nderivs >=2) {
    d2Py0 <- t(dPy) %*% dPy
    d2lPy0 <- t(dlPy) %*% dlPy

    for (i in 1:length(zval)){
      Pygz <- storePygz[,i]
      d2Pygz <- rep(0,nn)

      for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= j]
         Pjz <- storePjzsi[[i]][[j]]
         d2Pygz[ind] <- d2Pygz[ind] + as.vector(Pjz$d2fy)/Pjz$fy

         if (j < Jmax) for (k in (j+1):Jmax){ 
    	    indk <- (1:nn)[Jis >= k]
            Pkz <- storePjzsi[[i]][[k]]
 	    jrowuse <- (1:length(ind))[table(factor(indk, levels=ind))==1]
            temp <- (as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*(as.vector(Pkz$dfy)/Pkz$fy)
            d2Pygz[indk] <- d2Pygz[indk] + temp + t(temp)
         }
      }
      d2Py <- d2Py + Pygz*intwts[i]*d2Pygz
      d2lPy <- d2Py/Py      
    }

    d2Py <- d2Py - d2Py0
    d2lPy <- d2lPy - d2lPy0
}


# ------------------  Second Model part  -----------------------
# (same link is used here as for random intercept model)


NULL -> logfy -> dlogfy -> d2logfy -> fy -> dfy -> d2fy

out1 <- modelfn(y2,eta2,nderivs=nderivs,report="vals")
out2 <- modelfn(y2,eta2,nderivs=nderivs,report="logs")

fy <- out1$fy
logfy <- out2$logfy

if (nderivs >= 1) {
   dfy <-  out1$dfy
   dlogfy <- out2$dlogfy }

if (nderivs >= 2) {
   d2fy <- out1$d2fy
   d2logfy <- out2$d2logfy   }


# ------------------  Combined outputs  -----------------------

fy <- Py*fy
logfy <- lPy + logfy

if (nderivs >= 1) {
  dfy <- dPy*dfy
  dlogfy <- dlogfy + dlPy }

if (nderivs >= 2) {
  d2fy <- d2Py*d2fy
  d2logfy <- d2logfy + d2lPy }


list(fy=fy, dfy=dfy, d2fy=d2fy, logfy=logfy, dlogfy=dlogfy, d2logfy=d2logfy, error=NULL)
}


###############################################################################################
Clusbinlogis2 <- function(y, eta, X1, X2, nderivs=2, modelfn, sigma, Jis=NULL, nzval0=20, ...){
###############################################################################################

### Calculate joint density and its derivatives with respect to theta;
### with dimensions of nn, nn*ntheta and nn*ntheta*ntheta for fy, dfy and d2fy;
### served for MEtaStratModInf2 function only;

### "y" is the response of interest with dimension of (nn, Jmax+1) with the last
### column corresponding to the case-control response as a binary variable

### "eta" has a dimention of (nn, Jmax+1) with the last column corresponding to
###  the linear predictor for binary model part (z|y,x; theta2);
### sigma is not part of eta1 matrix so input separately
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

#zval <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
#nzval <- length(zval)
#intwts <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
#intwts[3] <- 8/15

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)

ntheta1 <- dim(X1)[2]+1
ntheta2 <- dim(X2)[2]
ntheta <- ntheta1+ntheta2

Py <- Pz <- rep(0,nn)
if (nderivs >=1) {
    dPy <- matrix(0,nrow=nn,ncol=ntheta1)
    dPz <- matrix(0,nrow=nn,ncol=ntheta2)
}
if (nderivs >= 1) {
    storePjzsi <- as.list(1:nzval)
    storePygz <- matrix(0,nrow=nn,ncol=nzval)
}

for (i in 1:length(zval)){
    Pygz <- rep(1,nn)
    if (nderivs >= 1) storePjzsi[[i]] <- as.list(1:Jmax)
    for (j in 1:Jmax){
      ind <- (1:nn)[Jis >= j]
      Pjz <- modelfn(y1[ind,j],sigma*zval[i]+eta1[ind,j],nderivs=nderivs,report="vals")
      if (nderivs >= 1) storePjzsi[[i]][[j]] <- Pjz
      Pygz[ind] <- Pygz[ind]*Pjz$fy
    }
    if (nderivs >= 1) storePygz[,i] <- Pygz
    Py <- Py + Pygz*intwts[i]
    if (nderivs >=1) {
      for (j in 1:Jmax){
        ind <- (1:nn)[Jis >= j]
        Xj <- if (length(ind)>1) cbind(X1[ind, ,j],sigma*rep(zval[i],length(ind)))
              else matrix(c(X1[ind, ,j],sigma*zval[i]),nrow=1)
        Pjz <- storePjzsi[[i]][[j]]
        dPy[ind, ] <- dPy[ind, ] + Pygz[ind]*intwts[i]*(as.vector(Pjz$dfy)/Pjz$fy)*Xj 
      }
    }
}

## TEST CODES;
#a<-matrix(1:12,4,3)
#a1<-array(t(a),c(3,1,4))
#b<-apply(a1,3, function(x) x %*% t(x))
#array(b,c(3,3,4))
##

if (nderivs >=2) {
    d2Py_a <- d2Py_b <- array(0,c(nn, ntheta1, ntheta1))

    for (i in 1:length(zval)){
      Pygz <- storePygz[,i]

      for (j in 1:Jmax){
         ind <- (1:nn)[Jis >= j]
         Xj <- if (length(ind)>1) cbind(X1[ind, ,j],sigma*rep(zval[i],length(ind))) 
               else  matrix(c(X1[ind, ,j],sigma*zval[i]),nrow=1)
         Pjz <- storePjzsi[[i]][[j]]

         Xjnew <- array(t(Xj),c(ntheta1,1,length(ind)))
         tmp <- apply(Xjnew,3,function(x) x %*% t(x))
         Xj2 <- array(tmp, c(ntheta1, ntheta1, length(ind)))
         Xj2 <- aperm(Xj2,c(3,1,2))
         d2Py_a[ind,,] <- d2Py_a[ind,,] + intwts[i]*Pygz[ind]*as.vector(Pjz$d2fy)/Pjz$fy*Xj2 

         if (j < Jmax) for (k in (j+1):Jmax){ 
            indk <- (1:nn)[Jis >= k]
            Xk <- if (length(indk)>1)  cbind(X1[indk, ,k],sigma*rep(zval[i],length(indk)))
                  else matrix(c(X1[indk, ,k],sigma*zval[i]),nrow=1)
            Pkz <- storePjzsi[[i]][[k]]
            jrowuse <- (1:length(ind))[table(factor(indk,levels=ind))==1]
            Xjjrowuse <- if (length(jrowuse)>1) Xj[jrowuse,] 
                         else matrix(Xj[jrowuse,], nrow=1)

	    tmp1 <- tmp2 <- NULL
 	    for (t in 1:ntheta1) {
		tmp1 <- cbind(tmp1,Xjjrowuse[,t]*Xk)
		tmp2 <- cbind(tmp2,Xk[,t]*Xjjrowuse)
	    }
	    Xjj2 <- array(tmp1, c(length(indk),ntheta1,ntheta1))
	    Xk2 <- array(tmp2, c(length(indk),ntheta1,ntheta1))
	 
	    part1 <- intwts[i]*Pygz[indk]*(as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*(as.vector(Pkz$dfy)/Pkz$fy)*Xjj2
	    part2 <- intwts[i]*Pygz[indk]*(as.vector(Pjz$dfy)/Pjz$fy)[jrowuse]*(as.vector(Pkz$dfy)/Pkz$fy)*Xk2
            d2Py_b[indk,,] <- d2Py_b[indk,,] + part1 + part2  
          }
       }
     }

    d2Py <- d2Py_a + d2Py_b
    d2Py[,ntheta1,ntheta1] <- d2Py[,ntheta1,ntheta1] + dPy[,ntheta1]
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



###################################################################################	
MEtaProspModInf2 <- function(theta, nderivs=2, y, x, wts=1, modelfn, inxStrat,
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


part1 <- RClusProspModInf(theta=theta1, nderivs=nderivs, y=y1, x=x1, wts=wts, modelfn=modelfn, 
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


###########################################################################################
MEtaStratModInf2<-function(theta,nderivs=2,ymat,xlist,Bcounts,ptildes,modelfn,inxStrat,
			   Jis=NULL, npar, nzval0=20, ...)
###########################################################################################
{ 

theta <- as.vector(theta)
ntheta <- length(theta)

theta1 <- theta[1:npar[1]]
theta2 <- theta[(npar[1]+1):sum(npar)]
ntheta1 <- length(theta1)
ntheta2 <- length(theta2)

betas  <- theta1[-ntheta1]
logsig <- theta1[ntheta1]
sigma  <- exp(logsig)

if (length(xlist) != 4) stop("xlist has wrong dimension!")
else if (is.null(names(xlist))) names(xlist) <- c("orig","1","2","xstrata")

strata <- xlist[[match("xstrata",names(xlist))]]
X1  <- (xlist[[match("orig",names(xlist))]])[strata==inxStrat,2:npar[1],1, drop=FALSE]  # reduce to cluster size of 1
X2 <- (xlist[[match("orig",names(xlist))]])[strata==inxStrat,(npar[1]+1):sum(npar),1]
X2_y1 <- (xlist[[match("1",names(xlist))]])[strata==inxStrat,(npar[1]+1):sum(npar),1]
X2_y2 <- (xlist[[match("2",names(xlist))]])[strata==inxStrat,(npar[1]+1):sum(npar),1] 

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

#zval <- sqrt(5 + (c(1,-1,0,-1,1)*sqrt(10)))*c(-1,-1,0,1,1)
#nzval <- length(zval)
#intwts <- (7 + 2*c(-1,1,0,1,-1)*sqrt(10))/60
#intwts[3] <- 8/15

nzval <- nzval0
ghq20 <-ghq(nzval,modified=FALSE)  # need R-library glmmML
zval <- ghq20$zeros*sqrt(2)
intwts <- ghq20$weights/sqrt(pi)

w11 <- Clusbinlogis2(y11,eta=etarray[,,2], X1=X1, X2=X2_y1, nderivs=nderivs, modelfn=modelfn, sigma=sigma, Jis=Jis, nzval0=nzval0, ...)
w21 <- Clusbinlogis2(y21,eta=etarray[,,3], X1=X1, X2=X2_y2, nderivs=nderivs, modelfn=modelfn, sigma=sigma, Jis=Jis, nzval0=nzval0, ...)
w12 <- Clusbinlogis2(y12,eta=etarray[,,2], X1=X1, X2=X2_y1, nderivs=nderivs, modelfn=modelfn, sigma=sigma, Jis=Jis, nzval0=nzval0, ...)
w22 <- Clusbinlogis2(y22,eta=etarray[,,3], X1=X1, X2=X2_y2, nderivs=nderivs, modelfn=modelfn, sigma=sigma, Jis=Jis, nzval0=nzval0, ...)

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


