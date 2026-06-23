#####################################################################################
rclusbin <- function(formula, data, weights=NULL, ClusInd=NULL, IntraClus=NULL,  
                xstrata=NULL, ystrata=NULL, obstype.name="obstype", NMat=NULL, xs.includes=FALSE, 
                MaxInClus=NULL, rmsingletons=FALSE, retrosamp="proband", gamma=NULL, nzval0=20, 
                fit=TRUE, devcheck=FALSE, linkname="logit", start=NULL, Qstart=NULL, sigma=NULL, 
                control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
#####################################################################################
{
# -- When we have clusters of size > 1, the "weights" give within-cluster counts which
#    is the total number of individuals with identical variable values in a cluster;
# -- When all clusters are of size 1, the "weights" give cluster counts which is the 
#    total number of clusters with identical variable values. 
# -- When retrosamp="gamma", the program only deals with 2-ystrata-2-gamma situation.

# Process formula -------------------------------------------------------------------

mf <- call <- match.call()
mf[[1]] <- as.name("model.frame")
mf$na.action <- as.name("na.keep")
NULL -> mf$weights -> mf$ClusInd -> mf$IntraClus
NULL -> mf$xstrata -> mf$ystrata -> mf$obstype.name -> mf$NMat -> mf$xs.includes 
NULL -> mf$MaxInClus -> mf$rmsingletons -> mf$retrosamp -> mf$gamma -> mf$nzval0
NULL-> mf$fit -> mf$devcheck -> mf$linkname -> mf$start -> mf$Qstart -> mf$sigma 
NULL -> mf$control -> mf$control.inner

mf <- eval(mf, sys.frame(sys.parent()))
y <- model.extract(mf,"response")
n <- if (is.matrix(y)) dim(y)[1] else length(y)
Terms <- attr(mf,"terms")
X <- model.matrix(Terms,mf)
if(!is.matrix(X)) X <- matrix(X)
XOrig <- X

if (is.null(weights)) w <- rep(1,n)
  else w <- data[,match(weights,names(data))]
if (is.null(ClusInd)) ClusInd <- 1:n
  else ClusInd <- data[,match(ClusInd,names(data))]
if (is.null(IntraClus)) IntraClus <- NULL
  else IntraClus <- data[,match(IntraClus,names(data))]

terms1 <- attr(terms(formula),"term.labels")
order1 <- attr(terms(formula),"order")
assign1 <- attr(X,"assign"); fnames<-names(attr(X,"contrasts"))


# Input checking ----------------------------------------------------------

if (is.na(match(obstype.name,names(data)))) 
   stop(paste("Dataframe did not have a column named",obstype.name))
obstype <- data[,match(obstype.name,names(data))]
if (any(is.na(data[,obstype.name])))
   stop(paste("Missing values not permitted in obstype variable",obstype.name))

#prosp <- all(obstype %in% c("uncond","y|x"))
prosp <- all(!is.na(match(obstype,c("uncond","y|x"))))
havestrata <- length(obstype[obstype=="strata"])>0
if (is.null(NMat)) haveNMat <- FALSE
  else haveNMat <- TRUE

if (havestrata && haveNMat)
  stop("NMat should only be provided when there are no strata observations!")

lname <- c("logit","probit","cloglog")
if (is.na(match(linkname, lname)))
   stop("This link is not implemented !")

haveretro <- any(!is.na(match(obstype,"retro")))
if (!haveretro) retrosamp <- NULL
if (haveretro && is.null(retrosamp))
  stop("A retrosamp scheme should be provided !")

if (is.null(retrosamp) || (!is.null(retrosamp) && retrosamp!="gamma")) 
    gamma <- NULL
retroname <- c("proband","allcontrol","gamma")
if (!is.null(retrosamp) && is.na(match(retrosamp, retroname)))
   stop("This retrosamp scheme is not implemented !")


if (any(is.na(w))) stop("Weights should not contain NAs !")
if (!is.null(sigma) && (length(sigma)>1 || sigma<=0))
      stop("Illegal sigma value !")


# Form the x-strata --------------------------------------------------------

if (!is.null(xstrata)) { ## Set xstrata=NULL if there is only one stratum exists !
  vx <- as.matrix(data[,match(xstrata,names(data))])
  vx1 <- vx[which(apply(!is.na(vx),1,sum)==dim(vx)[2]),,drop=FALSE]
  if (dim(unique(vx1))[1]==1) xstrata <- NULL }

xStrat <- if (is.null(xstrata)) rep(1,n) else {
    xstratform <- as.formula(paste("~","-1 + ",paste(xstrata,sep="",collapse=":")))
    data1 <- as.data.frame(lapply(data, factor))
    xstratmat <- model.matrix(xstratform, model.frame(xstratform,data1,
                              na.action=function(x)x))
    for(i in 1:ncol(xstratmat)) xstratmat[,i] <- i * xstratmat[,i]
    apply(xstratmat,1,sum)
}
names(xStrat) <- NULL
nxStrat <- max(xStrat)

xkey <- NULL
if (!is.null(xstrata)) {
   xkey <- dimnames(xstratmat)[[2]]
   names(xkey) <- 1:length(xkey) 
}

if (!is.null(gamma)) {
  if (is.vector(gamma)) gamma <- matrix(gamma, nrow=length(gamma), ncol=nxStrat)
  else  gamma <- as.matrix(gamma)
      
  if (any(dim(gamma) != c(2, nxStrat))) stop("gamma is in wrong dimention!")
}


# Form the y-strata --------------------------------------------------------

if (!is.null(ystrata)) {
    ystratform <- as.formula(paste("~","-1 + ",paste(ystrata,sep="",collapse=":")))
    data1 <- as.data.frame(lapply(data, factor))
    ystratmat <- model.matrix(ystratform, model.frame(ystratform,data1,
                              na.action=function(y)y))
    for(i in 1:ncol(ystratmat)) ystratmat[,i] <- i * ystratmat[,i]

    yStrat <- apply(ystratmat,1,sum)
    names(yStrat) <- NULL
}
else if (!is.null(retrosamp) && retrosamp=="gamma") 
     stop("ystrata must be provided when retrosamp=='gamma' !")


# Form the y-vector -----------------------------------------------------------------

y <- as.vector(y)
y <- 1 + 1*((yfac<-factor(y))==levels(yfac)[1])
if (max(y)> 2) stop("y must be binary")
y1name <- if(length(levels(yfac))==2) as.character(levels(yfac))[2]
          else paste("Not.", as.character(levels(yfac))[1],sep="")
ykey <- c(as.character(Terms[[2]]),y1name)


# Remove random missing values -------------------------------------------------------

mat <- cbind(y,X,xStrat,ClusInd)
if (!is.null(IntraClus)) mat <- cbind(mat,IntraClus)
miss1 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="uncond"]
miss2 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="retro"]
miss3 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="y|x"]

mat <- cbind(X,xStrat,ClusInd)
if (!is.null(IntraClus)) mat <- cbind(mat,IntraClus)
miss4 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="xonly"]

mat <- cbind(y,xStrat,ClusInd)
if (!is.null(IntraClus)) mat <- cbind(mat,IntraClus)
miss5 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="strata"]

toremove <- c(miss1,miss2,miss3,miss4,miss5)
missReport <- NULL
if(length(toremove)>0){
   #cat("\nObservations deleted due to missing data:\n")
   if ((length(miss1)>0))
       missReport <- rbind(missReport, paste("uncond:",length(miss1),"rows relating to", 
			   sum(w[miss1]),"observations"))
   if ((length(miss2)>0))
       missReport <- rbind(missReport, paste("retro:",length(miss2),"rows relating to",
			   sum(w[miss2]),"observations"))
   if ((length(miss3)>0))
       missReport <- rbind(missReport, paste("y|x:",length(miss3),"rows relating to",
                           sum(w[miss3]),"observations"))
   if ((length(miss4)>0))
       missReport <- rbind(missReport,paste("xonly:",length(miss4),"rows relating to",
			   sum(w[miss4]),"observations"))
   if ((length(miss5)>0))
       missReport <- rbind(missReport,paste("strata:",length(miss5),"rows relating to",
			   sum(w[miss5]),"observations"))
   dimnames(missReport) <- list(rep("",nrow(missReport)),rep("",ncol(missReport)))
   	
   w <- w[-toremove]
   ClusInd <- ClusInd[-toremove]
   if (!is.null(IntraClus)) IntraClus <- IntraClus[-toremove]
   obstype <- factor(as.character(obstype[-toremove]))
   if (!is.null(ystrata)) yStrat <- yStrat[-toremove]
   y <- y[-toremove]
   xStrat <- xStrat[-toremove]
   X <- X[-toremove,,drop=FALSE]
   n <- length(y)
}

clus <- as.numeric(factor(ClusInd))
newclus <- as.numeric(factor(interaction(xStrat, ClusInd)))
if (any(newclus!=clus)) stop("Subjects in same cluster must be in same xstrata !") 


# Stratum Counts Report -------------------------------------------------------------

## when we have clusters of size more than one, expand the data if each row
## indicates more than one individual (with weights>1)
if (length(unique(ClusInd))<n && any(w>1)) {
   ClusInd <- rep(ClusInd, w)
   if (!is.null(IntraClus)) IntraClus <- rep(IntraClus, w)
   obstype <- rep(obstype, w)
   y <- rep(y, w)
   if (!is.null(ystrata)) yStrat <- rep(yStrat,w)
   xStrat <- rep(xStrat, w)
   mw <- matrix(w, dim(X)[1], dim(X)[2])
   colnamesX <- colnames(X)
   X <- matrix(rep(X,mw), sum(w), dim(X)[2])
   colnames(X) <- colnamesX
   n <- length(y)
   w <- rep(1,n)
}

#cat("\nStrata Counts Report:\n")
obstypech <- as.character(obstype)
if (is.null(ystrata)) {
  z1 <- array3form(cbind(y,xStrat,obstypech), ClusInd=ClusInd, IntraClus=IntraClus, 
                 MaxInClus=MaxInClus, rmsingletons=rmsingletons, freqformat=TRUE)

  z1$y <- matrix(as.numeric(z1$XArray[,1,]), dim(z1$XArray)[1], dim(z1$XArray)[3])
  z1$xStrat <- as.numeric(z1$XArray[,-c(1,dim(z1$XArray)[2]),1])
  z1$obstype <- z1$XArray[,dim(z1$XArray)[2],1]

  if (!haveretro) z1$yStrat <- rep(1,dim(z1$y)[1])
  else if (retrosamp=="proband") z1$yStrat <- z1$y[,1]
       else { ## the only choice now is retrosamp=="allcontrol" 
            z1$yStrat <- rep(1, dim(z1$y)[1])
            z1$yStrat[apply(2-z1$y,1,sum,na.rm=TRUE) == 0] <- 2  }
}
else {
  z1 <- array3form(cbind(y,yStrat,xStrat,obstypech), ClusInd=ClusInd, IntraClus=IntraClus, 
                 MaxInClus=MaxInClus, rmsingletons=rmsingletons, freqformat=TRUE)

  z1$y <- matrix(as.numeric(z1$XArray[,1,]), dim(z1$XArray)[1], dim(z1$XArray)[3])
  z1$yStrat <- as.numeric(z1$XArray[,2,1])
  z1$xStrat <- as.numeric(z1$XArray[,-c(1,2,dim(z1$XArray)[2]),1])
  z1$obstype <- z1$XArray[,dim(z1$XArray)[2],1]
}

if (all(z1$Jis==1) && any(w>1)) {
    if (is.null(ystrata)) temp <- Unidata(cbind(obstype,xStrat,y),wts=w)
                else temp <- Unidata(cbind(obstype,xStrat,yStrat,y),wts=w)
    tempwts <- temp[temp[,dim(temp)[2]]>0, dim(temp)[2]]
    z1$counts <- tempwts 
}
z1$counts <- as.numeric(z1$counts)

nyStrat <- length(unique(z1$yStrat))
yStratfac <- factor(z1$yStrat)
xStratfac <- factor(z1$xStrat)
obstypefac <- factor(as.character(z1$obstype))
Clusize <- factor(z1$Jis)

StrReport <- xStrReport <- ClusReport <- NULL
sub <- (1:dim(z1$XArray)[1])[z1$obstype=="xonly"]
if(length(sub)>0){
   df <- data.frame(wts=z1$counts, obstype=obstypefac, yStrat=yStratfac,
                    xStrat=xStratfac)[-sub,]
   df1 <- data.frame(wts=z1$counts, obstype=obstypefac, yStrat=yStratfac,
                    xStrat=xStratfac)[sub,]   
}
else df <- data.frame(wts=z1$counts, obstype=obstypefac, yStrat=yStratfac,
                    xStrat=xStratfac)
cdf <- data.frame(wts=z1$counts, obstype=obstypefac, Clusize=Clusize, xStrat=xStratfac)

if (!prosp) StrReport <- ftable(xtabs(wts~obstype+yStrat+xStrat,data=df))
   else StrReport <- ftable(xtabs(wts~obstype+xStrat,data=df))
if(exists("df1")) xStrReport <-  xtabs(wts~obstype+xStrat,data=df1)

ClusReport <- ftable(xtabs(wts~obstype+Clusize+xStrat,data=cdf))

if (xs.includes) {
    oblevel0 <- levels(obstypefac)
    oblevel <- rep(oblevel0, rep(length(unique(yStratfac)),length(oblevel0)))  
    coblevel <- rep(oblevel0, rep(length(unique(Clusize)), length(oblevel0)))

    report <- ftable(xtabs(wts~obstype+yStrat+xStrat,data=df[df$obstype=="strata",]))
    strep <- xtabs(wts~yStrat+xStrat,data=df[df$obstype!="strata",])
    report[oblevel=="strata"] <- strep
    StrReport <- StrReport-report

    report <- ftable(xtabs(wts~obstype+Clusize+xStrat,data=cdf[df$obstype=="strata",]))
    strep <- xtabs(wts~Clusize+xStrat,data=cdf[df$obstype!="strata",])
    report[coblevel=="strata"] <- strep
    ClusReport <- ClusReport-report

    if (exists("df1")) {
        report1 <- xtabs(wts~obstype+xStrat,data=df1[df1$obstype=="strata",])
        strep1 <- xtabs(wts~xStrat,data=df1[df1$obstype!="strata",])
        report1[oblevel0=="strata"] <- strep1
        xStrReport <- xStrReport-report1
    }
}


# RETURN AT THIS POINT IF FIT=FALSE ------------------------------------------------------

if (!fit) {
    ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, 
                ClusReport=ClusReport, xkey=xkey, ykey=ykey, fit=fit, call=call, 
                assign1=assign1, fnames=fnames, terms1=terms1, order1=order1)
    class(ans) <- "rclusbin"
    return(ans) 
}


# Data reformating ---------------------------------------------------------

X <- as.matrix(X)
if (length(y) != dim(X)[1]) stop("y and X dimensions do not match")
inter2 <- rep(1,dim(X)[1]) ##match dim[2] of XArray to length(theta)
probandonly <- FALSE

id <- (1:n)[obstype != "strata"]
if (is.null(ystrata)) {
   z <- array3form(cbind(y,inter2,X,xStrat)[id,,drop=FALSE], ClusInd=ClusInd[id], 
          IntraClus=IntraClus[id], MaxInClus=MaxInClus, rmsingletons=rmsingletons, 
          freqformat=TRUE)

   z$y <- matrix(z$XArray[,1,], dim(z$XArray)[1], dim(z$XArray)[3])
   z$xStrat <- z$XArray[,dim(z$XArray)[2],1]
   z$XArray <- z$XArray[ ,-c(1,dim(z$XArray)[2]), , drop=FALSE]

   if (!haveretro) z$yStrat <- rep(1,dim(z$y)[1])
   else if (retrosamp=="proband") z$yStrat <- z$y[,1]
       else { ## the only choice now is retrosamp=="allcontrol" 
            z$yStrat <- rep(1, dim(z$y)[1])
            z$yStrat[apply(2-z$y,1,sum,na.rm=TRUE) == 0] <- 2  }
}
else  {
   z <- array3form(cbind(y,inter2,X,xStrat,yStrat)[id,,drop=FALSE], ClusInd=ClusInd[id], 
                IntraClus=IntraClus[id], MaxInClus=MaxInClus, rmsingletons=rmsingletons, 
                freqformat=TRUE)

   z$y <- matrix(z$XArray[,1,], dim(z$XArray)[1], dim(z$XArray)[3])
   z$xStrat <- z$XArray[,dim(z$XArray)[2]-1,1]
   z$yStrat <- z$XArray[,dim(z$XArray)[2],1]
   z$XArray <- z$XArray[ ,2:(dim(z$XArray)[2]-2), , drop=FALSE]
}

if(all(z$Jis==1)) {
  probandonly <- TRUE
  if(any(w[id]>1)) {
   if (is.null(ystrata))
       temp <- Unidata(cbind(xStrat,X[,dim(X)[2]:1],y)[id,,drop=FALSE],wts=w[id])
   else temp <- Unidata(cbind(yStrat,xStrat,X[,dim(X)[2]:1],y)[id,,drop=FALSE], 
                        wts=w[id])
   tempwts <- temp[temp[,dim(temp)[2]]>0, dim(temp)[2]]
   z$counts <- tempwts
  } 
}

z$Jis <- cbind(z$Jis, z$xStrat)
z$yStratOrg <- z$yStrat
z$yStrat <- cbind(z$yStrat, z$xStrat)

 
# Calculate rmat and Qmat ------------------------------------------

### Check provided NMat ###
if (haveNMat) {
   NMat <- as.matrix(NMat)
   if (any(dim(NMat) != c(nyStrat,nxStrat)) || min(NMat) <= 0)
     stop("Wrong NMat input !")
}

if (!prosp) {
   ### Calculate rmat ###
   if (haveNMat) {
     ind <- (1:dim(z1$XArray)[1])[z1$obstype=="uncond" | z1$obstype=="retro" ]
     nMat0 <- xtabs(wts~yStrat+xStrat,data=data.frame(wts=z1$counts, 
                    yStrat=yStratfac, xStrat=xStratfac)[ind,])
     rmat <- NMat - nMat0
     dimnames(rmat) <- dimnames(nMat0)[1:2]
   }
   else {
     rk <- rep(0,dim(z1$XArray)[1])
     if ((havestrata & xs.includes) | !havestrata) 
       rk[z1$obstype=="retro"] <- -1 * z1$counts[z1$obstype=="retro"]
     if  (havestrata & xs.includes) 
       rk[z1$obstype=="uncond"] <- -1 * z1$counts[z1$obstype=="uncond"]
     rk[z1$obstype=="strata"] <- 1 * z1$counts[z1$obstype=="strata"]

     rmat <- xtabs(rk~yStrat+xStrat, data=data.frame(rk=rk, yStrat=yStratfac,
                   xStrat=xStratfac))
     rmat <- matrix(rmat,dim(rmat)[1],dim(rmat)[2],dimnames=dimnames(rmat)[1:2])
   }

   if (dim(rmat)[1] == 1) { # only one level of y represented in the data
     if (dimnames(rmat)[[1]]=="1") rmat <- rbind(rmat,rep(0,dim(rmat)[2]))
     else rmat <- rbind(rep(0,dim(rmat)[2]),rmat)
   }

   ### Calculate Qmat ###
   if (is.null(Qstart)){
     if (!haveNMat) {
       if (havestrata) {
         ind <- if (xs.includes) (1:dim(z1$XArray)[1])[z1$obstype=="strata"]
           else (1:dim(z1$XArray)[1])[z1$obstype=="uncond" | z1$obstype=="retro" | 
                                      z1$obstype=="strata"]
         NMat <- xtabs(wts~yStrat+xStrat,data=data.frame(wts=z1$counts, 
                      yStrat=yStratfac, xStrat=xStratfac)[ind,])
       }
       else if (retrosamp=="proband" && any(z1$Jis>1)) {
         ind <- (1:n)[obstype=="uncond" | obstype=="retro"]
         NMat <- xtabs(wts~y+xStrat,data=data.frame(wts=rep(1,n),  
                      y=factor(y), xStrat=factor(xStrat))[ind,])
       }
       else stop("Qstart needs to be provided !")
     }
     
     if (is.matrix(NMat) && all(dim(NMat) == c(nyStrat,nxStrat)) && min(NMat)>0)
       	   Qmat <- sweep(NMat,2,apply(NMat,2,sum),"/")
     else  stop("Attempt to construct missing Qstart failed")
   }
   else { # Qstart is present
        if (is.vector(Qstart)) 
            Qstart <- matrix(Qstart, length(Qstart), nxStrat)
        if (is.matrix(Qstart) && (dim(Qstart)==c(nyStrat,nxStrat)) &&
           (round(apply(Qstart,2,sum),5)==rep(1,dim(Qstart)[2])))
              Qmat <- Qstart
        else stop("Illegal Qstart value !")  }
}


# Creat start if needed ----------------------------------------------

if (is.null(start)){
  haveretro <- length(obstype[obstype=="retro"])>0

  if (haveretro) {
    nMat <- xtabs(wts~yStrat+xStrat,data=data.frame(wts=z1$counts, 
                  yStrat=yStratfac, xStrat=xStratfac)[z1$obstype=="retro",])

    if (!haveNMat) {
       if (havestrata) {
         ind <- if (xs.includes) (1:dim(z1$XArray)[1])[z1$obstype=="strata"]
           else (1:dim(z1$XArray)[1])[z1$obstype=="uncond" | z1$obstype=="retro" | 
                                      z1$obstype=="strata"]
         NMat <- xtabs(wts~yStrat+xStrat,data=data.frame(wts=z1$counts, 
                     yStrat=yStratfac, xStrat=xStratfac)[ind,])
       }
       else if (retrosamp=="proband" && any(z1$Jis>1)) {
         ind <- (1:n)[obstype=="uncond" | obstype=="retro"]
         NMat <- xtabs(wts~y+xStrat,data=data.frame(wts=rep(1,n),  
                     y=factor(y), xStrat=factor(xStrat))[ind,])
       }
       else stop("Cannot construct starting values from this data !")
    }

    ## calculate weights using the n's and N's
    Nn <- (NMat+1)/(nMat+1)
    Nn <- sweep(Nn,2,apply(Nn,2,FUN=function(x)mean(x,na.rm=TRUE)),"/")
    Nnvec <- as.vector(Nn)
    xystrat <- (z$xStrat-1)*2 + z$yStratOrg
    wmult <- Nnvec[xystrat] #multiply by corresponding values in Nn
    if (any(is.na(wmult))) stop("Missing values in constructed wmult vector")
  }
  else wmult <- rep(1,dim(z$y)[1])

  # Calculate starting values from probands only
  # start <- glm(z$y[,1]==1 ~ z$XArray[ , -1, 1]-1, weights=z$counts*wmult, 
  #              family=eval(parse(text=paste("binomial(link=",linkname,
  #              ")",sep=""))))$coefficients
  # start <- c(start,0)

  # Calculate starting values from all individuals
  # An isproband term is added to account for higher effective sampling rate of probands
  #   if retrosamp="proband" and probandonly=FALSE
  # Use log(0.5) to start the log-scale parameter as necessary 

  ## Turn z$XArray back into a long matrix (1 row per individual)
  dimXA <- dim(z$XArray[,-1,,drop=FALSE])
  XMatLong <- matrix(aperm(z$XArray[,-1,,drop=FALSE],c(3,1,2)), nrow=dimXA[1]*dimXA[3], 
                     ncol=dimXA[2],byrow=FALSE)
  isproband <- matrix(0,nrow=dimXA[1],ncol=dimXA[3])
  isproband[,1] <- 1
  temcounts <- (z$counts*wmult)*matrix(1,nrow=dimXA[1],ncol=dimXA[3])
  XMatLong <- cbind(as.vector(t(temcounts)), as.vector(t(z$y)),
                    as.vector(t(isproband)), XMatLong)

  ## Drop off the rows of NAs
  XMatLong <- matrix(na.omit(XMatLong),ncol=ncol(XMatLong))
  colnames(XMatLong) <- c("counts", "y-var", "isproband", colnames(X))
   
  ## Check whether the model contains a dummy variable for proband status
  probin <- FALSE
  probcheck <- (XMatLong[,3]==XMatLong[,-c(1:3)])
  if (any(apply(probcheck,2,sum)==dim(XMatLong)[1])) {
      XMatLong <- XMatLong[,-3]
      probin <- TRUE  }

  ## Delete isproband if 'retrosamp!=proband' or 'probandonly==TRUE'
  if(!probin) {
    if (is.null(retrosamp) || (!is.null(retrosamp) && retrosamp!="proband") || 
        probandonly)  {XMatLong <- XMatLong[,-3]; probin <- TRUE }
  }
    
  start <- glm(XMatLong[,2]==1 ~ XMatLong[,-c(1,2)] - 1, weights=XMatLong[,1], 
               family=eval(parse(text=paste("binomial(link=",linkname,")",
               sep=""))))$coefficients
  if (!probin) start <- start[-1] ## exclude value for isproband

  ## Calculate starting value of sigma 
  sigma0 <- ifelse(!is.null(sigma), sigma, 0.5)
  
  if (is.null(sigma) && !probandonly) {
    current <- FALSE # we don't consider the following at current!
    if (current) {
      scounts <- (1:length(z$counts))[z$Jis[,1]>=2] #choose clusters with size>=2
      sn <- length(scounts)
      if(sn>=2) {#choose the first two individuals in each cluster & expand with counts 
         sy <- (z$y[scounts,])[rep(1:sn,z$counts[scounts]),1:2]
         sx0 <- z$XArray[scounts,-1,1:2,drop=FALSE]
         sx <- sx0[rep(1:sn,z$counts[scounts]),,]

         E1 <- mean(matrix(sx[,,1],ncol=dim(sx)[2]) %*% start)
         E2 <- mean(matrix(sx[,,2],ncol=dim(sx)[2]) %*% start)
         pq1 <- exp(E1)/(1+exp(E1))^2
         pq2 <- exp(E2)/(1+exp(E2))^2
         sigma0 <- sqrt(cov(sy[,1],sy[,2])/(pq1*pq2))
     }
    }
  }    

  start <- c(start, log(sigma0))
  #print(start)
}
else if (length(start) != (dim(X)[2]+1)){
    ## Match up any values in start with names in X, fill in holes with 0's.
    newstart <- rep(0,dim(X)[2]+1)
    posns <- match(names(start),c(dimnames(X)[[2]],"logscale"))   
    newstart[posns] <- start
    start <- newstart
}


# --------------- ESTIMATE GAMMA[1] (gamma0)  FROM DATA ----------------------
gamma0 <- NULL

if (!is.null(retrosamp) && retrosamp=="gamma"){
     if (havestrata) {## there is no NMat provided in the call!
         ind <- if (xs.includes) (1:dim(z1$XArray)[1])[z1$obstype=="strata"]
           else (1:dim(z1$XArray)[1])[z1$obstype=="uncond" | z1$obstype=="retro" | 
                                      z1$obstype=="strata"] 
     }
     else ind <- 1:dim(z1$XArray)[1]
     
     y1Strat <- rep(2, dim(z1$y)[1])
     y1Strat[apply(2-z1$y,1,sum,na.rm=TRUE) == 1] <- 1
     y1Stratfac <- factor(y1Strat)
     Mat1 <- xtabs(wts~y1Strat+xStrat,data=data.frame(wts=z1$counts, 
                   y1Strat=y1Stratfac, xStrat=xStratfac)[ind,])

     y2Strat <- rep(2, dim(z1$y)[1])
     y2Strat[apply(2-z1$y,1,sum,na.rm=TRUE)==1 & z1$yStrat==1] <- 1
     y2Stratfac <- factor(y2Strat)
     Mat2 <- xtabs(wts~y2Strat+xStrat,data=data.frame(wts=z1$counts, 
                   y2Strat=y2Stratfac, xStrat=xStratfac)[ind,])
 
     if (havestrata) gamma0 <- Mat2[1,]/Mat1[1,]
     else if (haveNMat) {
         Mat <- xtabs(wts~yStrat+xStrat,data=data.frame(wts=z1$counts, 
                      yStrat=yStratfac, xStrat=xStratfac))

         nom <- Mat2[1,]*NMat[1,]/Mat[1,]
         dem <- nom + (Mat1[1,]-Mat2[1,])*NMat[2,]/Mat[2,]
    	 gamma0 <- nom/dem
     }
     #else warning("Cannot construct gamma0 from this data !")
}


# ------------------------- Function call  ---------------------------

if (linkname=="logit") modelfn <- binlogistic 
  else if (linkname=="probit") modelfn <- binprobit
       else modelfn <- bincloglog

if (haveretro) {
   if (retrosamp=="proband") RStratModInf <- RProbandStratModInf 
   else if (retrosamp=="allcontrol") RStratModInf <- R0StratModInf 
        else RStratModInf <- R1StratModInf
}

if (!prosp) {
  if (devcheck) {
    divchk(theta=start, loglkfn=MLInf, nderivs=2, y=z$y, x=z$XArray, nzval0=nzval0,
         ProspModInf=RClusProspModInf, StratModInf=RStratModInf, modelfn=modelfn,
         Aposn=1:nrow(z$y), Acounts=z$counts, Bposn=1:nrow(z$y), 
         Bcounts=z$counts, rmat=rmat, Qmat=Qmat, xStrat=as.numeric(z$xStrat), 
         Jis=z$Jis, gamma=gamma, yStrat=z$yStrat)
    stop("Derivatives check on MLInf")
  }

 res <-  mlefn(theta=start, loglkfn=MLInf, y=z$y, x=z$XArray, nzval0=nzval0,
         ProspModInf=RClusProspModInf, StratModInf=RStratModInf, modelfn=modelfn,
         Aposn=1:nrow(z$y), Acounts=z$counts, Bposn=1:nrow(z$y), 
         Bcounts=z$counts, rmat=rmat, Qmat=Qmat, xStrat=as.numeric(z$xStrat), 
         Jis=z$Jis, gamma=gamma, yStrat=z$yStrat,
         control=control, control.inner=control.inner)  
}
else { 
 #res <- mlefn(theta=start, loglk=ProspInf, y=z$y, x=z$XArray,
 #             ProspModInf=RClusProspModInf, modelfn=modelfn, Aposn=1:nrow(z$y),
 #             Acounts=z$counts, xStrat=as.numeric(z$xStrat), Jis=z$Jis,
 #             control=control, control.inner=control.inner)

 ## the following gives identical results as above here:
 res <- mlefn(theta=start, loglkfn=RClusProspModInf, y=z$y, x=z$XArray, 
              wts=z$counts, modelfn=modelfn, Jis=z$Jis, gamma=NULL, yStrat=NULL, 
              control=control, control.inner=control.inner)  
}

names(res$theta) <- names(res$score) <- c(dimnames(X)[[2]],"logsigma")
dimnames(res$inf) <- list(names(res$theta),names(res$theta))
covmat <- solve(res$inf); dimnames(covmat) <- dimnames(res$inf)
dd <- sqrt(diag(covmat))
correlation <- covmat/outer(dd, dd)
#eta <- XOrig %*% res$theta[-length(res$theta)]
#sigma <- exp(res$theta[length(res$theta)])

# ------------------------------------------------------------------------------------
ans<-list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, 
          ClusReport=ClusReport, xkey=xkey, ykey=ykey, fit=fit, error=res$error, 
          coefficients=res$theta, loglk=res$loglk, score=res$score, inf=res$inf, 
          cov=covmat, cor=correlation, Qmat=res$extra$Qmat, gamma0=gamma0, 
          call=call, assign1=assign1, fnames=fnames, terms1=terms1, order1=order1)
class(ans)<-"rclusbin"
ans
}


#####################################################################################
print.rclusbin<-function(x,digits=max(3, getOption("digits") - 3), ...){
#####################################################################################
cat("\nCall:\n", deparse(x$call), "\n\n", sep=" ")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}

cat("\nCluster Size Report:\n")
print(x$ClusReport,quote = FALSE, row.names=FALSE)
cat("\nStratum Counts Report:\n")
print(x$StrReport, row.names=FALSE)
if (!is.null(x$xStrReport)) {
    cat("\nObservations of obstype==xonly\n")
    print(x$xStrReport)
}

if (!is.null(x$xkey)) {
   cat("\nKey to x-Strat:\n")
   print.default(x$xkey, quote=FALSE, row.names=FALSE)
}

cat("\n\nModel for prob of ", x$ykey[1],"=",x$ykey[2]," (y=1) given covariates\n",sep="")

if (x$fit) {
  if (x$error != 0) {cat("\n"); print("WARNING: FIT UNSUCCESSFUL") }
  cat("\nCoefficients:\n")
  print.default(format(coef(x), digits=digits), print.gap=2, quote = FALSE)
  cat("\nloglikelihood =",x$loglk," using",length(coef(x)),"parameters\n\n")

  #cat("\n")
  #cat("Covariance:\n")
  #print.default(format(x$cov, digits = digits), print.gap = 2, quote = FALSE)
  #cat("\n")
  #cat("Correlation:\n")
  #print.default(format(x$cor, digits = digits), print.gap = 2, quote = FALSE)
  #cat("\n")
}

invisible(x)
}


#####################################################################################
summary.rclusbin<-function(object,...){
#####################################################################################
z <- object
if (!z$fit) stop("No summary provided when fit=FALSE !")

assign1 <- z$assign1;
terms1 <- z$terms1;
order1 <- z$order1;
fnames <- z$fnames;

# find dfs for each variable
num <- -1
pos <- waldchis <- pvals <- coef.table1 <- NULL

for (i in 1:length(assign1)){
   if (assign1[i] > num) {
      pos<-c(pos, i)
      num <- assign1[i]
   }
}
pos1 <- c(pos[-1], length(assign1) +1)
values <- assign1[pos]

dfs <- lengths <- pos1 - pos

# do wald stats
if (length(pos) > 1) {
   for (i in 2: length(pos)){
      thetai <- z$coefficients[pos[i]:(pos[i] + lengths[i] -1)]
      covi <- z$cov[pos[i]:(pos[i] + lengths[i] - 1), pos[i]:(pos[i] + lengths[i] -1)]
      waldchi <- t(thetai) %*% solve(covi, thetai)
      waldchis <- c(waldchis, waldchi)
      pval <- (1 - pchisq(waldchi, df = dfs[i]))
      pvals <- c(pvals, pval)
   }
   coef.table1 <- cbind(dfs[-1], waldchis, pvals)
   dimnames(coef.table1)[[1]]<-terms1[values]
   dimnames(coef.table1)[[2]]<-c("Df", "Chi", "Pr(>Chi)")
}
#calculate t statistics for non factors
if(is.null(dim(z$cov))) s.err<-sqrt(z$cov) else s.err <- sqrt(diag(z$cov))
zstats <- z$coefficients/s.err
pvals <- (1-pnorm(abs(zstats)))*2
coef.table2<-cbind(z$coefficients, s.err, zstats, pvals)
dimnames(coef.table2)[[2]]<- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
if (assign1[1]==0) dimnames(coef.table2)[[1]]<-c("(Intercept)", terms1[values],"logsigma") 
else dimnames(coef.table2)[[1]]<-terms1[values]

ans<-list(missReprot=z$missReport, StrReport=z$StrReport, xStrReport=z$xStrReport,
          ClusReport=z$ClusReport, xkey=z$xkey, ykey=z$ykey, fit=z$fit, 
          error=z$error, coefficients=z$coefficients, loglk=z$loglk, call=z$call, 
          coef.table1=coef.table1, coef.table2=coef.table2)
class(ans)<-"summary.rclusbin"
ans
}


#####################################################################################
print.summary.rclusbin<-function(x, digits = max(3, getOption("digits") - 3), ...){
#####################################################################################
cat("\nCall:\n", deparse(x$call), "\n\n", sep=" ")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}

cat("\nCluster Size Report:\n")
print(x$ClusReport, row.names=FALSE)
cat("\nStratum Counts Report:\n")
print(x$StrReport, row.names=FALSE)
if (!is.null(x$xStrReport)) {
    cat("\nObservations of obstype==xonly\n")
    print(x$xStrReport)
}

if (!is.null(x$xkey)) {
   cat("\nKey to x-Strat:\n")
   print.default(x$xkey, quote=FALSE, row.names=FALSE)
}

cat("\n\nModel for prob of ", x$ykey[1],"=",x$ykey[2]," (y=1) given covariates\n",sep="")
if (x$error != 0) {cat("\n"); print("WARNING: FIT UNSUCCESSFUL") }
cat("\nloglikelihood =",x$loglk," using",length(coef(x)),"parameters\n\n")
if(!is.null(x$coef.table1)){
    cat("\n")
    cat("Wald Tests:\n")
    print(x$coef.table1, digits=4)
}
if(!is.null(x$coef.table2)){
    cat("\n")
    cat("Coefficients:\n")
    print(x$coef.table2, digits=4)
}

invisible(x)
}

