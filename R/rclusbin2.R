#############################################################################################
rclusbin2 <- function(formula1, formula2, weights=NULL, ClusInd.name=NULL, IntraClus.name=NULL, 
		yname, xstrata=NULL, ystrata, obstype.name="obstype", data, NMat=NULL, 
		xs.includes=FALSE, MaxInClus=NULL, rmsingletons=FALSE, retrosamp=TRUE, nzval0=20,
		fit=TRUE, devcheck=FALSE, linkname="logit", start=NULL, Qstart=NULL, sigma=NULL, 
                paruse="xis", control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
#############################################################################################
{
# -- When we have clusters of size > 1, the "weights" give within-cluster counts which
#    is the total number of individuals with identical variable values in a cluster;
# -- When all clusters are of size 1, the "weights" give cluster counts which is the 
#    total number of clusters with identical variable values. 
# -- retrosamp=TRUE should be forced here;
# -- ystrata should be the response variable for formula2 but in different dimentions;
# -- formula2 variables ends with sequential numbers indicating their within Cluster ID;

# -- devcheck: a temporary input to check derivatives


# Process formula -------------------------------------------------------------------

mf <- call <- match.call()
mf[[1]] <- as.name("model.frame")
names(mf)[1] <- "model"
mf$na.action <- as.name("na.keep")
NULL -> mf$ClusInd.name -> mf$IntraClus.name -> mf$yname
NULL -> mf$xstrata -> mf$ystrata -> mf$obstype.name -> mf$NMat -> mf$xs.includes 
NULL -> mf$MaxInClus -> mf$rmsingletons -> mf$retrosamp -> mf$nzval0
NULL-> mf$fit -> mf$linkname -> mf$start -> mf$Qstart -> mf$sigma 
NULL -> mf$control -> mf$control.inner -> mf$devcheck 


### First formula processing

mf1 <- mf[c("model","formula1","data","weights")]
mf1$na.action <- as.name("na.keep")
#resp1 <- as.character(mf1$formula[2])
resp1 <- yname
names(mf1)[2] <- "formula"
mf1 <- eval(mf1, sys.frame(sys.parent()))
y <- model.extract(mf1,"response")
names(y) <- NULL
n <- if (is.matrix(y)) dim(y)[1] else length(y)
Terms <- attr(mf1,"terms")
X <- model.matrix(Terms,mf1)
if(!is.matrix(X)) X <- matrix(X)
XOrig <- X
w <- model.extract(mf1,"weights")
if (is.null(w)) w <- rep(1,n)

terms1 <- attr(terms(formula1),"term.labels")
order1 <- attr(terms(formula1),"order")
assign1 <- attr(X,"assign"); fnames1 <- names(attr(X,"contrasts"))


if (is.null(ClusInd.name)) ClusInd <- 1:n
  else ClusInd <- data[,match(ClusInd.name,names(data))]
if (is.null(IntraClus.name)) IntraClus <- NULL
  else IntraClus <- data[,match(IntraClus.name,names(data))]



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
if (!haveretro) stop("A retrosamp scheme should be provided !")

if (any(is.na(w))) stop("Weights should not contain NAs !")
if (!is.null(sigma) && (length(sigma)>1 || sigma<=0))
      stop("Illegal sigma value !")

datnmiss <- data


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
else stop("ystrata should be provided !")

#print(ftable(yStrat~data$wave))  ### check y-strata


# Form the y-vector -----------------------------------------------------------------

y <- as.vector(y)
y <- 1 + 1*((yfac<-factor(y))==levels(yfac)[1])
if (max(y,na.rm=TRUE)> 2) stop("y must be binary")
y1name <- if(length(levels(yfac))==2) as.character(levels(yfac))[2]
          else paste("Not.", as.character(levels(yfac))[1],sep="")
ykey <- c(as.character(Terms[[2]]),y1name)

#print(ftable(y~data$wave))  ### check y-vector


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
   datnmiss <- data[-toremove,]  ###
   n <- length(y)
}

clus <- as.numeric(factor(ClusInd))
newclus <- as.numeric(factor(interaction(xStrat, ClusInd)))
if (any(newclus!=clus)) stop("Subjects in same cluster must be in same xstrata !") 


### Second formula processing

# out1 includes (XArray, Jis, counts) as outputs
out1 <- array3form2(datnmiss, ClusInd.name=ClusInd.name, IntraClus.name=IntraClus.name, 
                    MaxInClus=MaxInClus, rmsingletons=rmsingletons, freqformat=FALSE)

out2 <- aperm(out1$XArray,c(1,3,2))
name1 <- dimnames(out2)[[1]]
name2 <- dimnames(out2)[[2]]
name3 <- dimnames(out2)[[3]]

data0 <- matrix(out2,dim(out2)[1],dim(out2)[2]*dim(out2)[3])
data2 <- data.frame(data0)

vclass <- rep(NA, dim(out1$newdat)[2])
for (i in 1:length(vclass)) vclass[i] <- class(out1$newdat[,i])
vclass2 <- rep(vclass,rep(length(name2),length(name3)))

for (j in 1:dim(data2)[2])
  if ((vclass2[j] %in% c("character","factor"))==FALSE) data2[,j] <- as.numeric(data0[,j])

dimnames(data2) <- list(name1, paste(rep(name3,rep(length(name2),length(name3))), rep(name2,length(name3)),sep="."))
data2$counts <- out1$counts
obstype2 <- data2[,paste("obstype",name2[1],sep=".")]

model2 <- model.apply(formula2,data=data2,weights=counts)
y2 <- model2$y
X2 <- model2$X
#w2 <- model2$weights
#x2array <- array(X2,c(nrow(X2),ncol(X2),1))

assign2 <- model2$assign1
fnames2 <- model2$fnames1 
terms2 <- model2$terms1 
order2 <- model2$order1

data2_y1 <- data2_y2 <- data2
yid <- match(paste(resp1,name2[1],sep="."), names(data2))
data2_y1[,yid] <- rep(1, dim(data2)[1])
data2_y2[,yid] <- rep(0, dim(data2)[1])

X2_y1 <- model.apply(formula2,data=data2_y1,weights=counts)$X
X2_y2 <- model.apply(formula2,data=data2_y2,weights=counts)$X

#print(summary(factor(y2))) ## check y2
#print(names(data2))
#print(X2[1:10,])
#print(X2_y1[1:10,])
#print(X2_y2[1:10,])


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
  else if (retrosamp) z1$yStrat <- z1$y[,1]
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
                assign1=assign1, fnames1=fnames1, terms1=terms1, order1=order1,
		assign2=assign2, fnames2=fnames2, terms2=terms2, order2=order2)
    class(ans) <- "rclusbin2"
    return(ans) 
}


# Data reformating ---------------------------------------------------------

X <- as.matrix(X)
if (length(y) != dim(X)[1]) stop("y and X dimensions do not match")
inter2 <- rep(1,dim(X)[1]) ##match dim[2] of XArray to length(theta)

id <- (1:n)[obstype != "strata"]
z <- array3form(cbind(y,inter2,X,xStrat,yStrat)[id,,drop=FALSE], 
		ClusInd=ClusInd[id], IntraClus=IntraClus[id],
		MaxInClus=MaxInClus, rmsingletons=rmsingletons, freqformat=FALSE)

z$y <- matrix(z$XArray[,1,], dim(z$XArray)[1], dim(z$XArray)[3])
z$xStrat <- z$XArray[,dim(z$XArray)[2]-1,1]
z$yStrat <- z$XArray[,dim(z$XArray)[2],1]
z$XArray <- z$XArray[ ,2:(dim(z$XArray)[2]-2), , drop=FALSE]

probandonly <- FALSE
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
       else if (retrosamp && any(z1$Jis>1)) {
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
           (round(apply(Qstart,2,sum),0)==rep(1,dim(Qstart)[2])))
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
       else if (retrosamp && any(z1$Jis>1)) {
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
  #   if retrosamp and probandonly=FALSE
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

  ## Delete isproband if 'retrosamp' or 'probandonly==TRUE'
  if(!probin) {
    if (probandonly)  {XMatLong <- XMatLong[,-3]; probin <- TRUE }
  }
    
  start1 <- glm(XMatLong[,2]==1 ~ XMatLong[,-c(1,2)] - 1, weights=XMatLong[,1], 
               family=eval(parse(text=paste("binomial(link=",linkname,")",
               sep=""))))$coefficients
  if (!probin) start1 <- start1[-1] ## exclude value for isproband

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

  start1 <- c(start1, log(sigma0))
  #print(start1)

     offset <- rep(0,length(y2))
     if (any(obstype2 %in% "retro")) {
        	nMat <- xtabs(w2~y2+xStrat2,data=data.frame(y2,xStrat2=z$xStrat,w2=z$counts,obstype2)[obstype2=="retro",]) 
        	nMat <- matrix(nMat,dim(nMat)[1],dim(nMat)[2],dimnames=dimnames(nMat)[1:2])
        	offsets0 <- log( (nMat[1,]/nMat[2,])/(Qmat[1,]/Qmat[2,]))
        	offset <- ifelse((obstype2=="retro"), offsets0[z$xStrat], 0)
     }
        
     if (linkname=="probit") offset <- 0.43*offset
     #WHAT HERE??????  if (linkname=="cloglog") offset <- ??????
     
     ind2 <- (1:length(y2))[obstype2=="uncond" | obstype2=="retro" | obstype2=="y|x"]
     xarray2 <- array(X2,c(dim(X2),1))
     if (length(levels(factor(y2[ind2])))<2 | sum(z$counts[ind2]) < 50)
     	        stop("Cannot construct starting values from this data")
     templink <- if (linkname=="cloglog") "logit" else linkname
     # THERE IS AN R BUG here with glm with cloglog so use logit instead
     start2 <- glm((y2==1)[ind2]~xarray2[ind2, ,1]-1,weights=z$counts[ind2],
                   family=eval(parse(text=paste("binomial(link=", templink,")",sep=""))),
                   offset=offset[ind2])$coefficients

     start <- c(start1, start2)
     names(start) <- c(dimnames(X)[[2]], "log(sigma)", dimnames(X2)[[2]])

}
else { # Match up any values in start with names in Xs, fill in holes with 0's.
   newstart <- rep(0,1+dim(X)[2]+dim(X2)[2])
   if (length(start) != length(newstart)) {
 	vname <- c(dimnames(X)[[2]],dimnames(X2)[[2]])
   	posns <- match(names(start),vname)
        newstart[posns] <- start
   	start <- newstart
   }
}


# Form X array ----------------------------------------------

X2 <- X2[obstype2 != "strata",,drop=FALSE] 
X2_y1 <- X2_y1[obstype2 != "strata",,drop=FALSE] 
X2_y2 <- X2_y2[obstype2 != "strata",,drop=FALSE] 

#print(dim(z$XArray))
#print(dim(X2))

npar <- c(dim(z$XArray)[2], dim(X2)[2])
xarray <- xarray_y1 <- xarray_y2 <- array(0,c(dim(z$XArray)[1],sum(npar),dim(z$XArray)[3]))
dimnames(xarray) <- list(dimnames(z$XArray)[[1]], c(dimnames(z$XArray)[[2]], dimnames(X2)[[2]]),
			 dimnames(z$XArray)[[3]])
dimnames(xarray_y1) <- dimnames(xarray_y2) <- dimnames(xarray)

xarray[,1:npar[1],] <- z$XArray
xarray[,(npar[1]+1):sum(npar),1] <- X2
xarray_y1[,(npar[1]+1):sum(npar),1] <- X2_y1
xarray_y2[,(npar[1]+1):sum(npar),1] <- X2_y2
xlist <- list(xarray, xarray_y1, xarray_y2, as.numeric(z$xStrat))
names(xlist) <- c("orig","1","2","xstrata")

#print(z$XArray[1:6,,])
#print(head(X2))
#print(xarray[1:6,,])

#cat("\n Dimention of xarray: \n")
#print(dim(xarray))
#print(xarray[1:24,,])
#cat("\n Response of interest: \n")
#print(dim(z$y))
#print(z$y[1:24,])
#cat("\n Cluster size: \n")
#print(z$Jis[1:24,])
#cat("\n Case-control response: \n")
#print(length(y2))
#print(y2)
#cat("\n\n")


# ------------------------- Function call  ---------------------------

if (linkname=="logit") modelfn<-binlogistic
else if (linkname=="probit") modelfn<-binprobit
else if (linkname=="cloglog") modelfn<-bincloglog
else stop(paste("linkname",linkname,"not implemented"))


### NOTE: this code can only be tested without xStrat;
#divchk(theta=start, loglkfn=MEtaProspModInf2, nderivs=2, y=cbind(z$y,y2), 
#	x=xarray, npar=npar, modelfn=modelfn, Jis=z$Jis, yStrat=z$yStrat, inxStrat=1)
#stop("Yannan: test MEtaProspModInf2 (no xStrat)")

#cat("starting values --\n")
#print(start)
#cat("\nrmat --\n")
#print(rmat)
#cat("\nQmat--\n")
#print(Qmat)
#cat("\n\n")


if (devcheck) {
  divchk(theta=start, loglkfn=MLInf, nderivs=2, ProspModInf=MEtaProspModInf2, modelfn=modelfn, 
	 StratModInf=MEtaStratModInf2, y=cbind(z$y,y2), x=xarray,  npar=npar, nzval0=nzval0,
	 Aposn=1:nrow(z$y), Acounts=z$counts, Bposn=1:nrow(z$y), Bcounts=z$counts, 
         rmat=rmat, Qmat=Qmat, xStrat=as.numeric(z$xStrat), ymat=cbind(z$y,y2), xlist=xlist,
         Jis=z$Jis, yStrat=z$yStrat, paruse=paruse)
  stop("Derivatives check on MLInf")
}

res <-  mlefn(theta=start, loglkfn=MLInf, ProspModInf=MEtaProspModInf2, modelfn=modelfn, 
	 StratModInf=MEtaStratModInf2, y=cbind(z$y,y2), x=xarray,  npar=npar, nzval0=nzval0,
	 Aposn=1:nrow(z$y), Acounts=z$counts, Bposn=1:nrow(z$y), Bcounts=z$counts, 
         rmat=rmat, Qmat=Qmat, xStrat=as.numeric(z$xStrat), ymat=cbind(z$y,y2), xlist=xlist,
         Jis=z$Jis, yStrat=z$yStrat, paruse=paruse, control=control, control.inner=control.inner)  


names(res$theta) <- names(res$score) <- c(dimnames(X)[[2]], "log(sigma)", dimnames(X2)[[2]])
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
          cov=covmat, cor=correlation, Qmat=res$extra$Qmat, call=call, 
	  assign1=assign1, fnames1=fnames1, terms1=terms1, order1=order1,
	  assign2=assign2, fnames2=fnames2, terms2=terms2, order2=order2)
class(ans)<-"rclusbin2"
ans
}


#####################################################################################
print.rclusbin2<-function(x,digits=max(3, getOption("digits") - 3), ...){
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
summary.rclusbin2 <- function(object,...) {
#####################################################################################

# Works for R only so far

z <- object
if (!z$fit) stop("No summary provided when fit=FALSE !")

sigma <- z$cov
thetas <- z$coefficients    
assign1 <- z$assign1; assign2 <- z$assign2
terms1 <- z$terms1; terms2 <- z$terms2
order1 <- z$order1; order2 <- z$order2
fnames1 <- z$fnames1;  fnames2 <- z$fnames2
 

# calculate wald statistics ---------------------------------------------------

dowald <- function(assign, terms, thetas, sigma) {
    # find dfs for each variable
    num <- -1
    pos <- waldchis <- pvals <- coef.table <- NULL

    for (i in 1:length(assign)) {
    if (assign[i] > num) { 
       pos <- c(pos,i) 
       num <- assign[i] 
       }
    }
    pos1 <- c(pos[-1],length(assign)+1)
    values <- assign[pos]    
    dfs <- lengths <- pos1 - pos

    # do wald calculations
    if(length(pos)>1){
    for (i in 2:length(pos)) {
        thetai <- thetas[pos[i]:(pos[i]+lengths[i]-1)]
        covi <- sigma[pos[i]:(pos[i]+lengths[i]-1),pos[i]:(pos[i]+lengths[i]-1)]
        waldchi <- t(thetai) %*% solve(covi,thetai)
        waldchis <- c(waldchis,waldchi)
        pval <- (1-pchisq(waldchi,df=dfs[i]))
        pvals <- c(pvals,pval) 
    }
    coef.table <- cbind(dfs[-1], waldchis, pvals)
    dimnames(coef.table)[[1]] <- terms[values]
    dimnames(coef.table)[[2]] <- c("Df","Chi", "Pr(>Chi)")
    }
    coef.table
}

l1 <- length(assign1)+1; l2 <- length(assign2)
lt <- length(thetas)
ct1 <- dowald(assign1,terms1,thetas[1:l1],sigma[1:l1, 1:l1])
ct3 <- dowald(assign2,terms2,thetas[(l1+1):(l1+l2)], sigma[(l1+1):(l1+l2), (l1+1):(l1+l2)])


 # calculate t statistics for non factors ---------------------------------------------
if(is.null(dim(sigma))) s.err <- sqrt(sigma) else s.err <- sqrt(diag(sigma))

zstats <- thetas/s.err
pvals <- (1 - pnorm(abs(zstats)))*2
coef.table2 <- cbind(thetas, s.err, zstats, pvals)
dimnames(coef.table2)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

ct2 <- ct4 <- NULL

if (length(assign1) > 0) ct2 <- coef.table2[1:l1, ]
if (length(assign2) > 0) ct4 <- coef.table2[(l1+1):(l1+l2), ]

ans <- list(missReport=z$missReport, call=z$call, StrReport=z$StrReport, xStrReport=z$xStrReport,
       ClusReport=z$ClusReport, key=z$key, ykey=z$ykey, fit=z$fit, error=z$error, 
       loglk=z$loglk, coefficients=z$coefficients, coef.table1=ct1, coef.table2=ct2, 
       coef.table3=ct3, coef.table4=ct4)
class(ans) <- "summary.rclusbin2"
ans
}    


#####################################################################################
print.summary.rclusbin2 <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#####################################################################################    
cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}

cat("\nCluster Size Report:\n")
print(x$ClusReport, row.names=FALSE)
cat("\nStratum Counts Report:\n")

cat("\nStratum Counts Report:\n")
print(x$StrReport, row.names=FALSE)
if (!is.null(x$xStrReport)) {
    cat("\nObservations of obstype==xonly\n")
    print(x$xStrReport)
}

if (!is.null(x$key)) {
    cat("\nKey to x-Strat:\n")
    print.default(x$key,quote=FALSE, row.names=FALSE)
}

cat("\n\nModel of Interest (Random Intercept): ", x$ykey[1],"=",x$ykey[2]," given covariates\n",sep="")

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n\n")  
  
    if(!is.null(x$coef.table1)){
    cat("\n\"Y-Model\"\n") 
    cat("Wald Tests:\n")
    print(x$coef.table1,digits=4)
    cat("\n")
    }
    if(!is.null(x$coef.table2)){
    cat("Coefficients:\n")
    print(x$coef.table2,digits=4)
    }   

    if(!is.null(x$coef.table3)){
    cat("\n\"Z-Model\"\n") 
    cat("Wald Tests:\n")
    print(x$coef.table3,digits=4)
    cat("\n")
    }   
    if(!is.null(x$coef.table4)){
    cat("Coefficients:\n")
    print(x$coef.table4,digits=4)
    }   

}
else cat("\nCall requested that model not be fitted.\n")

invisible(x)
}

utils::globalVariables(c("counts","w"))
