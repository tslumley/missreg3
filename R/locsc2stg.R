#*****************************************
#**********  locsc2stg.R  ****************
#*****************************************

#####################################################################################
locsc2stg <- function(formula1, formula2, yCuts=NULL, weights=NULL, xstrata=NULL,
                     data=list(), obstype.name="obstype",  method="direct", fit=TRUE,
                     errdistrn="logistic", errmodpars=6, xs.includes=FALSE,
                     compactX=FALSE, compactY=TRUE, straty.maxnvals=20,
                     start=NULL, Qstart=NULL, deltastart=NULL, int.rescale=TRUE,
                     control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
#####################################################################################
{
# Process formula -------------------------------------------------------------------
# formula1 of the form y ~ x1 + x2..., formula2 of the form ~ x3 + x4...

mf <- call <- match.call()
mf[[1]] <- as.name("model.frame")
names(mf)[1] <- "model"
mf$na.action <- as.name("na.keep")
NULL -> mf$yCuts -> mf$xstrata -> mf$obstype.name -> mf$method -> mf$fit -> mf$errdistrn
NULL -> mf$compactX -> mf$compactY -> mf$straty.maxnvals
NULL -> mf$errmodpars -> mf$xs.includes -> mf$start -> mf$Qstart -> mf$deltastart
NULL -> mf$int.rescale -> mf$control -> mf$control.inner

# First formula processing ----------------------------------------------------------

mf1 <- mf[c("model","formula1","data","weights")]
mf1$na.action <- as.name("na.keep")
resp <- mf1$formula[2]
names(mf1)[2]<- "formula"
mf1 <- eval(mf1, sys.frame(sys.parent()))
Terms <- attr(mf1,"terms")
X1 <- model.matrix(Terms,mf1)
if(!is.matrix(X1)) X1 <- matrix(X1)
X1Orig <- X1
w <- model.extract(mf1,"weights")
terms1 <- attr(terms(formula1),"term.labels")
order1 <- attr(terms(formula1),"order")
assign1 <- attr(X1,"assign")
fnames1<-names(attr(X1,"contrasts"))
y <- model.extract(mf1,"response")

# Second formula processing ---------------------------------------------------------

mf2 <- mf[c("model","formula2","data")]
mf2$na.action <- as.name("na.keep")
names(mf2)[2] <- "formula"
mf2$formula <- as.formula(paste(resp,"~",mf2$formula[2]))
mf2 <- eval(mf2, sys.frame(sys.parent()))

Terms <- attr(mf2,"terms")
X2 <- model.matrix(Terms,mf2)
if(!is.matrix(X2)) X2 <- matrix(X2)
terms2 <- attr(terms(formula2),"term.labels")
order2 <- attr(terms(formula2),"order")
assign2 <- attr(X2,"assign")
fnames2<-names(attr(X2,"contrasts"))
colnames(X2) <- paste(colnames(X2), "2", sep = ".") # So names distinct
#------------------------------------------------------------------------------------


n <- if (is.matrix(y)) dim(y)[1] else length(y)
if(is.null(w)) w <- rep(1,n)

# Input checking ---------------------------------------------------------------------
#print("Input checking")

if (is.na(match(obstype.name,names(data)))) stop(paste("Dataframe did not have a column named",obstype.name))
obstype <- data[,match(obstype.name,names(data))] 
if (any(is.na(data[,obstype.name])))
   stop(paste("Missing values not permitted in obstype variable",obstype.name))

prospective <- all(obstype %in% c("uncond","y|x"))
# for S may need   prospective <- all(!is.na(match(obstype,c("uncond","y|x"))))

if (any(is.na(w))) stop("Weights should not contain NAs")

if (method == "ycutmeth") {
     if(!prospective & fit & is.null(yCuts)) stop("yCuts must be supplied if using method=ycutmeth")
}
else if (method=="direct"){
    if (!prospective & fit & is.null(yCuts)){
        if (is.null(start)) stop("Cannot calculate starting values for coefficients when yCuts is null")
#        if (is.null(deltastart)) print("Will use uniform deltastart starting values since yCuts is null")
    }
}
else
    stop(paste("method =",method, "not implemented"))

# Form the x-strata -----------------------------------------------------------------
#print("Form the x-strata")

xStrat <- if (is.null(xstrata)) rep(1,n) else {
    	xstratform <- as.formula(paste("~","-1 + ",paste(xstrata,sep="",collapse=":")))
	data1 <- as.data.frame(lapply(data, factor))
  	xstratmat <- model.matrix(xstratform, model.frame(xstratform,data1,na.action=function(x)x))
	for(i in 1:ncol(xstratmat)) xstratmat[,i] <- i * xstratmat[,i]
		apply(xstratmat,1,sum)
}
nStrat <- max(xStrat)

# Form the ystrata -----------------------------------------------------------------
yCutsKey <- NULL
if (is.null(yCuts)) {yStrat <- rep(1,n); nyStrat <- 1}
else { # Use yCuts to form y-strata
    #print("Form the ystrata")
    
    # Note that this is rather slow
    # To add later : match up the names with the xStrat internal ordering so we don't have to 
    # assume that the user knows the correct ordering of the yCut columns
    
    rngy <- range(y, na.rm=TRUE)
    rngdiff <- rngy[2] - rngy[1]
    yCutsrng <- range(yCuts,na.rm=TRUE)
    if (yCutsrng[1] <rngy[1] | yCutsrng[2]>rngy[2]) stop("yCuts must be inside range of y")

    if (!is.matrix(yCuts)) yCuts <- outer(as.vector(yCuts),rep(1,nStrat))
    if (ncol(yCuts) != nStrat) stop("No. of cols of yCuts matrix must equal number of x-strata")
    yStrat <- rep(NA,n)
    yCutsKey <- matrix(NA,ncol=nStrat,nrow=(nrow(yCuts)+1)) 
    yCutsNew <- matrix(NA,ncol=nStrat,nrow=(nrow(yCuts)+2)) 
    for(i in 1:nStrat){
        yCutsi <- c(rngy[1]-10*rngdiff, sort( yCuts[,i][!is.na(yCuts[,i])]), rngy[1]+10*rngdiff)
    	   temp <- cut(as.vector(y)[xStrat==i],yCutsi) 
    	   yStrat[xStrat==i] <- temp
        yCutsNew[1:length(yCutsi),i] <- yCutsi 
        theselevs <- levels(temp)
       	yCutsKey[1:length(theselevs),i] <- theselevs
    }
    yCuts <- yCutsNew

    yStrat <- as.numeric(yStrat)
    nyStrat <- max(yStrat)

}

# Remove random missing values --------------------------------------
#print("Remove random missing values")

mat <- cbind(y,X1,X2,xStrat)
miss1 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="uncond"]
miss2 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="retro"]
miss3 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="y|x"]
mat <- cbind(X1,X2,xStrat)
miss4 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="xonly"]
mat <- cbind(y,xStrat)
miss5 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="strata"]

toremove <- c(miss1,miss2,miss3,miss4,miss5)
missReport <- NULL
if(length(toremove)>0){	
#	   cat("\nObservations deleted due to missing data:\n")
   	if ((length(miss1)>0))
	  	   missReport <- rbind(missReport,paste("uncond:",length(miss1),"rows relating to",sum(w[miss1]),"observations"))
   	if ((length(miss2)>0))
		     missReport <- rbind(missReport,paste("retro:",length(miss2),"rows relating to",sum(w[miss2]),"observations")) 
	   if ((length(miss3)>0)) 
		     missReport <- rbind(missReport,paste("y|x:",length(miss3),"rows relating to",sum(w[miss3]),"observations")) 
	   if ((length(miss4)>0))
	     	missReport <- rbind(missReport,paste("xonly:",length(miss4),"rows relating to",sum(w[miss4]),"observations"))
	   if ((length(miss5)>0)) 
	     	missReport <- rbind(missReport,paste("strata:",length(miss5),"rows relating to",sum(w[miss5]),"observations")) 
    if (!is.null(missReport)) {colnames(missReport)<-"";rownames(missReport)<-rep("",nrow(missReport))}
    dimnames(missReport) <- list(rep("",nrow(missReport)),rep("",ncol(missReport)))
   	w <- w[-toremove]
	   obstype <- factor(as.character(obstype[-toremove]))
	   y <- y[-toremove]
	   if (!is.null(yCuts)) yStrat <- yStrat[-toremove]
	   xStrat <- xStrat[-toremove]
	   X1 <- X1[-toremove, ,drop=FALSE]	
	   X2 <- X2[-toremove, ,drop=FALSE]	
	   n <- length(y)
}
yStratfac <- factor(yStrat)
xStratfac <- factor(xStrat)

# Stratum Counts Report -------------------------------------------------------------

StrReport <- xStrReport <- key <- NULL

if (!is.null(yCuts)){
    sub <- (1:length(y))[obstype=="xonly"]
    if(length(sub)>0){
        df <- data.frame(w, obstype=factor(as.character(obstype)), yStrat=yStratfac, xStrat=xStratfac)[-sub,]
        df1 <- data.frame(w, obstype=factor(as.character(obstype)), yStrat=yStratfac, xStrat=xStratfac)[sub,]
    }
    else df <- data.frame(w,obstype=factor(as.character(obstype)),yStrat=yStratfac, xStrat=xStratfac)

    StrReport <- ftable(xtabs(w~obstype+yStrat+xStrat,data=df))
    if(exists("df1")) xStrReport <-  xtabs(w~obstype+xStrat,data=df1)

    if (xs.includes == TRUE) {
        oblevel0 <- levels(factor(as.character(obstype)))
        oblevel <- rep(oblevel0, rep(nyStrat,length(oblevel0)))  

        report <- ftable(xtabs(w~obstype+yStrat+xStrat,data=df[df$obstype=="strata",]))
        strep <- xtabs(w~yStrat+xStrat,data=df[df$obstype!="strata",])
        report[oblevel=="strata"] <- strep
        StrReport <- StrReport-report

        if (exists("df1")) {
            report1 <- xtabs(w~obstype+xStrat,data=df1[df1$obstype=="strata",])
            strep1 <- xtabs(w~xStrat,data=df1[df1$obstype!="strata",])
            report1[oblevel0=="strata"] <- strep1
            xStrReport <- xStrReport-report1
        }
    }
}
else { # yCuts is NULL so there are no y-strata
    df <- data.frame(w,obstype=factor(as.character(obstype)),xStrat=xStratfac)
    StrReport <- xtabs(w~obstype+xStrat,data=df)

    if (xs.includes == TRUE) {
        oblevel <- levels(factor(as.character(obstype)))
        report <- xtabs(w~obstype+xStrat,data=df[df$obstype=="strata",])
        strep <- xtabs(w~xStrat,data=df[df$obstype!="strata",])
        report[oblevel=="strata"] <- strep
        StrReport <- StrReport-report
    }
}

if (!(is.null(xstrata))) { # Form key to xStrat
    key <- dimnames(xstratmat)[[2]]
    names(key) <- 1:length(key)
}

if (!fit) {
     ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, yCutsKey=yCutsKey, fit=fit,
            call=call, assign1=assign1, assign2=assign2,
            fnames1=fnames1, fnames2=fnames2,terms1=terms1, terms2=terms2,
            order1=order1, order2=order2, n1=length(terms1)+1, n2=length(terms2)+1)
     class(ans) <- "locsc2stg"
     return(ans) # RETURN AT THIS POINT IF fit IS FALSE
}

# ------------------------------------------------------------------------------------
rmat <- Qmat <- NULL
if (!prospective & !is.null(yCuts)){
#    print("Construct Qmat and, if ycutmeth, also construct rmat")
    
    havestrata <- (length(obstype[obstype=="strata"])>0)
    if (method == "ycutmeth") {
        rk <- rep(0,n)
        if ((havestrata & xs.includes) | !havestrata) rk[obstype=="retro"] <- -1 * w[obstype=="retro"]
        if  (havestrata & xs.includes) rk[obstype=="uncond"] <- -1 * w[obstype=="uncond"]
        rk[obstype=="strata"] <- 1 * w[obstype=="strata"]
    
    
        rmat <- xtabs(rk~yStratfac+xStratfac)
        rmat <- matrix(rmat,nrow(rmat),ncol(rmat),dimnames=dimnames(rmat))
        rmat[is.na(yCutsKey)] <- NA
    
        if (nrow(rmat) == 1) stop("there is only 1 y-stratum present")
    }
    
    if (is.null(Qstart)){
        ind <- if (havestrata & !xs.includes) (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="strata"]
                              else (1:n)[obstype=="strata"]
        temp <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,w,obstype)[ind,])
    
        if (is.matrix(temp) && all(dim(temp) == c(nyStrat,nStrat)))	{
                  temp[is.na(yCutsKey)] <- NA
                  Qmat <- sweep(temp+1,2,apply(temp+1,2,FUN=function(x)sum(x,na.rm=TRUE)),"/")
        }
        else  stop("Attempt to construct missing Qstart failed")
    }
    else { # Qstart is present
        if(is.matrix(Qstart) && all(dim(Qstart) == c(nyStrat,nStrat)) &&  # matrix with correct dimensions
                all(round(apply(Qstart,2,sum),5) == rep(1,dim(Qstart)[2]))) # cols sum to 1
               Qmat <- Qstart
        else stop("illegal Qstart value")
    }
}
    
# Create start if needed ------------------------------------------------------------
# print("Create start if needed")
if (is.null(start)){ # Use lm to calculate starting values for the coefficients
    wmult <- rep(1,n)  
   ind <- 1:n
   if (!prospective) {
       if (all(!is.na(match(obstype,c("retro","xonly")))))  # all retro or xonly
 	          stop("Cannot construct starting values from this data")
   	
 	     ind <- (1:n)[obstype=="retro" | obstype=="y|x" | obstype=="uncond"]
 	     if (!is.na(match("retro",obstype))){ # adjust values of wmult if have retro obsns
           nMat1 <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,w,obstype)[obstype=="uncond" | obstype=="retro",])
           nMat <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,w,obstype)[obstype=="retro",])
   	       NMat <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,w,obstype)[obstype=="strata",])
   	       if (!xs.includes) NMat <- NMat + nMat1
           NMat[is.na(yCutsKey)] <- NA
     	    	Nn <- (NMat+1)/(nMat+1)
        	 	Nn <- sweep(Nn,2,apply(Nn,2,FUN=function(x)mean(x,na.rm=TRUE)),"/")
   		      ind1 <- (1:n)[obstype=="retro"] 
   		      Nnvec <- as.vector(Nn) 
   		      xystrat <- (xStrat-1)*(nyStrat)+yStrat
   		      wmult <- Nnvec[xystrat] # multiply by corresponding values in Nn
           if (any(is.na(wmult))) stop("Missing values in constructed wmult vector")
   	    }
	   }
	   gg <- lm(y[ind]~X1[ind,]-1,weights=w[ind]*wmult[ind])
	   start1 <- gg$coefficients
    scalef <- 1
    if (errdistrn=="logistic") scalef <- 1.8
    if (errdistrn=="t") {dfs <- errmodpars[1]; scalef <- sqrt(dfs/(dfs-2))}
	   start2 <- c(log(summary(gg)$sigma/scalef),rep(0,ncol(X2)-1))
	   start <- c(start1,start2)
}
else if (length(start) != (ncol(X1)+ncol(X2))){	

	# Match up any values in start with names in X, fill in holes with 0's.
    newstart <- rep(0,(dim(X1)[2]+dim(X2)[2]))
	   posns <- match(names(start),c(dimnames(X1)[[2]],dimnames(X2)[[2]]))   
    newstart[posns] <- start
	   start <- newstart
}

# Form xarray --------------------------------------------------------------

if (method=="direct")
    if (!prospective & method=="direct"){
    #    cat("Using direct method","\n\n")
        hvalue <- y[obstype=="strata"]
        hxStrat <- xStrat[obstype=="strata"]
        Cmult <- w[obstype=="strata"]
        if (compactY){
            # limit hvalues to straty.maxnvals  equally spaced possible values
            aaa <- min(hvalue); bbb <- max(hvalue);
            hvalue <- aaa + (round(straty.maxnvals*(hvalue-aaa)/(bbb-aaa)))*(bbb-aaa)/straty.maxnvals
            # compact down to distinct vales and multiplicities
            temp <- weightform(data.frame(hvalue,hxStrat),counts=Cmult)
            hvalue <- hvalue[temp$ind]; hxStrat <- hxStrat[temp$ind]; Cmult <- temp$counts
        }
    }
    else hvalue <- Cmult <- hxStrat <- numeric(0)
if (int.rescale){  # (standardise the variables
    muY <- mean(y,na.rm=TRUE)
    sdY <- sd(y,na.rm=TRUE)
    y <- (y-muY)/sdY
    if (!is.null(yCuts)) yCuts <- (yCuts-muY)/sdY
    if (!prospective & method=="direct") hvalue <- (hvalue-muY)/sdY

    scale1 <- apply(X1,2,FUN=function(x)sd(x,na.rm=TRUE))
    constant1 <- (scale1==0)
    loc1   <- apply(X1,2,FUN=function(x)mean(x,na.rm=TRUE))
    loc1[constant1] <- 0
    scale1[constant1] <- 1
    nrowsx <- nrow(X1)
    ncx1 <- ncol(X1)
    X1 <- (X1 - outer(rep(1,nrowsx),loc1)) / outer(rep(1,nrowsx),scale1)
    start1 <- start[1:ncx1]
    start1 <- (start1*scale1 + constant1*(sum(start1*loc1)-muY))/sdY

    scale2 <- apply(X2,2,FUN=function(x)sd(x,na.rm=TRUE))
    constant2 <- (scale2==0)
    if (sum(constant1)>1 | sum(constant2)>1) stop("More than 1 contant X-variable in the same model")
    loc2   <- apply(X2,2,FUN=function(x)mean(x,na.rm=TRUE))
    loc2[constant2] <- 0
    scale2[constant2] <- 1
    ncx2 <- ncol(X2)
    X2 <- (X2 - outer(rep(1,nrowsx),loc2)) / outer(rep(1,nrowsx),scale2)
    start2 <- start[(ncx1+1):(ncx1+ncx2)]
    start2 <- start2 * scale2 + constant2*(sum(start2*loc2)-log(sdY))

    start <- c(start1, start2)
}


ncolX1 <- ncol(X1)
X1 <- insertColumns(X1,matrix(0,n,ncol(X2),dimnames=list(NULL,dimnames(X2)[[2]])),ncol(X1)+1)
X2 <- insertColumns(X2,matrix(0,n,ncolX1),1)

X1 <- X1[obstype!="strata",]
X2 <- X2[obstype!="strata",]
w <- w[obstype!="strata"]
y <- y[obstype!="strata"]
xStrat <- xStrat[obstype!="strata"]
yStrat <- yStrat[obstype!="strata"]
obstype <- obstype[obstype!="strata"]
n <- length(y)

xarray <- array(0,c(n,ncol(X1),2))
xarray[, ,1] <- X1
xarray[, ,2] <- X2

Aposn <- (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="y|x"]
Acounts <- w[Aposn]
if (!prospective){
    Bposn <- (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="xonly"]
    Bcounts <- w[Bposn]
}
else Bposn <- Bcounts <- numeric(0)


deltawts <- NULL
if (!prospective & method=="direct" & is.null(deltastart)){
    # Preliminary calculations for starting values for the deltas in MLdirectInf
    deltawts <- NULL
    if (any(obstype=="retro" & !is.null(yCuts))) {
        if (is.null(Qmat)) print("SCREAM")
        indQ <- yStrat + (xStrat-1)*nyStrat
        deltawts <- as.vector(Qmat)[indQ]
        deltawts[is.na(deltawts) | obstype != "retro"] <- 1
        deltawts <- w * deltawts
    }
}
if (!prospective & compactX){  #   Compact down x-variable value information to distinct values accumulating counts
    temp <- weightform(data.frame(X1,X2,xStrat,obstype),counts=w,counts2=deltawts)
    Bposn <- temp$ind
    Bcounts <- temp$counts
    if (!is.null(deltawts)) deltawts <- temp$counts2
}
if (!prospective & any(obstype=="y|x", na.rm = TRUE)){ # Remove any y|x obsns from Bposn etc
     ind2 <- (1:n)[obstype=="y|x"]
     ToRemove <- (1:length(Bposn))[!is.na(match(Bposn,intersect(Bposn,ind2)))]
     Bposn <- Bposn[-ToRemove]; Bcounts <- Bcounts[-ToRemove]
     if (is.null(deltawts)) deltawts <- deltawts[-ToRemove]
}
deltamat <- NULL
if (!prospective & method=="direct"){
    lBmax <- 0
    for (j in 1:nStrat) lBmax <- max(length(Bposn[xStrat[Bposn]==j]),lBmax)
    if (!is.null(deltastart)) {
        deltamat <- deltastart
        if (any(dim(deltamat) != c(lBmax,nStrat))) stop(paste(
                  "Dimension problems with deltastart,dim=", dim(deltamat)))
    }
    else if (!is.null(yCuts)){ # continue calculations of starting values for deltas
        deltamat <- matrix(NA,lBmax,nStrat)
        for (j in 1:nStrat){
             deltamat[1:length(Bposn[xStrat[Bposn]==j]),j] <- deltawts[xStrat[Bposn]==j]
        }
        deltamat <- sweep(deltamat,2,apply(deltamat,2,FUN=function(x)sum(x,na.rm=TRUE)),"/")
    }
}

# Function call ------------------------------------------------------------
#print("Function call")

extra <- NULL

if (method=="direct"){
   if (errdistrn=="logistic") errdist <- logisterr
   else if (errdistrn=="normal") errdist <- normerr
   else if (errdistrn=="t") errdist <- terr
   else stop(paste("error distribution",errdistrn,"not implemented"))

   res <- mlefn(theta=start,MLdirectInf,deltamat=deltamat,
                    modelfn=locscale,hmodelfn=hylocscale,errdist=errdist,errmodpars=errmodpars,
                    x=xarray,y=as.matrix(y),xStrat=as.numeric(xStrat),
                    Aposn=Aposn,Acounts=Acounts,Bposn=Bposn,Bcounts=Bcounts,
                    hvalue=hvalue,Cmult=Cmult,hxStrat=hxStrat,
                    control=control, control.inner=control.inner)
}
else if (method=="ycutmeth"){
   if (errdistrn=="logistic") {errdist <- logisterr; errdistcdf <- logistcdf}
   else if (errdistrn=="normal") {errdist <- normerr; errdistcdf <- normcdf}
   else if (errdistrn=="t") {errdist <- terr; errdistcdf <- tcdf}
   else stop(paste("error distribution",errdistrn,"not implemented"))

   res <- mlefn(theta=start,MLInf,y=as.matrix(y),x=xarray,ProspModInf=MEtaProspModInf,StratModInf=MEtaStratModInf,
   modelfn=locscale,errdist=errdist,errmodpars=errmodpars,stratfn=locscstrat1,errdistcdf=errdistcdf,
             Aposn=Aposn,Acounts=Acounts,Bposn=Bposn,Bcounts=Bcounts,rmat=rmat,Qmat=Qmat,
             xStrat=as.numeric(xStrat),yCuts=yCuts,
             control=control, control.inner=control.inner)
}
else stop("method not implemented")

if (res$error != 1) {######

if (int.rescale){  # Translate back to the scale of the original variables
    theta1 <- res$theta[1:ncx1]
    theta1 <- (theta1*sdY/scale1) + constant1*(muY - sum((theta1*sdY/scale1)*loc1))
    theta2 <- res$theta[(ncx1+1):(ncx1+ncx2)]
    theta2 <- (theta2 / scale2) + constant2*((log(sdY)-sum((theta2/scale2)*loc2)))
    res$theta <- c(theta1, theta2)
    jacbn <- matrix(0,ncx1+ncx2,ncx1+ncx2)
    diag(jacbn)[1:ncx1] <- scale1/sdY
    jacbn[1:ncx1,1:ncx1][constant1,] <- jacbn[1:ncx1,1:ncx1][constant1,] + loc1/sdY
    diag(jacbn)[(ncx1+1): (ncx1+ncx2)] <- scale2
    if (any(constant2) & length(constant2)>1)
       jacbn[(ncx1+1): (ncx1+ncx2),(ncx1+1): (ncx1+ncx2)][constant2,] <-
                           jacbn[(ncx1+1): (ncx1+ncx2),(ncx1+1): (ncx1+ncx2)][constant2,] + loc2
    res$loglk <- res$loglk - sum(Acounts)*log(sdY)
    if (method=="direct") res$loglk <- res$loglk - sum(Cmult)*log(sdY)
    res$score <- as.vector(t(jacbn)%*%res$score)
    res$inf <- t(jacbn)%*%res$inf%*%jacbn
}
names(res$theta) <- names(res$score) <- dimnames(X1)[[2]]
dimnames(res$inf) <- list(dimnames(X1)[[2]],dimnames(X1)[[2]])
#covmat <- solve(res$inf) # change to avoid scale problems with inverse
temp <- 1/sqrt(diag(res$inf))
covmat <- diag(temp)%*%solve(diag(temp)%*%res$inf%*%diag(temp))%*%diag(temp)
dimnames(covmat) <- dimnames(res$inf)
dd <- sqrt(diag(covmat))
correlation <- covmat/outer(dd, dd)
pred <- X1Orig %*% res$theta[1:ncol(X1Orig)]

} ######
else pred <- covmat <- correlation <- NULL

# ------------------------------------------------------------------------------------
ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, yCutsKey=yCutsKey, fit=fit,
            error=res$error, coefficients=res$theta, loglk=res$loglk, score=res$score, inf=res$inf,
            fitted=pred, cov=covmat, cor=correlation,
            Qmat = res$extra$Qmat,deltamat=res$extra$deltamat,
            call=call, assign1=assign1, assign2=assign2,
            fnames1=fnames1, fnames2=fnames2,terms1=terms1, terms2=terms2,
            order1=order1, order2=order2, n1=length(terms1)+1, n2=length(terms2)+1)
class(ans) <- "locsc2stg"
ans
}


#####################################################################################
print.locsc2stg <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#####################################################################################
#n1 <- x$n1; n2 <- x$n2    
cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}
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
if (!is.null(x$yCutsKey)) {
    cat("\nKey to the y-Strat:\n")
    print.default(x$yCutsKey,quote=FALSE,digits=digits)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n\n")

#    cat("Covariance:\n")
#    print.default(format(x$cov, digits = digits), print.gap = 2, quote = FALSE)
#    cat("\n")
}
invisible(x)
}

#####################################################################################
summary.locsc2stg <- function(object,...) {
#####################################################################################

    z <- object
    sigma <- z$cov
    thetas <- z$coefficients    
    assign1 <- z$assign1; assign2 <- z$assign2
    terms1 <- z$terms1; terms2 <- z$terms2
    order1 <- z$order1; order2 <- z$order2
    fnames1 <- z$fnames1; fnames2 <- z$fnames2
    
    # calculate wald statistics

    dowald <- function(assign,terms,thetas,sigma){

    # find dfs for each variable
    num <- -1
    pos <- waldchis <- pvals <- coef.table <- NULL
    for (i in 1:length(assign)) {
   
        if (assign[i] > num) {
            pos <- c(pos, i)
            num <- assign[i]
        }
    }
    pos1 <- c(pos[-1], length(assign) + 1)
    values <- assign[pos]
    dfs <- lengths <- pos1 - pos
    if (length(pos) > 1) {
        for (i in 2:length(pos)) {
            thetai <- thetas[pos[i]:(pos[i] + lengths[i] - 1)]
            covi <- sigma[pos[i]:(pos[i] + lengths[i] - 1), pos[i]:(pos[i] + 
                lengths[i] - 1)]
            waldchi <- t(thetai) %*% solve(covi, thetai)
            waldchis <- c(waldchis, waldchi)
            pval <- (1 - pchisq(waldchi, df = dfs[i]))
            pvals <- c(pvals, pval)
        }
        coef.table <- cbind(dfs[-1], waldchis, pvals)
        dimnames(coef.table)[[1]] <- terms[values]
        dimnames(coef.table)[[2]] <- c("Df", "Chi", "Pr(>Chi)")
    }
coef.table
}

l1 <- length(assign1)
lt <- length(thetas)
ct1 <- dowald(assign1,terms1,thetas[1:l1],sigma[1:l1,1:l1])
ct3 <- dowald(assign2,paste(terms2,"2",sep="."),thetas[(l1+1):lt],sigma[(l1+1):lt,(l1+1):lt])

    # calculate t statistics for non factors ####################
    if(is.null(dim(sigma))) s.err <- sqrt(sigma) else s.err <- sqrt(diag(sigma))

    zstats <- thetas/s.err
    pvals <- (1 - pnorm(abs(zstats)))*2
    coef.table2 <- cbind(thetas, s.err, zstats, pvals)
    dimnames(coef.table2)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    ct2 <- ct4 <- NULL

    if (length(assign1) > 0) {
        ct2 <- coef.table2[1:length(assign1), ]
    }
    if (length(assign2) > 0) {
        ct4 <- coef.table2[(length(assign1) + 1):length(c(assign1,assign2)), ]
    }
    	
ans <- list(missReport=z$missReport, call=z$call,StrReport=z$StrReport, xStrReport=z$xStrReport,
       key=z$key, yCutsKey=z$yCutsKey, fit=z$fit, error=z$error, loglk=z$loglk, 
       coefficients=z$coefficients,coef.table1=ct1,coef.table2=ct2,coef.table3=ct3,coef.table4=ct4)
class(ans) <- "summary.locsc2stg"
ans
}    

#####################################################################################
print.summary.locsc2stg <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#####################################################################################    
cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}
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
if (!is.null(x$yCutsKey)) {
    cat("\nKey to the y-Strat:\n")
    print.default(x$yCutsKey,quote=FALSE,digits=digits)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n")
    if(!is.null(x$coef.table1)){
        cat("\nLocation Model:\n")
        cat("Wald Tests:\n")
        print(x$coef.table1,digits=4)
        cat("\n")
    }
    if(!is.null(x$coef.table2)){
    cat("Coefficients:\n")
    print(x$coef.table2,digits=4)
    } 
    if(!is.null(x$coef.table3) | !is.null(x$coef.table4))
    	   cat("\nScale Model:\n")
    if(!is.null(x$coef.table3)){
    cat("Wald Tests:\n")
        print(x$coef.table3,digits=4)
        cat("\n")
    }
    if(!is.null(x$coef.table4)){
        cat("Coefficients:\n")
        print(x$coef.table4,digits=4)
    }
}
else cat("\nCall requested that model not be fitted\n") 
invisible(x)
}


#####################################################################################
locsc2stgwtd <- function(formula1, formula2, yCuts=NULL, weights=NULL, xstrata=NULL,
                     data=list(), obstype.name="obstype", fit=TRUE, errdistrn="logistic", 
                     errmodpars=6, xs.includes=FALSE, start=NULL, celln=10, sampling="ycon",
                     control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
#####################################################################################
{
# Process formula -------------------------------------------------------------------
# formula1 of the form y ~ x1 + x2..., formula2 of the form ~ x3 + x4...

mf <- call <- match.call()
mf[[1]] <- as.name("model.frame")
names(mf)[1] <- "model"
mf$na.action <- as.name("na.keep")
NULL -> mf$yCuts -> mf$xstrata -> mf$obstype.name -> mf$fit -> mf$errdistrn
NULL -> mf$errmodpars -> mf$xs.includes -> mf$start -> mf$celln -> mf$sampling
NULL -> mf$control -> mf$control.inner

# First formula processing ----------------------------------------------------------

mf1 <- mf[c("model","formula1","data")]
mf1$na.action <- as.name("na.keep")
resp <- mf1$formula[2]
names(mf1)[2]<- "formula"
mf1 <- eval(mf1, sys.frame(sys.parent()))
Terms <- attr(mf1,"terms")
X1 <- model.matrix(Terms,mf1)
if(!is.matrix(X1)) X1 <- matrix(X1)
X1Orig <- X1
terms1 <- attr(terms(formula1),"term.labels")
order1 <- attr(terms(formula1),"order")
assign1 <- attr(X1,"assign")
fnames1<-names(attr(X1,"contrasts"))
y <- model.extract(mf1,"response")

# Second formula processing ---------------------------------------------------------

mf2 <- mf[c("model","formula2","data")]
mf2$na.action <- as.name("na.keep")
names(mf2)[2] <- "formula"
mf2$formula <- as.formula(paste(resp,"~",mf2$formula[2]))
mf2 <- eval(mf2, sys.frame(sys.parent()))

Terms <- attr(mf2,"terms")
X2 <- model.matrix(Terms,mf2)
if(!is.matrix(X2)) X2 <- matrix(X2)
terms2 <- attr(terms(formula2),"term.labels")
order2 <- attr(terms(formula2),"order")
assign2 <- attr(X2,"assign")
fnames2<-names(attr(X2,"contrasts"))
colnames(X2) <- paste(colnames(X2), "2", sep = ".") # So names distinct
#------------------------------------------------------------------------------------


n <- if (is.matrix(y)) dim(y)[1] else length(y)
w <- weights
if(is.null(w)) w <- rep(1,n)

# Input checking ---------------------------------------------------------------------
#print("Input checking")

if (is.na(match(obstype.name,names(data)))) stop(paste("Dataframe did not have a column named",obstype.name))
obstype <- data[,match(obstype.name,names(data))] 
if (any(is.na(data[,obstype.name])))
   stop(paste("Missing values not permitted in obstype variable",obstype.name))

prospective <- all(obstype %in% c("uncond","y|x"))
# for S may need   prospective <- all(!is.na(match(obstype,c("uncond","y|x"))))

if (any(is.na(w))) stop("Weights should not contain NAs")

# Form the x-strata -----------------------------------------------------------------
#print("Form the x-strata")

xStrat <- if (is.null(xstrata)) rep(1,n) else {
    	xstratform <- as.formula(paste("~","-1 + ",paste(xstrata,sep="",collapse=":")))
	data1 <- as.data.frame(lapply(data, factor))
  	xstratmat <- model.matrix(xstratform, model.frame(xstratform,data1,na.action=function(x)x))
	for(i in 1:ncol(xstratmat)) xstratmat[,i] <- i * xstratmat[,i]
		apply(xstratmat,1,sum)
}
nStrat <- max(xStrat)

# Form the ystrata -----------------------------------------------------------------
yCutsKey <- NULL
if (is.null(yCuts)) {yStrat <- rep(1,n); nyStrat <- 1}
else { # Use yCuts to form y-strata
    #print("Form the ystrata")
    
    # Note that this is rather slow
    # To add later : match up the names with the xStrat internal ordering so we don't have to 
    # assume that the user knows the correct ordering of the yCut columns
    
    rngy <- range(y, na.rm=TRUE)
    rngdiff <- rngy[2] - rngy[1]
    yCutsrng <- range(yCuts,na.rm=TRUE)
    if (yCutsrng[1] <rngy[1] | yCutsrng[2]>rngy[2]) stop("yCuts must be inside range of y")

    if (!is.matrix(yCuts)) yCuts <- outer(as.vector(yCuts),rep(1,nStrat))
    if (ncol(yCuts) != nStrat) stop("No. of cols of yCuts matrix must equal number of x-strata")
    yStrat <- rep(NA,n)
    yCutsKey <- matrix(NA,ncol=nStrat,nrow=(nrow(yCuts)+1)) 
    yCutsNew <- matrix(NA,ncol=nStrat,nrow=(nrow(yCuts)+2)) 
    for(i in 1:nStrat){
        yCutsi <- c(rngy[1]-10*rngdiff, sort( yCuts[,i][!is.na(yCuts[,i])]), rngy[1]+10*rngdiff)
    	   temp <- cut(as.vector(y)[xStrat==i],yCutsi) 
    	   yStrat[xStrat==i] <- temp
        yCutsNew[1:length(yCutsi),i] <- yCutsi 
        theselevs <- levels(temp)
       	yCutsKey[1:length(theselevs),i] <- theselevs
    }
    yCuts <- yCutsNew

    yStrat <- as.numeric(yStrat)
    nyStrat <- max(yStrat)

}

# Remove random missing values --------------------------------------
#print("Remove random missing values")

mat <- cbind(y,X1,X2,xStrat)
miss1 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="uncond"]
miss2 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="retro"]
miss3 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="y|x"]
mat <- cbind(X1,X2,xStrat)
miss4 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="xonly"]
mat <- cbind(y,xStrat)
miss5 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="strata"]

toremove <- c(miss1,miss2,miss3,miss4,miss5)
missReport <- NULL
if(length(toremove)>0){	
#	   cat("\nObservations deleted due to missing data:\n")
   	if ((length(miss1)>0))
	  	   missReport <- rbind(missReport,paste("uncond:",length(miss1),"rows relating to",sum(w[miss1]),"observations"))
   	if ((length(miss2)>0))
		     missReport <- rbind(missReport,paste("retro:",length(miss2),"rows relating to",sum(w[miss2]),"observations")) 
	   if ((length(miss3)>0)) 
		     missReport <- rbind(missReport,paste("y|x:",length(miss3),"rows relating to",sum(w[miss3]),"observations")) 
	   if ((length(miss4)>0))
	     	missReport <- rbind(missReport,paste("xonly:",length(miss4),"rows relating to",sum(w[miss4]),"observations"))
	   if ((length(miss5)>0)) 
	     	missReport <- rbind(missReport,paste("strata:",length(miss5),"rows relating to",sum(w[miss5]),"observations")) 
    if (!is.null(missReport)) {colnames(missReport)<-"";rownames(missReport)<-rep("",nrow(missReport))}
   	w <- w[-toremove]
	   obstype <- factor(as.character(obstype[-toremove]))
	   y <- y[-toremove]
	   if (!is.null(yCuts)) yStrat <- yStrat[-toremove]
	   xStrat <- xStrat[-toremove]
	   X1 <- X1[-toremove, ,drop=FALSE]	
	   X2 <- X2[-toremove, ,drop=FALSE]	
	   n <- length(y)
}
yStratfac <- factor(yStrat)
xStratfac <- factor(xStrat)

# Stratum Counts Report -------------------------------------------------------------

StrReport <- xStrReport <- key <- NULL

if (!is.null(yCuts)){
    sub <- (1:length(y))[obstype=="xonly"]
    if(length(sub)>0){
        df <- data.frame(w, obstype=factor(as.character(obstype)), yStrat=yStratfac, xStrat=xStratfac)[-sub,]
        df1 <- data.frame(w, obstype=factor(as.character(obstype)), yStrat=yStratfac, xStrat=xStratfac)[sub,]
    }
    else df <- data.frame(w,obstype=factor(as.character(obstype)),yStrat=yStratfac, xStrat=xStratfac)

    StrReport <- ftable(xtabs(w~obstype+yStrat+xStrat,data=df))
    if(exists("df1")) xStrReport <-  xtabs(w~obstype+xStrat,data=df1)

    if (xs.includes == TRUE) {
        oblevel0 <- levels(factor(as.character(obstype)))
        oblevel <- rep(oblevel0, rep(nyStrat,length(oblevel0)))  

        report <- ftable(xtabs(w~obstype+yStrat+xStrat,data=df[df$obstype=="strata",]))
        strep <- xtabs(w~yStrat+xStrat,data=df[df$obstype!="strata",])
        report[oblevel=="strata"] <- strep
        StrReport <- StrReport-report

        if (exists("df1")) {
            report1 <- xtabs(w~obstype+xStrat,data=df1[df1$obstype=="strata",])
            strep1 <- xtabs(w~xStrat,data=df1[df1$obstype!="strata",])
            report1[oblevel0=="strata"] <- strep1
            xStrReport <- xStrReport-report1
        }
    }
}
else { # yCuts is NULL so there are no y-strata
    df <- data.frame(w,obstype=factor(as.character(obstype)),xStrat=xStratfac)
    StrReport <- xtabs(w~obstype+xStrat,data=df)

    if (xs.includes == TRUE) {
        oblevel <- levels(factor(as.character(obstype)))
        report <- xtabs(w~obstype+xStrat,data=df[df$obstype=="strata",])
        strep <- xtabs(w~xStrat,data=df[df$obstype!="strata",])
        report[oblevel=="strata"] <- strep
        StrReport <- StrReport-report
    }
}

if (!(is.null(xstrata))) { # Form key to xStrat
    key <- dimnames(xstratmat)[[2]]
    names(key) <- 1:length(key)
}

if (!fit) {
     ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, yCutsKey=yCutsKey, fit=fit,
            call=call, assign1=assign1, assign2=assign2,
            fnames1=fnames1, fnames2=fnames2,terms1=terms1, terms2=terms2,
            order1=order1, order2=order2, n1=length(terms1)+1, n2=length(terms2)+1)
     class(ans) <- "locsc2stg"
     return(ans) # RETURN AT THIS POINT IF fit IS FALSE
}

# Construct required weights -------------------------------------------------------
wmult <- rep(1,n)
ind <- 1:n

if (!is.null(yCuts)) { #when have ystrata as real strati or post-strati, adjust for weights
  ind <- (1:n)[obstype=="retro" | obstype=="y|x" | obstype=="uncond"]
  nMat <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,w,obstype)[obstype=="uncond" | obstype=="retro",])
  NMat <- xtabs(w~yStrat+xStrat,data=data.frame(yStrat=yStratfac,xStrat=xStratfac,w,obstype)[obstype=="strata",])	
  if (!xs.includes) NMat <- NMat + nMat
  
  if (sum(NMat)==0) stop("Y is not fully observed, out of current consideration!")
  temp <- NMat/sum(NMat)
  if (sampling == "ycon") qvec <- t(temp)[t(nMat)!=0]
    else qvec <- temp[nMat!=0]

  ny <- nrow(nMat); nx <- ncol(nMat)	
  NMat0 <- NMat;  nMat0 <- nMat
  newyStrat <- yStrat;  newxStrat <- xStrat
  nomerge <- TRUE
	
  if (any(nMat < celln)) { # merge cells that have less than 'celln' observations 
	if (sampling == "xycon") warning(paste("Less than",celln,"observations in some cells!"))
	 
	if (sampling == "xcon" | sampling == "uncon") {	 
	  nomerge <- FALSE 
	  qvec <- NULL; total <- rep(NA, nx)
      for (j in 1:nx)  {
       idj <- 1:ny
	   for (i in 1:ny) { #merge by column!
	     count1 <- count2 <- 0; bward <- FALSE
	     new.n <- nMat[i,j]; new.N <- NMat[i,j]
	     if (new.n < celln) {
	       while (new.n < celln) {
	         if ((i+count1) < ny) {
	             count1 <- count1+1
	             new.n <- sum(nMat0[i:(i+count1),j])
	             new.N <- sum(NMat0[i:(i+count1),j])
	         }
	         else  {
	             if ((i-count2) > 1) {
	               count2 <- count2+1
	               new.n <- sum(nMat0[(i-count2):ny,j])
	               new.N <- sum(NMat0[(i-count2):ny,j])
                   bward <- TRUE 
	             }
	             else { 
	               warning(paste("Less than",celln,"observations in Xstrata",j,"!"))
	               break() }
	         }	          
	       }
	      
	       if (bward) {
	         nMat[(i-count2):ny,j] <- rep(new.n,ny-i+count2+1)
	         NMat[(i-count2):ny,j] <- rep(new.N,ny-i+count2+1)
	         newyStrat[xStrat==j][!is.na(match(yStrat[xStrat==j], (i-count2):ny))] <- i-count2
	         idj[(i-count2):ny] <- i-count2
	       }
	       else {
	         nMat[i:(i+count1),j] <- rep(new.n,count1+1)
	         NMat[i:(i+count1),j] <- rep(new.N,count1+1)
	         newyStrat[xStrat==j][!is.na(match(yStrat[xStrat==j], i:(i+count1)))] <- i
	         idj[i:(i+count1)] <- i
	       }
	    } 
	  }
      nj <- nMat[unique(idj),j]; Nj <- NMat[unique(idj),j]
	  total[j] <- sum(Nj/nj)    
	  qvec <- c(qvec, Nj/sum(NMat0))  #ordered by Y then X
	 }  
    }
    
    if (sampling == "ycon") {
	  qvec <- NULL 
      for (i in 1:ny)  {
       idi <- 1:nx
	   for (j in 1:nx) { #merge by row!
	     count1 <- count2 <- 0; bward <- FALSE
	     new.n <- nMat[i,j]; new.N <- NMat[i,j]
	     if (new.n < celln) {
	       while (new.n < celln) {
	         if ((j+count1) < nx) {
	             count1 <- count1+1
	             new.n <- sum(nMat0[i,j:(j+count1)])
	             new.N <- sum(NMat0[i,j:(j+count1)])
	         }
	         else  {
	             if ((j-count2) > 1) {
	               count2 <- count2+1
	               new.n <- sum(nMat0[i,(j-count2):nx])
	               new.N <- sum(NMat0[i,(j-count2):nx])
                   bward <- TRUE 
	             }
	             else { 
	               warning(paste("Less than",celln,"observations in Ystrata",i,"!"))
	               break() }
	         }	          
	       }
	      
	       if (bward) {
	         nMat[i,(j-count2):nx] <- rep(new.n,nx-j+count2+1)
	         NMat[i,(j-count2):nx] <- rep(new.N,nx-j+count2+1)
	         newxStrat[yStrat==i][!is.na(match(xStrat[yStrat==i], (j-count2):nx))] <- j-count2
	         idi[(j-count2):nx] <- j-count2
	       }
	       else {
	         nMat[i,j:(j+count1)] <- rep(new.n,count1+1)
	         NMat[i,j:(j+count1)] <- rep(new.N,count1+1)
	         newxStrat[yStrat==i][!is.na(match(xStrat[yStrat==i], j:(j+count1)))] <- j
	         idi[j:(j+count1)] <- j
	       }
	    } 
	  }
      Ni <- NMat[i,unique(idi)]
	  qvec <- c(qvec, Ni/sum(NMat0))  #ordered by X then Y
	 }  
    }
    
  }

  if (any(nMat==0)) Nn <- (NMat+1)/(nMat+1)
    else Nn <- NMat/nMat
  if (sampling == "xcon" | sampling == "xycon") { #make weights sum to 1 in each xstrata 
    if (nomerge) total <- apply(Nn,2,FUN=function(x)sum(x,na.rm=TRUE))
    Nn <- sweep(Nn,2,total,"/") 
  }
  Nnvec <- as.vector(Nn)
  xystrat <- (xStrat-1)*(nyStrat)+yStrat
  wmult <- Nnvec[xystrat] # multiply by corresponding values in Nn
  if (any(is.na(wmult))) stop("Missing values in constructed 'wmult' vector for weights")
}
w <- w*wmult
 	
# Create start if needed ------------------------------------------------------------
# print("Create start if needed")

if (is.null(start)){ # Use lm to calculate starting values for the coefficients
	gg <- lm(y[ind]~X1[ind,]-1,weights=w[ind])
	start1 <- gg$coefficients
    scalef <- 1
    if (errdistrn=="logistic") scalef <- 1.8
    if (errdistrn=="t") {dfs <- errmodpars[1]; scalef <- sqrt(dfs/(dfs-2))}
	   start2 <- c(log(summary(gg)$sigma/scalef),rep(0,ncol(X2)-1))
	   start <- c(start1,start2)
}
else if (length(start) != (ncol(X1)+ncol(X2))){	

	# Match up any values in start with names in X, fill in holes with 0's.
	newstart <- rep(0,(dim(X1)[2]+dim(X2)[2]))
	posns <- match(names(start),c(dimnames(X1)[[2]],dimnames(X2)[[2]]))   
    newstart[posns] <- start
	start <- newstart
}

# Form xarray --------------------------------------------------------------

ncolX1 <- ncol(X1)
X1 <- insertColumns(X1,matrix(0,n,ncol(X2),dimnames=list(NULL,dimnames(X2)[[2]])),ncol(X1)+1)
X2 <- insertColumns(X2,matrix(0,n,ncolX1),1)

X1 <- X1[obstype!="strata",]
X2 <- X2[obstype!="strata",]
w <- w[obstype!="strata"]
y <- y[obstype!="strata"]
if (!is.null(yCuts)) {
  newyStrat <- newyStrat[obstype!="strata"]
  newxStrat <- newxStrat[obstype!="strata"]
}
obstype <- obstype[obstype!="strata"]
n <- length(y)

xarray <- array(0,c(n,ncol(X1),2))
xarray[, ,1] <- X1
xarray[, ,2] <- X2

# Function call ------------------------------------------------------------
#print("Function call")

if (errdistrn=="logistic") errdist <- logisterr
   else if (errdistrn=="normal") errdist <- normerr
   else if (errdistrn=="t") errdist <- terr
   else stop(paste("error distribution",errdistrn,"not implemented"))

#cat("Check score and inf:\n")
#check <- divchk(theta=start, loglkfn=locscale.log,modelfn=locscale,errdist=errdist,y=y,x=xarray,wts=w)
#cat("\n")

### obtain estimates ###
res <- mlefn(theta=start, loglkfn=locscale.log, modelfn=locscale, errdist=errdist, errmodpars=errmodpars, 
             y=y, x=xarray, wts=w, control=control, control.inner=control.inner)

if (!is.null(yCuts)) {
  ### obtain covariance matrix with adjustment ###
  mx <- rbind(xarray[,,1],xarray[,,2])
  eta <- matrix(mx %*% as.vector(res$theta),ncol=2)	
  z <- locscale(y,eta,nderivs=2,errdist=errdist, ...)

  ntheta <-length(start)
  dldtheta <- matrix(0,n,ntheta)
  for (j in 1:2)  #for no. of etas
    dldtheta <- dldtheta + diag(w*z$dlogfy[,j]) %*% xarray[,,j]
    
  d2ldthetaP <- NULL; ii <- 0
  vmid <- matrix(0,ntheta,ntheta)
  if (sampling == "ycon") stratum <- as.numeric(interaction(newxStrat,newyStrat))
    else  stratum <- as.numeric(interaction(newyStrat,newxStrat))
  
  for (i in sort(unique(stratum))) {
    indi <- (1:n)[stratum==i]; ii <- ii + 1
    if (length(indi) < 2) qvec[ii] <- NA
     else {
       vmid <- vmid + length(indi)*empvar(dldtheta[indi,])    
       Pi <- -1/qvec[ii]*apply(dldtheta[indi,],2,sum)
       d2ldthetaP <- cbind(d2ldthetaP, Pi)
     }
  }
  
  qvec <- qvec[!is.na(qvec)]
  qmat <- diag(qvec) - qvec %*% t(qvec)
  A <- 1/sum(NMat0)*(d2ldthetaP %*% qmat %*% t(d2ldthetaP))
  
  temp <- solve(res$inf,vmid+A)
  covmat <- solve(res$inf,t(temp))
}
else covmat <- solve(res$inf)

names(res$theta) <- names(res$score) <- dimnames(X1)[[2]]
dimnames(covmat) <- list(dimnames(X1)[[2]],dimnames(X1)[[2]])
dd <- sqrt(diag(covmat))
correlation <- covmat/outer(dd,dd)
pred <- X1Orig %*% res$theta[1:ncol(X1Orig)]
# ------------------------------------------------------------------------------------
ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, yCutsKey=yCutsKey, 
			fit=fit, error=res$error, coefficients=res$theta, loglk=res$loglk, score=res$score, inf=res$inf,
            fitted=pred, cov=covmat, cor=correlation,call=call, assign1=assign1, assign2=assign2,
            fnames1=fnames1, fnames2=fnames2,terms1=terms1, terms2=terms2,
            order1=order1, order2=order2, n1=length(terms1)+1, n2=length(terms2)+1)
class(ans) <- "locsc2stgwtd"
ans
}

#####################################################################################
print.locsc2stgwtd <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#####################################################################################
#n1 <- x$n1; n2 <- x$n2    
cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}
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
if (!is.null(x$yCutsKey)) {
    cat("\nKey to the y-Strat:\n")
    print.default(x$yCutsKey,quote=FALSE,digits=digits)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n\n")

#    cat("Covariance:\n")
#    print.default(format(x$cov, digits = digits), print.gap = 2, quote = FALSE)
#    cat("\n")
}
invisible(x)
}

#####################################################################################
summary.locsc2stgwtd <- function(object,...) {
#####################################################################################

    z <- object
    sigma <- z$cov
    thetas <- z$coefficients    
    assign1 <- z$assign1; assign2 <- z$assign2
    terms1 <- z$terms1; terms2 <- z$terms2
    order1 <- z$order1; order2 <- z$order2
    fnames1 <- z$fnames1; fnames2 <- z$fnames2
    
    # calculate wald statistics

    dowald <- function(assign,terms,thetas,sigma){

    # find dfs for each variable
    num <- -1
    pos <- waldchis <- pvals <- coef.table <- NULL
    for (i in 1:length(assign)) {
   
        if (assign[i] > num) {
            pos <- c(pos, i)
            num <- assign[i]
        }
    }
    pos1 <- c(pos[-1], length(assign) + 1)
    values <- assign[pos]
    dfs <- lengths <- pos1 - pos
    if (length(pos) > 1) {
        for (i in 2:length(pos)) {
            thetai <- thetas[pos[i]:(pos[i] + lengths[i] - 1)]
            covi <- sigma[pos[i]:(pos[i] + lengths[i] - 1), pos[i]:(pos[i] + 
                lengths[i] - 1)]
            waldchi <- t(thetai) %*% solve(covi, thetai)
            waldchis <- c(waldchis, waldchi)
            pval <- (1 - pchisq(waldchi, df = dfs[i]))
            pvals <- c(pvals, pval)
        }
        coef.table <- cbind(dfs[-1], waldchis, pvals)
        dimnames(coef.table)[[1]] <- terms[values]
        dimnames(coef.table)[[2]] <- c("Df", "Chi", "Pr(>Chi)")
    }
coef.table
}

l1 <- length(assign1)
lt <- length(thetas)
ct1 <- dowald(assign1,terms1,thetas[1:l1],sigma[1:l1,1:l1])
ct3 <- dowald(assign2,paste(terms2,"2",sep="."),thetas[(l1+1):lt],sigma[(l1+1):lt,(l1+1):lt])

    # calculate t statistics for non factors ####################
    if(is.null(dim(sigma))) s.err <- sqrt(sigma) else s.err <- sqrt(diag(sigma))

    zstats <- thetas/s.err
    pvals <- (1 - pnorm(abs(zstats)))*2
    coef.table2 <- cbind(thetas, s.err, zstats, pvals)
    dimnames(coef.table2)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

    ct2 <- ct4 <- NULL

    if (length(assign1) > 0) {
        ct2 <- coef.table2[1:length(assign1), ]
    }
    if (length(assign2) > 0) {
        ct4 <- coef.table2[(length(assign1) + 1):length(c(assign1,assign2)), ]
    }
    	
ans <- list(missReport=z$missReport, call=z$call,StrReport=z$StrReport, xStrReport=z$xStrReport,
       key=z$key, yCutsKey=z$yCutsKey, fit=z$fit, error=z$error, loglk=z$loglk, coefficients=z$coefficients,
      coef.table1=ct1,coef.table2=ct2, coef.table3=ct3,coef.table4=ct4)
class(ans) <- "summary.locsc2stgwtd"
ans
}    

#####################################################################################
print.summary.locsc2stgwtd <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#####################################################################################    
cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}
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
if (!is.null(x$yCutsKey)) {
    cat("\nKey to the y-Strat:\n")
    print.default(x$yCutsKey,quote=FALSE,digits=digits)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n")
    if(!is.null(x$coef.table1)){
        cat("\nLocation Model:\n")
        cat("Wald Tests:\n")
        print(x$coef.table1,digits=4)
        cat("\n")
    }
    if(!is.null(x$coef.table2)){
    cat("Coefficients:\n")
    print(x$coef.table2,digits=4)
    } 
    if(!is.null(x$coef.table3) | !is.null(x$coef.table4))
    	   cat("\nScale Model:\n")
    if(!is.null(x$coef.table3)){
    cat("Wald Tests:\n")
        print(x$coef.table3,digits=4)
        cat("\n")
    }
    if(!is.null(x$coef.table4)){
        cat("Coefficients:\n")
        print(x$coef.table4,digits=4)
    }
}
else cat("\nCall requested that model not be fitted\n") 
invisible(x)
}

