############################################################################################
bin2stg <- function(formula, weights=NULL, xstrata=NULL, obstype.name="obstype", data,
               fit=TRUE, xs.includes=FALSE, linkname="logit", start=NULL, Qstart=NULL,
               int.rescale=TRUE, off.set=NULL,
               control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
############################################################################################
{
# Process formula -------------------------------------------------------------------

mf <- call <- match.call()
mf[[1]] <- as.name("model.frame")
mf$na.action <- as.name("na.keep")
NULL -> mf$xstrata -> mf$obstype.name -> mf$fit -> mf$xs.includes -> mf$linkname 
NULL -> mf$start -> mf$Qstart -> mf$int.rescale -> mf$off.set -> mf$control -> mf$control.inner
 
mf <- eval(mf, sys.frame(sys.parent()))
y <- model.extract(mf,"response")
n <- if (is.matrix(y)) dim(y)[1] else length(y)
Terms <- attr(mf,"terms")
X <- model.matrix(Terms,mf)
if(!is.matrix(X)) X <- matrix(X)
XOrig <- X
w <- model.extract(mf,"weights")
if(is.null(w)) w <- rep(1,n)
terms1 <- attr(terms(formula),"term.labels")
order1 <- attr(terms(formula),"order")
assign1 <- attr(X,"assign"); fnames<-names(attr(X,"contrasts"))

# Input checking ---------------------------------------------------------------------

if (is.na(match(obstype.name,names(data)))) stop(paste("Dataframe did not have a column named",obstype.name))
obstype <- data[,match(obstype.name,names(data))]
if (any(is.na(data[,obstype.name])))
    stop(paste("Missing values not permitted in obstype variable",obstype.name))
    
prospective <- all(obstype %in% c("uncond","y|x"))
# for S may need prospective <- all(!is.na(match(obstype,c("uncond","y|x"))))

if (any(is.na(w))) stop("Weights should not contain NAs")


# Form the off.set -----------------------------------------------------------------

if (is.null(off.set)) off.set <- rep(0, n)
    else off.set <- as.vector(off.set)


# Form the x-strata -----------------------------------------------------------------

xStrat <- if (is.null(xstrata)) rep(1,n) else {
    xstratform <- as.formula(paste("~","-1 + ",paste(xstrata,sep="",collapse=":")))
    data1 <- as.data.frame(lapply(data, factor))
   	xstratmat <- model.matrix(xstratform, model.frame(xstratform,data1,na.action=function(x)x))
	   for(i in 1:ncol(xstratmat)) xstratmat[,i] <- i * xstratmat[,i]
		  apply(xstratmat,1,sum)
}
nStrat <- max(xStrat)


# Form the y-vector -----------------------------------------------------------------

matrixY <- FALSE

if ((is.vector(y)|is.factor(y))&!is.list(y)) {
    yfac<-factor(y)
    if (length(levels(yfac))>1){
        y <- 1 + 1*((yfac)==levels(yfac)[1])
        y1name <- if (length(levels(yfac))==2)  as.character(levels(yfac))[2]
                  else paste("Not.",as.character(levels(yfac))[1],sep="")
    }
    else if (length(levels(yfac))==1){
        y <- 1*((yfac)==levels(yfac)[1])
        y1name <- c(as.character(levels(yfac)),NA)
    }
    else stop("Illegal y-data")
}
else if (is.matrix(y)) {# expand out into vector form
     matrixY <- TRUE
     y1name <- colnames(y)[1]
     cnames <- dimnames(X)[2]
     X <- kronecker(X,c(1,1))
     dimnames(X)[2] <- cnames
     w <- as.vector(t(y*w))
     obstype <- factor(as.vector(t(as.matrix(data.frame(obstype,obstype)))))
     off.set <- as.vector(rbind(off.set, off.set))
     xStrat <- as.vector(rbind(xStrat,xStrat))
     y <- rep(c(1,2),n)
     n <- 2*n
}
else stop("Illegal y-data")

off.set <- matrix(off.set, n, 1)


# Remove random missing values -------------------------------------------------------

mat <- cbind(y,X,xStrat)
miss1 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="uncond"]
miss2 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="retro"]
miss3 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="y|x"]
mat <- cbind(X,xStrat)
miss4 <- (1:n)[apply(mat,1,function(x) any(is.na(x))) & obstype=="xonly"]
mat <- cbind(y,xStrat)
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
	   obstype <- factor(as.character(obstype[-toremove]))
	   y <- y[-toremove]
	   xStrat <- xStrat[-toremove]
	   X <- X[-toremove,,drop=FALSE]
           off.set <- off.set[-toremove,,drop=FALSE]
	   n <- length(y)
}

# Stratum Counts Report -------------------------------------------------------------
#cat("\nStrata Counts Report:\n")

StrReport <- xStrReport <- key <- NULL
if (matrixY) yKey <- c("y",y1name)
   else yKey <- c(as.character(Terms[[2]]), y1name)

sub <- (1:length(y))[obstype=="xonly"]
if(length(sub)>0){
	   df <- data.frame(w, obstype=factor(as.character(obstype)), y, xStrat)[-sub,]
   	df1 <- data.frame(w, obstype=factor(as.character(obstype)), y, xStrat)[sub,]
}
else df <- data.frame(w,obstype=factor(as.character(obstype)),y,xStrat) 

StrReport <- ftable(xtabs(w~obstype+y+xStrat,data=df))
if(exists("df1")) xStrReport <-  xtabs(w~obstype+xStrat,data=df1)

if (xs.includes == TRUE) {
    oblevel0 <- levels(factor(as.character(obstype)))
    oblevel <- rep(oblevel0, rep(2,length(oblevel0)))  #y is binary

    report <- ftable(xtabs(w~obstype+y+xStrat,data=df[df$obstype=="strata",]))
    strep <- xtabs(w~y+xStrat,data=df[df$obstype!="strata",])
    report[oblevel=="strata"] <- strep
    StrReport <- StrReport-report

    if (exists("df1")) {
        report1 <- xtabs(w~obstype+xStrat,data=df1[df1$obstype=="strata",])
        strep1 <- xtabs(w~xStrat,data=df1[df1$obstype!="strata",])
        report1[oblevel0=="strata"] <- strep1
        xStrReport <- xStrReport-report1
    }
}

if (!(is.null(xstrata))) {
    key <- dimnames(xstratmat)[[2]]
    names(key) <- 1:length(key)
}

if (!fit) {
    ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, 
     		fit=fit, call=call, assign1=assign1, fnames=fnames, terms1=terms1, order1=order1)
    class(ans) <- "bin2stg"
    return(ans) # RETURN AT THIS POINT IF fit IS FALSE
}

# ------------------------------------------------------------------------------------
rmat <- Qmat <- NULL

if (!prospective) {
   	havestrata <- length(obstype[obstype=="strata"])>0
   	rk <- rep(0,n)
   	if ((havestrata & xs.includes) | !havestrata) rk[obstype=="retro"] <- -1 * w[obstype=="retro"]
   	if  (havestrata & xs.includes) rk[obstype=="uncond"] <- -1 * w[obstype=="uncond"]
   	rk[obstype=="strata"] <- 1 * w[obstype=="strata"]


   	rmat <- xtabs(rk~y+xStrat)
   	rmat <- matrix(rmat,dim(rmat)[1],dim(rmat)[2],dimnames=dimnames(rmat)[1:2])
   	if (dim(rmat)[1] == 1){ # only one level of y represented in the data
      	  if (dimnames(rmat)[[1]]=="1") rmat <- rbind(rmat,rep(0,dim(rmat)[2]))
      		  else rmat <- rbind(rep(0,dim(rmat)[2]),rmat)
	   }
   	if (is.null(Qstart)){
   	#    ind <- if (havestrata & !xs.includes) (1:n)[obstype=="uncond"|obstype=="strata"]
   	     ind <- if (havestrata & !xs.includes) (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="strata"]
                              else (1:n)[obstype=="strata"]
       	 temp <- xtabs(w~y+xStrat,data=data.frame(y,xStrat,w,obstype)[ind,])

      	  if (is.matrix(temp) && all(dim(temp) == c(2,nStrat)) && min(temp)>0)
       	                 Qmat <- sweep(temp,2,apply(temp,2,sum),"/")
      		 else  stop("Attempt to construct missing Qstart failed")
   	}
   	else { # Qstart is present
      	  	if(is.matrix(Qstart) && all(dim(Qstart) == c(2,nStrat)) &&  # matrix with correct dimensions
         	  (round(apply(Qstart,2,sum),5) == rep(1,dim(Qstart)[2]))) # cols sum to 1
              		Qmat <- Qstart
        	else stop("illegal Qstart value")
   	}
}
# Create start if needed ------------------------------------------------------------
if (is.null(start)){ # use glm() to obtain starting values
     offset <- rep(0,n)
     if (any(obstype %in% "retro")) {
        	nMat <- xtabs(w~y+xStrat,data=data.frame(y,xStrat,w,obstype)[obstype=="retro",]) 
        	nMat <- matrix(nMat,dim(nMat)[1],dim(nMat)[2],dimnames=dimnames(nMat)[1:2])
        	offsets0 <- log( (nMat[1,]/nMat[2,])/(Qmat[1,]/Qmat[2,]))
        	offset <- ifelse((obstype=="retro"), offsets0[xStrat], 0)
     }
        
     if (linkname=="probit") offset <- 0.43*offset
     #WHAT HERE??????  if (linkname=="cloglog") offset <- ??????
     
     ind <- (1:n)[obstype=="uncond" | obstype=="retro" | obstype=="y|x"]
     xarray <- array(X,c(dim(X),1))
     if (length(levels(factor(y[ind])))<2 | sum(w[ind]) < 50)
     	        stop("Cannot construct starting values from this data")
	  # Note that Splus 3 is ignoring the offset, needs to be fixed
     templink <- if (linkname=="cloglog") "logit" else linkname
     # THERE IS AN R BUG here with glm with cloglog so use logit instead
     start <- glm((y==1)[ind]~xarray[ind, ,1]-1,weights=w[ind],
                   family=eval(parse(text=paste("binomial(link=", templink,")",sep=""))),
                   offset=offset[ind])$coefficients
}
else if (length(start) != dim(X)[2]){
   	# Match up any values in start with names in X, fill in holes with 0's.
   	newstart <- rep(0,dim(X)[2])
   	posns <- match(names(start),dimnames(X)[[2]])
    newstart[posns] <- start
   	start <- newstart
}


# Rescale variables ----------------------------------------------------------
if (int.rescale){  # (standardise the variables
    scaleval <- apply(X,2,FUN=function(x)sd(x,na.rm=TRUE))
    constant <- (scaleval==0)
    loc   <- apply(X,2,FUN=function(x)mean(x,na.rm=TRUE))
    loc[constant] <- 0
    scaleval[constant] <- 1
    nrowsx <- nrow(X)
    ncX <- ncol(X)
    X <- (X - outer(rep(1,nrowsx),loc)) / outer(rep(1,nrowsx),scaleval)
    start <- start*scaleval + constant*sum(start*loc)
}
# Form xarray ----------------------------------------------------------------
X <- X[obstype != "strata",,drop=FALSE]
w <- w[obstype != "strata"]
y <- y[obstype != "strata"]
off.set <- off.set[obstype != "strata",,drop=FALSE]
xStrat <- xStrat[obstype != "strata"]
obstype <- obstype[obstype != "strata"]
n <- length(y)
xarray <- array(X,c(n,ncol(X),1))

Aposn <- (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="y|x"]
Acounts <- w[Aposn]
if (!prospective) {
	Bposn <- (1:n)[obstype=="uncond"|obstype=="retro"|obstype=="xonly"]
	Bcounts <- w[Bposn]
}
else Bposn <- Bcounts <- numeric(0)

# ----------------------FUNCTION CALL----------------------------------------------
if (linkname=="logit") {modelfn<-binlogistic; stratfn<-binlogisticstrat}
else if (linkname=="probit") {modelfn<-binprobit; stratfn<-binprobitstrat}
else if (linkname=="cloglog") {modelfn<-bincloglog; stratfn<-bincloglogstrat}
else stop(paste("linkname",linkname,"not implemented"))

res <- mlefn(theta=start,MLInf,y=as.matrix(y),x=xarray,ProspModInf=MEtaProspModInf,
	     StratModInf=MEtaStratModInf, modelfn=modelfn, stratfn=stratfn,
 	     Aposn=Aposn,Acounts=Acounts,Bposn=Bposn,Bcounts=Bcounts,rmat=rmat,Qmat=Qmat,
             xStrat=as.numeric(xStrat), off.set=off.set, 
	     control=control, control.inner=control.inner)   

if (int.rescale){  # Translate back to the scale of the original variables
    res$theta <- res$theta/scaleval - constant*sum((res$theta/scaleval)*loc)
    jacbn <- diag(scaleval)
    jacbn[constant,] <- jacbn[1:ncX,1:ncX][constant,] + loc
    res$score <- as.vector(t(jacbn)%*%res$score)
    res$inf <- t(jacbn)%*%res$inf%*%jacbn
}

names(res$theta) <- names(res$score) <- dimnames(X)[[2]]
dimnames(res$inf) <- list(dimnames(X)[[2]],dimnames(X)[[2]])
#covmat <- solve(res$inf) # change to avoid scale problems with inverse
temp <- 1/sqrt(diag(res$inf))
covmat <- diag(temp)%*%solve(diag(temp)%*%res$inf%*%diag(temp))%*%diag(temp)
dimnames(covmat) <- dimnames(res$inf)
dd <- sqrt(diag(covmat))
correlation <- covmat/outer(dd, dd)
eta <- XOrig %*% res$theta
pred <- NULL
if (linkname=="logit") pred <- plogis(eta)
else if (linkname=="probit") pred <- pnorm(eta)
else if (linkname=="cloglog") pred <- 1 -exp(-exp(eta))
else stop("Link not implemented") # this should have already stopped much earlier
# ------------------------------------------------------------------------------------
ans <- list(missReport=missReport, StrReport=StrReport, xStrReport=xStrReport, key=key, yKey=yKey, fit=fit, 
		     	error=res$error,coefficients=res$theta, loglk=res$loglk, score=res$score, inf=res$inf,
            fitted=pred, cov=covmat, cor=correlation,
            Qmat = res$extra$Qmat,
            call=call, assign1=assign1, fnames=fnames, terms1=terms1, order1=order1)
class(ans) <- "bin2stg"
ans
}

#####################################################################################
print.bin2stg <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#####################################################################################    

cat("\nCall:\n", deparse(x$call), "\n", sep = "")

if (!is.null(x$missReport)){
	   cat("\nObservations deleted due to missing data:")
    print(x$missReport,quote = FALSE, row.names=FALSE)
}
cat("\nStratum Counts Report:\n")
print(x$StrReport,quote = FALSE, row.names=FALSE)
if (!is.null(x$xStrReport)) {
    cat("\nObservations of obstype==xonly\n")
    print(x$xStrReport)
}

cat("\nModel for prob of ", x$yKey[1],"=",x$yKey[2]," (y=1)  given covariates\n",sep="")
if (!is.null(x$key)) {
    cat("\nKey to x-Strat:\n")
    print.default(x$key,quote=FALSE, row.names=FALSE)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
        quote = FALSE)
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n\n")
#   cat("\n")
#   cat("Covariance:\n")
#   print.default(format(x$cov, digits = digits), print.gap = 2, quote = FALSE)
}
invisible(x)
}

#####################################################################################
summary.bin2stg <- function(object,...) {
#####################################################################################

# Works for R only so far

    z <- object
    sigma <- z$cov
    thetas <- z$coefficients    
    assign1 <- z$assign1
    terms1 <- z$terms1
    order1 <- z$order1
    fnames <- z$fnames
 
    # calculate wald statistics

    # find dfs for each variable
    num <- -1
    pos <- waldchis <- pvals <- coef.table1 <- NULL

    for (i in 1:length(assign1)) {
    if (assign1[i] > num) { 
       pos <- c(pos,i) 
       num <- assign1[i] 
       }
    }
    pos1 <- c(pos[-1],length(assign1)+1)
    values <- assign1[pos]
    
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

    coef.table1 <- cbind(dfs[-1],waldchis,pvals)
    dimnames(coef.table1)[[1]] <- terms1[values]
    dimnames(coef.table1)[[2]] <- c("Df","Chi", "Pr(>Chi)")
    }

    # calculate t statistics for non factors ####################
    if(is.null(dim(sigma))) s.err <- sqrt(sigma) else s.err <- sqrt(diag(sigma))

    zstats <- thetas/s.err
    pvals <- (1 - pnorm(abs(zstats)))*2
    coef.table2 <- cbind(thetas, s.err, zstats, pvals)
    dimnames(coef.table2)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

ans <- list(missReport=z$missReport, call=z$call,StrReport=z$StrReport, xStrReport=z$xStrReport,
       key=z$key, yKey=z$yKey, fit=z$fit, error=z$error, loglk=z$loglk, coefficients=z$coefficients,
       coef.table1=coef.table1,coef.table2=coef.table2)
class(ans) <- "summary.bin2stg"
ans
}    

#####################################################################################
print.summary.bin2stg <- function (x, digits = max(3, getOption("digits") - 3), ...) {
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

cat("\nModel for prob of ", x$yKey[1],"=",x$yKey[2]," (y=1)  given covariates\n",sep="")
if (!is.null(x$key)) {
    cat("\nKey to x-Strat:\n")
    print.default(x$key, quote=FALSE, row.names=FALSE)
}

if (x$fit) {
   if (x$error != 0) print("WARNING, FIT UNSUCCESSFUL")
    cat("\nloglikelihood =",x$loglk, "  using", length(coef(x)), "parameters\n\n")    
    if(!is.null(x$coef.table1)){
    cat("Wald Tests:\n")
    print(x$coef.table1,digits=4)
    cat("\n")
    }
    if(!is.null(x$coef.table2)){
    cat("Coefficients:\n")
    print(x$coef.table2,digits=4)
    }   
}
invisible(x)
}
