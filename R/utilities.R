###############################################################################
printlist <- function(x){
###############################################################################

    namesx <- names(x)
    for (i in 1:length(namesx)){
        print(paste("$",namesx[i],sep=""),quote=FALSE)
        cat("\n")
        print(x[[i]],quote=FALSE)
        if (i<length(namesx)) cat("\n\n")
    }
}

###############################################################################
na.keep <- function(X){X}
###############################################################################

##############################
ArbRound <- function(x,k=10){
##############################
# "Rounds" elements of a vector so that there are (at most)
#  k equally spaced distinct values between the minimun and maximum
#  The minimum and maximum are preserved
   a <- min(x); b <- max(x);
   x <- a + (round(k*(x-a)/(b-a)))*(b-a)/k
   x
}


#######################################################
insertColumns <- function(X,Y,pos){
#######################################################
# inserts columns Y into matrix X, starting at pos

if (nrow(X)!=nrow(Y))
      stop(paste("Row dimensions do not match"))
if ((pos>ncol(X)+1)|(pos<1))
      stop(paste("pos is not in correct range"))
if(pos==ncol(X)+1)  
      cbind(X,Y)	
else {  leftX <- X[,1:pos-1]
	rightX <- X[,pos:ncol(X)]
	cbind(leftX,Y,rightX)}
}


#############################################################################
bind.3darray <- function(a,b)
#############################################################################
{
#   rbind  each b[, ,j] matrix under a[, ,j]
   dima <- dim(a)
   if (is.null(a)) result <- b
   else{
     result <- array(0,dima+c(dim(b)[1],0,0))
     if (length(dima) != 3 | any(dima[2:3] != dim(b)[2:3])) 
         stop("bind.3darray: Objects not compatible 3d arrays")
     for (j in (1:dima[3])) result[, ,j] <- rbind(a[, ,j],b[, ,j])
   }
   result
}



###########################################################################################
array3form <- function(XMat, ClusInd, IntraClus=NULL, MaxInClus=NULL, rmsingletons=FALSE, freqformat=TRUE)
###########################################################################################
{

# Takes cluster data in which rows are individuals and clustermembership is given by the
# column vector ClusInd and makes a 3-D array in which the 1st dim indexes clusters, the 2nd
# dim is the cols of XMat and the 3rd indexes individuals within the cluster.
# The 3rd index runs from 1 to the maximum number of individuals in a cluster
# (or MaxInClus if supplied).  Array positions for non-existent individuals are
# packed out with missing values (NA)

# If MaxInClus is supplied, we take only the first MaxInClus individuals in any cluster


if (!is.matrix(XMat)) stop("array3form: XMat is not a matrix")
XMatcnames <- colnames(XMat)
oo <- if (is.null(IntraClus)) order(ClusInd) else order(ClusInd, IntraClus)
ClusInd <- ClusInd[oo]
XMat <- XMat[oo,]
Jis <- table(ClusInd)
nClus <- length(Jis)
Jmax <- max(Jis)
ncols <- ncol(XMat)

# For clusters with Jis <- Jmax, we will pack out XMat with rows of missing values
# so that each cluster contributes Jmax rows.
# We will then turn the result into a nClus x ncols x Jmax array called XArray


# Construct Addressing vectors iVec and JVec
JMat <- matrix(NA,nrow=nClus,ncol=Jmax,byrow=TRUE)
for (j in 1:Jmax){
     tempvecj <- rep(NA,Jmax)
     tempvecj[1:j] <- 1:j
     if (any(Jis==j)) {
        JMat[Jis==j,] <- outer(rep(1,sum(Jis==j)),tempvecj)
     }
}

JVec <- as.vector(t(JMat)) ; JVec <- JVec[!is.na(JVec)]
iVec <- factor(ClusInd); levels(iVec) <- 1:nClus; iVec <- as.integer(iVec)
# First XArray is a matrix packed out with rows of NA's
XArray <- matrix(NA,nrow=nClus*Jmax,ncol=ncols)
XArray[(iVec-1)*Jmax+JVec,] <- XMat   # XMat expanded to equal no. of rows per cluster
# Now turn it into the desired array
XArray <- array(t(XArray),c(ncols,Jmax,nClus))
XArray <- aperm(XArray,c(3,1,2))

ClusLab <- unique(ClusInd)  # Start to store Cluster labels from ClusInd
if (rmsingletons) {
    XArray <- XArray[Jis>1, , ]
    ClusLab <- ClusLab[Jis>1]
    Jis <- Jis[Jis>1]
    nClus <- length(Jis)
}
if (!is.null(MaxInClus)){
    Jmax <- min(Jmax,MaxInClus)
    Jis[Jis>MaxInClus] <- MaxInClus
    XArray <- XArray[, ,1:Jmax]
}
counts <- rep(1,nClus)
if (freqformat) {
    # Turn XArray into nClus rows so we can condense it into a frequency format
    XMatWide <- matrix(aperm(XArray,c(2,3,1)),nrow=nClus,ncol=Jmax*ncols,byrow=TRUE)
    z <- weightform(XMatWide)
    Jis <- Jis[z$ind]; counts <- z$counts; nClus <- length(z$counts); ClusLab <- 1:nClus
    XArray <- array(XMatWide[z$ind,],c(nClus,ncols,Jmax))
}
dimnames(XArray) <- list(ClusLab,XMatcnames,paste("eta",1:Jmax,sep="."))
list(XArray=XArray,Jis=Jis,counts=counts)

} # Endfn array3form


############################################
replaceNA <- function(x,with="miss"){
############################################
# Value is x with any missing values replaced by value of with
# Works for vectors, whether numeric, character or factor
# Work with any matrices I've tried it on too (CW)

  if (is.factor(x)){
      levsx <- levels(x)
      x <- as.character(x)
      x[is.na(x)] <- with
      x <- factor(x,levels=c(levsx,with))
  }
  else x[is.na(x)] <- with
  x
}


##############################################################################	
weightform <- function(X, counts=NULL, counts2=NULL){
##############################################################################	
# take dataframe or matrix X
#       and return ids of unique rows in "ind" and the number of repetitions in "wts"
# counts is a separate counts vector going along with the rows of the matrix
#  (if you have a counts vector in X, remove it before sending it to weightform()

    nr <- nrow(X)
    nc <- ncol(X)
    if (is.data.frame(X) | is.character(X)) {
         X <- as.data.frame(X)
         for (j in 1:nc) X[,j] <- as.numeric(X[,j])
         X <- as.matrix(X)
    }
    missing.replace <- max(X, na.rm=TRUE) + 9
    X[is.na(X)] <- missing.replace
    if (is.null(counts)) counts <- rep(1,nr)
    sumcounts <- sum(counts)
    arg <- "X[,1]"
    if (nc > 1) for (j in (2:nc)) arg <- paste(arg,",X[,",j,"]",sep="")
    o <- eval(parse(text=paste("order(",arg,")",sep="")))

    X <- X[o,]
    counts <- counts[o]
    if (!is.null(counts2)) counts2 <- counts2[o]
     orig <- (1:nr)[o]

    wts2 <- NULL
    ind <- 1
    wts <- counts[1]
    if (!is.null(counts2)) wts2 <- counts2[1]
    for (j in 2:nr) {
        if (sum(abs(X[j,]-X[j-1,])) != 0) {
             ind <- c(ind,j)
             wts <- c(wts,counts[j])
             if (!is.null(counts2)) wts2 <- c(wts2,counts2[j])
        }
        else {
              wts[length(wts)] <- wts[length(wts)] + counts[j]
              if (!is.null(counts2)) wts2[length(wts2)] <- wts2[length(wts2)] + counts2[j]
        }
    }

    list(o=o,ind=orig[ind],counts=wts,counts2=wts2)
}

##############################################################################	
dfcompress <- function(X,counts.name=NULL){
##############################################################################	
# Compresses a dataframe or matrix to the unique rows plus frequencies
# If a non-null counts.name is given, a column within X with that name
# is treated as already containing frequencies. These are then accumulated.
# If counts.name is null, a final column called "counts" is added to X

    nc <- ncol(X)
    replacecol <- FALSE
    Xsend <- X
    if (!is.null(counts.name)){
        counts.col <- match(counts.name,dimnames(X)[[2]])
        if (is.na(match(counts.col,1:nc))) stop("no column found with that name")
        if (!is.numeric(X[,counts.col])) stop("Non-numeric column specified for counts")
        counts <- X[,counts.col]
        Xsend <- X[,-counts.col]
        replacecol <- TRUE
    }
    else{
        counts.name <- "counts"
        counts <- rep(1,nrow(X))
    }
    z <- weightform(Xsend,counts=counts)
    counts <- z$counts
    if (replacecol){
        X <- X[z$ind,]
        X[,counts.col] <- counts
    }
    else X <- data.frame(X[z$ind,],counts)
    X
}

###################################
concat.dfs <- function(df1,df2){
###################################

cnames <- unique(c(names(df1),names(df2)))
nr1 <- dim(df1)[1]
nr2 <- dim(df2)[1]

merged <- NULL
icol <- 0
for (i in 1:length(cnames)){
     i1 <- match(cnames[i],names(df1))
     i2 <- match(cnames[i],names(df2))
     if (!is.na(i1) & !is.na(i2)){
         icol <- icol + 1
         coli <- if (is.numeric(df1[,i1]) & is.numeric(df2[,i2]))
c(df1[,i1],df2[,i2])
                  else
factor(c(as.character(df1[,i1]),as.character(df2[,i2])))
     }
     else if (is.na(i1)){
         icol <- icol + 1
         coli <- if (is.numeric(df2[,i2]))   c(rep(NA,nr1),df2[,i2])
                  else factor(c(rep(NA,nr1),as.character(df2[,i2])))
     }
     else {
         icol <- icol + 1
         coli <- if (is.numeric(df1[,i1]))   c(df1[,i1],rep(NA,nr2))
                  else factor(c(as.character(df1[,i1]),rep(NA,nr2)))
     }
     merged <- if (i==1) data.frame(coli) else data.frame(merged,coli)
     dimnames(merged)[[2]][icol] <- cnames[i]
}
merged
}


##############################################################################	
weightformvec <- function(y){
##############################################################################	
# take vector y and return unique values and the number of repetitions
   wts <- table(y)
   yuse <- as.numeric(names(wts))
   names(yuse) <- NULL
   names(wts) <- 1:length(wts)
   list(y=yuse,wts=wts)
}


##############################################################################	
weightformmat <- function(X,missing.replace= -999){
##############################################################################	
# take matrix X and return ids of unique rows in "ind" and the number of repetitions in "wts"
# Since numerical comparisons of subsequent rows are done, we replace NA's by missing.replace

   nr <- nrow(X)
   nc <- ncol(X)
   X[is.na(X)] <- missing.replace
    arg <- "X[,1]"
   if (nc > 1) for (j in (2:nc)) arg <- paste(arg,",X[,",j,"]",sep="")
   o <- eval(parse(text=paste("order(",arg,")",sep="")))

   ind <- o[1]
   wts <- 1
   for (j in 2:length(o))
        if (sum(abs(X[o[j],]-X[o[j-1],])) != 0) {
                   ind <- c(ind,o[j])
                   wts <- c(wts,1)
                   }
        else wts[length(wts)] <- wts[length(wts)] + 1

   if (sum(wts) != nr) print("***weightformmat:  ERROR ***")

   list(o=o,ind=ind,wts=wts)
}

#############################################################################
matincltotals <- function(mat){
#############################################################################
   if (!is.matrix(mat)) stop("matincltotals: Argument not a matrix")
   if (length(dimnames(mat))<2) dimnames(mat) <- 
                  list(as.character(1:nrow(mat)),as.character(1:ncol(mat)))
   mat2 <- cbind(mat,apply(mat,1,sum))
   mat2 <- rbind(mat2,apply(mat2,2,sum))
   dimnames(mat2) <- list(c(dimnames(mat)[[1]],"Total"),c(dimnames(mat)[[2]],"Total"))
   mat2
   }



##############################################################################	
rbivnorm<- function(n,rho,mux=0,sdx=1,muy=0,sdy=1)
##############################################################################	
{
x<-rnorm(n)
y<- rho*x + sqrt(1-rho^2)*rnorm(n)
x <- mux + sdx*x
y <- muy + sdy*y
cbind(x,y)
}

############################################################
simsummary <- function(title.pr,thetas,vars,theta,nround=5){
############################################################
  cat(title.pr,"\n")
  meanest <- apply(thetas,2,mean,na.rm=TRUE)
  nthetas <- apply(1*!is.na(thetas),2,sum)
  cat("mean   =",round(meanest,nround),"\n")
  varest <- apply(thetas,2,var,na.rm=TRUE)
  cat("var    =",round(varest,nround),"\n")
  cat("AvVar  =",round(apply(vars,2,mean,na.rm=TRUE),nround),"\n")
  nvars <- apply(1*!is.na(vars),2,sum)
  sds <- sqrt(vars)
  cat("AvSd   =",round(apply(sds,2,mean,na.rm=TRUE),nround),"\n")
  cat("Bias z =", round((meanest-theta)/sqrt(varest/nvars),2),"\n")
  cat("nthetas=", nthetas,"\n")
  cat("nvars  =", nvars,"\n")

  cin <- rep(NA,length(theta))
  for (j in (1:length(theta))) {
       diffs <- (thetas[,j]-theta[j])/sds[,j]
       cin[j] <- length(diffs[abs(diffs) < 1.96])/nrow(thetas)
       }
  print(c("coverage",round(100*cin,1)))
  cat("....................................................................")
  cat("\n")
  varest
}


####################################################################	
oddsratio <- function(y,wts=NULL){
####################################################################

# calculate the oddsratio of two 0/1 binary variables 
# need to make some improvement
   
   if (!is.matrix(y)) stop("y should be a matrix !")
   if (dim(y)[2] != 2) stop("y should have two columns indicating two variables!")

   y1 <- as.numeric(factor(y[,1]))
   y2 <- as.numeric(factor(y[,2]))
   if (length(unique(y1))>2 || length(unique(y2))>2) 
	stop("both variables should be binary !")
   y <- cbind(2-y1,2-y2)   
   if (is.null(wts)) wts <- rep(1,nrow(y))
	
   y11 <- y12 <- y21 <- y22 <- 0
   for (i in 1:nrow(y)){
      if (sum(y[i,]==c(0,0))==2)  y11 <- y11 + wts[i]
      else if (sum(y[i,]==c(0,1))==2)  y12 <- y12 + wts[i]
           else if (sum(y[i,]==c(1,0))==2)  y21 <- y21 + wts[i]
                else y22 <- y22 + wts[i]
   }

   oddsR <- (y11 * y22)/(y12 * y21)
}


######################################
Unidata <- function(vmat,wts=NULL) {
######################################

# a function used to set up a data matrix with each row representing a 
# unique set of values, with an extra last column giving the total number 
# of observations falling into each row category (wts);

vmat <- as.matrix(vmat)
if (any(is.na(vmat)))
  stop("There are NA in vmat!")
  
if (is.null(wts)) wts <- rep(1,nrow(vmat))
else {
   wts <- as.vector(wts)
   if (any(is.na(wts))) stop("There are NA in wts!")
     else if (length(wts) != nrow(vmat)) stop("wts has different length from vmat!")
}

nvar <- ncol(vmat)
nvalue <- rep(NA,nvar)
vvalue <- NULL
for (i in 1:nvar) {
  vvalue <- c(vvalue,list(sort(unique(vmat[,i]))))
  nvalue[i] <- length(vvalue[[i]])
}

torow <- prod(nvalue)
nvmat <- matrix(NA,nrow=torow,ncol=nvar)
nvmat[,1] <- rep(vvalue[[1]],torow/nvalue[1])
if (nvar >1) {
  for (j in 2:nvar)
    nvmat[,j] <- rep(rep(vvalue[[j]],rep(prod(nvalue[1:(j-1)]),nvalue[j])),
                   torow/prod(nvalue[1:j]))
}
w <- rep(0,torow)
for (k in 1:torow)
  for (i in 1:nrow(vmat))
    if (all(vmat[i,]==nvmat[k,])) w[k] <- w[k] + wts[i]
nvmat <- cbind(nvmat,w)
if (!is.null(dimnames(vmat)[[2]]))
  dimnames(nvmat) <- list(NULL,c(dimnames(vmat)[[2]],"wts"))
nvmat 
}


###########################################################
divchk <- function(theta,loglkfn,nderivs=2, ...){
###########################################################

# this function is used to check derivatives

   ntheta <- length(theta)
   if (nderivs < 1) stop("divchk: asked for no derivatives")

   w <- loglkfn(theta, nderivs, ...)

   score.est <- rep(0,ntheta)

   for (ii in 1:ntheta) {
       ei <- rep(0,ntheta)
       ei[ii] <- 1
       hi <- max(1,abs(theta[ii]))*1e-5

       score.est[ii] <- (loglkfn(theta+hi*ei, nderivs=0, ...)$loglk -
             loglkfn(theta-hi*ei, nderivs=0, ...)$loglk)/(2*hi)
       }
   print("score")
   print(w$score)
   print("est. score")
   print(score.est)
   print("diff")
   print(round(w$score-score.est,8))

   if (nderivs < 2) return(invisible(w))

   print("")

   inf.est <- matrix(0,ntheta,ntheta)

   for (ii in 1:ntheta) {
       ei <- rep(0,ntheta)
       ei[ii] <- 1
       hi <- max(1,abs(theta[ii]))*1e-5

       inf.est[ii,] <- -(loglkfn(theta+hi*ei, nderivs=1, ...)$score -
             loglkfn(theta-hi*ei, nderivs=1, ...)$score)/(2*hi)
       }

   print("inf")
   print(w$inf)
   print("est. inf")
   print(inf.est)
   print("diff")
   print(round(w$inf-inf.est,8))
   print("")
  invisible(w)
}


######################################################
empvar <- function(x,nreps=rep(1,nrow(x))){
######################################################
# variance of an n*p matirx allowing for replicated rows
# x should be a column matrix
# applied for weighted method in variance adjustment

  n <- sum(nreps)
  diff <- (x-outer(rep(1,nrow(x)),apply(nreps*x,2,sum))/n)
  (t(diff)%*%diag(nreps)%*%diff)/(n-1)
}

