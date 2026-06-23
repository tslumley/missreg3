##########################################################################################
mlefn.control <-
    function(messg="", niter=20, tol=1e-8, guide="uphill",print.progress=2,
            max.eigenrat = 0.05, n.earlyit=0,
            constrain="no",fixed=NULL,Aconstrain=NULL,cconstrain=NULL, ...)
##########################################################################################
{
list(messg=messg, niter=niter, tol=tol, guide=guide, print.progress=print.progress,
           max.eigenrat=max.eigenrat, n.earlyit=n.earlyit, constrain=constrain,
            fixed=fixed, Aconstrain=Aconstrain, cconstrain=cconstrain  )
}
##########################################################################################
mlefn.control.inner <-
    function(messg="Inner:", niter=20, tol=1e-8, guide="auto",print.progress=0,
            max.eigenrat = 0.05, n.earlyit=0,
            constrain="no",fixed=NULL,Aconstrain=NULL,cconstrain=NULL, ...)
##########################################################################################
{
list(messg=messg, niter=niter, tol=tol, guide=guide, print.progress=print.progress,
           max.eigenrat=max.eigenrat, n.earlyit=n.earlyit, constrain=constrain,
            fixed=fixed, Aconstrain=Aconstrain, cconstrain=cconstrain  )
}
#######################################################################################
qr.solve.cw <- function(a, b, tol = 1e-7, messg=""){
#######################################################################################
#This is a hack at qr.solve, which is the engine of solve, to report errors
    if( !is.qr(a) )	a <- qr(a, tol = tol)
    nc <- ncol(a$qr)
    errorcode <- 0
    if( a$rank != nc ) errorcode <- 1
#                     	stop("singular matrix `a' in solve")
    if( missing(b) ) {
	              if( nc != nrow(a$qr) ) errorcode <- 2
#	                  stop("only square matrices can be inverted")#Doesn't apply here
	              b <- diag(1,nc)
        }
    soln <- NULL
    if (errorcode == 0) soln <- qr.coef(a,b)
    else cat(messg, "singular matrix `a' in solve","\n")
    list(soln=soln,err=errorcode)
}

#######################################################################################
mlefn <- function(theta,loglkfn, control=mlefn.control(...), ...){
#######################################################################################
      greenstadt.step <- function(inf,score,max.eigenrat,iter,notuphill=FALSE,messg=messg,print.progress=1){
          # Take a step of the form described on p 601 of Wild and Seber, 1989
          # All eigenvalues of inf that are less than max(abs(eigs$values))*max.eigenrat
          #    are reset to that value
              eigs <- eigen(inf,symmetric=TRUE)
              if (notuphill & print.progress>=1)
                  cat(messg,"mlefn: eigenvalue ajustment on iteration",iter, "\n",
                  "eigenvals = ", eigs$values, "\n")
              min.evaluse <- max(abs(eigs$values))*max.eigenrat
              mod.evals <- ifelse(eigs$values>min.evaluse,eigs$values,min.evaluse)
              step <- as.vector((eigs$vector%*%((1/mod.evals)*t(eigs$vector)))%*% score)
              step
             }

messg <- control$messg
niter <- control$niter; tol <- control$tol; guide <- control$guide; print.progress <-
control$print.progress; max.eigenrat <- control$max.eigenrat; n.earlyit <- control$n.earlyit; constrain <-
control$constrain; fixed <- control$fixed; Aconstrain <- control$Aconstrain; cconstrain <- control$cconstrain


if (guide != "no") reverse <- ifelse(guide=="uphill",1,-1) # -1 for downhill
theta0 <- theta
ntheta <- length(theta)

Z <- NULL
if(constrain != "no"){
       if (constrain=="fix"){
             fixed <- sort(fixed)
             Aconstrain <- (diag(1,ntheta))[fixed,]
             if (length(fixed) == 1)  Aconstrain <- matrix(Aconstrain,1,ntheta)
             cconstrain <- theta[fixed]
             }
       nrA <- nrow(Aconstrain)
       ncA <- ncol(Aconstrain)
       if (nrA >= ntheta) stop("mlefn: too many constraints")
       if (ncA != ntheta) stop("mlefn: number of constraint cols .ne. ntheta")
       qrA <- qr(t(Aconstrain))
       if (qrA$rank<nrA) stop("mlefn: constraint matrix not full row rank")
       Q <- qr.Q(qrA,complete=TRUE)
       R <- qr.R(qrA,complete=FALSE)
       S <- Q[,1:nrA,drop=FALSE]%*%t(backsolve(R,diag(rep(1,nrA))))
       Z <- Q[,(nrA+1):ncA,drop=FALSE]
       theta0 <- S%*%cconstrain + (Z%*%t(Z))%*%(theta0-S%*%cconstrain)
                       #  Force theta0 to satisfy the constraints
       }

conv <- FALSE

extra <- NULL  # extra is a device that allows loglkfn to pass information it
               # constructed last iteration back to itself for the next iteration

for (j in 1:niter) {
    counter <- j; 
    ww <- loglkfn(theta0, nderivs=2, extra=extra, ...)
    if (!is.null(ww$error)) if (ww$error != 0){
           if (print.progress>=1)
                    cat(messg,"mlefn: loglk evaluation failed at iteration",j,"-- error =1","\n")
           return(list(theta=theta0,error=1))
           }
    extra <- ww$extra    #  capture info passed extra to pass back in subsequent calls
    ww$loglk -> old.loglk -> loglk
    score <- if (constrain == "no") ww$score   else t(Z) %*% ww$score
    inf <- if (constrain == "no")  ww$inf   else t(Z) %*% ww$inf %*% Z

    if (j<= n.earlyit && length(theta0)>1) step <- 
          greenstadt.step(reverse*inf,reverse*score,max.eigenrat,iter=j,
notuphill=FALSE, messg=messg, print.progress=print.progress)
    else {   #
        soln <- qr.solve.cw(inf,score,messg=messg)
        if (soln$err != 0){
             if (print.progress>=1)
                    cat(messg,"mlefn: qr.solve.cw failed on iteration",j," -- error =2","\n")
             return(list(theta=theta0,loglk=ww$loglk,score=ww$score,
                   inf=ww$inf, constrscore=score, constrinf=inf, Z=Z, error=2))
                  }
         step <- soln$soln
         }
    phistep <- step
    if (constrain != "no") step <- Z %*% step
 
    if (guide != "no"){
        if (is.na(reverse*sum(score*phistep))) {
           if (print.progress>=1)
              cat(messg,"mlefn: NA produced in sum(score*phistep) on iteration",j," -- error =3","\n")
           return(list(theta=theta0,loglk=ww$loglk,score=ww$score,
                  inf=ww$inf, constrscore=score, constrinf=inf, Z=Z, error=3))
        }
 
#       ENSURE STEP DIRECTION UPHILL and not too long
        if ((reverse*sum(score*phistep))<0) {
           if (length(theta0)>1) step <- greenstadt.step(reverse*inf,reverse*score,
                                           max.eigenrat,iter=j,notuphill=TRUE,
                                           messg=messg, print.progress=print.progress)
           else step <- -step
        }

#        ENSURE STEP INCREASES LIKELIHOOD
         lambda <- 2
         step.reduces <- FALSE
         while (!step.reduces) {
             lambda <- lambda/2
             theta1 <- theta0 + lambda*step
             w <- loglkfn(theta1, nderivs=0,extra=extra, ...)
             loglk <- w$loglk
             success.eval <- TRUE
             if (!is.null(w$error)) if (w$error != 0)  success.eval <- FALSE
             if (is.null(loglk))  success.eval <- FALSE
             if (!is.null(loglk))  if (is.na(loglk))  success.eval <- FALSE
             if (!success.eval)
                  if (print.progress>=1) cat(messg,"mlefn: loglk evaluation failed on halving at iteration",j,"\n")
                       # carry on and hope next step halving brings us back
                               # into an area where the likelhood can be calculated
             if (success.eval){
                  if (print.progress>=2 & lambda<1) cat(messg,"iter=",counter,"  lambda=",
                         lambda,"  old=",old.loglk,"   new=",loglk,"\n")
                  delta.loglk <- loglk - old.loglk
                  if(!is.na(delta.loglk) & ((reverse*delta.loglk) > -abs(old.loglk)*tol)){
                          step.reduces <- TRUE
                          extra <- ww$extra #capture any info passed in extra
                          }
                  }
             if (!step.reduces & lambda < 0.001){
                  if (print.progress>=1)
                          cat(messg,"mlefn: Too many halvings on iteration",j," -- error =4","\n")
                  return(list(theta=theta0,loglk=ww$loglk,score=ww$score,
                         inf=ww$inf, constrscore=score, constrinf=inf, Z=Z, error=4))
                 #  return inf at the point where we could not progress from
                  }
             }
         }
     else  # no guidance, straight Newton or Fisher Scoring
         {
         lambda <- 1
         theta1 <- theta0 + step
         }
    
#   CHECK FOR CONVERGENCE

    if (print.progress>=2)
         cat(messg,"iteration", counter,"  loglk=",ww$loglk,"\n","  theta1 =",theta1,"\n")
    theta0 <- theta1
    conv <- FALSE
    if ((lambda >= 1) & (max(abs(step)) < tol*max(max(abs(theta)),1)))  conv <- TRUE
#             (Do not allow to converge while halving)

    if (conv){
       constrinf <- if (constrain=="no") NULL else inf
       constrscore <- if (constrain=="no") NULL else score
       if (print.progress>=1) cat(messg,"mlefn took",counter,"iterations\n")
       return(list(theta=as.vector(theta1),loglk=ww$loglk,score=ww$score,inf=ww$inf,
          constrscore=constrscore, constrinf=constrinf,Z=Z,counter=counter,extra=extra, error=0))
       }
} # End of iterations loop

# if (print.progress>=2) cat(messg,"WARNING mlefn: Too many iterations -- error=5","\n")
return(list(theta=theta1,loglk=ww$loglk,score=ww$score,
                   inf=ww$inf, constrscore=score, constrinf=inf, Z=Z, error=5,messg=messg))
}

