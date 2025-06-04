censre4d <- function(id,id2,C,Q,X,Z,alpha,D,S,eps=.001,eps2=1e-5,burnin=25,nsize=500,
                   silent=F,maxloops=100,iseed) {
#
# Purpose: Implements Monte Carlo EM estimator for mixed effects model with 
#          left and/or right censored data. Assumes normal errors and 
#          normal distributions for random effects.
# 
# Required:
#
#  id = unique identifier for each individual/unit
#
#  id2 = replicate number for each individual/unit
#
#  C = censoring indicator; C = -1 means left censored, C = 0 means
#      uncensored, C = 1 means right censored
#
#  Q = response variable; if C = 0, Q is the uncensored observation; if
#      C != 0, Q is the censoring level
#
#  X = design matrix for fixed effects; generate using model.matrix or
#      by hand
#
#  Z = design matrix for random effects; generate using model.matrix or  
#      by hand
#
#  alpha = initial values for fixed effects (optional)
#
#  D = initial guess for the covariance matrix of the random effects
#
#  S = initial guess for the within-subject residual variance
#
# Optional:
#
#  eps = max. relative change in parameters between successive iterates
#        to achieve convergence; typically this will be larger than one
#        might use in a deterministic algorithm
#  
#  eps2 = small number added to the absolute parameter to get around the 
#         denominator close to 0 case in determining convergence defined as
#         max(abs(theta_new - theta_old)/(abs(theta_old)+eps2)) < eps        
#
#  burnin = number of Monte Carlo samples to discard during Gibbs 
#           sampling; the sampler typically burns in quite rapidly
#           in this application
#
#  nsize = starting number of Monte Carlo samples to compute expected 
#          sufficient statistics on each person during E-step; this
#          number is doubled if the absolute change in any parameter is
#          less than the estimated Monte Carlo standard error for that
#          parameter
#
#  silent = if T, do not display intermediate results
#
#  maxloops = maximun number of EM loops
#
#  iseed = random number seed (a negative integer)
#
# Outputs:
#
#   A list with convergence indicator/criteria, parameter estimates,
#   variance estimate for fixed effects, expected residual E(y-X*alpha|Q,C)
#   and expected random effects E(b|Q,C).
#
# Reference: Hughes JP. Mixed effects models with censored data with 
#            application to HIV RNA levels. Biometrics 55:625-629, 1999.
#
  nsubj <- length(unique(id))
  n <- rle(id)$lengths
  replen <- table(id,id2)[1,1]
  p <- ncol(X)
  q <- ncol(Z)
  sumn <- sum(n)
  if (missing(alpha)) alpha <- solve(t(X)%*%X)%*%t(X)%*%Q
  info <- matrix(0,p,p)
  if (length(S)==1) { S <- diag(rep(S,replen)) }
  critd <- rep(0, p+(q*(q+1))/2+(replen*(replen+1))/2)
  storage.mode(X) <- "double"
  storage.mode(Z) <- "double"
  storage.mode(D) <- "double"
  storage.mode(S) <- "double"
  storage.mode(info) <- "double"
  Er <- rep(0, sumn)
  storage.mode(Er) <- "double"
  set.seed(iseed)
  zzz <- .Fortran("censre4d",
                  as.integer(burnin),
                  nsize=as.integer(nsize),
                  as.integer(nsubj),
                  as.integer(sumn),
                  as.integer(n),
                  as.integer(p),
                  as.integer(q),
                  as.integer(replen),
                  X,
                  Z,
                  as.double(Q),
                  as.integer(C),
                  alpha=as.double(alpha),
                  D=D,
                  S=S,
                  info=info,
                  as.integer(silent),
                  as.double(eps),
                  as.double(eps2),
                  as.integer(maxloops),
                  niter=integer(1),
                  critmx=double(1),
                  critmc=double(1),
                  critd=as.double(critd),
                  Er=Er)
  
  var.alpha <- solve(zzz$info)
  Eb  <- matrix(NA, nsubj, q)
  D1  <- zzz$D
  S1  <- zzz$S
  cnt <- 0
  for (uid in unique(id)) {
    ind <- which(id==uid)
    ni  <- length(ind)
    Zi  <- Z[ind,]
    Vi  <- Zi%*%D1%*%t(Zi) + kronecker(diag(ni/replen),S1)
    Eb[cnt+1,] <- as.vector(D1%*%t(Zi)%*%solve(Vi)%*%zzz$Er[ind])
    cnt <- cnt + 1
  }
  
  list(conv=as.numeric(abs(max(zzz$critmx))<=eps), alpha=zzz$alpha, D=zzz$D, S=zzz$S, 
       var=var.alpha, nsize=zzz$nsize, niter=zzz$niter, critmx=zzz$critmx, 
       critmc=zzz$critmc, critd=zzz$critd, Er=zzz$Er, Eb=Eb)
  
}
