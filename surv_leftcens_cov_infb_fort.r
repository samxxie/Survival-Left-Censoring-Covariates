
surv.leftcens.cov.inf.fort <- function(T, delta, id1, Q, D, id2, repl, V, X, Z, sample_interval,
    objFit, eb=objFit$eb, eyc=objFit$eyc, nsize=objFit$nsize, iseed=300, maxiter2=50, nburnin=1000)
{ ##
  ## Estimate the parameters in the generalized linear mixed model with covariates subject to measurement errors
  ## and left-censoring due to detection limits.
  ## 
  ## y: vector of binary responses
  ## Q: matrix of observed and left-censored biomarkers, each column corresponds to one biomarker
  ## C: matrix of censoring indicators, corresponds to Q with the same dimension, 0=observed, 1=left-censored
  ## D: matrix of covariates for the binary response and biomarkers, the response and biomarkers each correspond 
  ##    to one block of columns in D, with the span length of the columns given in p
  ## p: vector of integers, corresponds to the span lengths of columns in D for the binary response and biomarkers
  ## id: vector of ids, same length as the binary responses y
  ## exc: matrix of expected value of "centered" X (i.e., after subtracting the fixed effects on the biomarkers) 
  ## NIter: the maximal number of iterations of the MCNR algorithm
  ## objFit: return list from survival left censoring on covariate model fitting, including the following components: 
  ##          betas  - vector of parameters: beta_y, beta_b, vec(beta_x)
  ##          sigb11 - variance of the between subject random variation b_{i1} on the binary response
  ##          sigb2  - variance matrix of the between subject random variation b_{i,-1} on the biomarkers
  ##          sige   - variance matrix of the measurement errors on the biomarkers
  ## tol: convergence criterion of the MCNR algorithm
  ## delta1: a small number for preventing the denominator from being too close to 0 in evaluating the convergence
  ## iseed: random seed of the MCNR algorithm
  ## maxiter2: the maximal number of iterations used in finding the mode of b_i with the Newton-Raphson algorithm
  ## nburnin: the burn-in period in the MCMC sampling
  ## nsize0: the initial monte-carlo sample size, adaptively changed in the MCNR algorithm
  ## maxnsize: the maximal monte-carlo sample size for the MCNR algorithm
  ## nbatches: the number of batches used in calculating the monte-carlo standard errors of the parameters.
  ## 
  ## Return a list with the following components:
  ##   conv: the convergence criterion achieved
  ##   convd: the detailed convergence criteria for all parameters - a vector
  ##   mcdel: the monte-carlo standard error criterion - the minimum of the ratios of absolute parameter changes
  ##          to the monte-carlo standard errors of the parameters beta_y and beta_b
  ##   nIter: number of iterations used
  ##   nsize: the monte-carlo sample size at the final iteration
  ##   estimate: the parameter estimate at the final iteration - a list of the same components as objFit$estimate
  ##   exc: the updated expected values of X after subtracting the fixed effects on the biomarkers.
  ##
  
  ## Get the dimensions
  m   <- length(T)
  n   <- as.vector(table(id2))
  sumn <- sum(n)
  px  <- ncol(X)
  q   <- ncol(Q)
  qb  <- 2*q
  pv  <- ifelse(is.null(dim(V)), 1, ncol(V))
  dim(V) <- c(m, pv)
  etimes0 <- sort(T[delta==1])
  etimes1 <- unique(etimes0)
  ne0 <- length(etimes0)
  ne1 <- length(etimes1)
  ne1p1 <- ne1 + 1
  ne12u <- ne1*(ne1+1)/2
  td <- sapply(T, function(t) sum(etimes1<=t))
  sumtd <- sum(td)
  sumtdpm <- sumtd + m
  
  maxni <- max(n)
  maxrepl <- max(repl)
  pxtq <- px*q
  pvpq <- pv + q
  nfixpar <- pxtq + pvpq
  qbu <- qb*(qb+1)/2
  qeu <- q*(q+1)/2
  ntotpar <- ne1 + nfixpar + qbu + qeu
  
  
  ## Sort the subject level data and visit level data
  order.subj <- order(id1)
  id1 <- id1[order.subj]
  T <- T[order.subj] 
  delta <- delta[order.subj]
  V <- V[order.subj,]
  eb <- eb[order.subj,]
  
  order.svis <- order(id2, repl)
  id2 <- id2[order.svis]
  repl <- repl[order.svis]
  Q <- Q[order.svis,]
  D <- D[order.svis,]
  X <- X[order.svis,]
  Z <- Z[order.svis,]
  eyc <- eyc[order.svis,]
  
  
  ## Make sure the data are in correct modes
  T <- as.vector(T)
  delta <- as.vector(delta)
  id1 <- as.vector(id1)
  Q <- as.matrix(Q)
  D <- as.matrix(D)
  id2 <- as.vector(id2)
  repl <- as.vector(repl) 
  V <- as.matrix(V)
  X <- as.matrix(X) 
  Z <- as.matrix(Z)
  
  eb0 <- eb
  eyc0 <- eyc
  eb <- as.matrix(eb)
  eyc  <- as.matrix(eyc)
  
  
  ## Assign the initial values
  if (is.null(objFit)) {
    stop("The parameter objFit must be specified!\n")
  } else {
    h0t <- as.vector(objFit$estimate$h0t)
    alphamat  <- matrix(objFit$estimate$alpha, px, q, byrow=F)
    beta <- as.vector(objFit$estimate$beta)
    gamma <- as.vector(objFit$estimate$gamma)
    sigb <- as.matrix(objFit$estimate$sigb)
    sige <- as.matrix(objFit$estimate$sige)
  }
  
  
  ## Create the work arrays
  lwka  <- m*(1+qb) + maxni*(2+px+q*3) + maxrepl*3 + (ne0+1) + ne1*9 + ne12u + ne1p1*(2+px+q) +
           ntotpar*(1+ntotpar*5) + ntotpar^3 + q*(maxni*2+ne1p1) + q*maxni*ne1 + q*qb^2*(maxni+ne1) +
		   q^2*(maxni*3+maxni^2+ne1*2+ne12u) + q^2*maxni*ne1*2 + q^3*(maxni+maxni^2) + q^3*maxni*ne1 + 
		   q^4*maxni^2 + qb*(ne1+ne12u) + qb^2*ne1 + qb^2*q^2*maxni + (qb+q*maxni)*nsize + sumn*q +
		   sumtd*4 + sumtdpm*(2+px)
  
  liwka <- m*3 + maxni*q + sumn
  #browser()
  
  
  set.seed(iseed)
  #browser()
  retList <- .Fortran("survlinf",
                      time=as.double(T),
					  event=as.integer(delta),
					  id1=as.integer(id1),
					  m=as.integer(m),
					  Q=as.double(Q),
					  ldq=as.integer(sumn),
					  sumn=as.integer(sumn),
					  qe=as.integer(q),
					  D=as.integer(D),
					  ldd=as.integer(sumn),
					  id2=as.integer(id2),
					  repl=as.integer(repl),
                      V=as.double(V),
					  ldv=as.integer(m),
					  pv=as.integer(pv),
					  X=as.double(X),
					  ldx=as.integer(sumn),
					  px=as.integer(px),
					  Z=as.double(Z),
					  ldz=as.integer(sumn),
					  samp_pd=as.double(sample_interval),
					  eb=as.double(eb),
					  ldeb=as.integer(m),
					  eyc=as.double(eyc),
					  ldeyc=as.integer(sumn),
					  h0t=as.double(h0t),
                      alphamat=as.double(alphamat),
					  ldam=as.integer(px),
					  beta=as.double(beta),
					  gamma=as.double(gamma),
					  sigb=as.double(sigb),
					  ldsb=as.integer(qb),
					  sige=as.double(sige),
					  ldse=as.integer(q),
                      maxiter2=as.integer(maxiter2),
					  nburnin=as.integer(nburnin),
					  nsize=as.integer(nsize),
					  origParm=double(ntotpar),
					  origGrad=double(ntotpar),
					  origInfo=double(ntotpar^2),
					  ldoi=as.integer(ntotpar),
                      natrParm=double(ntotpar),
					  natrGrad=double(ntotpar),
					  natrInfo=double(ntotpar^2),
					  ldni=as.integer(ntotpar),
					  wkarr=double(lwka),
					  iwkarr=integer(liwka),
                      DUP=F)
  
  list(origParm=retList$origParm, gradient=retList$origGrad, info=matrix(retList$origInfo,ntotpar,ntotpar),
       natrParm=retList$natrParm, gradientt=retList$natrGrad, infot=matrix(retList$natrInfo,ntotpar,ntotpar))
  
}
