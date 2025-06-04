
surv.leftcens.cov.est.fort <- function(T, delta, id1, Q, D, id2, repl, V, X, Z, sample_interval, eb, eyc, 
    NIter=50, initVal=NULL, tol=1e-2, delta1=1e-3, iseed=300, maxiter2=50, nburnin=1000, nsize0=500, 
    maxnsize=320000, nbatchs=25)
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
  ## initVal: list of initial values, including the following components: 
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
  ##   estimate: the parameter estimate at the final iteration - a list of the same components as initVal
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
  td <- sapply(T, function(t) sum(etimes1<=t))
  sumtd <- sum(td)
  sumtdpm <- sumtd + m
  
  maxni <- max(n)
  maxrepl <- max(repl)
  pxtq <- px*q
  pvpq <- pv + q
  nfixpar <- pxtq + pvpq
  qeu <- q*(q+1)/2
  ntotpar2 <- ne1 + nfixpar + qeu
  
  
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
  if (is.null(initVal)) {
    stop("The parameter initVal must be specified!\n")
  } else {
    h0t <- as.vector(initVal$h0t)
    alphamat  <- matrix(initVal$alpha, px, q, byrow=F)
    beta <- as.vector(initVal$beta)
    gamma <- as.vector(initVal$gamma)
    sigb <- as.matrix(initVal$sigb)
    sige <- as.matrix(initVal$sige)
  }
  nsize <- nsize0
  
  
  ## Create the work arrays
  #lwka  <- m*qb + maxni*(px+4*q+2) + maxnsize*qb + maxrepl*3 + (ne0 + 1) + ne1*(pxtq+2*qb+2*q+14) +
  #         ne1^3 + ne1p1*(px+2) + ntotpar2*4 + ntotpar2^2*5 + ntotpar2^3 + pxtq*q + pxtq^2 + sumn*q +
  #         sumtd*4 + sumtdpm*(px+2)
  lwka  <- m*(1+qb) + max(ne1,qeu)^3 + maxni*(2+px+4*q) + maxrepl*3 + (ne0 + 1) + ne1*(14+pxtq+2*q+2*qb) +
           ne1p1*(2+px) + ntotpar2*(4+5*ntotpar2) + ntotpar2^3 + pxtq^2 + q*(maxni+pxtq) + qb*maxnsize + 
		   sumn*q + sumtd*4 + sumtdpm*(2+px)
  #lwka <- floor(1.15*lwka)
  
  #liwka <- maxnsize + sumn + 3*m + maxni*q + nbatchs
  liwka <- m*3 + maxni*q + maxnsize + nbatchs + sumn
  #liwka <- floor(1.15*liwka)
  #browser()
  
  
  set.seed(iseed)
  retList <- .Fortran("survlccv",
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
                      maxiter=as.integer(NIter),
                      eps1=as.double(tol),
                      eps2=as.double(delta1),
                      maxiter2=as.integer(maxiter2),
                      nburnin=as.integer(nburnin),
                      nsize=as.integer(nsize),
                      maxnsize=as.integer(maxnsize),
                      nbatchs=as.integer(nbatchs),
                      hconv=double(1),
                      sconv=double(ne1),
                      pconv=double(1),
                      dconv=double(nfixpar+qb*(qb+1)/2+qeu),
                      mmcdel=double(1),
                      dmcdel=double(pvpq),
                      iter=integer(1),
                      wkarr=double(lwka),
                      iwkarr=integer(liwka),
                      DUP=F)
  #browser()
  
  list(hconv=retList$hconv, sconv=retList$sconv, pconv=retList$pconv, dconv=retList$dconv, mmcdel=retList$mmcdel, 
       dmcdel=retList$dmcdel, nIter=retList$iter, nsize=retList$nsize, estimate=list(h0t=retList$h0t,
         alpha=retList$alphamat,beta=retList$beta,gamma=retList$gamma,sigb=matrix(retList$sigb,qb,qb,byrow=F),
         sige=matrix(retList$sige,q,q,byrow=F)), 
       eb=matrix(retList$eb,m,qb,byrow=F), eyc=matrix(retList$eyc,sumn,q,byrow=F))
  
}
