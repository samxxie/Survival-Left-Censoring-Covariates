
gen_int_dat_from_surv <- function(survdat, int_len) {
  
  ##
  ## Generate interval style data (multiple rows per id) from survival data
  ##   (one row per id)
  ##
  
  outdat <- matrix(NA, nrow=sum(ceiling(survdat$time/int_len)), ncol=5)
  colnames(outdat) <- c("id", "repl", "tstart", "tend", "delta")
  intervaldat <- data.frame(outdat)
  
  cnt <- 0
  for (id in sort(unique(survdat$id))) {
    
	nint <- ceiling(survdat$time[survdat$id==id]/int_len)
	
	if (nint > 1) {
	  intervaldat[cnt+(1:(nint-1)),] <- cbind(id=rep(id,nint-1), repl=1:(nint-1),
	    tstart=int_len*(0:(nint-2)), tend=int_len*(1:(nint-1)), delta=rep(0,nint-1))
    }
    
    intervaldat[cnt+nint,] <- c(id=id, repl=nint, tstart=int_len*(nint-1),
	  tend=survdat$time[survdat$id==id], delta=survdat$delta[survdat$id==id])
	
	cnt <- cnt + nint
  }
  
  intervaldat
  
}
