
##
## R code for applying surv_leftcens_cov to example data
##

#datcodedir <- "C:\\Users\\xxie\\EinsteinMed Dropbox\\Xianhong (Sam) Xie\\ODrivePartial\\Xianhong Xie\\Misc\\Survival_Left_Censoring_Covariate\\GitHub"
datcodedir <- "\\\\tsclient\\C\\Users\\xxie\\EinsteinMed Dropbox\\Xianhong (Sam) Xie\\ODrivePartial\\Xianhong Xie\\Misc\\Survival_Left_Censoring_Covariate\\GitHub"

setwd(datcodedir)

library(nlme)
library(survival)
#dyn.load("censre4d_i386.dll")
#dyn.load("surv_leftcens_cov_v1_i386.dll")
dyn.load("censre4d_x64.dll")
dyn.load("surv_leftcens_cov_v1_x64.dll")
source("censre4d.r")
source("surv_leftcens_cov_estc_fort.r")
source("surv_leftcens_cov_infb_fort.r")
source("gen_int_dat_from_surv.r")


### Load the example data

longidat <- read.csv(file="example_dat1_longi.csv")
survdat  <- read.csv(file="example_dat1_surv.csv")

intervaldat  <- gen_int_dat_from_surv(survdat, 1/2)

intervaldat2 <- merge(intervaldat, longidat, by=c("id","repl")) # adding covariates to the interval data
intervaldat2 <- intervaldat2[order(intervaldat2$id,intervaldat2$repl),]

intervaldat2$q1h <- intervaldat2$q1
intervaldat2$q1h[intervaldat2$c1==1] <- intervaldat2$q1h[intervaldat2$c1==1] - log10(2)


### Call the functions

## Naive method

# Fit the model
coxfit1 <- coxph(Surv(tstart,tend,delta) ~ x1bl + q1h + q2, data=intervaldat2)
summary(coxfit1)


## Two-step approach

dat2     <- longidat
dat2$q1h <- dat2$q1
dat2$q1h[dat2$c1==1] <- dat2$q1[dat2$c1==1] - log10(2)
idname  <- "id"
id2name <- "repl"
dat2    <- dat2[order(dat2[,idname],dat2[,id2name]),] #make sure the data is ordered by id then id2
id      <- dat2[,idname]
id      <- as.numeric(factor(id)) # if id is too big for the integer type, transform it to sequential numbers starting from 1
id2     <- as.vector(unlist(sapply(table(id),function(i) 1:i))) # get the sequential replicate numbers starting from 1
visit2  <- dat2$repl # sequential visit numbers allowing for intermittent missing (starting from 1)
Q  <- cbind(dat2$q1, dat2$q2)
nq <- ncol(Q)
Q2 <- cbind(dat2$q1h, dat2$q2) # HDL initial response value
C  <- cbind(dat2$c1, 0)
covnames <- c("x1", "x2")
D1 <- cbind(Intercept=1, dat2[,covnames]) # 3 covs including the intercept
colnames(D1)[-1] <- covnames
np <- ncol(D1)
Z1 <- cbind(Intercept=1, slope=visit2)

idr   <- rep(id,  rep(nq,length(id)))
rrdat <- rep(id2, rep(nq,length(id)))
visit2r <- rep(visit2, rep(nq,length(id)))
vtype <- factor(rep(1:nq,length(id)))
Xrdat <- as.vector(t(as.matrix(Q2)))
D1r   <- D1[rep(1:nrow(D1),rep(nq,length(id))),]

lmedat <- data.frame(idr=idr, rrdat=rrdat, vtype=vtype, visit2=visit2r, Xrdat=Xrdat, D1r)
  
options(contrasts=c(factor="contr.treatment", ordered="contr.poly"))
if (length(covnames)==1) {
  lmeform <- paste("Xrdat", paste("vtype*",covnames,"-1",sep=""), sep=" ~ ")
} else { #length(covnames) > 1
  lmeform <- paste("Xrdat", paste("vtype*(",paste(covnames,collapse="+"),")-1",sep=""), sep=" ~ ")
}
fit1   <- lme(as.formula(lmeform), random=~vtype*visit2-1|idr, weights=varIdent(form=~1|vtype), 
              corr=corSymm(form=~1|idr/rrdat), data=lmedat, method="ML", control=lmeControl(returnObject=TRUE))

alpha1 <- as.vector(fit1$coef$fixed)
G1     <- fit1$sigma^2 * pdMatrix(fit1$modelStruct$reStruct)[[1]]
Zg     <- diag(2*nq) # random intercepts and slopes by vtype (hence 2*nq)
Zg[-(1:(nq+1)),nq+1] <- 1
G1b    <- Zg%*%G1%*%t(Zg)
R1     <- fit1$sigma^2 * diag(c(1,coef(fit1$modelStruct$varStruct,unconstr=F))) %*%
            corMatrix(fit1$modelStruct$corStruct)[[1]] %*% diag(c(1,coef(fit1$modelStruct$varStruct,unconstr=F)))

TM     <- diag(nq)
TM[1,] <- 1
alpmat <- rbind(alpha1[1:nq],
                cbind(alpha1[nq+(1:(np-1))],matrix(alpha1[-(1:(nq+np-1))],np-1,nq-1,byrow=T))%*%TM)
  
Vtyper <- t(sapply(vtype, function(x) (x==rep(1:nq,rep(np,nq)))*1))
Xr     <- D1r[,rep(1:np,nq)]*Vtyper
Zr     <- t(sapply(vtype,function(x) (x==(1:nq))*1))
Zrv    <- cbind(Zr, kronecker(visit2,diag(nq)))
  
exb  <- matrix(predict(fit1,level=1), nrow(Q), nq, byrow=T) # predicted values excluding measurement errors

outdat <- data.frame(cbind(id=id, repl=id2, exb))
names(outdat)[-(1:2)] <- c("q1ne", "q2ne")

intervaldat3 <- merge(intervaldat2, outdat, by=c("id","repl")) # adding predicted covariates to the interval data
intervaldat3 <- intervaldat3[order(intervaldat3$id,intervaldat3$repl),]

coxfit2 <- coxph(Surv(tstart,tend,delta) ~ x1bl + q1ne + q2ne, data=intervaldat3)
summary(coxfit2)


## Joint model

seed <- 100

m    <- length(unique(longidat$id))
x1bl <- longidat$x1bl[cumsum(c(1,rle(longidat$id)$lengths[-m]))]

ebinit  <- as.matrix(ranef(fit1))%*%t(Zg)
excinit <- Q2 - as.matrix(D1)%*%alpmat

bhaz2 <- basehaz(coxfit2, centered=F)[,"hazard"]
h0t2  <- diff(unique(sort(c(0,bhaz2))))

sjmfit <- surv.leftcens.cov.est.fort(survdat$time, survdat$delta, survdat$id, Q, C, longidat$id, longidat$repl, 
	         x1bl, D1, Z1, 1/2, ebinit, excinit, initVal=list(h0t=h0t2,alpha=as.vector(alpmat),beta=coef(coxfit2)[1],
	         gamma=coef(coxfit2)[-1],sigb=G1b,sige=R1), tol=3e-2)

sjminf <- surv.leftcens.cov.inf.fort(survdat$time, survdat$delta, survdat$id, Q, C, longidat$id, longidat$repl, 
  	         x1bl, D1, Z1, 1/2, objFit=sjmfit)


peRes3 <- c(sjmfit$estimate$beta, sjmfit$estimate$gamma)

apVar  <- solve(sjminf$info)
separ  <- sqrt(diag(apVar))

aptVar <- solve(sjminf$infot)
setpar <- sqrt(diag(aptVar))

seRes  <- setpar[length(sjmfit$estimate$h0t)+nq*np+(1:3)]


### format the inference outputs (Joint Method)

nhp <- length(sjmfit$estimate$h0t)
nfa <- np * nq # alpha matrix
nfs <- 3 # beta & gamma
nvr <- 4 * 5 / 2 # sigb
nve <- 2 * 3 / 2 # sige
ntotpar <- nhp + nfa + nfs + nvr + nve 
# nparms=40, order of parameters: h0(t), np*nq (column major), beta, gamma, sigb (lower triange columns), sige (lower triangle columns)


lcl95 <- rep(NA, ntotpar)
lcl95[1:nhp] <- exp(sjminf$natrParm[1:nhp] - 1.96*setpar[1:nhp]) # h0t
lcl95[nhp+(1:(nfa+nfs))] <- sjminf$natrParm[nhp+(1:(nfa+nfs))] - 1.96*setpar[nhp+(1:(nfa+nfs))] # alpha/beta/gamma

ucl95 <- rep(NA, ntotpar)
ucl95[1:nhp] <- exp(sjminf$natrParm[1:nhp] + 1.96*setpar[1:nhp]) # h0t
ucl95[nhp+(1:(nfa+nfs))] <- sjminf$natrParm[nhp+(1:(nfa+nfs))] + 1.96*setpar[nhp+(1:(nfa+nfs))] # alpha/beta/gamma

zval <- rep(NA, ntotpar)
zval[nhp+(1:(nfa+nfs))] <- sjminf$natrParm[nhp+(1:(nfa+nfs))] / setpar[nhp+(1:(nfa+nfs))]

pval <- rep(NA, ntotpar)
pval[nhp+(1:(nfa+nfs))] <- 2 * (1 - pnorm(abs(zval[nhp+(1:(nfa+nfs))])))


# inference on the first biomarker
mreg_bmk1_jm <- cbind(Estimate=sjminf$natrParm[nhp+(1:np)], 
    SE=setpar[nhp+(1:np)], 
    LCL=lcl95[nhp+(1:np)], UCL=ucl95[nhp+(1:np)],
    Z=zval[nhp+(1:np)], P=pval[nhp+(1:np)])
rownames(mreg_bmk1_jm) <- c("Intercept", covnames)
mreg_bmk1_jm


# inference on the second biomarker
mreg_bmk2_jm <- cbind(Estimate=sjminf$natrParm[nhp+np+(1:np)], 
    SE=setpar[nhp+np+(1:np)], 
    LCL=lcl95[nhp+np+(1:np)], UCL=ucl95[nhp+np+(1:np)],
    Z=zval[nhp+np+(1:np)], P=pval[nhp+np+(1:np)])
rownames(mreg_bmk2_jm) <- c("Intercept", covnames)
mreg_bmk2_jm


# inference on the variance-covariance parameters
mreg_vcov_inf_jm <- cbind(Estimate=sjminf$origParm[nhp+nfa+nfs+(1:(nvr+nve))],
                      SE=separ[nhp+nfa+nfs+(1:(nvr+nve))])
rownames(mreg_vcov_inf_jm) <- c("sigb11","sigb21","sigb31","sigb41",
                             "sigb22","sigb32","sigb42","sigb33",
							 "sigb43","sigb44","sige11","sige21",
							 "sige22") # random intercept followed by solope by biomarker
mreg_vcov_inf_jm


# inference on baseline harzard
h0t_inf3 <- cbind(Time=sort(unique(survdat$time[survdat$event==1])),
  Estimate=sjminf$origParm[1:nhp], SE=separ[1:nhp])
h0t_inf3


# inference on survival coefs
sreg_coefs3 <- cbind(Estimate=sjminf$natrParm[nhp+nfa+(1:nfs)], 
    SE=setpar[nhp+nfa+(1:nfs)], 
    Z=zval[nhp+nfa+(1:nfs)], P=pval[nhp+nfa+(1:nfs)],
	HR=exp(sjminf$natrParm[nhp+nfa+(1:nfs)]),
    LCL=exp(lcl95[nhp+nfa+(1:nfs)]), UCL=exp(ucl95[nhp+nfa+(1:nfs)]))
rownames(sreg_coefs3) <- c("x1bl", "q1ne", "q2ne")
sreg_coefs3


##########

#save.image("example_apply_results.RData")

#############################################################################
