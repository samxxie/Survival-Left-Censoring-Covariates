c
c     Generalize the code for more than two biomarkers and optimize the code
c
      
      subroutine censre4d(burnin,nsize,nsubj,sumn,n,p,q,rs,X,Z,
     *                    YQ,YC,alpha,D,S,info,silent,eps,eps2,maxloops,
     *                    niter,critmx,critmc,critd,Er)
c     
c Routine to implement EM algorithm to estimate parameters of mixed
c effects linear model in the presence of right and/or left censored
c data.
c
c Arguments:
c
c burnin (I) = number of Monte Carlo samples to discard during Gibbs
c           sampling; the sampler typically burns in quite rapidly
c           in this application
c nsize (I) = starting number of Monte Carlo samples to compute expected
c          sufficient statistics on each person during E-step; this
c          number is doubled if the absolute change in any parameter
c          is less than the estimated Monte Carlo standard error of that 
c          parameter
c nsubj (I) = number of subjects/clusters in the data set
c sumn (I) = total number of observations
c n (I) = vector of length nsubj giving the number of observations
c      per subject/cluster
c p (I) = number of fixed effects (columns of design matrix)
c q (I) = number of random effects (columns of design matrix)
c rs(I) = number of measurement types per replicate
c X (I) = design matrix for fixed effects
c Z (I) = design matrix for random effects
c YQ (I) = response (see censre2.s)
c YC (I) = censoring indicator (see censre2.s)
c alpha (IO) = on input, initial values for fixed effects; on output,
c           the final fitted values
c D (IO) = on input, initial values for covariance matrix of random effects;
c       on output, the final fitted values
c S (IO) = on input, initial value for within-replicate residual variance;
c       on output, the final fitted value
c info (O) = information matrix for the fixed effects
c silent (I) = if 0, print out intermediate results; otherwise, don't print
c           intermediate results
c eps (I) = max. relative change in parameters between successive iterates
c        to achieve convergence; typically this will be larger than one
c        might use in a deterministic algorithm
c maxloops (I) = maximum number of EM loops
c
c Author:
c   Jim Hughes (jphughes@u.washington.edu)
c   Dept. of Biostat.
c   Univ. of WA
c
c Revision history
c 4/28/99 - Fixed computation of the likelihood; incorrect computation
c           caused the gibbs sample size to increase beyond necessity
c           in some situations.
c
c 1/30/07 - revised censre2.f to censre3.f to call R functions pnorm, qnorm and 
c           random # routines using c wrappers
c
c 2/6/07 - changed the name of sumX orutine to mysumX to avoid conflicts in Splus
c
      Parameter(MAXP=100,MAXQ=20,MAXR=10,MAXN=200)
      Parameter(MAXPQ=MAXP+(MAXQ*(MAXQ+1))/2+(MAXR*(MAXR+1))/2)
      implicit double precision(A-H,O-Z)
      integer burnin,nsize,nsubj,n(*),p,q,rs,YC(*),sumn,ind,silent
      integer npar,loop,maxloops,nrepl,sumn2,niter,i,ierr
      double precision X(sumn,*),Z(sumn,*),YQ(*),critmx,critmc,critd(*),
     $     Er(*)
      double precision alpha(*),D(q,*),S(rs,*),info(p,*),eps,eps2
      double precision Eee(MAXR,MAXR),Ebb(MAXQ,MAXQ),XWY(MAXP),ui(MAXN)
      double precision Vcc(MAXP,MAXP)
      double precision Wi(MAXN,MAXN),DZW(MAXQ,MAXN),DZ(MAXQ,MAXN),
     $     XWX(MAXP,MAXP)
      double precision XW(MAXP,MAXN),r(MAXN),rr(MAXN,MAXN),
     $     Vmc(MAXP,MAXP)
      double precision kR(MAXN,MAXN),RW(MAXN,MAXN),RWR(MAXN,MAXN),
     $     RWRR(MAXN,MAXN)
      double precision y(MAXN),vy(MAXN,MAXN),det(2),va(MAXP),
     $     trRW(MAXR,MAXR)
      double precision newa(MAXP),newD(MAXQ,MAXQ),newS(MAXR,MAXR),mcdel,
     $     delmx
      double precision delta(MAXPQ),rlogl,llogl,ln10,ldetV,ldetD,junk(1)
      data ln10/2.302585/
      
c      call random(1,junk,iseed)
c      call dblepr("junk",4,junk,1)
c      call intpr("burnin",-1,burnin,1)
c	  call intpr("nsize",-1,nsize,1)
c      call intpr("nsubj",-1,nsubj,1)
c	  call intpr("sumn",-1,sumn,1)
c	  call intpr("n",-1,n,nsubj)
c	  call intpr("p",-1,p,1)
c	  call intpr("q",-1,q,1)
c	  call intpr("rs",-1,rs,1)
c	  call dblepr("X(,1)",-1,X(1,1),sumn)
c	  call dblepr("X(,30)",-1,X(1,30),sumn)
c	  call dblepr("Z(,1)",-1,Z(1,1),sumn)
c	  call dblepr("Z(,3)",-1,Z(1,3),sumn)
c	  call dblepr("YQ",-1,YQ,sumn)
c	  call intpr("YC",-1,YC,sumn)
c	  call dblepr("alpha",-1,alpha,30)
c	  call dblepr("D",-1,D,9)
c	  call dblepr("S",-1,S,9)
c	  call dblepr("info(1,1)",-1,info,30)
c	  call dblepr("info(1,30)",-1,info,30)
c	  call intpr("silent",-1,silent,1)
c	  call dblepr("eps",-1,eps,1)
c	  call dblepr("eps2",-1,eps2,1)
c	  call intpr("maxloops",-1,maxloops,1)
c	  call intpr("niter",-1,niter,1)
c	  call dblepr("critmx",-1,critmx,1)
c	  call dblepr("critmc",-1,critmc,1)
c	  call dblepr("critd",-1,critd,42)
c	  call dblepr("Er",-1,Er,sumn)
c      return

      npar = p + (q*(q+1))/2 + (rs*(rs+1))/2
      sumn2 = sumn / rs
      loop = 0
 10   loop = loop+1
      ind = 0
c     
c Following loop does the matrix operations which constitute the
c E step of the EM algorithm
c     
      call veczero(p,XWY)
      call matzero(MAXQ,q,q,Ebb)
      call matzero(MAXR,rs,rs,Eee)
      call matzero(MAXP,p,p,XWX)
      call matzero(MAXP,p,p,Vcc)
      call matzero(MAXP,p,p,Vmc)
c      call matzero(MAXN,n(i),n(i),kR)
      call matzero(MAXN,MAXN,MAXN,kR)
      
c
      do 20 i = 1,nsubj
        nrepl = n(i) / rs
c        call intpr("i",1,i,1)
c		 return
        call makeui(sumn,p,ind,X,n(i),alpha,ui)
c      call dblepr("ui",2,ui,n(i))
        call kron(MAXN,S,rs,nrepl,kR)
        call makeVi(MAXN,sumn,q,ind,Z,D,kR,n(i),Wi)
c      call dblepr("Wi",2,Wi,n(i))
        call gibbs(MAXN,burnin,nsize,ui,Wi,n(i),YQ(ind+1),YC(ind+1),
     *             r,rr,y,vy)
c       call dblepr("r",1,r,n(i))
c       call dblepr("rr",2,rr,n(i))
c       call dblepr("y",1,y,n(i))
c r = E(Y - X%*%alpha | YQ,YC,alpha)
c rr = E((Y-X%*%alpha)(Y-X%*%alpha)^T | YQ,YC,alpha)
c y = E(Y | YQ,YC,alpha)
        call saver(n(i),r,ind,Er)
        call invert(MAXN,n(i),Wi,det,ierr)
c      call intpr("ierr",4,ierr,1)
c        ldetV = det(2)*ln10 + log(det(1))
        call makeDZW(MAXN,MAXQ,sumn,q,n(i),ind,D,Z,Wi,DZW,DZ)
c        call intpr('MAXQ',-1,MAXQ,1)
c        call intpr('q',-1,q,1)
c        call intpr('ni',-1,n(i),1)
c        do 15 k=1,n(i)
c 15     call dblepr('DZW',-1,DZW(1,k),q)
c        call dblepr('r',-1,r,n(i))
c        call predRe(MAXQ,q,n(i),DZW,r,Eba(q*ind/rs+1))
c        call dblepr('Eb',-1,Eba(q*ind/rs+1),q)
c        if (i.eq.100) return
        call makeXW(MAXN,MAXP,sumn,p,n(i),ind,X,Wi,XW)
        call mysumX(MAXN,MAXP,sumn,p,n(i),ind,X,XW,y,vy,r,rr,XWY,XWX,
     *            Vcc,Vmc)
        call makeRW(MAXN,kR,Wi,n(i),rr,RW,RWR,RWRR)
        call blktrace(MAXN,MAXR,RWRR,n(i),rs,trRW)
        call sumEbb(MAXN,MAXQ,q,n(i),DZW,rr,D,DZ,Ebb)
        call sumEee(MAXR,S,trRW,rs,nrepl,Eee)
c        call dblepr("Ebb",3,Ebb,q)
c        call dblepr("Eee",3,Eee,1)
c        return
        ind = ind + n(i)
 20   continue
      
c      return
c
      call makeinfo(MAXP,p,XWX,Vcc,info)
      call invert(MAXP,p,XWX,det,ierr)
c
c M-step
c
      call mstep(MAXP,MAXQ,MAXR,nsubj,p,q,rs,sumn2,XWX,XWY,Ebb,Eee,
     *           newa,newD,newS)
c      newS = 1.01436
c
c compute monte Carlo uncertainty in estimate of alpha
c
      call varmc(MAXP,p,XWX,Vmc,va)
c
c check for convergence
c
      call pctchg(MAXQ,MAXR,p,q,rs,alpha,D,S,newa,newD,newS,va,delta,
     *            mcdel,eps2)
      if (silent.eq.0) then
        call DBLEPR("newa",-1,newa,p)
        do 30 i = 1,q
 30     call DBLEPR("newD",-1,newD(1,i),q)
        do 40 i = 1,rs
 40     call DBLEPR("newS",-1,newS(1,i),rs)
        call DBLEPR("Relative change in parameters",-1,delta,npar)
      endif
c     
      delmx = delta(1)
      do 50 i = 1,npar
 50   if (delta(i).gt.delmx) delmx=delta(i)
      if (delmx.le.eps) goto 60
      if (mcdel.lt.1.0) then
        nsize=nsize*2
        call intpr("Monte Carlo sample size increased to ",-1,nsize,1)
      else if (mcdel.ge.20) then
        nsize=nsize/2
        call intpr("Monte Carlo sample size decreased to ",-1,nsize,1)
      endif
c     
      call copy(MAXQ,MAXR,p,q,rs,newa,newD,newS,alpha,D,S)
      if (loop.gt.maxloops) goto 60
      goto 10
 60   continue
      
      niter = loop
      critmx = delmx
      critmc = mcdel
      do 70 i = 1,npar
        critd(i)=delta(i)
 70   continue
      
      return
      end
      
      
      subroutine makeui(lda,p,ind,X,n,alpha,ui)
      implicit double precision(A-H,O-Z)
      integer lda,p,ind,n,i,j
      double precision X(lda,*),alpha(*),ui(*)
      
      do 5 i=1,n
        ui(i) = 0.0d0
 5    continue
      
      do 10 j = 1,p
        do 10 i = 1,n
          ui(i) = ui(i) + X(ind+i,j)*alpha(j)
 10   continue
      return
      end
      
      subroutine kron(ldkr,S,rs,nrepl,kR)
      implicit double precision(A-H,O-Z)
      integer ldkr,rs,nrepl,i,j,l
      double precision S(rs,*),kR(ldkr,*)
c
      do 10 i = 1,nrepl
        do 10 l = 1,rs
          do 10 j = 1,rs
            kR((i-1)*rs+j,(i-1)*rs+l) = S(j,l)
 10   continue
      return
      end
      
      subroutine makeVi(ldn,lda,q,ind,Z,D,kR,n,Vi)
      implicit double precision(A-H,O-Z)
      integer ldn,lda,q,ind,n,i,j,k,l
      double precision Z(lda,*),D(q,*),kR(ldn,*),Vi(ldn,*)
c     
      do 10 l = 1,n
        do 10 i = 1,n
          Vi(i,l) = kR(i,l)
          do 10 j = 1,q
            do 10 k = 1,q
              Vi(i,l) = Vi(i,l) + Z(ind+i,k)*D(k,j)*Z(ind+l,j)
 10   continue
      return
      end
      
      subroutine saver(ni,r,ind,Er)
      implicit double precision(A-H,O-Z)
      integer ni,ind,j
      double precision r(*),Er(*)
      
      do 10 j=1,ni
        Er(ind+j) = r(j)
 10   continue
      return
      end
      
      subroutine makeDZW(ldn,ldq,lda,q,n,ind,D,Z,Wi,DZW,DZ)
      implicit double precision(A-H,O-Z)
      integer ldn,ldq,lda,q,ind,n,i,j,k,l
      double precision Wi(ldn,*),DZW(ldq,*),DZ(ldq,*)
      double precision Z(lda,*),D(q,*)
      
      do 10 j = 1,n
      do 10 i = 1,q
        DZW(i,j) = 0.0D0
        DZ(i,j) = 0.0D0
        do 10 k = 1,q
          DZ(i,j) = DZ(i,j) + D(i,k)*Z(ind+j,k)
          do 10 l = 1,n
            DZW(i,j) = DZW(i,j) + D(i,k)*Z(ind+l,k)*Wi(l,j)
 10   continue
      return
      end
      
      subroutine predRe(ldq,q,n,DZW,r,Eb)
      implicit double precision(A-H,O-Z)
      integer ldq,q,n,j,k
      double precision DZW(ldq,*),r(*),Eb(*)
      
      do 10 k=1,q
        Eb(k) = 0.0d0
        do 10 j=1,n
          Eb(k)=Eb(k)+DZW(k,j)*r(j)
 10   continue
      return
      end
      
      subroutine makeXW(ldn,ldp,lda,p,n,ind,X,Wi,XW)
      implicit double precision(A-H,O-Z)
      integer ldn,ldp,lda,p,n,ind,i,j,k
      double precision X(lda,*),Wi(ldn,*),XW(ldp,*)
      
      do 10 j = 1,n
      do 10 i = 1,p
        XW(i,j) = 0.0D0
        do 10 k = 1,n
          XW(i,j) = XW(i,j) + X(ind+k,i)*Wi(k,j)
 10   continue
      return
      end
      
      subroutine mysumX(ldn,ldp,lda,p,n,ind,X,XW,y,vy,r,rr,XWY,XWX,
     *                  Vcc,Vmc)
      implicit double precision(A-H,O-Z)
      integer ldn,ldp,lda,p,n,ind,i,j,k,l
      double precision X(lda,*),XWX(ldp,*),XW(ldp,*),XWY(*)
      double precision Vcc(ldp,*),y(*),r(*),rr(ldn,*),vy(ldn,*),
     $     Vmc(ldp,*)
      
      do 10 j = 1,n
      do 10 i = 1,p
        XWY(i) = XWY(i) + XW(i,j)*y(j)
 10   continue
      do 20 j = 1,p
      do 20 i = 1,p
        do 20 k = 1,n
          XWX(i,j) = XWX(i,j) + XW(i,k)*X(ind+k,j)
          do 20 l = 1,n
            Vcc(i,j) = Vcc(i,j) + XW(i,k)*(rr(k,l)-r(k)*r(l))*XW(j,l)
            Vmc(i,j) = Vmc(i,j) + XW(i,k)*vy(k,l)*XW(j,l)
 20   continue
      return
      end
      
      subroutine makeRW(ldkr,kR,Wi,n,rr,RW,RWR,RWRR)
      implicit double precision(A-H,O-Z)
      integer ldkr,n,i,j,k,l
      double precision kR(ldkr,*),Wi(ldkr,*),rr(ldkr,*)
      double precision RW(ldkr,*),RWR(ldkr,*),RWRR(ldkr,*)
      
      do 10 j = 1,n
        do 10 i = 1,n
          RW(i,j) = 0.0D0
          do 10 k = 1,n
            RW(i,j) = RW(i,j) + kR(i,k)*Wi(k,j)
 10   continue      
      
      do 20 j = 1,n
        do 20 i = 1,n
          RWR(i,j) = 0.0D0
          do 20 k = 1,n
            RWR(i,j) = RWR(i,j) + RW(i,k)*kR(k,j)
 20   continue      
      
      do 40 j = 1,n
        do 40 i = 1,n
          RWRR(i,j) = 0.0D0
          do 30 l = 1,n
            do 30 k = 1,n
              RWRR(i,j) = RWRR(i,j) + RW(i,k)*rr(k,l)*RW(j,l)
 30       continue
          RWRR(i,j) = RWRR(i,j) - RWR(i,j)
 40   continue
      return
      end
      
      subroutine blktrace(ldn,ldtr,RWRR,n,rs,trRW)
      implicit double precision(A-H,O-Z)
      integer ldn,ldtr,n,rs,nrepl,i,j,k
      double precision RWRR(ldn,*),trRW(ldtr,*)
      
      nrepl = n / rs
      do 20 j = 1,rs
        do 20 i = 1,rs
          trRW(i,j) = 0.0D0
          do 20 k = 1,nrepl
            trRW(i,j) = trRW(i,j) + RWRR((k-1)*rs+i,(k-1)*rs+j)
 20   continue
      return
      end
      
      subroutine sumEbb(ldn,ldq,q,n,DZW,rr,D,DZ,Ebb)
      implicit double precision(A-H,O-Z)
      integer ldn,ldq,q,n,i,j,k,l
      double precision DZW(ldq,*),rr(ldn,*),D(q,*),DZ(ldq,*),Ebb(ldq,*)
      
      do 20 j = 1,q
      do 20 i = 1,q
        Ebb(i,j) = Ebb(i,j) + D(i,j)
        do 20 k = 1,n
        do 10 l = 1,n
          Ebb(i,j) = Ebb(i,j) + DZW(i,k)*rr(k,l)*DZW(j,l)
 10     continue
        Ebb(i,j) = Ebb(i,j) - DZW(i,k)*DZ(j,k)
 20   continue
      return
      end
      
      subroutine sumEee(ldtr,S,trRW,rs,nrepl,Eee)
      implicit double precision(A-H,O-Z)
      integer ldtr,rs,nrepl,i,j
      double precision S(rs,*),trRW(ldtr,*),Eee(ldtr,*)
      
      do 10 j = 1,rs
        do 10 i = 1,rs
          Eee(i,j) = Eee(i,j) + nrepl*S(i,j) + trRW(i,j)
 10   continue
      return
      end
      
      subroutine makeinfo(ldp,p,XWX,Vcc,info)
      implicit double precision(A-H,O-Z)
      integer ldp,p,i,j
      double precision XWX(ldp,*),Vcc(ldp,*),info(p,*)
      
      do 10 j = 1,p
      do 10 i = 1,p
 10   info(i,j) = XWX(i,j)-Vcc(i,j)
      return
      end
      
      subroutine gibbs(ldn,burnin,nsize,mean,W,m,YQ,YC,r,rr,Z,vy)
      Parameter(MAXP=100,MAXQ=20,MAXN=200)
      implicit double precision(A-H,O-Z)
      integer ldn,burnin,nsize,m,i,j,l,ierr
      integer i1,i2,one,YC(*)
      double precision mean(*),W(ldn,*),YQ(*),r(*),rr(ldn,*),vy(ldn,*)
      double precision mwt(MAXN,MAXN),cvar(MAXN),u(MAXN),p(MAXN),x(MAXN)
      double precision Win(MAXN,MAXN),Z(MAXN),d(2),sdev1
      double precision zero,cpnorm,cqnorm
      data one/1/,zero/0.0D0/
C     
      call random(m,u)
c      call dblepr("u",1,u,m)
      do 30 i = 1,m
        Z(i) = YQ(i)-mean(i)
        if (YC(i).eq.0) then
          x(i) = Z(i)
        else
          sdev1 = dsqrt(W(i,i))
          p(i) = cpnorm(Z(i),zero,sdev1)
c         call dblepr("p",1,p(i),1)
          if (YC(i).lt.0) then
            u(i) = p(i)*u(i)
            x(i) = cqnorm(u(i),zero,sdev1)
          else
            u(i) = p(i)+(1.0-p(i))*u(i)
            x(i) = cqnorm(u(i),zero,sdev1)
          endif
          
          call dropi(MAXN,m,W,i,Win)
          
          if (m.gt.1) call invert(MAXN,m-1,Win,d,ierr)
          cvar(i) = W(i,i)
          i1 = 0
          do 20 j = 1,m
            mwt(i,j) = 0.0D0
            if (j.ne.i) then
              i1 = i1+1
              i2 = 0
              do 10 l = 1,m
                if (l.ne.i) then
                  i2 = i2+1
                  mwt(i,j) = mwt(i,j) + W(i,l)*Win(i2,i1)
                endif
 10           continue
              cvar(i) = cvar(i) - mwt(i,j)*W(j,i)
            endif
 20       continue
          cvar(i) = dsqrt(cvar(i))
        endif
 30   continue
c     
c       call DBLEPR("mwt",3,mwt,m)
c       call DBLEPR("cvar",4,cvar,m)
c       call DBLEPR("Z",1,Z,m)
c       call dblepr("p",1,p,m)
c       call DBLEPR("x",1,x,m)
      
      call gsample(ldn,burnin,m,mwt,cvar,Z,YC,x,r,rr,vy)
      call veczero(m,r)
      call matzero(MAXN,m,m,rr)
      call gsample(ldn,nsize,m,mwt,cvar,Z,YC,x,r,rr,vy)
      
      do 40 i = 1,m
        Z(i) = r(i)+mean(i)
 40   continue
      return
      end
      
      subroutine dropi(ldn,n,W,i,V)
      implicit double precision(A-H,O-Z)
      integer ldn,n,i,j,l,i1,i2
      double precision W(ldn,*),V(ldn,*)
      i1 = 0
      do 20 j=1,n
        if (j.eq.i) goto 20
        i1 = i1+1
        i2 = 0
        do 10 l=1,n
          if (l.eq.i) goto 10
          i2 = i2+1
          V(i1,i2) = W(j,l)
 10     continue
 20   continue
      return
      end
      
      subroutine invert(lda,n,a,d,info)
      implicit double precision(A-H,O-Z)
      PARAMETER(MAXN=200)
      integer lda,n,ipvt(MAXN)
      integer job,info
      double precision a(lda,*),d(2),work(MAXN)
      data job/11/
c     
      call DGEFA(a,lda,n,ipvt,info)
      if (info.gt.0) call intpr("info = ",7,info,1)
      call DGEDI(a,lda,n,ipvt,d,work,job)
      return
      end
      
      subroutine mstep(ldp,ldq,lds,nsubj,p,q,rs,sumn,XWX,XWY,Ebb,Eee,
     *                 newa,newD,newS)
      implicit double precision(A-H,O-Z)
      integer ldp,ldq,lds,nsubj,p,q,rs,sumn,i,j
      double precision XWX(ldp,*),XWY(*),Ebb(ldq,*),Eee(lds,*)
      double precision newa(*),newD(ldq,*),newS(lds,*)
c     
c update u
      do 10 i = 1,p
         newa(i) = 0.0
         do 10 j = 1,p
           newa(i) = newa(i) + XWX(i,j)*XWY(j)
 10   continue
c update D and S
      do 20 i = 1,q
      do 20 j = i,q
         newD(i,j) = Ebb(i,j)/nsubj
         newD(j,i) = newD(i,j)
 20   continue
      do 30 i=1,rs
      do 30 j=i,rs
        newS(i,j) = Eee(i,j)/sumn
        newS(j,i) = newS(i,j)
 30   continue
      return
      end
      
      subroutine varmc(ldp,p,XWX,Vmc,a)
      implicit double precision(A-H,O-Z)
      integer ldp,p,i,j,k
      double precision XWX(ldp,*),Vmc(ldp,*),a(*)
c     
      do 10 i = 1,p
        a(i) = 0.0D0
        do 10 j = 1,p
        do 10 k = 1,p
          a(i) = a(i) + XWX(i,j)*Vmc(j,k)*XWX(k,i) 
 10   continue
      return
      end
      
      subroutine copy(ldq,lds,p,q,rs,newa,newD,newS,a,D,S)
      implicit double precision(A-H,O-Z)
      integer ldq,lds,p,q,rs,i,j
      double precision a(*),D(q,*),S(rs,*)
      double precision newa(*),newD(ldq,*),newS(lds,*)
c     
      do 10 i = 1,p
 10   a(i) = newa(i)
      do 20 j = 1,q
      do 20 i = 1,q
        D(i,j) = newD(i,j)
 20   continue
      do 30 j = 1,rs
      do 30 i = 1,rs
        S(i,j) = newS(i,j)
 30   continue
      return
      end
      
      subroutine matzero(lda,p,q,A)
      implicit double precision(A-H,O-Z)
      integer lda,p,q,i,j
      double precision A(lda,*)
      do 10 j = 1,q
      do 10 i = 1,p
 10   A(i,j) = 0.0D0
      return
      end
      
      subroutine veczero(n,A)
      implicit double precision(A-H,O-Z)
      integer n,i
      double precision A(*)
      do 10 i = 1,n
 10   A(i) = 0.0D0
      return
      end
      
      subroutine matcopy(lda,a,b,from,to)
      integer lda,a,b,i,j
      double precision from(lda,*),to(lda,*)
c
      do 10 j = 1,b
      do 10 i = 1,a
 10   to(i,j) = from(i,j)
      return
      end
      
      subroutine pctchg(lpq,lds,p,q,rs,a,D,S,newa,newD,newS,va,delta,
     *                  mcdel,eps2)
      implicit double precision(A-H,O-Z)
      integer lpq,lds,p,q,rs,ijk,i,j
      double precision a(*),D(q,*),S(rs,*),delta(*),va(*),mcdel,temp
      double precision newa(*),newD(lpq,*),newS(lds,*)
      double precision ulim,eps2
      data ulim/40.0D0/
c     
      ijk = 0
      mcdel = ulim
      do 10 i = 1,p
      ijk = ijk+1
      if (va(i).gt.0.0D0) then
        temp = abs(a(i)-newa(i))/sqrt(va(i))
      else
        temp = ulim
      endif
      if (temp.lt.mcdel) mcdel=temp
c      if (abs(a(i)).lt.eps2) then 
c        delta(ijk) = abs(a(i)-newa(i))
c      else 
c        delta(ijk) = abs((a(i)-newa(i))/a(i))
c      endif
      delta(ijk)=abs(newa(i)-a(i))/(abs(a(i))+eps2)
 10   continue
      do 20 j = 1,q
      do 20 i = 1,j
      ijk=ijk+1
c      if (abs(D(i,j)).lt.eps2) then 
c        delta(ijk) = abs(D(i,j)-newD(i,j))
c      else 
c        delta(ijk) = abs((D(i,j)-newD(i,j))/D(i,j))
c      endif
      delta(ijk)=abs(newD(i,j)-D(i,j))/(abs(D(i,j))+eps2)
 20   continue
      do 30 j = 1,rs
      do 30 i = 1,j
      ijk=ijk+1
c      if (abs(S(i,j)).lt.eps2) then 
c        delta(ijk) = abs(S(i,j)-newS(i,j))
c      else 
c        delta(ijk) = abs((S(i,j)-newS(i,j))/S(i,j))
c      endif
      delta(ijk)=abs(newS(i,j)-S(i,j))/(abs(S(i,j))+eps2)
 30   continue
      return
      end
      
      subroutine gsample(ldn,n,m,mwt,sdev,q,cens,x,r,rr,vy)
      implicit double precision(A-H,O-Z)
      PARAMETER(MAXN=200)
      integer ldn,n,m,cens(*),i,j,l,i1,ione
      double precision u(MAXN),cpnorm,cqnorm,p,cmu,rmin,rmax,one
      double precision mwt(ldn,*),sdev(*),q(*),r(*),x(*),rr(ldn,*),
     $     vy(ldn,*),eps
      data eps/1D-8/,one/1.0D0/,ione/1/
c     
      do 5 i = 1,m
        if (cens(i).eq.0) then
          x(i) = q(i)
        endif
 5    continue       
      
      do 40 i1 = 1,n
        call random(m,u)
        do 20 i = 1,m
          if (cens(i).ne.0) then
            cmu = 0.0D0
            do 10 l = 1,m
c Note that mwt(i,i) = 0
 10         cmu = cmu + mwt(i,l)*x(l)
            p = cpnorm(q(i),cmu,sdev(i))
            rmin=0.0D0
            rmax=p
            if (cens(i).gt.0) then 
              rmin=p
              rmax=1.0D0
            endif
            u(i) = rmin + u(i)*(rmax-rmin)
c          if (u(i).lt.eps.or.u(i).gt.one-eps) then
c            x(i) = q(i)
c          else
            x(i) = cqnorm(u(i),cmu,sdev(i))
c          endif
          endif
 20     continue
c        call dblepr("rmin",4,rmin,1)
c        call dblepr("rmax",4,rmax,1)
c        call dblepr("u",1,u,m)
c        call dblepr("cmu",3,cmu,1)
c        call dblepr("sdev",3,sdev,m)
c        call dblepr("x",1,x,m)
        do 30 j = 1,m
          r(j) = r(j) + x(j)
          do 30 i = 1,j
            rr(i,j) = rr(i,j) + x(i)*x(j)
 30     continue
 40   continue
      do 50 j = 1,m
        r(j) = r(j)/n
        do 50 i = 1,j
          rr(i,j) = rr(i,j)/n
          rr(j,i) = rr(i,j)
c vy is meant to be an estimate of the Monte Carlo variance of E(r) = E(Y);
c It is not quite right since the generated x's were not independent; fix this
c in a later version
          vy(i,j) = (rr(i,j) - r(i)*r(j))/n
          vy(j,i) = vy(i,j)
 50   continue
      return
      end
      
      
      subroutine gibbs2(ldn,burnin,nsize,mean,W,m,YQ,YC,r,rr,r3,r4,ldn2,
     *                  Z,Win,Vy)
      Parameter(MAXP=100,MAXQ=20,MAXN=200)
      implicit double precision(A-H,O-Z)
      integer ldn,burnin,nsize,m,YC(*),ldn2
      integer i1,i2,one,m2,i,j,l,ierr
      double precision mean(*),W(ldn,*),YQ(*),r(*),rr(ldn,*),r3(ldn,*),
     $     r4(ldn2,*)
      double precision Z(*),Win(ldn,*),Vy(ldn,*)
      double precision mwt(MAXN,MAXN),cvar(MAXN),u(MAXN),p(MAXN),x(MAXN)
      double precision d(2),sdev1
      double precision zero,cpnorm,cqnorm
      data one/1/,zero/0.0D0/
c     
      m2 = m*m
      call random(m,u)
      
      do 30 i = 1,m
        Z(i) = YQ(i)-mean(i)
        if (YC(i).eq.0) then
          x(i) = Z(i)
        else
          sdev1 = dsqrt(W(i,i))
          p(i) = cpnorm(Z(i),zero,sdev1)
          if (YC(i).lt.0) then
            u(i) = p(i)*u(i)
            x(i) = cqnorm(u(i),zero,sdev1)
          else
            u(i) = p(i)+(1.0d0-p(i))*u(i)
            x(i) = cqnorm(u(i),zero,sdev1)
          endif
          
          call dropi(ldn,m,W,i,Win)
          
          if (m.gt.1) call invert(ldn,m-1,Win,d,ierr)
          cvar(i) = W(i,i)
          i1 = 0
          do 20 j = 1,m
            mwt(i,j) = 0.0D0
            if (j.ne.i) then
              i1 = i1+1
              i2 = 0
              do 10 l = 1,m
                if (l.ne.i) then
                  i2 = i2+1
                  mwt(i,j) = mwt(i,j) + W(i,l)*Win(i2,i1)
                endif
 10           continue
              cvar(i) = cvar(i) - mwt(i,j)*W(j,i)
            endif
 20       continue
          cvar(i) = dsqrt(cvar(i))
        endif
 30   continue
      
      call gsample2(ldn,burnin,m,mwt,MAXN,cvar,Z,YC,x,r,rr,r3,r4,ldn2,
     *     Vy)
      call veczero(m,r)
      call matzero(ldn,m,m,rr)
      call matzero(ldn,m,m2,r3)
      call matzero(ldn2,m2,m2,r4)
      call gsample2(ldn,nsize,m,mwt,MAXN,cvar,Z,YC,x,r,rr,r3,r4,ldn2,
     *     Vy)
      
      do 40 i = 1,m
        Z(i) = r(i)+mean(i)
 40   continue
      
      return
      end
      
      
      subroutine gsample2(ldn,n,m,mwt,ldmw,sdev,q,cens,x,r,rr,r3,r4,
     *                    ldn2,Vy)
      Parameter(MAXP=100,MAXQ=20,MAXN=200)
      implicit double precision(A-H,O-Z)
      integer ldn,n,m,ldmw,cens(*),ldn2,ione,m2,i,j,l,i1,j1,j2,l1,l2
      double precision u(MAXN),cpnorm,cqnorm,p,cmu,rmin,rmax
      double precision mwt(ldmw,*),sdev(*),q(*),x(*),r(*),rr(ldn,*)
      double precision r3(ldn,*),r4(ldn2,*),Vy(ldn,*),eps,one
      data eps/1D-8/,one/1.0D0/,ione/1/
c     
      m2 = m*m
      do 5 i = 1,m
        if (cens(i).eq.0) then
           x(i) = q(i)
        end if
 5    continue       
      
      do 70 i1 = 1,n
        call random(m,u)
        do 20 i = 1,m
          if (cens(i).ne.0) then
            cmu = 0.0D0
            do 10 l = 1,m
c Note that mwt(i,i) = 0
 10         cmu = cmu + mwt(i,l)*x(l)
            p = cpnorm(q(i),cmu,sdev(i))
            rmin=0.0D0
            rmax=p
            if (cens(i).gt.0) then 
              rmin=p
              rmax=1.0D0
            endif
            u(i) = rmin + u(i)*(rmax-rmin)
c          if (u(i).lt.eps.or.u(i).gt.one-eps) then
c            x(i) = q(i)
c          else
            x(i) = cqnorm(u(i),cmu,sdev(i))
c          endif
          endif
 20     continue
        
        do 30 j = 1,m
          r(j) = r(j) + x(j)
          do 30 i = 1,j
            rr(i,j) = rr(i,j) + x(i)*x(j)
 30     continue
        
        do 60 j = 1,m2
          j1 = (j-1)/m + 1
          j2 = j - (j1-1)*m
          do 40 i=1,m
             r3(i,j) = r3(i,j) + x(i)*x(j1)*x(j2)
 40       continue
          do 50 l=1,j
            l1 = (l-1)/m + 1
            l2 = l - (l1-1)*m
            r4(l,j) = r4(l,j) + x(l1)*x(l2)*x(j1)*x(j2)
 50       continue
 60     continue
 70   continue
      
      do 80 j = 1,m
        r(j) = r(j)/n
        do 80 i = 1,j
           rr(i,j) = rr(i,j)/n
           rr(j,i) = rr(i,j)
           Vy(i,j) = (rr(i,j) - r(i)*r(j))/n
           Vy(j,i) = Vy(i,j)
 80   continue
      
      do 110 j=1,m2
        j1 = (j-1)/m + 1
        j2 = j - (j1-1)*m
        do 90 i=1,m
          r3(i,j) = r3(i,j)/n
 90     continue
        do 100 l=1,j
          r4(l,j) = r4(l,j)/n
          r4(j,l) = r4(l,j)
 100    continue
 110  continue
      
      return
      end
      
      
      subroutine random(n,x)
      integer i, n
      double precision crunif, x(*)
c     
      call rndstart()
      do 10 i=1,n
        x(i) = crunif()
 10   continue
      call rndend()
      return
      end
      
      
ccccccc Routines from CMLIB ccccccccccccccc
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(*)
      DOUBLE PRECISION A(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
C
C     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),DET(2),WORK(*)
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     DOUBLE PRECISION(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN DABS,MOD
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
C        OR DGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
