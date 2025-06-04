
C
C     Subroutines for survival analysis with left censoring covariate(s). 
C     This part of codes estimates the parameters.
C

C      
C   ChangeLog:
C
      
C     Copy a vector to another vector
      subroutine veccopy(a,la,b)
      
      integer la,i
      double precision a(*),b(*)
      
      do 10 i=1,la
        b(i)=a(i)
 10   continue
      
      return
      end
      
      
C     Copy a matrix to another matrix
      subroutine mtxcopy(a,lda,m,n,b,ldb)
      
      integer lda,m,n,ldb,i,j
      double precision a(lda,*),b(ldb,*)
      
      do 10 j=1,n
        do 10 i=1,m
          b(i,j)=a(i,j)
 10   continue
      
      return
      end
      
      
C     Generate the components (alphamat,beta,gamma) out of the vector betas
      subroutine genbtcom(betas,alphamat,ldam,px,q,beta,pv,gamma)
      
      integer ldam,px,q,pv
      double precision betas(*),alphamat(ldam,*),beta(*),gamma(*)
      
      integer k,l
      
      do 10 k=1,q
        do 10 l=1,px
          alphamat(l,k) = betas((k-1)*px+l)
 10   continue
      
      l = px*q
      do 20 k=1,pv
        beta(k) = betas(l+k)
 20   continue
      
      l = l + pv
      do 30 k=1,q
        gamma(k) = betas(l+k)
 30   continue
      
      return
      end
      
      
C     Convert an integer to string without leading white space
      subroutine int2str(k,kstr,klen)
      
      integer k,klen
      character kstr*(*)
      
      integer i,j,l,mod
      character char
      
c     calculate the length of input integer k
      klen = 1
      i = k
      
 10   i = i/10
      if (i.eq.0) go to 20
      klen = klen + 1
      go to 10
 20   continue
c      call intpr('klen',-1,klen,1)
c      return
      
c     output the individual digits to kstr
      i = k
      do 30 j=1,klen
        l = 10**(klen-j)
        kstr(j:j) = char(48+i/l)
        i = mod(i,l)
 30   continue
c      call intpr(kstr,klen,0,0)
      
      return
      end
      
      
C     Compute the negative log posterior density of bi given yi and etc
      double precision function blogden(bi,qb,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse)
      
      parameter (MAXQ=4)
      integer qb,ldsb,deltai,tid,ldxt,px,ldzt,ldam,q,ldyc,ni,ldzi,ldse
      double precision bi(*),sigbinv(ldsb,*),vibeta,h0t(*),xit(ldxt,*),
     $     zit(ldzt,*),alphamat(ldam,*),gamma(*),yic(ldyc,*),zi(ldzi,*),
     $     sigeinv(ldse,*)
      
      integer j,k,l,s,idx
      double precision lli,bim(2,MAXQ),rij(MAXQ),tmpvec1(MAXQ),tmpsum1,
     $     tmpsum2,dexp
      
c      call intpr('qb',-1,qb,1)
c      call dblepr('bi',-1,bi,qb)
c      call intpr('ldsb',-1,ldsb,1)
c      call dblepr('sigbinv(,1)',-1,sigbinv(1,1),qb)
c      call dblepr('sigbinv(,2)',-1,sigbinv(1,2),qb)
c      call dblepr('sigbinv(,3)',-1,sigbinv(1,3),qb)
c      call dblepr('sigbinv(,4)',-1,sigbinv(1,4),qb)
c      return
      
      lli = 0.0d0
      
c     contribution to the log-likelihood due to the normality of bi
      do 10 k=1,qb
c       add the diagonal terms first b/c of the double subtractions right after
        lli = lli + sigbinv(k,k)*bi(k)**2/2.0d0
        do 10 l=1,k
          lli = lli - bi(l)*sigbinv(l,k)*bi(k)
 10   continue
c      call dblepr('lli',-1,lli,1)
c      return
      
c     transform bi from the vector form to matrix form (dim of 2 by q)
      do 20 k=1,q
        do 20 l=1,2
          bim(l,k) = bi((k-1)*2+l)
 20   continue
c      call dblepr('bim(,1)',-1,bim(1,1),2)
c      call dblepr('bim(,q)',-1,bim(1,q),2)
c      return
      
c     contribution to the log-likelihood due to the normality of ri_s
      do 60 j=1,ni
        do 40 k=1,q
          tmpsum1 = 0.0d0
          do 30 l=1,2
            tmpsum1 = tmpsum1 + zi(j,l)*bim(l,k)
 30       continue
          rij(k) = yic(j,k) - tmpsum1
 40     continue
        do 50 k=1,q
          lli = lli + sigeinv(k,k)*rij(k)**2/2.0d0
          do 50 l=1,k
            lli = lli - rij(l)*sigeinv(l,k)*rij(k)
 50     continue                
 60   continue
c      call dblepr('lli',-1,lli,1)
c      return
      
c     contribution to the log-likelihood due to being an event: log(ht)
      if (deltai .eq. 1) then
        do 70 k=1,q
          do 70 l=1,2
            lli = lli + zit(tid+1,l)*bim(l,k)*gamma(k)
 70     continue
      end if
c      call dblepr('lli',-1,lli,1)
c      return
      
c     contribution to the log-likelihood due to log(St)=-Ht
      if (tid .gt. 0) then
        tmpsum1 = 0.0d0
        do 110 k=1,tid
          tmpsum2 = 0.0d0
          
          do 100 l=1,q
            tmpvec1(l) = 0.0d0
            do 80 s=1,px
              tmpvec1(l) = tmpvec1(l) + xit(k,s)*alphamat(s,l)
 80         continue
            do 90 s=1,2
              tmpvec1(l) = tmpvec1(l) + zit(k,s)*bim(s,l)
 90         continue
            tmpsum2 = tmpsum2 + tmpvec1(l)*gamma(l)
 100      continue
          
          tmpsum1 = tmpsum1 + h0t(k)*dexp(tmpsum2)
 110    continue
        lli = lli - dexp(vibeta)*tmpsum1
      end if
c      call dblepr('lli',-1,lli,1)
c      return
      
      blogden = -lli
      
      return
      end
      
      
C     Compute the negative log posterior density of bi given yi (ver 2)
      double precision function blogden2(bi,qb,dt,s,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse)
      
      parameter (MAXQ=4)
      integer qb,ldsb,deltai,tid,ldxt,px,ldzt,ldam,q,ldyc,ni,ldzi,ldse
      double precision bi(*),dt(*),s,sigbinv(ldsb,*),vibeta,h0t(*),
     $     xit(ldxt,*),zit(ldzt,*),alphamat(ldam,*),gamma(*),
     $     yic(ldyc,*),zi(ldzi,*),sigeinv(ldse,*),blogden
      
      integer k
      double precision bwk(2*MAXQ)
      
      do 10 k=1,qb
        bwk(k)=bi(k)+s*dt(k)
 10   continue
      
      blogden2=blogden(bwk,qb,sigbinv,ldsb,deltai,vibeta,h0t,tid,
     $     xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse)
      
      return
      end
      
      
C     Compute the gradient of negative log posterior density of bi given yi
      subroutine getgrdb(bi,qb,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse,grad)
      
      parameter (MAXQ=4)
      integer qb,ldsb,deltai,tid,ldxt,px,ldzt,ldam,q,ldyc,ni,ldzi,ldse
      double precision bi(*),sigbinv(ldsb,*),vibeta,h0t(*),xit(ldxt,*),
     $     zit(ldzt,*),alphamat(ldam,*),gamma(*),yic(ldyc,*),zi(ldzi,*),
     $     sigeinv(ldse,*),grad(*)
      
      integer j,k,l,s,idx
      double precision bim(2,MAXQ),tmpmat1(MAXQ,2*MAXQ),rij(MAXQ),
     $     tmpsum1,tmpsum2,tmpval1,dexp
      
c      call intpr('q',-1,q,1)
c      call dblepr('tmpmat1(,1)',-1,tmpmat1(1,1),q)
c      call dblepr('tmpmat1(,4)',-1,tmpmat1(1,4),q)
c      return
      
      do 10 k=1,qb
        grad(k) = 0.0d0
        do 10 l=1,qb
          grad(k) = grad(k) + sigbinv(l,k)*bi(l)
 10   continue
c      call dblepr('grad',-1,grad,qb)
c      return
      
      
c     transform bi from the vector form to matrix form (dim of 2 by q)
      do 20 k=1,q
        do 20 l=1,2
          bim(l,k) = bi((k-1)*2+l)
 20   continue
c      call dblepr('bim(,1)',-1,bim(1,1),2)
c      call dblepr('bim(,q)',-1,bim(1,q),2)
c      return
      
      call matzero(MAXQ,q,qb,tmpmat1)      
c      call dblepr('tmpmat1(,1)',-1,tmpmat1(1,1),q)
c      call dblepr('tmpmat1(,4)',-1,tmpmat1(1,4),q)
c      return
      do 70 j=1,ni
        do 50 k=1,q
          tmpsum1 = 0.0d0
          do 30 l=1,2
            tmpsum1 = tmpsum1 + zi(j,l)*bim(l,k)
 30       continue                    
          rij(k) = yic(j,k) - tmpsum1
          do 40 l=1,2
            tmpmat1(k,(k-1)*2+l) = zi(j,l)
 40       continue
 50     continue
        
        do 60 s=1,qb
          do 60 l=1,q
            do 60 k=1,q
              grad(s) = grad(s) - rij(k)*sigeinv(k,l)*tmpmat1(l,s)
 60     continue
 70   continue
c      call dblepr('grad',-1,grad,qb)
c      return
      
      
      if (deltai .eq. 1) then
        do 80 k=1,q
          do 80 l=1,2
            idx = (k-1)*2 + l
            grad(idx) = grad(idx) - gamma(k)*zit(tid+1,l)
 80     continue
      end if
c      call dblepr('grad',-1,grad,qb)
c      return
      
      
      if (tid .gt. 0) then
        do 130 j=1,tid
          tmpsum1 = 0.0d0
          do 110 k=1,q
            tmpsum2 = 0.0d0
            do 90 l=1,px
              tmpsum2 = tmpsum2 + xit(j,l)*alphamat(l,k)
 90         continue
            do 100 l=1,2
              tmpsum2 = tmpsum2 + zit(j,l)*bim(l,k)
 100        continue            
            tmpsum1 = tmpsum1 + tmpsum2*gamma(k)
 110      continue
          tmpval1 = h0t(j)*dexp(vibeta+tmpsum1)
          do 120 k=1,q
            do 120 l=1,2
              idx = (k-1)*2 + l
              grad(idx) = grad(idx) + tmpval1*gamma(k)*zit(j,l)
 120      continue
 130    continue        
      end if
c      call dblepr('grad',-1,grad,qb)
      
      return
      end
      
      
C     Compute the hessian of negative log posterior density of bi given yi (i.e., 
C     information matrix)
      subroutine gethessb(bi,qb,sigbinv,ldsb,vibeta,h0t,tid,xit,ldxt,px,
     $     zit,ldzt,alphamat,ldam,q,gamma,zi,ldzi,ni,sigeinv,ldse,
     $     hess,ldhs)
      
      parameter (MAXQ=4)
      integer qb,ldsb,tid,ldxt,px,ldzt,ldam,q,ldzi,ni,ldse,ldhs
      double precision bi(*),sigbinv(ldsb,*),vibeta,h0t(*),xit(ldxt,*),
     $     zit(ldzt,*),alphamat(ldam,*),gamma(*),zi(ldzi,*),
     $     sigeinv(ldse,*),hess(ldhs,*)
      
      integer j,k,l,s,t,idx
      double precision bim(2,MAXQ),tmpmat1(MAXQ,2*MAXQ),tmpvec1(2*MAXQ),
     $     tmpsum1,tmpsum2,tmpval1,dexp
      
      do 10 k=1,qb
        do 10 l=1,k
          hess(l,k) = sigbinv(l,k)
 10   continue
c      call dblepr('hess(,1)',-1,hess(1,1),qb)
c      call dblepr('hess(,4)',-1,hess(1,4),qb)
c      return
      
      call matzero(MAXQ,q,qb,tmpmat1)      
      do 40 j=1,ni
        do 20 k=1,q
          do 20 l=1,2
            tmpmat1(k,(k-1)*2+l) = zi(j,l)
 20     continue
        do 30 k=1,qb
          do 30 l=1,k
            do 30 s=1,q
              do 30 t=1,q
                hess(l,k) = hess(l,k) + 
     $                tmpmat1(t,l)*sigeinv(t,s)*tmpmat1(s,k)               
 30     continue
 40   continue
c      call dblepr('hess(,1)',-1,hess(1,1),qb)
c      call dblepr('hess(,4)',-1,hess(1,4),qb)
c      return
      
      if (tid .gt. 0) then
c       transform bi from the vector form to matrix form (dim of 2 by q)
        do 45 k=1,q
          do 45 l=1,2
            bim(l,k) = bi((k-1)*2+l)
 45     continue
        
        do 100 j=1,tid
          do 50 k=1,q
            do 50 l=1,2
              tmpvec1((k-1)*2+l) = gamma(k)*zit(j,l)
 50       continue
          
          tmpsum1 = 0.0d0
          do 80 k=1,q
            tmpsum2 = 0.0d0
            do 60 l=1,px
              tmpsum2 = tmpsum2 + xit(j,l)*alphamat(l,k)
 60         continue
            do 70 l=1,2
              tmpsum2 = tmpsum2 + zit(j,l)*bim(l,k)
 70         continue            
            tmpsum1 = tmpsum1 + tmpsum2*gamma(k)
 80       continue
          tmpval1 = h0t(j)*dexp(vibeta+tmpsum1)
          
          do 90 k=1,qb
            do 90 l=1,k
              hess(l,k) = hess(l,k) + tmpval1*tmpvec1(l)*tmpvec1(k)
 90       continue
 100    continue
      end if
c      call dblepr('hess(,1)',-1,hess(1,1),qb)
c      call dblepr('hess(,4)',-1,hess(1,4),qb)
c      return
      
      do 110 k=1,qb
        do 110 l=1,k
          if (l .ne. k) hess(k,l) = hess(l,k)
 110  continue
      
      return
      end
      
      
C     Golden section search
c      subroutine gssch(bc,delta,sigbinv,ldsb,qb,yi,ni,etaif,betar2,xic,
c     $          ldxc,px,sigeinv,ldse,maxiter3,eps3,sopt,bwk)
c
c      integer ldsb,qb,yi(*),ni,ldxc,px,ldse,maxiter3
c      double precision bc(*),delta(*),sigbinv(ldsb,*),etaif(*),
c     $     betar2(*),xic(ldxc,*),sigeinv(ldse,*),eps3,sopt,bwk(*)
c
c      integer k
c      double precision a,b,c,dsqrt,x1,f1,x2,f2,dabs,blogden2
c
c      c = (dsqrt(5.0d0)-1.0d0)/2.0d0
c      a = 0.0d0
c      b = 1.0d0
c      x1 = c*a + (1.0d0-c)*b
c      f1 = blogden2(bc,delta,x1,sigbinv,ldsb,qb,yi,ni,
c     $     etaif,betar2,xic,ldxc,px,sigeinv,ldse,bwk)
c      x2 = (1.0d0-c)*a + c*b
c      f2 = blogden2(bc,delta,x2,sigbinv,ldsb,qb,yi,ni,
c     $     etaif,betar2,xic,ldxc,px,sigeinv,ldse,bwk)
c      do 10 k=1,maxiter3
c        if (f1.lt.f2) then
c          b = x2
c          x2 = x1
c          f2 = f1
c          x1 = c*a + (1.0d0-c)*b
c          f1 = blogden2(bc,delta,x1,sigbinv,ldsb,qb,yi,ni,
c     $         etaif,betar2,xic,ldxc,px,sigeinv,ldse,bwk)
c        else
c          a = x1
c          x1 = x2
c          f1 = f2
c          x2 = (1.0d0-c)*a + c*b
c          f2 = blogden2(bc,delta,x2,sigbinv,ldsb,qb,yi,ni,
c     $         etaif,betar2,xic,ldxc,px,sigeinv,ldse,bwk)
c        end if
c        if (dabs(x2-x1).lt.eps3) goto 20
c 10   continue
c
c 20   sopt = (x1+x2)/2.0d0
c
c      return
c      end
      
      
C     Line search with "Newton" direction: Armijo rule and backtracking method
      subroutine linesch(bi,qb,f0,grad,inc,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse,maxiter,lambda,info,iter)
      
      integer qb,ldsb,deltai,tid,ldxt,px,ldzt,ldam,q,ldyc,ni,ldzi,ldse,
     $     maxiter,info,iter
      double precision bi(*),f0,grad(*),inc(*),sigbinv(ldsb,*),vibeta,
     $     h0t(*),xit(ldxt,*),zit(ldzt,*),alphamat(ldam,*),gamma(*),
     $     yic(ldyc,*),zi(ldzi,*),sigeinv(ldse,*),lambda
      
      integer k
      double precision alpha,beta,fc,tmpval1,blogden2
      
      alpha = 1.0d-2
      beta  = 5.0d-1
      
      iter = 0
      lambda = 1.0d0
      fc = blogden2(bi,qb,inc,lambda,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse)
c      call dblepr('fc',-1,fc,1)
c      return
      
      tmpval1 = 0.0d0
      do 10 k=1,qb
        tmpval1 = tmpval1 + inc(k)*grad(k)
 10   continue
      
 20   continue
      if ((iter .lt. maxiter) .and. (fc .gt. (f0+alpha*lambda*tmpval1)))
     $     then
        iter = iter + 1
        lambda = beta*lambda
        fc = blogden2(bi,qb,inc,lambda,sigbinv,ldsb,deltai,
     $       vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,
     $       gamma,yic,ldyc,ni,zi,ldzi,sigeinv,ldse)           
        go to 20           
      end if
c      call intpr('iter',-1,iter,1)
c      return
      
      if (fc .le. (f0+alpha*lambda*tmpval1)) then 
        info = 0
      else 
        info = 1
      end if
c      call intpr('info',-1,info,1)
      
      return
      end
      
      
C     N-step Newton-Raphson with line searching
      subroutine nsnrls(bi0,qb,nstep,eps1,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse,maxiter3,bi1,h1,ldhs,info,
     $     iter)
      
      parameter (MAXQ=4)
      parameter (MAXQB=2*MAXQ)
      integer qb,nstep,ldsb,deltai,tid,ldxt,px,ldzt,ldam,q,ldyc,ni,ldzi,
     $     ldse,maxiter3,ldhs,info,iter
      double precision bi0(*),eps1,sigbinv(ldsb,*),vibeta,
     $     h0t(*),xit(ldxt,*),zit(ldzt,*),alphamat(ldam,*),gamma(*),
     $     yic(ldyc,*),zi(ldzi,*),sigeinv(ldse,*),bi1(*),h1(ldhs,*)
      
      integer k,hcnt,ldhc,lddm,ierr,lwork,ipvt(MAXQB),job,iter2
      double precision eps2,f1,blogden,g1(MAXQB),tmpval1,dabs,
     $     hc(MAXQB,MAXQB),hc2(MAXQB,MAXQB),eigvals(MAXQB),
     $     dummyv(MAXQB),dummym(MAXQB,MAXQB),work(6*MAXQB),inc(MAXQB),
     $     tmpscale,lambda
      character*1 jobvl,jobvr
      
c      call intpr('qb',-1,qb,1)
c      call dblepr('bi0',-1,bi0,qb)
c      call intpr('nstep',-1,nstep,1)
c      call dblepr('eps1',-1,eps1,1)
c      call intpr('ldsb',-1,ldsb,1)
c      call dblepr('sigbinv(,1)',-1,sigbinv(1,1),qb)
c      call dblepr('sigbinv(,4)',-1,sigbinv(1,4),qb)
c      call intpr('deltai',-1,deltai,1)
c      call dblepr('vibeta',-1,vibeta,1)
c      call intpr('tid',-1,tid,1)
c      call dblepr('h0t',-1,h0t,tid)
c      call intpr('ldxt',-1,ldxt,1)
c      call intpr('px',-1,px,1)
c      call dblepr('xit(,1)',-1,xit(1,1),tid+1)
c      call dblepr('xit(,3)',-1,xit(1,3),tid+1)
c      call intpr('ldzt',-1,ldzt,1)
c      call dblepr('zit(,1)',-1,zit(1,1),tid+1)
c      call dblepr('zit(,2)',-1,zit(1,2),tid+1)
c      call intpr('ldam',-1,ldam,1)
c      call intpr('q',-1,q,1)
c      call dblepr('alphamat(,1)',-1,alphamat(1,1),px)
c      call dblepr('alphamat(,2)',-1,alphamat(1,2),px)
c      call dblepr('gamma',-1,gamma,q)
c      call intpr('ldyc',-1,ldyc,1)
c      call intpr('ni',-1,ni,1)
c      call dblepr('yic(,1)',-1,yic(1,1),ni)
c      call dblepr('yic(,2)',-1,yic(1,2),ni)
c      call intpr('ldzi',-1,ldzi,1)
c      call dblepr('zi(,1)',-1,zi(1,1),ni)
c      call dblepr('zi(,2)',-1,zi(1,2),ni)
c      call intpr('ldse',-1,ldse,1)
c      call dblepr('sigeinv(,1)',-1,sigeinv(1,1),q)
c      call dblepr('sigeinv(,2)',-1,sigeinv(1,2),q)
c      call intpr('maxiter3',-1,maxiter3,1)
c      call dblepr('bi1',-1,bi1,qb)
c      call intpr('ldhs',-1,ldhs,1)
c      call dblepr('h1(,1)',-1,h1(1,1),qb)
c      call dblepr('h1(,4)',-1,h1(1,4),qb)
c      call intpr('info',-1,info,1)
c      call intpr('iter',-1,iter,1)
c      return
      
      eps2 = 1.0d-2
      ldhc = MAXQB
      lddm = MAXQB
      lwork = 6*MAXQB
      jobvl = 'N'
      jobvr = 'N'
      job = 0
      
      iter = 0
      call veccopy(bi0,qb,bi1)
c      call dblepr('bi1',-1,bi1,qb)
c      return
      
      f1 = blogden(bi1,qb,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse)
c      call dblepr('f1',-1,f1,1)
c      return
      call getgrdb(bi1,qb,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse,g1)
c      call dblepr('g1',-1,g1,qb)
c      return
      
 10   tmpval1 = dabs(g1(1))
      do 20 k=2,qb
        if (dabs(g1(k)) .gt. tmpval1) tmpval1 = dabs(g1(k))
 20   continue
c      call dblepr('tmpval1',-1,tmpval1,1)
c      return
      
c      call intpr('iter',-1,iter,1)
c      call dblepr('eps1',-1,eps1,1)      
      if ((iter .lt. nstep) .and. (tmpval1 .gt. eps1)) then
        call gethessb(bi1,qb,sigbinv,ldsb,vibeta,h0t,tid,xit,ldxt,px,
     $        zit,ldzt,alphamat,ldam,q,gamma,zi,ldzi,ni,sigeinv,ldse,
     $        h1,ldhs)
c        call dblepr('h1(,1)',-1,h1(1,1),qb)
c        call dblepr('h1(,2)',-1,h1(1,2),qb)
c        call dblepr('h1(,3)',-1,h1(1,3),qb)
c        call dblepr('h1(,4)',-1,h1(1,4),qb)
c        return
        call mtxcopy(h1,ldhs,qb,qb,hc,ldhc)
        hcnt = 0
c        call dblepr('hc(,1)',-1,hc(1,1),qb)
c        call dblepr('hc(,4)',-1,hc(1,4),qb)
c        return
        
 30     call mtxcopy(hc,ldhc,qb,qb,hc2,ldhc)
c        call dblepr('hc2(,1)',-1,hc2(1,1),qb)
c        call dblepr('hc2(,4)',-1,hc2(1,4),qb)
c        return
        call dgeev(jobvl,jobvr,qb,hc2,ldhc,eigvals,dummyv,dummym,lddm,
     $       dummym,lddm,work,lwork,ierr)  
c        call intpr('ierr',-1,ierr,1)
c        call dblepr('eigvals',-1,eigvals,qb)
c        return
        
        tmpval1 = eigvals(1)
        do 40 k=2,qb
          if (eigvals(k) .lt. tmpval1) tmpval1 = eigvals(k)
 40     continue
c        call dblepr('tmpval1',-1,tmpval1,1)
c        return
        
        if (tmpval1 .gt. 0.0d0) then
          call mtxcopy(hc,ldhc,qb,qb,hc2,ldhc)
c        call dblepr('hc2(,1)',-1,hc2(1,1),qb)
c        call dblepr('hc2(,4)',-1,hc2(1,4),qb)
c        return
          call dgefa(hc2,ldhc,qb,ipvt,ierr)
c        call intpr('ierr',-1,ierr,1)
c        return
          call veccopy(g1,qb,inc)
c        call dblepr('inc',-1,inc,qb)
c        return
          call dgesl(hc2,ldhc,qb,ipvt,inc,job)
c        call dblepr('inc',-1,inc,qb)
c        return
          call dscal(qb,-1.0d0,inc,1)
c        call dblepr('inc',-1,inc,qb)
c        return
          go to 60
        else
          tmpscale = dabs(tmpval1) + eps2
          do 50 k=1,qb
            hc(k,k) = hc(k,k) + tmpscale
 50       continue
          hcnt = hcnt + 1
          go to 30
        end if
        
 60     call linesch(bi1,qb,f1,g1,inc,sigbinv,ldsb,deltai,
     $       vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $       yic,ldyc,ni,zi,ldzi,sigeinv,ldse,maxiter3,lambda,ierr,
     $       iter2)
c        call dblepr('lambda',-1,lambda,1)
c        return
        iter = iter + 1
        do 70 k=1,qb
          bi1(k) = bi1(k) + lambda*inc(k)
 70     continue
c      call dblepr('bi1',-1,bi1,qb)
c      return
        
        f1 = blogden(bi1,qb,sigbinv,ldsb,deltai,
     $       vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $       yic,ldyc,ni,zi,ldzi,sigeinv,ldse)
c      call dblepr('f1',-1,f1,1)
c      return
        call getgrdb(bi1,qb,sigbinv,ldsb,deltai,
     $       vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $       yic,ldyc,ni,zi,ldzi,sigeinv,ldse,g1)        
c      call dblepr('g1',-1,g1,qb)
c      return
        go to 10           
      end if
      
      tmpval1 = dabs(g1(1))
      do 80 k=2,qb
        if (dabs(g1(k)) .gt. tmpval1) tmpval1 = dabs(g1(k))
 80   continue
c      call dblepr('tmpval1',-1,tmpval1,1)
c      return
      
      if (tmpval1 .le. eps1) then
        info = 0
      else 
        info = 1
      end if
c      call intpr('info',-1,info,1)
      
      call gethessb(bi1,qb,sigbinv,ldsb,vibeta,h0t,tid,xit,ldxt,px,
     $        zit,ldzt,alphamat,ldam,q,gamma,zi,ldzi,ni,sigeinv,ldse,
     $        h1,ldhs)
c      call dblepr('h1(,1)',-1,h1(1,1),qb)
c      call dblepr('h1(,2)',-1,h1(1,2),qb)
c      call dblepr('h1(,3)',-1,h1(1,3),qb)
c      call dblepr('h1(,4)',-1,h1(1,4),qb)
      
      return
      end
      
      
C     Set up Gibbs sampler for drawing xi given bi and yi
      subroutine setgibbs(sige,ldse,qe,gmat,ldgm,cstd)
      
      parameter(MAXQ=4)
      integer ldse,qe,ldgm
      double precision sige(ldse,*),gmat(ldgm,*),cstd(*)
      
      integer k,l,s,ldtm,qe2,cnt1,cnt2,ierr
      double precision tmpmat(MAXQ-1,MAXQ-1),det(2),dsqrt
      
      call matzero(ldgm,qe,qe,gmat)
      
c     save sige[-k,-k] to tmpmat for each k from 1 to qe
      ldtm = MAXQ - 1
      qe2  = qe - 1
      do 60 k=1,qe
        cnt1 = 0
        do 20 l=1,qe
          if (l.ne.k) then
            cnt1=cnt1+1
            cnt2 = 0
            do 10 s=1,qe
              if (s.ne.k) then
                cnt2=cnt2+1
                tmpmat(cnt2,cnt1)=sige(s,l)
              end if
 10         continue
          end if
 20     continue
        call invert(ldtm,qe2,tmpmat,det,ierr)
        
c       compute the conditional mean and variance for each k being conditioned
        cnt1 = 0
        do 40 l=1,qe
          if (l.ne.k) then
            cnt1 = cnt1+1
            cnt2 = 0
            do 30 s=1,qe
              if (s.ne.k) then
                 cnt2 = cnt2 + 1
                 gmat(l,k)=gmat(l,k)+sige(s,k)*tmpmat(cnt2,cnt1)
              end if
 30         continue
          end if
 40     continue
        
        cstd(k) = sige(k,k)
        do 50 l=1,qe
          cstd(k)=cstd(k)-gmat(l,k)*sige(l,k)
 50     continue
        cstd(k) = dsqrt(cstd(k))
 60   continue
      
      return
      end
      
      
C     Get Monte-Carlo initial value for bi and yic
      subroutine getycbm0(qic,ldqc,ni,q,dila,dil,di,lddi,bi0,qb,
     $     sigbinv,ldsb,deltai,vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,
     $     alphamat,ldam,gamma,yic,ldyc,zi,ldzi,sigeinv,ldse,
     $     maxiter2,eps1,maxiter3,gmat,ldgm,cstd,bimode,bicurv,ldbc,
     $     bhinv,ldbh,bim0,yic0,ldyc2)
      
c     Note: bi0 is the initial bi for Newton-Raphson with line search method,
c           yic is the input yic for the Newton-Raphson method, and (bim0, yic0) 
c           is the output (bi, yic) for initiating the MCMC algorithm on the pair.
      parameter (MAXQ=4)
      integer ldqc,ni,q,dila,dil(*),lddi,qb,ldsb,deltai,tid,ldxt,px,
     $     ldzt,ldam,ldyc,ldzi,ldse,maxiter2,maxiter3,ldgm,ldbc,ldbh,
     $     ldyc2
      integer di(lddi,*)
      double precision qic(ldqc,*),bi0(*),sigbinv(ldsb,*),vibeta,
     $     h0t(*),xit(ldxt,*),zit(ldzt,*),alphamat(ldam,*),gamma(*),
     $     yic(ldyc,*),zi(ldzi,*),sigeinv(ldse,*),eps1,gmat(ldgm,*),
     $     cstd(*),bimode(*),bicurv(ldbc,*),bhinv(ldbh,*),bim0(*),
     $     yic0(ldyc2,*)
      
      integer j,k,l,s,ierr,iter,jpvt(2*MAXQ)
      double precision det(2),r(2*MAXQ),bimm(2,MAXQ),tmpvec1(MAXQ),
     $     cmu,pp,cpnorm,cqnorm
      
c      call intpr('ldqc',-1,ldqc,1)
c      call intpr('ni',-1,ni,1)
c      call intpr('q',-1,q,1)	 
c      call dblepr('qic(,1)',-1,qic(1,1),ni)
c      call dblepr('qic(,2)',-1,qic(1,2),ni)
c      call intpr('dila',-1,dila,1)
c      call intpr('dil',-1,dil,ni)
c      call intpr('lddi',-1,lddi,1)
c      call dblepr('di(,1)',-1,di(1,1),ni)
c      call dblepr('di(,2)',-1,di(1,2),ni)
c      call intpr('ldddi',-1,ldddi,1)
c      call dblepr('ddi(,1)',-1,ddi(1,1),ni)
c      call dblepr('ddi(,2)',-1,ddi(1,2),ni)
c      call intpr('qb',-1,qb,1)
c      call dblepr('bi0',-1,bi0,qb)
c      call intpr('ldsb',-1,ldsb,1)
c      call dblepr('sigbinv(,1)',-1,sigbinv(1,1),qb)
c      call dblepr('sigbinv(,4)',-1,sigbinv(1,4),qb)     
c      call intpr('deltai',-1,deltai,1)
c      call dblepr('vibeta',-1,vibeta,1)
c      call intpr('tid',-1,tid,1)
c      call dblepr('h0t',-1,h0t,tid)
c      call intpr('ldxt',-1,ldxt,1)
c      call intpr('px',-1,px,1)
c      call dblepr('xit(,1)',-1,xit(1,1),tid+1)
c      call dblepr('xit(,3)',-1,xit(1,3),tid+1)
c      call intpr('ldzt',-1,ldzt,1)
c      call dblepr('zit(,1)',-1,zit(1,1),tid+1)
c      call dblepr('zit(,2)',-1,zit(1,2),tid+1)
c      call intpr('ldam',-1,ldam,1)
c      call dblepr('alphamat(,1)',-1,alphamat(1,1),px)
c      call dblepr('alphamat(,2)',-1,alphamat(1,2),px)
c      call dblepr('gamma',-1,gamma,q)
c      call intpr('ldyc',-1,ldyc,1)
c      call dblepr('yic(,1)',-1,yic(1,1),ni)
c      call dblepr('yic(,2)',-1,yic(1,2),ni)
c      call intpr('ldzi',-1,ldzi,1)
c      call dblepr('zi(,1)',-1,zi(1,1),ni)
c      call dblepr('zi(,2)',-1,zi(1,2),ni)
c      call intpr('ldse',-1,ldse,1)
c      call dblepr('sigeinv(,1)',-1,sigeinv(1,1),q)
c      call dblepr('sigeinv(,2)',-1,sigeinv(1,2),q)     
c      call intpr('maxiter2',-1,maxiter2,1)
c      call dblepr('eps1',-1,eps1,1)
c      call intpr('maxiter3',-1,maxiter3,1)
c      call intpr('ldgm',-1,ldgm,1)
c      call dblepr('gmat(,1)',-1,gmat(1,1),q)
c      call dblepr('gmat(,2)',-1,gmat(1,2),q)
c      call dblepr('cstd',-1,cstd,q)
c      return
      
c     make sure yic0 is the same as qic for those uncensored observations
      do 10 j=1,ni
        do 10 k=1,q
          if (di(j,k).eq.0) yic(j,k)=qic(j,k)
 10   continue
      
      call nsnrls(bi0,qb,maxiter2,eps1,sigbinv,ldsb,deltai,
     $     vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $     yic,ldyc,ni,zi,ldzi,sigeinv,ldse,maxiter3,bimode,bicurv,ldbc,
     $     ierr,iter)
c      call intpr('ierr',-1,ierr,1)
c      call intpr('iter',-1,iter,1)
c      call dblepr('bimode',-1,bimode,qb)
c      call dblepr('bicurv(,1)',-1,bicurv(1,1),qb)
c      call dblepr('bicurv(,2)',-1,bicurv(1,2),qb)
c      call dblepr('bicurv(,3)',-1,bicurv(1,3),qb)
c      call dblepr('bicurv(,4)',-1,bicurv(1,4),qb)
c      return
      
      call mtxcopy(bicurv,ldbc,qb,qb,bhinv,ldbh)
      call invert(ldbh,qb,bhinv,det,ierr)
      call dchdc(bhinv,ldbh,qb,r,jpvt,0,ierr)
c     note: ierr is the index of last positive diagonal element of the cholesky factor
c      call intpr('ierr',-1,ierr,1)
c      return
      
      do 20 k=1,qb-1
        do 20 l=k+1,qb
          bhinv(l,k) = 0.0d0
 20   continue
c      call dblepr('bhinv(,1)',-1,bhinv(1,1),qb)
c      call dblepr('bhinv(,2)',-1,bhinv(1,2),qb)
c      call dblepr('bhinv(,3)',-1,bhinv(1,3),qb)
c      call dblepr('bhinv(,4)',-1,bhinv(1,4),qb)
c      return
      
      call randomn(qb,r)
      do 30 k=1,qb
        bim0(k) = bimode(k)
        do 30 l=1,qb
          bim0(k)=bim0(k)+bhinv(l,k)*r(l)
 30   continue
c      call dblepr('bim0',-1,bim0,qb)
c      return
      
      call mtxcopy(yic,ldyc,ni,q,yic0,ldyc2)
c      call dblepr('yic0(,1)',-1,yic0(1,1),ni)
c      call dblepr('yic0(,2)',-1,yic0(1,2),ni)
c      return
      
      if (dila.eq.1) then
        do 40 k=1,q
          do 40 l=1,2
            bimm(l,k) = bim0((k-1)*2+l)
 40     continue
c      end if
        
c      call mtxcopy(yic,ldyc,ni,q,yic0,ldyc2)
        do 80 j=1,ni        
          if (dil(j).eq.1) then
            do 50 k=1,q
              tmpvec1(k) = 0.0d0
              do 50 l=1,2
                tmpvec1(k) = tmpvec1(k) + zi(j,l)*bimm(l,k)
 50         continue
            
            do 70 l=1,q
              if (di(j,l).eq.1) then
                cmu = tmpvec1(l)
                do 60 s=1,q
                  if (s.ne.l) cmu = cmu + 
     $                  gmat(s,l)*(yic0(j,s)-tmpvec1(s))
 60             continue
                call random(1,r)
                pp = cpnorm(qic(j,l),cmu,cstd(l))*r(1)
                yic0(j,l)=cqnorm(pp,cmu,cstd(l))
              end if
 70         continue
          end if
 80     continue
      end if
      
      return
      end
      
      
C     Metropolized Gibbs sampler
      subroutine mgibbs(qic,ldqc,ni,q,dila,dil,di,lddi,sigbinv,ldsb,qb,
     $     deltai,vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,
     $     gamma,zi,ldzi,sigeinv,ldse,gmat,ldgm,cstd,bimode,bicurv,
     $     ldbc,bhinv,ldbh,bim,yicm,ldycm,biarr,ldba,nburnin,nsize,
     $     eycbi,ldey,eycb2i,ldey2,wkzb,ldwz)
      
c     Note: bim and yicm the initial bim0 and yic0 at the beginning and
c        are overwritten with the last bim and yic at the end of MCMC.
c        wkzb is a work array of the dimension q*maxni corresp. to t(Zi*bimm).
      parameter (MAXQ=4)
      integer ldqc,ni,q,dila,dil(*),lddi,ldsb,qb,deltai,tid,ldxt,px,
     $     ldzt,ldam,ldzi,ldse,ldgm,ldbc,ldbh,ldycm,ldba,nburnin,nsize,
     $     ldey,ldey2,ldwz
      integer di(lddi,*)
      double precision qic(ldqc,*),sigbinv(ldsb,*),vibeta,
     $     h0t(*),xit(ldxt,*),zit(ldzt,*),alphamat(ldam,*),gamma(*),
     $     zi(ldzi,*),sigeinv(ldse,*),gmat(ldgm,*),cstd(*),bimode(*),
     $     bicurv(ldbc,*),bhinv(ldbh,*),bim(*),yicm(ldycm,*),
     $     biarr(ldba,*),eycbi(ldey,*),eycb2i(ldey2,*),wkzb(ldwz,*)
      
      integer j,k,l,s
      double precision r(2*MAXQ),bic(2*MAXQ),bimm(2,MAXQ),cmu,pp,ap,
     $     blogden,cpnorm,cqnorm
      
      call matzero(ldey,ni,q,eycbi)
      call matzero(ldey2,q,q,eycb2i)
c      call dblepr('eycbi(,1)',-1,eycbi(1,1),ni)
c      call dblepr('eycbi(,2)',-1,eycbi(1,2),ni)
c      call dblepr('eycb2i(,1)',-1,eycb2i(1,1),q)
c      call dblepr('eycb2i(,2)',-1,eycb2i(1,2),q)
c      return
      
c      call dblepr('bimode',-1,bimode,qb)
c      call dblepr('bhinv(,1)',-1,bhinv(1,1),qb)
c      call dblepr('bhinv(,2)',-1,bhinv(1,2),qb)
c      call dblepr('bhinv(,3)',-1,bhinv(1,3),qb)
c      call dblepr('bhinv(,4)',-1,bhinv(1,4),qb)      
c      call dblepr('bim',-1,bim,qb)
c      call intpr('nburnin',-1,nburnin,1)
c      call intpr('nsize',-1,nsize,1)
c      return
      do 105 k=1,nburnin+nsize
        call randomn(qb,r)
c        call dblepr('r',-1,r,qb)
c        return
        do 10 l=1,qb
          bic(l) = bimode(l)
          do 10 s=1,qb
            bic(l)=bic(l)+bhinv(s,l)*r(s)
 10     continue
c        call dblepr('bic',-1,bic,qb)
c        return
        
        ap = blogden(bim,qb,sigbinv,ldsb,deltai,vibeta,h0t,tid,xit,
     $       ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,yicm,ldycm,ni,
     $       zi,ldzi,sigeinv,ldse) - blogden(bic,qb,sigbinv,ldsb,deltai,
     $       vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,
     $       yicm,ldycm,ni,zi,ldzi,sigeinv,ldse)
c        call dblepr('ap',-1,ap,1)
c        return
        do 20 l=1,qb
          ap=ap + bicurv(l,l)*((bim(l)-bimode(l))**2 - 
     $          (bic(l)-bimode(l))**2)/2.0d0
          do 20 s=1,l
            ap=ap-(bim(s)-bimode(s))*bicurv(s,l)*(bim(l)-bimode(l))+
     $            (bic(s)-bimode(s))*bicurv(s,l)*(bic(l)-bimode(l))
 20     continue
c        call dblepr('ap',-1,ap,1)
c        return
        call random(1,r)
        if (r(1).lt.dexp(ap)) then
          do 30 l=1,qb
            bim(l)=bic(l)
 30       continue
        end if
c        call dblepr('bim',-1,bim,qb)
c        return
        
        
        if ((k.gt.nburnin).or.(dila.eq.1)) then
          if (k.gt.nburnin) then
            do 40 l=1,qb
              biarr(l,k-nburnin)=bim(l)
 40         continue
          end if 
          
          do 50 l=1,q
            do 50 s=1,2
              bimm(s,l) = bim((l-1)*2+s)
 50       continue
          
c          call intpr('ldwz',-1,ldwz,1)
c          return
          do 100 j=1,ni
            if ((k.gt.nburnin).or.(dil(j).eq.1)) then
              do 60 l=1,q
                wkzb(l,j) = 0.0d0
                do 60 s=1,2
                   wkzb(l,j) = wkzb(l,j) + zi(j,s)*bimm(s,l)
 60             continue
            end if
            
            if (dil(j).eq.1) then            
              do 80 l=1,q
                if (di(j,l).eq.1) then
                  cmu = wkzb(l,j)
                  do 70 s=1,q
                    if (s.ne.l) cmu = cmu + 
     $                    gmat(s,l)*(yicm(j,s)-wkzb(s,j))
 70               continue
                  call random(1,r)
                  pp = cpnorm(qic(j,l),cmu,cstd(l))*r(1)
                  yicm(j,l)=cqnorm(pp,cmu,cstd(l))
                end if
 80           continue
            end if
            
            if (k.gt.nburnin) then                      
              do 90 l=1,q
                eycbi(j,l)=eycbi(j,l)+yicm(j,l)-wkzb(l,j)
                do 90 s=1,l
                  eycb2i(s,l)=eycb2i(s,l)+(yicm(j,s)-wkzb(s,j))*
     $                  (yicm(j,l)-wkzb(l,j))
 90           continue
            end if
 100      continue
        end if
c        call dblepr('biarr(,1)',-1,biarr(1,1),qb)
c        call dblepr('biarr(,nsize)',-1,biarr(1,nsize),qb)
c        call dblepr('eycbi(,1)',-1,eycbi(1,1),qb)
c        call dblepr('eycbi(,2)',-1,eycbi(1,2),qb)
c        call dblepr('eycb2i(,1)',-1,eycb2i(1,1),q)
c        call dblepr('eycb2i(,2)',-1,eycb2i(1,2),q)
c        return
 105  continue
      
      do 110 j=1,ni
        do 110 l=1,q
          eycbi(j,l)=eycbi(j,l)/nsize
 110  continue
      
      do 120 l=1,q
        do 120 s=1,l
          eycb2i(s,l)=eycb2i(s,l)/nsize
          if (s.ne.l) eycb2i(l,s)=eycb2i(s,l)
 120  continue
      
      return
      end
      
      
C     Evaluate Monte-Carlo estimates
      subroutine evalwpb(biarr,ldba,qb,nsize,deltai,vibeta,h0t,tid,
     $     xit,ldxt,px,zit,ldzt,alphamat,ldam,q,gamma,vi,pv,zitz,ldzz,
     $     Eb,Ebb,ldeb,Evibpoi,gamhexpx,Ehexp,Ehexpxz,gamexpxm,ldgx,
     $     expxzm,ldem,gamhexp2,ldgh,hexpxzc,ldez,hexpxz2,ldez2,      
     $     fderivar,ldfv,batchs,szbatch,nbatchs,tmpvec1,tmpvec2,
     $     Ehvibpoi,Evibpob,ldvb,Ehvibpob,ldvb2,hzitz,ldhz,xitalp,ldxa)
      
      parameter (MAXQ=4,MAXPX=50,MAXPV=20)
      integer ldba,qb,nsize,deltai,tid,ldxt,px,ldzt,ldam,q,pv,ldzz,ldeb,
     $     ldgx,ldem,ldgh,ldez,ldez2,ldfv,batchs(*),szbatch(*),nbatchs,
     $     ldvb,ldvb2,ldhz,ldxa
      double precision biarr(ldba,*),vibeta,h0t(*),xit(ldxt,*),
     $     zit(ldzt,*),alphamat(ldam,*),gamma(*),vi(*),zitz(ldzz,*),
     $     Eb(*),Ebb(ldeb,*),Evibpoi(*),gamhexpx(*),Ehexp,Ehexpxz(*),
     $     gamexpxm(ldgx,*),expxzm(ldem,*),gamhexp2(ldgh,*),
     $     hexpxzc(ldez,*),hexpxz2(ldez2,*),fderivar(ldfv,*),tmpvec1(*),
     $     tmpvec2(*),Ehvibpoi(*),Evibpob(ldvb,*),Ehvibpob(ldvb2,*),
     $     hzitz(ldhz,*),xitalp(ldxa,*)
      
      integer j,k,l,s,t,idx1,ldbm,ldhz2,pvpq,bchcnt,offset,bchsize
      double precision tmpsum1,tmpsum2,dexp,bim(2,MAXQ),tmpvec3(MAXQ),
     $     tmpvec4(4),tmpvec5(MAXQ*MAXPX),tmpvec6(MAXQ),tmpvec7(MAXPX),
     $     hvibpobz(MAXQ,MAXQ),tmpmeani(MAXPV+MAXQ),tmpmeano(MAXPV+MAXQ)
      
      ldbm = MAXQ
      ldhz2 = MAXQ
      pvpq = pv + q
      
      bchsize = nsize/nbatchs
      do 10 j=1,nbatchs
        szbatch(j)=bchsize
        offset = (j-1)*bchsize
        do 10 k=1,bchsize
          batchs(offset+k) = j
 10   continue
      if (nsize.gt.(nbatchs*bchsize)) then
        offset = nbatchs*bchsize
        do 20 k=1,nsize-nbatchs*bchsize
          batchs(offset+k) = nbatchs
 20     continue
        szbatch(nbatchs)=szbatch(nbatchs)+nsize-nbatchs*bchsize
      end if
      
      call matzero(ldxa,tid,q,xitalp)
      do 40 k=1,tid
        tmpvec1(k) = vibeta
        do 30 l=1,q
          do 30 s=1,px
            tmpvec1(k) = tmpvec1(k) + xit(k,s)*alphamat(s,l)*gamma(l)
            xitalp(k,l) = xitalp(k,l) + xit(k,s)*alphamat(s,l)
 30     continue
        tmpvec1(k) = dexp(tmpvec1(k))
        tmpvec2(k) = h0t(k)*tmpvec1(k)
        
        do 40 l=1,4
          hzitz(k,l) = tmpvec2(k)*zitz(k,l)
 40   continue
c      call dblepr('xitalp(,1)',-1,xitalp(1,1),tid)
c      call dblepr('xitalp(,2)',-1,xitalp(1,2),tid)
c      call dblepr('tmpvec1',-1,tmpvec1,tid)
c      call dblepr('tmpvec2',-1,tmpvec2,tid)
c      call dblepr('hzitz(,1)',-1,hzitz(1,1),tid)
c      call dblepr('hzitz(,2)',-1,hzitz(1,2),tid)
c      call dblepr('hzitz(,3)',-1,hzitz(1,3),tid)
c      call dblepr('hzitz(,4)',-1,hzitz(1,4),tid)
c      return
      
      call veczero(qb,Eb)
      call matzero(ldeb,qb,qb,Ebb)
      call veczero(tid,Evibpoi)
      call matzero(ldvb,tid,qb,Evibpob)
      call matzero(ldhz2,q,q,hvibpobz)
      call veczero(pvpq,tmpmeani)
      call veczero(pvpq,tmpmeano)
      call matzero(ldfv,pvpq,pvpq,fderivar)
      
      bchcnt = 1
      do 190 k=1,nsize        
        if (batchs(k).ne.bchcnt) then
          do 50 l=1,pvpq
            tmpmeani(l) = tmpmeani(l)/szbatch(bchcnt)
            tmpmeano(l) = tmpmeano(l) + tmpmeani(l)
            do 50 s=1,l
              fderivar(s,l) = fderivar(s,l) + tmpmeani(l)*tmpmeani(s)
 50       continue
          
          bchcnt = batchs(k)
          call veczero(pvpq,tmpmeani)
        end if
        
        do 60 l=1,qb
          Eb(l) = Eb(l) + biarr(l,k)
          do 60 s=1,l
            Ebb(s,l) = Ebb(s,l) + biarr(s,k)*biarr(l,k)
 60     continue
        
        do 70 l=1,q
          do 70 s=1,2
            bim(s,l) = biarr((l-1)*2+s,k)
 70     continue
        
        tmpsum2 = 0.0d0
        call veczero(4,tmpvec4)
        do 140 l=1,tid
          tmpsum1 = 0.0d0
          do 90 s=1,q
            tmpvec3(s) = 0.0d0
            do 80 t=1,2
              tmpvec3(s) = tmpvec3(s) + zit(l,t)*bim(t,s)
 80         continue
            tmpsum1 = tmpsum1 + tmpvec3(s)*gamma(s)
 90       continue
          Evibpoi(l) = Evibpoi(l) + dexp(tmpsum1)
          
          do 100 s=1,qb
            Evibpob(l,s) = Evibpob(l,s) + dexp(tmpsum1)*biarr(s,k)
 100      continue
          
          do 110 s=1,4
            tmpvec4(s) = tmpvec4(s) + dexp(tmpsum1)*hzitz(l,s)
 110      continue
          
          tmpsum2 = tmpsum2 + tmpvec2(l)*dexp(tmpsum1)
          
          if (deltai.eq.0) then
            do 120 s=1,q
              tmpmeani(pv+s) = tmpmeani(pv+s) + 
     $              tmpvec2(l)*dexp(tmpsum1)*(xitalp(l,s)+tmpvec3(s))
 120        continue
          else
            do 130 s=1,q
              tmpmeani(pv+s) = tmpmeani(pv+s) - 
     $              tmpvec2(l)*dexp(tmpsum1)*(xitalp(l,s)+tmpvec3(s))
 130        continue            
          end if
 140    continue
c        call dblepr('bim(,1)',-1,bim(1,1),2)
c        call dblepr('bim(,q)',-1,bim(1,q),2)
c        call dblepr('tmpvec4',-1,tmpvec4,4)
c        return
        
        do 150 l=1,q
          do 150 s=1,l
            do 150 t=1,2
              do 150 j=1,2
                hvibpobz(s,l) = hvibpobz(s,l) + 
     $               bim(j,s)*tmpvec4((t-1)*2+j)*bim(t,l)
 150    continue
        
        if (deltai.eq.0) then
          do 160 l=1,pv
            tmpmeani(l) = tmpmeani(l) + tmpsum2*vi(l)
 160      continue
        else
          do 170 l=1,pv
            tmpmeani(l) = tmpmeani(l) - tmpsum2*vi(l)
 170      continue
          
          do 180 l=1,q
            do 180 s=1,2
              tmpmeani(pv+l) = tmpmeani(pv+l) + zit(tid+1,s)*bim(s,l)
 180      continue
        end if        
 190  continue
c      call dblepr('eb',-1,eb,qb)      
c      call dblepr('ebb(,1)',-1,ebb(1,1),qb)      
c      call dblepr('ebb(,4)',-1,ebb(1,4),qb)
c      call dblepr('tmpmeani',-1,tmpmeani,pvpq)
c      call dblepr('tmpmeano',-1,tmpmeano,pvpq)      
c      call dblepr('fderivar(,1)',-1,fderivar(1,1),pvpq)      
c      call dblepr('fderivar(,3)',-1,fderivar(1,3),pvpq)      
c      call dblepr('Evibpoi',-1,Evibpoi,tid)
c      call dblepr('Evibpob(,1)',-1,Evibpob(1,1),tid)      
c      call dblepr('Evibpob(,4)',-1,Evibpob(1,4),tid)      
c      call dblepr('hvibpobz(,1)',-1,hvibpobz(1,1),q)      
c      call dblepr('hvibpobz(,2)',-1,hvibpobz(1,2),q)      
c      return
      
      
c     use the computed quantities from the above simulations to compute the components needed
      do 210 k=1,qb
        Eb(k) = Eb(k) / nsize
        do 210 l=1,k
          Ebb(l,k) = Ebb(l,k) / nsize
          if (l.ne.k) Ebb(k,l) = Ebb(l,k)
 210  continue
c      call dblepr('eb',-1,eb,qb)      
c      call dblepr('ebb(,1)',-1,ebb(1,1),qb)      
c      call dblepr('ebb(,4)',-1,ebb(1,4),qb)
c      return
      
      
      Ehexp = 0.0d0
      call veczero(q,Ehexpxz)
      call matzero(ldgh,px*q,px*q,gamhexp2)
      call matzero(ldez,q,px*q,hexpxzc)
      call matzero(ldez2,q,q,hexpxz2)
      call veczero(px,tmpvec7)
      do 290 k=1,tid
        Evibpoi(k) = tmpvec1(k)*Evibpoi(k)/nsize
        Ehvibpoi(k) = h0t(k)*Evibpoi(k)
        Ehexp = Ehexp + Ehvibpoi(k)
        do 220 l=1,qb
          Evibpob(k,l) = tmpvec1(k)*Evibpob(k,l)/nsize
          Ehvibpob(k,l) = h0t(k)*Evibpob(k,l)
 220    continue
        
        do 270 l=1,q
          Ehexpxz(l) = Ehexpxz(l) + Ehvibpoi(k)*xitalp(k,l)
          expxzm(k,l) = Evibpoi(k)*xitalp(k,l)
          
          tmpsum1 = Ehvibpoi(k)*xitalp(k,l)
          do 230 s=1,2
            tmpsum1 = tmpsum1 + zit(k,s)*Ehvibpob(k,(l-1)*2+s)
 230      continue
          
          do 240 s=1,q
            do 240 t=1,px 
              idx1 = (s-1)*px + t
              hexpxzc(l,idx1)=hexpxzc(l,idx1)+tmpsum1*gamma(s)*xit(k,t)              
              if (s.eq.l) hexpxzc(l,idx1)=hexpxzc(l,idx1)+
     $             Ehvibpoi(k)*xit(k,t)
 240      continue
          
          do 250 s=1,px
            idx1 = (l-1)*px + s
            tmpvec5(idx1) = gamma(l)*xit(k,s)
            gamexpxm(k,idx1) = gamma(l)*Evibpoi(k)*xit(k,s)
 250      continue
          
          tmpvec6(l) = 0.0d0
          do 260 s=1,2
            idx1 = (l-1)*2 + s
            Ehexpxz(l) = Ehexpxz(l) + zit(k,s)*Ehvibpob(k,idx1)
            expxzm(k,l) = expxzm(k,l) + zit(k,s)*Evibpob(k,idx1)
            tmpvec6(l) = tmpvec6(l) + zit(k,s)*Ehvibpob(k,idx1)
 260      continue
          
          do 270 s=1,l
            hexpxz2(s,l)=hexpxz2(s,l) +
     $            xitalp(k,s)*Ehvibpoi(k)*xitalp(k,l) +
     $            xitalp(k,s)*tmpvec6(l) +
     $            xitalp(k,l)*tmpvec6(s)
 270    continue      
        
        do 280 l=1,px*q
          do 280 s=1,l
            gamhexp2(s,l) = gamhexp2(s,l) + 
     $            Ehvibpoi(k)*tmpvec5(s)*tmpvec5(l)
 280    continue
        
        do 290 l=1,px
          tmpvec7(l) = tmpvec7(l) + Ehvibpoi(k)*xit(k,l)
 290  continue
c      call dblepr('Evibpoi',-1,Evibpoi,tid)
c      call dblepr('Ehvibpoi',-1,Ehvibpoi,tid)
c      call dblepr('Evibpob(,1)',-1,Evibpob(1,1),tid)      
c      call dblepr('Evibpob(,4)',-1,Evibpob(1,4),tid)      
c      call dblepr('Ehvibpob(,1)',-1,Ehvibpob(1,1),tid)      
c      call dblepr('Ehvibpob(,4)',-1,Ehvibpob(1,4),tid)      
c      return
      
      
      do 310 k=1,q
        do 300 l=1,k
          hvibpobz(l,k) =  hvibpobz(l,k)/nsize
          hexpxz2(l,k) =  hexpxz2(l,k) + hvibpobz(l,k)
          if (l.ne.k) then 
            hvibpobz(k,l) = hvibpobz(l,k)
            hexpxz2(k,l) = hexpxz2(l,k)
          end if
 300    continue
        
        do 310 l=1,px
          gamhexpx((k-1)*px+l) = gamma(k)*tmpvec7(l)
 310  continue
      
      do 320 k=1,px*q
        do 320 l=1,k
          if (l.ne.k) gamhexp2(k,l)=gamhexp2(l,k)
 320  continue
      
      do 330 k=1,pvpq
        tmpmeani(k) = tmpmeani(k)/szbatch(bchcnt)
        tmpmeano(k) = tmpmeano(k) + tmpmeani(k)
        do 330 l=1,k
          fderivar(l,k) = (fderivar(l,k) + tmpmeani(l)*tmpmeani(k) - 
     $          tmpmeano(l)*tmpmeano(k)/nbatchs) / (nbatchs*(nbatchs-1))
          if (l.ne.k) fderivar(k,l) = fderivar(l,k)        
 330  continue      
c      call dblepr('Ehexp',-1,Ehexp,1)
c      call dblepr('gamhexpx',-1,gamhexpx,px*q)     
c      call dblepr('Ehexpxz',-1,Ehexpxz,q)     
c      call dblepr('gamexpxm(,1)',-1,gamexpxm(1,1),tid)     
c      call dblepr('gamexpxm(,6)',-1,gamexpxm(1,6),tid)     
c      call dblepr('expxzm(,1)',-1,expxzm(1,1),tid)     
c      call dblepr('expxzm(,2)',-1,expxzm(1,2),tid)     
c      call dblepr('gamhexp2(,1)',-1,gamhexp2(1,1),px*q)     
c      call dblepr('gamhexp2(,6)',-1,gamhexp2(1,6),px*q)     
c      call dblepr('hexpxzc(,1)',-1,hexpxzc(1,1),q)     
c      call dblepr('hexpxzc(,6)',-1,hexpxzc(1,6),q)     
c      call dblepr('fderivar(,1)',-1,fderivar(1,1),pvpq)     
c      call dblepr('fderivar(,3)',-1,fderivar(1,3),pvpq)     
c      call dblepr('hvibpobz(,1)',-1,hvibpobz(1,1),q)     
c      call dblepr('hvibpobz(,2)',-1,hvibpobz(1,2),q)     
c      call dblepr('hexpxz2(,1)',-1,hexpxz2(1,1),q)     
c      call dblepr('hexpxz2(,2)',-1,hexpxz2(1,2),q)     
      
      return
      end
            
      
C     Sum over the parts for iteration step A
      subroutine sumAll1(Eb,qb,Ebb,ldeb,Evibpoi,tid,gamhexpx,px,q,Ehexp,
     $     Ehexpxz,gamexpxm,ldgx,expxzm,ldem,gamhexp2,ldgh,hexpxzc,ldez,
     $     hexpxz2,ldez2,fderivar,ldfv,pv,deltai,h0t,vi,eycbi,ldey,
     $     eycb2i,ldey2,xi,ldxi,ni,xit,ldxt,zit,ldzt,alphamat,ldam,
     $     gamma,sigeinv,ldse,krwe2,ldke,dqe,lddq,tqeqe,ldte,sumebb,
     $     ldsb,sumdf,ne1,sumd2f,ldsd,summct,ldsmc)
      
      parameter (MAXQ=4,MAXPX=50)
      integer qb,ldeb,tid,px,q,ldgx,ldem,ldgh,ldez,ldez2,ldfv,pv,deltai,
     $     ldey,ldey2,ldxi,ni,ldxt,ldzt,ldam,ldse,ldke,lddq,ldte,ldsb,
     $     ne1,ldsd,ldsmc
      double precision Eb(*),Ebb(ldeb,*),Evibpoi(*),gamhexpx(*),Ehexp,
     $     Ehexpxz(*),gamexpxm(ldgx,*),expxzm(ldem,*),gamhexp2(ldgh,*),
     $     hexpxzc(ldez,*),hexpxz2(ldez2,*),fderivar(ldfv,*),h0t(*),
     $     vi(*),eycbi(ldey,*),eycb2i(ldey2,*),xi(ldxi,*),xit(ldxt,*),
     $     zit(ldzt,*),alphamat(ldam,*),gamma(*),sigeinv(ldse,*),
     $     krwe2(ldke,*),dqe(lddq,*),tqeqe(ldte,*),sumebb(ldsb,*),
     $     sumdf(*),sumd2f(ldsd,*),summct(ldsmc,*)
      
      integer j,k,l,s,t,idx1,idx2,idx3,pxtq,pvpq,nfixpar,qeu,qea,ldtm,
     $     ldtm3,offset
      double precision tmpvec1(MAXQ),tmpvec2(MAXQ**2),
     $     tmpmat1(MAXQ,MAXPX*MAXQ),tmpmat2(MAXQ,MAXPX*MAXQ),
     $     tmpmat3(MAXQ**2,MAXPX*MAXQ),tmpmat4(MAXQ,MAXQ),
     $     tmpmat5(MAXQ**2,MAXQ**2),tmpmat6(MAXQ**2,MAXQ**2)
      
c      call dblepr('gamhexp2(,1)',-1,gamhexp2(1,1),6)     
c      call dblepr('gamhexp2(,6)',-1,gamhexp2(1,6),6)     
c      return
      
      pxtq = px*q
      pvpq = pv + q
      nfixpar = pxtq + pvpq
      qeu = q*(q+1)/2
      qea = q**2
      ldtm  = MAXQ
c      ldtm2 = MAXQ
      ldtm3 = MAXQ**2
c      ldtm4 = MAXQ
c      ldtm5 = MAXQ**2
c      ldtm6 = MAXQ**2
      
c      call dblepr('sumEbb(,1)',-1,sumEbb(1,1),qb)     
      do 10 k=1,qb
        do 10 l=1,k
          sumEbb(l,k)=sumEbb(l,k)+Ebb(l,k)
c          if (l.ne.k) sumEbb(k,l)=sumEbb(l,k)
 10   continue
c      call dblepr('Ebb(,1)',-1,Ebb(1,1),qb)     
c      call dblepr('sumEbb(,1)',-1,sumEbb(1,1),qb)     
c      call dblepr('sumEbb(,4)',-1,sumEbb(1,4),qb)     
c      return
      
c      do 20 k=1,q
c        do 20 l=1,k
c          sumEee(l,k)=sumEee(l,k)+eycb2i(l,k)
cc          if (l.ne.k) sumEee(k,l)=sumEee(l,k)          
c 20   continue
      
      
      if (tid.gt.0) then
        do 50 k=1,tid
          sumdf(k) = sumdf(k) - Evibpoi(k)
          
          do 30 l=1,pxtq
            idx1 = ne1 + l
            sumd2f(k,idx1) = sumd2f(k,idx1) - gamexpxm(k,l)
 30       continue
          
          do 40 l=1,pv
            idx1 = ne1 + pxtq + l
            sumd2f(k,idx1) = sumd2f(k,idx1) - Evibpoi(k)*vi(l)
 40       continue   
          
          do 50 l=1,q
            idx1 = ne1 + pxtq + pv + l
            sumd2f(k,idx1) = sumd2f(k,idx1) - expxzm(k,l)
 50     continue            
        
        if (deltai.eq.1) then
          sumdf(tid) = sumdf(tid) + 1/h0t(tid)
          sumd2f(tid,tid) = sumd2f(tid,tid) - 1.0d0/h0t(tid)**2
        end if        
      end if
c      call dblepr('sumDf',-1,sumDf,ne1)     
c      call dblepr('sumD2f(,1)',-1,sumD2f(1,1),ne1)     
c      call dblepr('sumD2f(,tid)',-1,sumD2f(1,tid),ne1)     
c      return
      
      
      call matzero(ldtm,q,pxtq,tmpmat1)
      call matzero(ldtm3,qea,pxtq,tmpmat3)
      do 90 j=1,ni
        do 70 k=1,q
          do 60 l=1,px
            tmpmat1(k,(k-1)*px+l)=xi(j,l)          
 60       continue           
          
          tmpvec1(k) = 0.0d0
          do 70 l=1,q
            tmpvec1(k)=tmpvec1(k)+sigeinv(k,l)*eycbi(j,l)
 70     continue
        
        do 80 k=1,q
          do 80 l=1,px
            idx1 = (k-1)*px + l
            do 80 s=1,q
              tmpmat2(s,idx1) = 0.0d0
              do 80 t=1,q
              tmpmat2(s,idx1)=tmpmat2(s,idx1)+
     $                sigeinv(s,t)*tmpmat1(t,idx1)
 80     continue
        
        do 90 k=1,pxtq
          idx1 = ne1 + k
          do 90 l=1,q
            sumdf(idx1) = sumdf(idx1) + tmpvec1(l)*tmpmat1(l,k)
            do 90 s=1,q               
              idx2 = (l-1)*q + s
              tmpmat3(idx2,k)=tmpmat3(idx2,k)+tmpmat2(l,k)*tmpvec1(s)+
     $             tmpvec1(l)*tmpmat2(s,k)
              
              do 90 t=1,k
                idx3 = ne1 + t
                sumd2f(idx3,idx1) = sumd2f(idx3,idx1) -
     $               tmpmat1(s,t)*sigeinv(s,l)*tmpmat1(l,k)
 90   continue
c      call dblepr('sumD2f(57:62,57)',-1,sumD2f(57,57),6)     
c      call dblepr('sumD2f(57:62,62)',-1,sumD2f(57,62),6)     
c      return
      
      
      do 100 k=1,q
        do 100 l=1,px
          idx1 = (k-1)*px + l
          idx2 = ne1 + idx1
          sumdf(idx2)=sumdf(idx2)+deltai*gamma(k)*xit(tid+1,l)-
     $         gamhexpx(idx1)
 100  continue
      
      do 110 k=1,pv
        idx1 = ne1 + pxtq + k
        sumdf(idx1) = sumdf(idx1) + (deltai-Ehexp)*vi(k)
 110  continue
      
      if (deltai.eq.1) then
        do 130 k=1,q
          idx1 = ne1 + pxtq + pv + k
          do 120 l=1,px
            sumdf(idx1)=sumdf(idx1)+xit(tid+1,l)*alphamat(l,k)
            
            idx2 = ne1 + (k-1)*px + l
            sumd2f(idx2,idx1)=sumd2f(idx2,idx1)+xit(tid+1,l)
 120      continue
          
          do 130 l=1,2
c            idx2 = (k-1)*2 + l
            sumdf(idx1)=sumdf(idx1)+zit(tid+1,l)*Eb((k-1)*2+l)
 130    continue
      end if
      
      do 140 k=1,q
        idx1 = ne1 + pxtq + pv + k
        sumdf(idx1)=sumdf(idx1)-Ehexpxz(k)
        do 140 l=1,k
          idx2 = ne1 + pxtq + pv + l
          sumd2f(idx2,idx1)=sumd2f(idx2,idx1)-hexpxz2(l,k)
 140  continue
      
      
      do 170 k=1,pxtq
        idx1 = ne1 + k
        do 150 l=1,k
          idx2 = ne1 + l
          sumd2f(idx2,idx1)=sumd2f(idx2,idx1)-gamhexp2(l,k)
 150    continue
        
        do 160 l=1,pv
          idx2 = ne1 + pxtq + l
          sumd2f(idx1,idx2)=sumd2f(idx1,idx2)-gamhexpx(k)*vi(l)
 160    continue
        
        do 170 l=1,q
          idx2 = ne1 + pxtq + pv + l
          sumd2f(idx1,idx2)=sumd2f(idx1,idx2)-hexpxzc(l,k)
 170  continue
      
      
      do 190 k=1,pv
        idx1 = ne1 + pxtq + k
        do 180 l=1,k
          idx2 = ne1 + pxtq + l
          sumd2f(idx2,idx1)=sumd2f(idx2,idx1)-Ehexp*vi(l)*vi(k)
 180    continue
        
        do 190 l=1,q
          idx2 = ne1 + pxtq + pv + l
          sumd2f(idx1,idx2)=sumd2f(idx1,idx2)-vi(k)*Ehexpxz(l)
 190  continue
      
      
      do 200 k=1,q
        do 200 l=1,q
          idx1 = (k-1)*q + l
          tmpvec2(idx1) = ni*sigeinv(l,k)
          tmpmat4(l,k) = 0.0d0
          do 200 s=1,q
            do 200 t=1,q
              idx2 = (s-1)*q + t
              tmpvec2(idx1)=tmpvec2(idx1)-eycb2i(t,s)*krwe2(idx2,idx1)
              tmpmat4(l,k)=tmpmat4(l,k)+
     $             sigeinv(l,t)*eycb2i(t,s)*sigeinv(s,k)
 200  continue
      
      do 210 k=1,q
        do 210 l=1,q
          idx1 = (k-1)*q + l
          do 210 s=1,q
            do 210 t=1,q
              idx2 = (s-1)*q + t
              tmpmat5(idx2,idx1) = tmpmat4(s,k)*sigeinv(t,l)+
     $             sigeinv(s,k)*tmpmat4(t,l)
 210  continue
      
      do 220 k=1,q
        do 220 l=1,q
          idx1 = (k-1)*q + l
          do 220 s=1,q
            do 220 t=1,q
              idx2 = (s-1)*q + t
              tmpmat6(idx2,idx1) = 0.0d0
              do 220 j=1,qea
              tmpmat6(idx2,idx1)=tmpmat6(idx2,idx1) + 
     $                ni*tqeqe(idx2,j)*krwe2(j,idx1) -
     $                tmpmat5(idx2,j)*tqeqe(j,idx1)
 220  continue
      
      offset = ne1 + nfixpar
      do 250 k=1,qeu
        idx1 = offset + k
        do 230 l=1,q
          do 230 s=1,q
            idx2 = (l-1)*q + s
            sumdf(idx1)=sumdf(idx1)-tmpvec2(idx2)*dqe(idx2,k)/2.0d0
 230    continue
        
        do 240 l=1,q
          do 240 s=1,px
            idx2 = (l-1)*px + s
            idx3 = ne1 + idx2
            do 240 t=1,qea
              sumd2f(idx3,idx1)=sumd2f(idx3,idx1) -
     $              tmpmat3(t,idx2)*dqe(t,k)/2.0d0
 240    continue
        
        do 250 l=1,k
          idx2 = offset + l
          do 250 s=1,qea
            do 250 t=1,qea
              sumd2f(idx2,idx1)=sumd2f(idx2,idx1) + 
     $              dqe(t,l)*tmpmat6(t,s)*dqe(s,k)/2.0d0
 250  continue
      
      
      do 260 k=1,pvpq
        do 260 l=1,k
          summct(l,k) = summct(l,k) + fderivar(l,k)
 260  continue
c      call dblepr('sumDf',-1,sumDf,68)     
c      call dblepr('sumD2f(,1)',-1,sumD2f(1,1),68)
c      call dblepr('sumD2f(,tid)',-1,sumD2f(1,tid),68)
c      call dblepr('sumD2f(,ne1)',-1,sumD2f(1,ne1),68)
      
c      call dblepr('sumD2f(,ne1+1)',-1,sumD2f(1,ne1+1),68)
c      call dblepr('sumD2f(,ne1+6)',-1,sumD2f(1,ne1+6),68)
      
c      call dblepr('sumD2f(,ne1+7)',-1,sumD2f(1,ne1+7),68)
c      call dblepr('sumD2f(,ne1+8)',-1,sumD2f(1,ne1+8),68)
c      call dblepr('sumD2f(,ne1+9)',-1,sumD2f(1,ne1+9),68)
c      call dblepr('sumD2f(,ne1+10)',-1,sumD2f(1,ne1+10),68)
c      call dblepr('sumD2f(,ne1+11)',-1,sumD2f(1,ne1+11),68)
c      call dblepr('sumD2f(,ne1+12)',-1,sumD2f(1,ne1+12),68)
      
      return
      end
      
      
c     Get the Jacobian and derivative of Jacobian on log(h0t)
      subroutine getjcdlh(h0t,ne1,lh0t,Jc,ldj,dJ,lddj)
      
      integer ne1,ldj,lddj
      double precision h0t(*),lh0t(*),Jc(ldj,*),dJ(lddj,*)
      
      integer j,dlog
      
      call matzero(ldj,ne1,ne1,Jc)
      call matzero(lddj,ne1**2,ne1,dJ)
      do 10 j=1,ne1
        lh0t(j) = dlog(h0t(j))
        Jc(j,j) = h0t(j)
        dJ((j-1)*ne1+j,j) = h0t(j)
 10   continue
      
      return
      end
      
      
c     Get the Jacobian and derivative of Jacobian on log(chol) transformed
c     parameters      
      subroutine getjcdlc(vmat,ldvm,qe,vmat2,ldvm2,Jc,ldj,dJ,lddj)
      
      parameter (MAXQ=4)
      integer ldvm,qe,ldvm2,ldj,lddj
      double precision vmat(ldvm,*),vmat2(ldvm2,*),Jc(ldj,*),dJ(lddj,*)
      
      integer qeu,jpvt(MAXQ),ierr,i,j,k,idx1,dlog,dexp
      double precision tmp1(MAXQ)
      
      qeu = qe*(qe+1)/2
      
c     log(chol) transformation on the parameters
      call mtxcopy(vmat,ldvm,qe,qe,vmat2,ldvm2)
      call dchdc(vmat2,ldvm2,qe,tmp1,jpvt,0,ierr)
      
      do 10 i=1,qe
        vmat2(i,i)=dlog(vmat2(i,i))
 10   continue
      
c     get the first derivative (Jacobian)
      call matzero(ldj,qeu,qeu,Jc)
      
      do 40 j=1,qe
        do 40 i=1,j
          idx1=j*(j-1)/2+i
          if (i.lt.j) then
            do 20 k=1,i
              if (k.lt.i) then
                Jc(idx1,i*(i-1)/2+k)=vmat2(k,j)
                Jc(idx1,j*(j-1)/2+k)=vmat2(k,i)
              else
                Jc(idx1,i*(i-1)/2+i)=dexp(vmat2(i,i))*vmat2(i,j)
                Jc(idx1,j*(j-1)/2+i)=dexp(vmat2(i,i))
              end if
 20         continue
          else
            do 30 k=1,i
              if (k.lt.i) then
                Jc(idx1,i*(i-1)/2+k)=2.0d0*vmat2(k,i)
              else
                Jc(idx1,i*(i-1)/2+i)=2.0d0*dexp(2.0d0*vmat2(i,i))
              end if
 30         continue            
          end if
 40   continue
      
c      call dblepr('Jc(,1)',-1,Jc(1,1),6)
c      call dblepr('Jc(,2)',-1,Jc(1,2),6)
c      call dblepr('Jc(,3)',-1,Jc(1,3),6)
c      call dblepr('Jc(,4)',-1,Jc(1,4),6)
c      call dblepr('Jc(,5)',-1,Jc(1,5),6)
c      call dblepr('Jc(,6)',-1,Jc(1,6),6)
      
c     get the second derivative
      call matzero(lddj,qeu**2,qeu,dJ)
      
      do 70 j=1,qe
        do 70 i=1,j
          idx1=j*(j-1)/2+i
          if (i.lt.j) then
            do 50 k=1,i
              if (k.lt.i) then
                dJ((i*(i-1)/2+k-1)*qeu+idx1,j*(j-1)/2+k)=1.0d0
                dJ((j*(j-1)/2+k-1)*qeu+idx1,i*(i-1)/2+k)=1.0d0
              else
                dJ((i*(i-1)/2+i-1)*qeu+idx1,i*(i-1)/2+i)=
     $                dexp(vmat2(i,i))*vmat2(i,j)
                dJ((i*(i-1)/2+i-1)*qeu+idx1,j*(j-1)/2+i)=
     $                dexp(vmat2(i,i))
                dJ((j*(j-1)/2+i-1)*qeu+idx1,i*(i-1)/2+i)=
     $                dexp(vmat2(i,i))
              end if
 50         continue
          else
            do 60 k=1,i
              if (k.lt.i) then
                dJ((i*(i-1)/2+k-1)*qeu+idx1,i*(i-1)/2+k)=2.0d0
              else
                dJ((i*(i-1)/2+i-1)*qeu+idx1,i*(i-1)/2+i)=
     $                4.0d0*dexp(2.0d0*vmat2(i,i))
              end if
 60         continue
          end if
 70   continue
      
c      call dblepr('dJ(,1)',-1,dJ(1,1),36)
c      call dblepr('dJ(,2)',-1,dJ(1,2),36)
c      call dblepr('dJ(,3)',-1,dJ(1,3),36)
c      call dblepr('dJ(,4)',-1,dJ(1,4),36)
c      call dblepr('dJ(,5)',-1,dJ(1,5),36)
c      call dblepr('dJ(,6)',-1,dJ(1,6),36)
      
      return
      end

      
C     Update the parameters
      subroutine updatep(sumEbb,ldsb,qb,m,sumdf,ne1,px,q,pv,sumd2f,ldsd,
     $     summct,ldsmc,h0t,betas,sige,ldse,P2,ldp2,newh0t,newbetas,
     $     newsigb,ldsb2,newsige,ldse2,mcbetat,ldmt,mcvar,J2,ldj2,dJ2,
     $     lddj,natrparm,nnatparm,natrgrad,natrhess,ldnh,natrhinv,ldni,
     $     dJwk,lddw)
      
      parameter (MAXQ=4,MAXQEU=MAXQ*(MAXQ+1)/2)
      integer ldsb,qb,m,ne1,px,q,pv,ldsd,ldsmc,ldse,ldp2,ldsb2,ldse2,
     $     ldmt,ldj2,lddj,ldnh,ldni,lddw
      double precision sumEbb(ldsb,*),sumdf(*),sumd2f(ldsd,*),
     $     summct(ldsmc,*),h0t(*),betas(*),sige(ldse,*),P2(ldp2,*),
     $     newh0t(*),newbetas(*),newsigb(ldsb2,*),newsige(ldse2,*),
     $     mcbetat(ldmt,*),mcvar(*),J2(ldj2,*),dJ2(lddj,*),natrparm(*),
     $     nnatparm(*),natrgrad(*),natrhess(ldnh,*),natrhinv(ldni,*),
     $     dJwk(lddw,*)
      
      integer k,l,s,t,ierr,nfixpar,qeu,ntotpar2,pvpq,ldtm,ldjt,offset,
     $     idx1
      double precision det(2),tmpmat1(MAXQ,MAXQ),Jctmp(MAXQEU,MAXQEU),
     $     dexp    
      
      nfixpar = px*q + pv + q
      qeu = q*(q+1)/2
      ntotpar2 = ne1 + nfixpar + qeu
      pvpq = pv + q
      ldtm = MAXQ
      ldjt = MAXQEU
      
      do 10 k=1,qb
        do 10 l=1,k
          newsigb(l,k)=sumEbb(l,k)/m
          if (l.ne.k) then
            sumEbb(k,l)=sumEbb(l,k)
            newsigb(k,l)=newsigb(l,k)
          end if
 10   continue
c      call dblepr('newsigb(,1)',-1,newsigb(1,1),qb)
c      call dblepr('newsigb(,4)',-1,newsigb(1,4),qb)
c      return
      
      
      do 20 k=1,ntotpar2
        do 20 l=1,k
          if (l.ne.k) sumd2f(k,l)=sumd2f(l,k)
 20   continue
      
      do 30 k=1,pvpq
        do 30 l=1,k
          if (l.ne.k) summct(k,l)=summct(l,k)
 30   continue
      
      
      call matzero(ldj2,ntotpar2,ntotpar2,J2)
      call getjcdlh(h0t,ne1,newh0t,J2,ldj2,dJwk,lddw)
      
      call matzero(lddj,ntotpar2**2,ntotpar2,dJ2)
      do 40 k=1,ne1
        do 40 l=1,ne1
          do 40 s=1,ne1
            dJ2((l-1)*ntotpar2+s,k) = dJwk((l-1)*ne1+s,k)
 40   continue
      
      
      do 50 k=1,nfixpar
        idx1 = ne1 + k
        J2(idx1,idx1) = 1.0d0
 50   continue
      
      
      call getjcdlc(sige,ldse,q,tmpmat1,ldtm,Jctmp,ldjt,dJwk,lddw)
      
      offset = ne1 + nfixpar
      do 60 k=1,qeu
        do 60 l=1,qeu
          J2(offset+l,offset+k) = Jctmp(l,k)
 60   continue
      
      do 70 k=1,qeu
        do 70 l=1,qeu
          do 70 s=1,qeu
            dJ2((offset+l-1)*ntotpar2+offset+s,offset+k) = 
     $            dJwk((l-1)*qeu+s,k)
 70   continue
      
      
      call veccopy(newh0t,ne1,natrparm)
      call veccopy(betas,nfixpar,natrparm(ne1+1))
      do 80 k=1,q
        do 80 l=1,k
          idx1 = offset + q*(l-1) - (l-1)*(l-2)/2 + (k-l+1)
          natrparm(idx1) = tmpmat1(l,k)
 80   continue
c      call dblepr('natrparm',-1,natrparm,ntotpar2)
c      return
      
      call matzero(ldni,ntotpar2,ntotpar2,natrhinv)
      do 90 k=1,ntotpar2
        natrgrad(k) = 0.0d0
        do 90 l=1,ntotpar2
          do 90 s=1,ntotpar2
            natrgrad(k) = natrgrad(k) + sumdf(s)*J2(s,l)*P2(l,k)
            
            idx1 = (l-1)*ntotpar2 + s
            do 90 t=1,k
              natrhinv(t,k)=natrhinv(t,k)+J2(s,t)*sumd2f(s,l)*J2(l,k)
              
              if (t.eq.l) natrhinv(t,k)=natrhinv(t,k)+
     $             sumdf(s)*dJ2(idx1,k)
 90   continue
c      call dblepr('natrgrad',-1,natrgrad,ntotpar2)
      
      do 100 k=1,ntotpar2
        do 100 l=1,k
          if (l.ne.k) natrhinv(k,l) = natrhinv(l,k)
 100  continue
      
      do 120 k=1,ntotpar2
        do 120 l=1,k
          natrhess(l,k) = 0.0d0
          do 110 s=1,ntotpar2
            do 110 t=1,ntotpar2
              natrhess(l,k) = natrhess(l,k)+P2(t,l)*natrhinv(t,s)*
     $              P2(s,k)
 110      continue
          if (l.ne.k) natrhess(k,l) = natrhess(l,k)
 120  continue
c      call dblepr('natrhess(,1)',-1,natrhess(1,1),ntotpar2)
c      call dblepr('natrhess(,68)',-1,natrhess(1,68),ntotpar2)
c      return
      
      call mtxcopy(natrhess,ldnh,ntotpar2,ntotpar2,natrhinv,ldni)
      call invert(ldni,ntotpar2,natrhinv,det,ierr)
c      call dblepr('natrhinv(,1)',-1,natrhinv(1,1),ntotpar2)
c      call dblepr('natrhinv(,68)',-1,natrhinv(1,68),ntotpar2)
c      return
      
      do 140 k=1,ntotpar2
        nnatparm(k) = natrparm(k)
        do 130 l=1,ntotpar2
          nnatparm(k) = nnatparm(k) - natrhinv(k,l)*natrgrad(l)
 130    continue
        
        if (k.le.ne1) then
          newh0t(k) = dexp(nnatparm(k))
        else if (k.le.(ne1+nfixpar)) then
          newbetas(k-ne1)=nnatparm(k)
        end if
 140  continue
c      call dblepr('nnatparm',-1,nnatparm,ntotpar2)
c      call dblepr('newh0t',-1,newh0t,ne1)
c      call dblepr('newbetas',-1,newbetas,nfixpar)
c      return
      
      idx1 = offset
      do 150 k=1,q
        do 150 l=k,q
          idx1 = idx1 + 1
          tmpmat1(l,k) = nnatparm(idx1)
          if (l.eq.k) then
            tmpmat1(l,k) = dexp(tmpmat1(l,k))
          else 
            tmpmat1(k,l) = 0.0d0
          end if
 150  continue
c      call dblepr('tmpmat1(,1)',-1,tmpmat1(1,1),q)
c      call dblepr('tmpmat1(,2)',-1,tmpmat1(1,2),q)
c      return
      
      do 170 k=1,q
        do 170 l=1,k
          newsige(l,k) = 0.0d0
          do 160 s=1,q
            newsige(l,k) = newsige(l,k) + tmpmat1(l,s)*tmpmat1(k,s)
 160      continue
          if (l.ne.k) newsige(k,l) = newsige(l,k)
 170  continue
c      call dblepr('newsige(,1)',-1,newsige(1,1),q)
c      call dblepr('newsige(,2)',-1,newsige(1,2),q)
c      return
      
      
      offset = ne1 + px*q
      do 200 k=1,pvpq
        do 190 l=1,k
          mcbetat(l,k)=0.0d0
          do 180 s=1,pvpq
            do 180 t=1,pvpq
              mcbetat(l,k)=mcbetat(l,k)+natrhinv(offset+l,offset+t)*
     $              summct(t,s)*natrhinv(offset+s,offset+k)
 180      continue
          if (l.ne.k) mcbetat(k,l)=mcbetat(l,k)
 190    continue
        mcvar(k)=mcbetat(k,k)
 200  continue
c      call dblepr('mcbetat(,1)',-1,mcbetat(1,1),pvpq)
c      call dblepr('mcbetat(,3)',-1,mcbetat(1,3),pvpq)
c      call dblepr('mcvar',-1,mcvar,pvpq)
      
      return
      end
      
      
C     Generate block diagonal matrix
c      subroutine blkdiag(x,p,dx,blockmtx,ldbm)
c
c      integer p(*),dx,ldbm
c      double precision x(*),blockmtx(ldbm,*)
c
c      integer ps2,sump,k,l,offset
c
c      ps2 = sump(p,dx)
c      call matzero(ldbm,dx,ps2,blockmtx)
c
c      offset=0
c      do 20 k=1,dx
c        do 10 l=1,p(k)
c          blockmtx(k,offset+l)=x(offset+l)
c 10     continue
c        offset=offset+p(k)
c 20   continue
c
c      return
c      end
            
      
C     Assessing percent change of the parameters and monte-carlo SE
      subroutine pctchg(h0t,ne1,betas,nfixpar,px,q,sigb,ldsb,qb,sige,
     $     ldse,newh0t,newbetas,newsigb,ldsb2,newsige,ldse2,eps,mcvar,
     $     hconv,sconv,pconv,dconv,mmcdel,dmcdel)
      
      parameter (MAXQ=4,MAXQB=2*MAXQ)
      integer ne1,nfixpar,px,q,ldsb,qb,ldse,ldsb2,ldse2
      double precision h0t(*),betas(*),sigb(ldsb,*),sige(ldse,*),
     $     newh0t(*),newbetas(*),newsigb(ldsb2,*),newsige(ldse2,*),eps,
     $     mcvar(*),hconv,sconv(*),pconv,dconv(*),mmcdel,dmcdel(*)
      
      integer k,l,jpvt(MAXQB),ierr,pxtq,pvpq,qbu,offset,idx1
      double precision tmp1(MAXQB),cholb1(MAXQB,MAXQB),
     $     cholb2(MAXQB,MAXQB),chole1(MAXQ,MAXQ),chole2(MAXQ,MAXQ),
     $     dabs,dsqrt
      
      pxtq = px*q
      pvpq = nfixpar - px*q
      qbu = qb*(qb+1)/2
c      call intpr('nfixpar',-1,nfixpar,1)
c      call intpr('pxtq',-1,pxtq,1)
c      call intpr('pvpq',-1,pvpq,1)
c      call dblepr('mmcdel',-1,mmcdel,1)
c      call intpr('pvpq',-1,pvpq,1)
c      call dblepr('dmcdel',-1,dmcdel,pvpq)
c      return
      
      hconv = 0.0d0
      do 10 k=1,ne1
        sconv(k)=dabs(newh0t(k)-h0t(k))/h0t(k)
        if (sconv(k).gt.hconv) hconv=sconv(k)
 10   continue
c      call dblepr('sconv',-1,sconv,ne1)
c      call dblepr('hconv',-1,hconv,1)
c      call intpr('pvpq',-1,pvpq,1)
c      call dblepr('dmcdel',-1,dmcdel,pvpq)
c      return
      
      pconv = 0.0d0
      do 20 k=1,nfixpar
        dconv(k)=dabs(newbetas(k)-betas(k))/(dabs(betas(k))+eps)
        if (dconv(k).gt.pconv) pconv=dconv(k)
 20   continue
      
      call mtxcopy(sigb,ldsb,qb,qb,cholb1,MAXQB)
      call mtxcopy(newsigb,ldsb2,qb,qb,cholb2,MAXQB)
      call dchdc(cholb1,MAXQB,qb,tmp1,jpvt,0,ierr)
      call dchdc(cholb2,MAXQB,qb,tmp1,jpvt,0,ierr)
      
      offset = nfixpar
      do 30 k=1,qb
        do 30 l=1,k
          idx1 = offset + qb*(l-1) - (l-1)*(l-2)/2 + (k-l+1)
          dconv(idx1)=dabs(cholb2(l,k)-cholb1(l,k))/
     $          (dabs(cholb1(l,k))+eps)
          if (dconv(idx1).gt.pconv) pconv=dconv(idx1)
 30   continue
      
      call mtxcopy(sige,ldse,q,q,chole1,MAXQ)
      call mtxcopy(newsige,ldse2,q,q,chole2,MAXQ)
      call dchdc(chole1,MAXQ,q,tmp1,jpvt,0,ierr)
      call dchdc(chole2,MAXQ,q,tmp1,jpvt,0,ierr)
      
      offset = nfixpar + qbu
      do 40 k=1,q
        do 40 l=1,k
          idx1 = offset + q*(l-1) - (l-1)*(l-2)/2 + (k-l+1)
          dconv(idx1)=dabs(chole2(l,k)-chole1(l,k))/
     $         (dabs(chole1(l,k))+eps)
          if (dconv(idx1).gt.pconv) pconv=dconv(idx1)
 40   continue
c      call dblepr('dconv',-1,dconv,nfixpar+qbu+q*(q+1)/2)
c      call dblepr('pconv',-1,pconv,1)
c      call intpr('pvpq',-1,pvpq,1)
c      call dblepr('dmcdel',-1,dmcdel,pvpq)
c      return
      
      offset = pxtq
      idx1 = offset + 1
      dmcdel(1)=dabs(newbetas(idx1)-betas(idx1))/dsqrt(mcvar(1))
      mmcdel = dmcdel(1)
      do 50 k=2,pvpq
        idx1 = offset + k
        dmcdel(k)=dabs(newbetas(idx1)-betas(idx1))/dsqrt(mcvar(k))
        if ((dmcdel(k)).lt.mmcdel) mmcdel=dmcdel(k)
 50   continue
c      call dblepr('dmcdel',-1,dmcdel,pvpq)
c      call dblepr('mmcdel',-1,mmcdel,1)
      
      return
      end
      
      
C     Copy the parameters from new set to current set
      subroutine pcopy(newh0t,ne1,newbetas,nfixpar,newsigb,ldsb,qb,
     $     newsige,ldse,q,h0t,betas,sigb,ldsb2,sige,ldse2)
      
      integer ne1,nfixpar,ldsb,qb,ldse,q,ldsb2,ldse2
      double precision newh0t(*),newbetas(*),newsigb(ldsb,*),
     $     newsige(ldse,*),h0t(*),betas(*),sigb(ldsb2,*),sige(ldse2,*)
      
      integer k,l
      
      do 10 k=1,ne1
        h0t(k)=newh0t(k)
 10   continue
      
      do 20 k=1,nfixpar
        betas(k)=newbetas(k)
 20   continue
      
      do 30 k=1,qb
        do 30 l=1,k
          sigb(l,k)=newsigb(l,k)
          if (l.ne.k) sigb(k,l)=sigb(l,k)
 30   continue
      
      do 40 k=1,q
        do 40 l=1,k
          sige(l,k)=newsige(l,k)
          if (l.ne.k) sige(k,l)=sige(l,k)
 40   continue
      
      return
      end
      
      
C     Main subroutine for fitting survival outcome with left-censored covariates
      
      subroutine survlccv(time,event,id1,m,Q,ldq,sumn,qe,D,ldd,id2,repl,
     $     V,ldv,pv,X,ldx,px,Z,ldz,samp_pd,eb,ldeb,eyc,ldeyc,h0t,
     $     alphamat,ldam,beta,gamma,sigb,ldsb,sige,ldse,maxiter,eps1,
     $     eps2,maxiter2,nburnin,nsize,maxnsize,nbatchs,hconv,sconv,
     $     pconv,dconv,mmcdel,dmcdel,iter,wkarr,iwkarr)
      
      parameter (MAXQ=4,MAXPX=50,MAXPV=20)
      parameter (MAXQEA=MAXQ**2,MAXQEU=MAXQ*(MAXQ+1)/2,MAXQB=2*MAXQ,
     $     MAXPXTQ=MAXPX*MAXQ,MAXPVPQ=MAXPV+MAXQ)
      integer event(*),id1(*),m,ldq,sumn,qe,ldd,id2(*),repl(*),ldv,pv,
     $     ldx,px,ldz,ldeb,ldeyc,ldam,ldsb,ldse,maxiter,maxiter2,
     $     nburnin,nsize,maxnsize,nbatchs,iter,iwkarr(*)
      integer D(ldd,*)
      double precision time(*),Q(ldq,*),V(ldv,*),X(ldx,*),Z(ldz,*),
     $     samp_pd,eb(ldeb,*),eyc(ldeyc,*),h0t(*),alphamat(ldam,*),
     $     beta(*),gamma(*),sigb(ldsb,*),sige(ldse,*),eps1,eps2,hconv,
     $     sconv(*),pconv,dconv(*),mmcdel,dmcdel(*),wkarr(*)
      
      integer i,j,k,l,qb,pxtq,pvpq,nfixpar,qea,qeu,ne0,ne1,ntotpar2,
     $     maxiter3,cid,idx1,idx2,offset1,offset2,offset3,maxrepl,
     $     maxni,sumtd,sumtdpm,tid,tidp1,ni,dila,dilj,replidx,replidx2,
     $     sseq(MAXQ),ierr,klen,ks,kt,int,min,max,mod
      integer ina,ietimes0,ietimes1,itd,idla,idl,it2imap,ixt,izt,iztz,
     $     ip2,idi,ibatchs,iszbatch,ivbeta,isumdf,isumd2f,
     $     iebp,ieycp,ixi,izi,iqic,ieyic,iyic0,ixit,izit,izitz,
     $     ibiarr,ieycbi,iwkzb,iEvibpoi,igamexxm,
     $     iexpxzm,igamhex2,ihexpxzc,itmpvec1,itmpvec2,
     $     iEhvibpo,iEvibpob,iEhvibpb,ihzitz,ixitalp,inewh0t,
     $     iJ2,idJ2,inatrpar,innatpar,inatgrad,inathess,inathinv,idJwk
      integer ldxt,ldzt,ldztz,ldp2,lddi,ldsd,ldebp,ldeycp,ldxi,
     $     ldzi,ldqic,ldeyic,ldyic0,ldxit,ldzit,ldzitz,ldbiarr,ldeycbi,
     $     ldwkzb,ldgx,ldem,ldgh,ldez,ldvb,ldvb2,ldhz,ldxa,ldj2,
     $     lddj2,ldnh,ldni,lddjwk
      integer ldsb2,ldsb3,ldsb4,ldse2,ldse3,ldse4,ldtq,lddq,ldpq,ldke,
     $     ldgm,ldbc,ldbh,ldey2,ldebbi,ldfv,ldmcbet,ldsmc,ldez2
      double precision eps1b,tmp1,etimes1j,tqeqe(MAXQEA,MAXQEA),
     $     dqe(MAXQEA,MAXQEU),pqe(MAXQEU,MAXQEU),sigbinv(MAXQB,MAXQB),
     $     sigeinv(MAXQ,MAXQ),krwe2(MAXQEA,MAXQEA),gmat(MAXQ,MAXQ),
     $     cstd(MAXQ),sumEbb(MAXQB,MAXQB),sumEee(MAXQ,MAXQ),
     $     sumMct(MAXPVPQ,MAXPVPQ),ebi0(MAXQB),bimode(MAXQB),
     $     bicurv(MAXQB,MAXQB),bhinv(MAXQB,MAXQB),
     $     bim0(MAXQB),eycb2i(MAXQ,MAXQ),vi(MAXPV),ebi(MAXQB),
     $     ebbi(MAXQB,MAXQB),gamhexpx(MAXPXTQ),Ehexp,Ehexpxz(MAXQ),
     $     Ehexpxz2(MAXQ,MAXQ),fderivar(MAXPVPQ,MAXPVPQ),
     $     newsigb(MAXQB,MAXQB),newsige(MAXQ,MAXQ),
     $     mcbetat(MAXPVPQ,MAXPVPQ),mcvar(MAXPVPQ),
     $     betas(MAXPXTQ+MAXPVPQ),newbetas(MAXPXTQ+MAXPVPQ),det(2)
      character kstr*2,labstr*20
      
      qb = 2*qe
      pxtq = px*qe
      pvpq = pv + qe
      nfixpar = pxtq + pvpq
      qea = qe**2
      qeu = qe*(qe+1)/2
      
      ldsb2 = MAXQB
      ldsb3 = MAXQB
      ldsb4 = MAXQB
      ldse2 = MAXQ
      ldse3 = MAXQ
      ldse4 = MAXQ
      
      ldtq = MAXQEA
      lddq = MAXQEA
      ldpq = MAXQEU	  
      ldke = MAXQEA	  
      ldgm = MAXQ
      ldbc = MAXQB
      ldbh = MAXQB
      ldey2 = MAXQ
      ldebbi = MAXQB
      ldfv = MAXPVPQ
      ldmcbet = MAXPVPQ
      ldsmc = MAXPVPQ
      ldez2 = MAXQ
      
      eps1b = 1.0d-6
      maxiter3 = 50
      
c      call dblepr('time',-1,time,m)
c      call intpr('event',-1,event,m)
c      call intpr('id1',-1,id1,m)
c      call dblepr('Q(,1)',-1,Q(1,1),sumn)
c      call dblepr('Q(,2)',-1,Q(1,2),sumn)
c      call intpr('ldq',-1,ldq,1)
c      call intpr('qe',-1,qe,1)
c      call intpr('D(,1)',-1,D(1,1),sumn)
c      call intpr('D(,2)',-1,D(1,2),sumn)
c      call intpr('ldd',-1,ldd,1)
c      call intpr('id2',-1,id2,sumn)
c      call intpr('repl',-1,repl,sumn)
c      call dblepr('V(,1)',-1,V(1,1),m)
c      call intpr('ldv',-1,ldv,1)
c      call intpr('pv',-1,pv,1)
c      call dblepr('X(,1)',-1,X(1,1),sumn)
c      call dblepr('X(,2)',-1,X(1,2),sumn)
c      call dblepr('X(,3)',-1,X(1,3),sumn)
c      call intpr('ldx',-1,ldx,1)
c      call intpr('px',-1,px,1)
c      call dblepr('Z(,1)',-1,Z(1,1),sumn)
c      call dblepr('Z(,2)',-1,Z(1,2),sumn)
c      call intpr('ldz',-1,ldz,1)
c      call dblepr('eb(,1)',-1,eb(1,1),m)
c      call dblepr('eb(,2)',-1,eb(1,2),m)
c      call intpr('ldeb',-1,ldeb,1)
c      call dblepr('eyc(,1)',-1,eyc(1,1),sumn)
c      call dblepr('eyc(,2)',-1,eyc(1,2),sumn)
c      call intpr('ldeyc',-1,ldeyc,1)
c      call dblepr('h0t',-1,h0t,56)
c      call dblepr('alphamat(,1)',-1,alphamat(1,1),px)
c      call dblepr('alphamat(,2)',-1,alphamat(1,2),px)
c      call intpr('ldam',-1,ldam,1)
c      call dblepr('beta',-1,beta,pv)
c      call dblepr('gamma',-1,h0t,qe)
c      call dblepr('sigb(,1)',-1,sigb(1,1),4)
c      call dblepr('sigb(,2)',-1,sigb(1,2),4)
c      call dblepr('sigb(,3)',-1,sigb(1,3),4)
c      call dblepr('sigb(,4)',-1,sigb(1,4),4)
c      call intpr('ldsb',-1,ldsb,1)
c      call dblepr('sige(,1)',-1,sige(1,1),qe)
c      call dblepr('sige(,2)',-1,sige(1,2),qe)
c      call intpr('ldse',-1,ldse,1)
c      call intpr('maxiter',-1,maxiter,1)
c      call dblepr('eps1',-1,eps1,1)
c      call dblepr('eps2',-1,eps2,1)
c      call intpr('maxiter2',-1,maxiter2,1)
c      call intpr('nburnin',-1,nburnin,1)
c      call intpr('nsize',-1,nsize,1)
c      call intpr('maxnsize',-1,maxnsize,1)
c      call intpr('nbatchs',-1,nbatchs,1)
c      call dblepr('ddi(,1)',-1,ddi(1,1),MAXNIS)
c      call dblepr('ddi(,2)',-1,ddi(1,2),MAXNIS)
c      return
      
c     calculate the length of each id: assume id1 and id2 are sorted and share ids;
c       also assume repl_s take integers 1, 2, ..., and are the actual visit numbers,
c       i.e., (repl_s-1)*samp_pd are the actual times in years, every subject must have
c       baseline visit 1.
      ina = 1
      ietimes0 = 2
c     note: leave one extra double precision for possible memory access past lower bound
      ne0 = 0
      do 10 i=1,m
        iwkarr(ina+i-1) = 0
        if (event(i).eq.1) then
          ne0=ne0+1
          wkarr(ietimes0+ne0-1) = time(i)
        end if
 10   continue
c      call intpr('ne0',-1,ne0,1)
c      call dblepr('etimes0',-1,wkarr(ietimes0),ne0)
c      return
      
      cid = id2(1)
      idx1 = ina
      maxrepl = 0
      do 20 k=1,sumn
        if (.not.(id2(k).eq.cid)) then
          cid=id2(k)
          idx1=idx1 + 1
        end if
        
        iwkarr(idx1)=iwkarr(idx1)+1
        
        if (repl(k).gt.maxrepl) maxrepl=repl(k)
 20   continue
c      call intpr('ni_s',-1,iwkarr(ina),m)
c      call intpr('maxrepl',-1,maxrepl,1)
c      return
      
      
c     sort the event times with insertion sort (not the fastest way but the code is simple)
      do 40 k=2,ne0
        tmp1 = wkarr(ietimes0+k-1)
        l = k - 1
 30     continue
        if ((l.ge.1).and.(wkarr(ietimes0+l-1).gt.tmp1)) then
          wkarr(ietimes0+l)=wkarr(ietimes0+l-1)
          l = l - 1
          go to 30
        end if
        wkarr(ietimes0+l)=tmp1
 40   continue
c      call dblepr('etimes0',-1,wkarr(ietimes0),ne0)
c	  return
      
      ietimes1 = ietimes0 + ne0
      ne1 = 1
      wkarr(ietimes1+ne1-1) = wkarr(ietimes0)
      do 50 k=2,ne0
        if (wkarr(ietimes0+k-1).ne.wkarr(ietimes1+ne1-1)) then
          ne1 = ne1 + 1
          wkarr(ietimes1+ne1-1) = wkarr(ietimes0+k-1)
        end if
 50   continue
c      call intpr('ne0',-1,ne0,1)
c      call intpr('ne1',-1,ne1,1)
c      call dblepr('etimes1',-1,wkarr(ietimes1),ne1)
c      call intpr('nfixpar',-1,nfixpar,1)
c      call intpr('qeu',-1,qeu,1)
c	  return
      
      ntotpar2 = ne1 + nfixpar + qeu
c      call intpr('ntotpar2',-1,ntotpar2,1)
c	  return
      
      
c     calculate the number of event time at or before each observed or censored time
      itd = ina + m
      sumtd = 0
      do 60 i=1,m
        idx1 = itd + i - 1
        iwkarr(idx1) = 0
        tmp1 = time(i)
        do 60 k=1,ne1
          if (wkarr(ietimes1+k-1).le.tmp1) then
            iwkarr(idx1)=iwkarr(idx1) + 1
            sumtd=sumtd+1
          end if
 60   continue
c      call intpr('sumtd',-1,sumtd,1)
c      call intpr('td_s',-1,iwkarr(itd),m)
c	  return
      
      
c     calculate the summary left-censoring indicator by each subject and each subject/visit
      idla = itd + m
      idl = idla + m
      maxni = 0
      offset1 = 0
      do 90 i=1,m
        ni = iwkarr(ina+i-1)
        if (ni.gt.maxni) maxni=ni
        dila = 0
        do 80 j=1,ni
          dilj = 0
          idx1 = offset1 + j
          do 70 k=1,qe
            if (D(idx1,k).eq.1) then
              dila = 1
              dilj = 1
            end if
 70       continue
          iwkarr(idl+idx1-1) = dilj
 80     continue
        iwkarr(idla+i-1) = dila
        offset1 = offset1 + ni
 90   continue
c      call intpr('maxni',-1,maxni,1)
c      call intpr('dla',-1,iwkarr(idla),m)
c      call intpr('dl',-1,iwkarr(idl),sumn)
c	  return
      
      
c     generate the t2imap data matrix
      it2imap = ietimes1 + ne1
      do 100 k=1,maxrepl
        wkarr(it2imap+k-1) = k
        wkarr(it2imap+maxrepl+k-1) = (k-1)*samp_pd
        wkarr(it2imap+2*maxrepl+k-1) = k*samp_pd
 100  continue
c      call dblepr('t2imap(,1)',-1,wkarr(it2imap),maxrepl)
c      call dblepr('t2imap(,2)',-1,wkarr(it2imap+maxrepl),maxrepl)
c      call dblepr('t2imap(,3)',-1,wkarr(it2imap+2*maxrepl),maxrepl)	  
c	  return
      
      
c     generate the matrices xt, zt, ztz
      sumtdpm = sumtd + m
      ixt = it2imap + maxrepl*3
      ldxt = sumtdpm
      izt = ixt + ldxt*px
      ldzt = sumtdpm
      iztz = izt + ldzt*2
      ldztz = sumtd
      offset1 = 0
      offset2 = 0
      offset3 = 0
      do 200 i=1,m
        ni = iwkarr(ina+i-1)        
        tid = iwkarr(itd+i-1)
        tidp1 = tid + 1
        if (tid.gt.0) then      
          do 150 j=1,tid
            etimes1j = wkarr(ietimes1+j-1)
            do 110 k=1,maxrepl
              if ((wkarr(it2imap+maxrepl+k-1).lt.etimes1j).and.
     $              (etimes1j.le.wkarr(it2imap+2*maxrepl+k-1))) 
     $              replidx=int(wkarr(it2imap+k-1))
 110        continue
            
            do 120 k=1,ni
              if (repl(offset3+k).le.replidx) replidx2=k
 120        continue      
            idx1 = offset3 + replidx2
            
            do 130 k=1,px
              wkarr(ixt+(k-1)*sumtdpm+offset1+j-1)=X(idx1,k)
 130        continue
            
            do 140 k=1,2
              wkarr(izt+(k-1)*sumtdpm+offset1+j-1)=Z(idx1,k)
 140        continue
            
            do 150 k=1,2
              do 150 l=1,2
                idx2 = (k-1)*2 + l
                wkarr(iztz+(idx2-1)*sumtd+offset2+j-1)=Z(idx1,l)*
     $               Z(idx1,k)
 150      continue
        end if
        
        tmp1 = time(i)
        do 160 k=1,maxrepl
          if ((wkarr(it2imap+maxrepl+k-1).lt.tmp1).and.
     $          (tmp1.le.wkarr(it2imap+2*maxrepl+k-1))) 
     $          replidx=int(wkarr(it2imap+k-1))
 160    continue
        
        do 170 k=1,ni
          if (repl(offset3+k).le.replidx) replidx2=k
 170    continue      
        idx1 = offset3 + replidx2
        
        do 180 k=1,px
          wkarr(ixt+(k-1)*sumtdpm+offset1+tid)=X(idx1,k)
 180    continue
        
        do 190 k=1,2
          wkarr(izt+(k-1)*sumtdpm+offset1+tid)=Z(idx1,k)
 190    continue
        
        offset1 = offset1 + tidp1
        offset2 = offset2 + tid
        offset3 = offset3 + ni
 200  continue
c      call dblepr('xt(,1)',-1,wkarr(ixt),sumtdpm)
c      call dblepr('xt(,2)',-1,wkarr(ixt+sumtdpm),sumtdpm)
c      call dblepr('xt(,3)',-1,wkarr(ixt+2*sumtdpm),sumtdpm)	  
c      call dblepr('zt(,1)',-1,wkarr(izt),sumtdpm)
c      call dblepr('zt(,2)',-1,wkarr(izt+sumtdpm),sumtdpm)
c      call dblepr('ztz(,1)',-1,wkarr(iztz),sumtd)
c      call dblepr('ztz(,2)',-1,wkarr(iztz+sumtd),sumtd)
c      call dblepr('ztz(,3)',-1,wkarr(iztz+2*sumtd),sumtd)
c      call dblepr('ztz(,4)',-1,wkarr(iztz+3*sumtd),sumtd)
c	  return
      
      
c     generate the communication matrix for sige
      call matzero(ldtq,qea,qea,tqeqe)
      do 210 i=1,qe
        do 210 j=1,qe
          tqeqe((i-1)*qe+j,(j-1)*qe+i) = 1.0d0
 210  continue
c      call dblepr('tqeqe(,1)',-1,tqeqe(1,1),qea)
c      call dblepr('tqeqe(,2)',-1,tqeqe(1,2),qea)
c      call dblepr('tqeqe(,3)',-1,tqeqe(1,3),qea)
c      call dblepr('tqeqe(,4)',-1,tqeqe(1,4),qea)
c	  return
      
      
c     generate the duplication matrix for sige
      call matzero(lddq,qea,qeu,dqe)
      sseq(1)=0
      do 220 j=1,qe-1
        sseq(j+1)=sseq(j)+j
 220  continue
      
      do 230 j=1,qe
        do 230 i=1,qe
          idx1 = (j-1)*qe + i
          if (i.le.j) then
            dqe(idx1,sseq(j)+i) = 1.0d0
          else
            dqe(idx1,sseq(i)+j) = 1.0d0
          end if
 230  continue
c      call dblepr('dqe(,1)',-1,dqe(1,1),qea)
c      call dblepr('dqe(,2)',-1,dqe(1,2),qea)
c      call dblepr('dqe(,3)',-1,dqe(1,3),qea)
c	  return
      
      
c     generate the permutation matrix for sige
      call matzero(ldpq,qeu,qeu,pqe)
      do 240 j=1,qe
        do 240 i=1,j
          kt = (j-1)*j/2 + i
          ks = qe*(i-1) - (i-1)*(i-2)/2 + (j-i+1)
          pqe(kt,ks) = 1.0d0
 240  continue
c      call dblepr('pqe(,1)',-1,pqe(1,1),qeu)
c      call dblepr('pqe(,2)',-1,pqe(1,2),qeu)
c      call dblepr('pqe(,3)',-1,pqe(1,3),qeu)
c	  return
      
      
c     populate the combined permutation matrix
      ip2 = iztz + ldztz*4
      ldp2 = ntotpar2
      call matzero(ldp2,ntotpar2,ntotpar2,wkarr(ip2))
      do 250 j=1,ne1+nfixpar
        wkarr(ip2+(j-1)*ntotpar2+j-1) = 1.0d0
 250  continue
      
      offset1 = ne1+nfixpar
      do 260 j=1,qeu
        do 260 i=1,qeu
          wkarr(ip2+(offset1+j-1)*ntotpar2+offset1+i-1) = pqe(i,j)
 260  continue
c      call dblepr('p2(,1)',-1,wkarr(ip2),ntotpar2)
c      call dblepr('p2(,56)',-1,wkarr(ip2+55*ntotpar2),ntotpar2)
c      call dblepr('p2(,57)',-1,wkarr(ip2+56*ntotpar2),ntotpar2)
c      call dblepr('p2(,65)',-1,wkarr(ip2+64*ntotpar2),ntotpar2)
c      call dblepr('p2(,66)',-1,wkarr(ip2+65*ntotpar2),ntotpar2)
c      call dblepr('p2(,68)',-1,wkarr(ip2+67*ntotpar2),ntotpar2)
c	  return
            
c     Make sure yic0 is the same as qic for those uncensored observations
      offset1 = 0
      do 266 i=1,m
        ni = iwkarr(ina+i-1)
        do 264 k=1,qe
          do 264 j=1,ni
            if (D(offset1+j,k).eq.0) then
              idx1 = offset1 + j
              eyc(idx1,k) = Q(idx1,k)
              do 262 l=1,px
                eyc(idx1,k)=eyc(idx1,k)-X(idx1,l)*alphamat(l,k)
 262          continue
            end if
 264    continue
        offset1 = offset1 + ni
 266  continue
c      call dblepr('eyc(,1)',-1,eyc(1,1),sumn)
c      call dblepr('eyc(,2)',-1,eyc(1,2),sumn)
c      return
      
      
c     Start the iteration      
      idi = idl + sumn
      lddi = maxni
      ibatchs = idi + lddi*qe
      iszbatch = ibatchs + maxnsize
c     dimension of szbatch is nbatchs by 1 (a vector)
      
c      ibetas = ip2 + ldp2*ntotpar2
      ivbeta = ip2 + ldp2*ntotpar2
      isumdf = ivbeta + m
      isumd2f = isumdf + ntotpar2
      ldsd = ntotpar2
c      isummct = isumd2f + ldsd*ntotpar2
c      ldsmc = pvpq
      iebp = isumd2f + ldsd*ntotpar2
      ldebp = m
      ieycp = iebp + ldebp*qb
      ldeycp = sumn
      ixi = ieycp + ldeycp*qe
      ldxi = maxni
      izi = ixi + ldxi*px
      ldzi = maxni
      iqic = izi + ldzi*2
      ldqic = maxni
      ieyic = iqic + ldqic*qe
      ldeyic = maxni
      iyic0 = ieyic + ldeyic*qe
      ldyic0 = maxni      
      ixit = iyic0 + ldyic0*qe
      ldxit = ne1 + 1
      izit = ixit + ldxit*px
      ldzit = ne1 + 1
      izitz = izit + ldzit*2
      ldzitz = ne1
      ibiarr = izitz + ldzitz*4
      ldbiarr = qb
      ieycbi = ibiarr + ldbiarr*maxnsize
      ldeycbi = maxni
      iwkzb = ieycbi + ldeycbi*qe
      ldwkzb = qe
      iEvibpoi = iwkzb + ldwkzb*maxni
c      igamhex1 = iEvibpoi + ne1
c      iEhexpxz = igamhex1 + pxtq
      igamexxm = iEvibpoi + ne1
      ldgx = ne1
      iexpxzm = igamexxm + ldgx*pxtq
      ldem = ne1
      igamhex2 = iexpxzm + ldem*qe
      ldgh = pxtq
      ihexpxzc = igamhex2 + ldgh*pxtq
      ldez = qe
c      ihexpxz2 = ihexpxzc + ldez*pxtq
c      ldez2 = qe
      itmpvec1 = ihexpxzc + ldez*pxtq
      itmpvec2 = itmpvec1 + ne1
      iEhvibpo = itmpvec2 + ne1
      iEvibpob = iEhvibpo + ne1
      ldvb = ne1
      iEhvibpb = iEvibpob + ldvb*qb
      ldvb2 = ne1
      ihzitz = iEhvibpb + ldvb2*qb
      ldhz = ne1
      ixitalp = ihzitz + ldhz*4
      ldxa = ne1
      inewh0t = ixitalp + ldxa*qe
c      inewbets = inewh0t + ne1
      iJ2 = inewh0t + ne1
      ldj2 = ntotpar2
      idJ2 = iJ2 + ldj2*ntotpar2
      lddj2 = ntotpar2**2
      inatrpar = idJ2 + lddj2*ntotpar2
      innatpar = inatrpar + ntotpar2
      inatgrad = innatpar + ntotpar2
      inathess = inatgrad + ntotpar2
      ldnh = ntotpar2
      inathinv = inathess + ldnh*ntotpar2
      ldni = ntotpar2
      idJwk = inathinv + ldni*ntotpar2
      lddjwk = max(ne1,qeu)**2
c     dimension of dJwk is max(ne1,qeu)**2 by max(ne1,qeu)
      
      do 440 iter=1,maxiter
        call ifprintf(iter,0,1)
        
        offset1 = pxtq + pv
        do 270 k=1,qe
c          wkarr(ibetas+offset1+k-1)=gamma(k)
          betas(offset1+k)=gamma(k)
          do 270 l=1,px
            idx1 = (k-1)*px + l
c            wkarr(ibetas+idx1-1)=alphamat(l,k)
            betas(idx1)=alphamat(l,k)
 270    continue
        
        do 280 k=1,pv
c          wkarr(ibetas+pxtq+k-1)=beta(k)
          betas(pxtq+k)=beta(k)
 280    continue
c      call dblepr('betas',-1,betas,nfixpar)
c	  return
        
        do 290 i=1,m
          idx1 = ivbeta + i - 1 
          wkarr(idx1)=0.0d0
          do 290 k=1,pv
            wkarr(idx1)=wkarr(idx1)+V(i,k)*beta(k)
 290    continue
c      call dblepr('Vbeta',-1,wkarr(ivbeta),m)
c	  return
        
        call mtxcopy(sigb,ldsb,qb,qb,sigbinv,ldsb2)
        call invert(ldsb2,qb,sigbinv,det,ierr)
c      call dblepr('sigbinv(,1)',-1,sigbinv(1,1),qb)
c      call dblepr('sigbinv(,2)',-1,sigbinv(1,2),qb)
c      call dblepr('sigbinv(,3)',-1,sigbinv(1,3),qb)
c      call dblepr('sigbinv(,4)',-1,sigbinv(1,4),qb)
c	  return
        
        call mtxcopy(sige,ldse,qe,qe,sigeinv,ldse2)
        call invert(ldse2,qe,sigeinv,det,ierr)
c      call dblepr('sigeinv(,1)',-1,sigeinv(1,1),qe)
c      call dblepr('sigeinv(,2)',-1,sigeinv(1,2),qe)
c	  return
        
c        call matzero(ldke,qea,qea,krwe2)
        do 300 j=1,qe
          do 300 i=1,qe
            do 300 l=1,qe
              do 300 k=1,qe
                krwe2((i-1)*qe+k,(j-1)*qe+l) = sigeinv(k,l)*sigeinv(i,j)
 300    continue
c      call dblepr('krwe2(,1)',-1,krwe2(1,1),qea)
c      call dblepr('krwe2(,2)',-1,krwe2(1,2),qea)
c      call dblepr('krwe2(,3)',-1,krwe2(1,3),qea)
c      call dblepr('krwe2(,4)',-1,krwe2(1,4),qea)
c	  return
        
        call setgibbs(sige,ldse,qe,gmat,ldgm,cstd)
c      call dblepr('gmat(,1)',-1,gmat(1,1),qe)
c      call dblepr('gmat(,2)',-1,gmat(1,2),qe)
c      call dblepr('cstd',-1,cstd,qe)
c	  return
        
        call matzero(ldsb3,qb,qb,sumEbb)
        call matzero(ldse3,qe,qe,sumEee)
        call veczero(ntotpar2,wkarr(isumdf))
        call matzero(ldsd,ntotpar2,ntotpar2,wkarr(isumd2f))
        call matzero(ldsmc,pvpq,pvpq,sumMct)
        
        call mtxcopy(eb,ldeb,m,2,wkarr(iebp),ldebp)
        call mtxcopy(eyc,ldeyc,sumn,qe,wkarr(ieycp),ldeycp)
        
        offset1 = 0
        offset2 = 0
        offset3 = 0
        do 430 i=1,m
c          call ifprintf(i,1,0)
          if (mod(i,5).eq.0) call ifprintf(i,1,0)
          
          ni = iwkarr(ina+i-1)  
          
          do 310 k=1,qb
            ebi0(k) = eb(i,k)
 310      continue
c      call dblepr('ebi0',-1,ebi0,qb)
c	  return
          
          do 350 j=1,ni
            idx1 = offset3 + j
            
            do 320 k=1,px
              wkarr(ixi+(k-1)*ldxi+j-1) = X(idx1,k)
 320        continue
            
            do 330 k=1,2
              wkarr(izi+(k-1)*ldzi+j-1) = Z(idx1,k)
 330        continue
            
            do 350 k=1,qe
              idx2 = iqic + (k-1)*ldqic + j - 1
              wkarr(idx2) = Q(idx1,k)
              do 340 l=1,px
                wkarr(idx2)=wkarr(idx2)-X(idx1,l)*alphamat(l,k)
 340          continue
              
              wkarr(ieyic+(k-1)*ldeyic+j-1) = eyc(idx1,k)
              
              iwkarr(idi+(k-1)*lddi+j-1) = D(idx1,k)
c              ddi(j,k) = D(idx1,k)*1.0d0
 350      continue
c      call dblepr('xi(,1)',-1,wkarr(ixi),ldxi)
c      call dblepr('xi(,2)',-1,wkarr(ixi+ldxi),ldxi)
c      call dblepr('xi(,3)',-1,wkarr(ixi+2*ldxi),ldxi)
c      call dblepr('zi(,1)',-1,wkarr(izi),ldzi)
c      call dblepr('zi(,2)',-1,wkarr(izi+ldzi),ldzi)
c      call dblepr('qic(,1)',-1,wkarr(iqic),ldqic)
c      call dblepr('qic(,2)',-1,wkarr(iqic+ldqic),ldqic)
c      call dblepr('eyic(,1)',-1,wkarr(ieyic),ldeyic)
c      call dblepr('eyic(,2)',-1,wkarr(ieyic+ldeyic),ldeyic)
c      call intpr('di(,1)',-1,iwkarr(idi),lddi)
c      call intpr('di(,2)',-1,iwkarr(idi+lddi),lddi)
c      call dblepr('ddi(,1)',-1,ddi(1,1),ni)
c      call dblepr('ddi(,2)',-1,ddi(1,2),ni)
c	  return
          
          tid = iwkarr(itd+i-1)
          tidp1 = tid + 1
          do 380 j=1,tidp1
            idx1 = offset1 + j
            
            do 360 k=1,px
              wkarr(ixit+(k-1)*ldxit+j-1)=wkarr(ixt+(k-1)*ldxt+idx1-1)
 360        continue
            
            do 370 k=1,2
              wkarr(izit+(k-1)*ldzit+j-1)=wkarr(izt+(k-1)*ldzt+idx1-1)
 370        continue          
 380      continue
c      call dblepr('xit(,1)',-1,wkarr(ixit),ldxit)
c      call dblepr('xit(,2)',-1,wkarr(ixit+ldxit),ldxit)
c      call dblepr('xit(,3)',-1,wkarr(ixit+2*ldxit),ldxit)
c      call dblepr('zit(,1)',-1,wkarr(izit),ldzit)
c      call dblepr('zit(,2)',-1,wkarr(izit+ldzit),ldzit)
c	  return
          
          if (tid.gt.0) then
            do 390 j=1,tid
              idx2 = offset2 + j
              do 390 k=1,4
                 wkarr(izitz+(k-1)*ldzitz+j-1)=
     $                wkarr(iztz+(k-1)*ldztz+idx2-1)
 390        continue
          end if
c      call dblepr('zitz(,1)',-1,wkarr(izitz),ldzitz)
c      call dblepr('zitz(,2)',-1,wkarr(izitz+ldzitz),ldzitz)
c      call dblepr('zitz(,3)',-1,wkarr(izitz+2*ldzitz),ldzitz)
c      call dblepr('zitz(,4)',-1,wkarr(izitz+3*ldzitz),ldzitz)
c	  return
          
          do 400 k=1,pv
            vi(k) = V(i,k)
 400      continue
c      call dblepr('vi',-1,vi,pv)
c	  return
          
          
c      call intpr('ldqic',-1,ldqic,1)
c      call intpr('ni',-1,ni,1)
c      call intpr('qe',-1,qe,1)		  
c      call dblepr('qic(,1)',-1,wkarr(iqic),ldqic)
c      call dblepr('qic(,2)',-1,wkarr(iqic+ldqic),ldqic)		  
c      call intpr('i',-1,i,1)
c      call intpr('dila',-1,iwkarr(idla+i-1),1)
c      call intpr('offset3',-1,offset3,1)
c      call intpr('dil',-1,iwkarr(idl+offset3),ni)
c      call intpr('lddi',-1,lddi,1)
c      call intpr('di(,1)',-1,iwkarr(idi),lddi)
c      call intpr('di(,2)',-1,iwkarr(idi+lddi),lddi)
c      call intpr('ldddi',-1,ldddi,1)
c      call dblepr('ddi(,1)',-1,ddi(1,1),ldddi)
c      call dblepr('ddi(,2)',-1,ddi(1,2),ldddi)
c      call intpr('qb',-1,qb,1)
c      call dblepr('ebi0',-1,ebi0,qb)
c      call intpr('ldsb2',-1,ldsb2,1)
c      call dblepr('sigbinv(,1)',-1,sigbinv(1,1),ldsb2)
c      call dblepr('sigbinv(,4)',-1,sigbinv(1,4),ldsb2)
c      call intpr('deltai',-1,event(i),1)
c      call dblepr('vibeta',-1,wkarr(ivbeta+i-1),1)
c      call intpr('tid',-1,tid,1)
c      call dblepr('h0t',-1,h0t,tid)
c      call intpr('ldxit',-1,ldxit,1)
c      call intpr('px',-1,px,1)
c      call dblepr('xit(,1)',-1,wkarr(ixit),ldxit)
c      call dblepr('xit(,3)',-1,wkarr(ixit+2*ldxit),ldxit)
c      call intpr('ldzit',-1,ldzit,1)
c      call dblepr('zit(,1)',-1,wkarr(izit),ldzit)
c      call dblepr('zit(,2)',-1,wkarr(izit+ldzit),ldzit)
c      call intpr('ldam',-1,ldam,1)
c      call dblepr('alphamat(,1)',-1,alphamat(1,1),ldam)
c      call dblepr('alphamat(,2)',-1,alphamat(1,2),ldam)
c      call dblepr('gamma',-1,gamma,qe)
c      call intpr('ldeyic',-1,ldeyic,1)
c      call dblepr('eyic(,1)',-1,wkarr(ieyic),ldeyic)
c      call dblepr('eyic(,2)',-1,wkarr(ieyic+ldeyic),ldeyic)
c      call intpr('ldzi',-1,ldzi,1)
c      call dblepr('zi(,1)',-1,wkarr(izi),ldzi)
c      call dblepr('zi(,2)',-1,wkarr(izi+ldzi),ldzi)
c      call intpr('ldse2',-1,ldse2,1)
c      call dblepr('sigeinv(,1)',-1,sigeinv(1,1),ldse2)
c      call dblepr('sigeinv(,2)',-1,sigeinv(1,2),ldse2)
c      call intpr('maxiter2',-1,maxiter2,1)
c      call dblepr('eps1b',-1,eps1b,1)
c      call intpr('maxiter3',-1,maxiter3,1)
c      call intpr('ldgm',-1,ldgm,1)
c      call dblepr('gmat(,1)',-1,gmat(1,1),ldgm)
c      call dblepr('gmat(,2)',-1,gmat(1,2),ldgm)
c      call dblepr('cstd',-1,cstd,qe)
c      call dblepr('bimode',-1,bimode,qb)
c      call intpr('ldbc',-1,ldbc,1)
c      call dblepr('bicurv(,1)',-1,bicurv(1,1),ldbc)
c      call dblepr('bicurv(,4)',-1,bicurv(1,4),ldbc)
c      call intpr('ldbh',-1,ldbh,1)
c      call dblepr('bhinv(,1)',-1,bhinv(1,1),ldbh)
c      call dblepr('bhinv(,4)',-1,bhinv(1,4),ldbh)
c      call dblepr('bim0',-1,bim0,qb)
c      call intpr('ldyic0',-1,ldyic0,1)
c      call dblepr('yic0(,1)',-1,wkarr(iyic0),ldyic0)
c      call dblepr('yic0(,2)',-1,wkarr(iyic0+ldyic0),ldyic0)
c	  return
          call getycbm0(wkarr(iqic),ldqic,ni,qe,iwkarr(idla+i-1),
     $         iwkarr(idl+offset3),iwkarr(idi),lddi,ebi0,qb,sigbinv,
     $         ldsb2,event(i),wkarr(ivbeta+i-1),h0t,tid,wkarr(ixit),
     $         ldxit,px,wkarr(izit),ldzit,alphamat,ldam,gamma,
     $         wkarr(ieyic),ldeyic,wkarr(izi),ldzi,sigeinv,ldse2,
     $         maxiter2,eps1b,maxiter3,gmat,ldgm,cstd,bimode,bicurv,
     $         ldbc,bhinv,ldbh,bim0,wkarr(iyic0),ldyic0)
c          if (i.eq.2) then
c          call ifprintf(0,0,0)
c          call dblepr('bimode',-1,bimode,qb)
c          call dblepr('bicurv(,1)',-1,bicurv(1,1),qb)
c          call dblepr('bicurv(,4)',-1,bicurv(1,4),qb)
c          call dblepr('bim0',-1,bim0,qb)
c          call dblepr('yic0(,1)',-1,wkarr(iyic0),ldyic0)
c          call dblepr('yic0(,2)',-1,wkarr(iyic0+ldyic0),ldyic0)          
c          return
c          end if
          
          call mgibbs(wkarr(iqic),ldqic,ni,qe,iwkarr(idla+i-1),
     $         iwkarr(idl+offset3),iwkarr(idi),lddi,sigbinv,ldsb2,qb,
     $         event(i),wkarr(ivbeta+i-1),h0t,tid,wkarr(ixit),ldxit,px,
     $         wkarr(izit),ldzit,alphamat,ldam,gamma,wkarr(izi),ldzi,
     $         sigeinv,ldse2,gmat,ldgm,cstd,bimode,bicurv,ldbc,bhinv,
     $         ldbh,bim0,wkarr(iyic0),ldyic0,wkarr(ibiarr),ldbiarr,
     $         nburnin,nsize,wkarr(ieycbi),ldeycbi,eycb2i,ldey2,
     $         wkarr(iwkzb),ldwkzb)
c          if (i.eq.2) then
c          call dblepr('biarr(,1)',-1,wkarr(ibiarr),qb)
c          call dblepr('biarr(,4000)',-1,wkarr(ibiarr+(4000-1)*ldbiarr),
c     $         qb)
c          call dblepr('eycbi(,1)',-1,wkarr(ieycbi),qb)
c          call dblepr('eycbi(,2)',-1,wkarr(ieycbi+ldeycbi),qb)
c          call dblepr('eycb2i(,1)',-1,eycb2i(1,1),qe)
c          call dblepr('eycb2i(,2)',-1,eycb2i(1,2),qe)          
c	   return
c          end if
          
          call evalwpb(wkarr(ibiarr),ldbiarr,qb,nsize,event(i),
     $         wkarr(ivbeta+i-1),h0t,tid,wkarr(ixit),ldxit,px,
     $         wkarr(izit),ldzit,alphamat,ldam,qe,gamma,vi,pv,
     $         wkarr(izitz),ldzitz,Ebi,Ebbi,ldebbi,wkarr(iEvibpoi),
     $         gamhexpx,Ehexp,Ehexpxz,wkarr(igamexxm),
     $         ldgx,wkarr(iexpxzm),ldem,wkarr(igamhex2),ldgh,
     $         wkarr(ihexpxzc),ldez,Ehexpxz2,ldez2,      
     $         fderivar,ldfv,iwkarr(ibatchs),iwkarr(iszbatch),nbatchs,
     $         wkarr(itmpvec1),wkarr(itmpvec2),wkarr(iEhvibpo),
     $         wkarr(iEvibpob),ldvb,wkarr(iEhvibpb),ldvb2,wkarr(ihzitz),
     $         ldhz,wkarr(ixitalp),ldxa)
c          if (i.eq.2) then
c          call dblepr('Evibpoi',-1,wkarr(iEvibpoi),tid)
c          call dblepr('Ehexpxzc(,1)',-1,wkarr(ihexpxzc),qe)
c          call dblepr('Ehexpxzc(,6)',-1,wkarr(ihexpxzc+5*ldez),qe)
c          call dblepr('fderivar(,1)',-1,fderivar(1,1),pvpq)
c          call dblepr('fderivar(,3)',-1,fderivar(1,3),pvpq)
c	  return
c          end if
          
          do 410 k=1,qb
            eb(i,k) = ebi(k)
 410      continue
c          call dblepr('ebi',-1,ebi,qb)          
          
          do 420 k=1,qe
            do 420 j=1,ni
              idx1 = offset3 + j
              eyc(idx1,k) = wkarr(ieycbi+(k-1)*ldeycbi+j-1)
              do 420 l=1,2
                eyc(idx1,k) = eyc(idx1,k) + wkarr(izi+(l-1)*ldzi+j-1)*
     $                ebi((k-1)*2+l)
 420      continue
c          call dblepr('eyic(,1)',-1,eyc(offset3+1,1),ni)          
c          call dblepr('eyic(,2)',-1,eyc(offset3+1,2),ni)          
c          return
          
          call sumAll1(Ebi,qb,Ebbi,ldebbi,wkarr(iEvibpoi),tid,
     $         gamhexpx,px,qe,Ehexp,Ehexpxz,
     $         wkarr(igamexxm),ldgx,wkarr(iexpxzm),ldem,wkarr(igamhex2),
     $         ldgh,wkarr(ihexpxzc),ldez,Ehexpxz2,ldez2,fderivar,
     $         ldfv,pv,event(i),h0t,vi,wkarr(ieycbi),ldeycbi,eycb2i,
     $         ldey2,wkarr(ixi),ldxi,ni,wkarr(ixit),ldxit,wkarr(izit),
     $         ldzit,alphamat,ldam,gamma,sigeinv,ldse2,krwe2,ldke,dqe,
     $         lddq,tqeqe,ldtq,sumebb,ldsb3,wkarr(isumdf),ne1,
     $         wkarr(isumd2f),ldsd,summct,ldsmc)
          
c          if (i.eq.300) then
c          if ((iter.eq.2).and.(i.eq.300)) then
c          call ifprintf(0,0,0)
c          call dblepr('sumDf',-1,wkarr(isumdf),ntotpar2)          
c          call dblepr('sumD2f(,1)',-1,wkarr(isumd2f),ntotpar2)          
c          call dblepr('sumD2f(,56)',-1,wkarr(isumd2f+55*ldsd),ntotpar2)          
c          call dblepr('sumD2f(,57)',-1,wkarr(isumd2f+56*ldsd),ntotpar2)          
c          call dblepr('sumD2f(,62)',-1,wkarr(isumd2f+61*ldsd),ntotpar2)          
c          call dblepr('sumD2f(,63)',-1,wkarr(isumd2f+62*ldsd),ntotpar2)          
c          call dblepr('sumD2f(,64)',-1,wkarr(isumd2f+63*ldsd),ntotpar2)          
c          call dblepr('sumD2f(,65)',-1,wkarr(isumd2f+64*ldsd),ntotpar2)          
c          call dblepr('sumD2f(,66)',-1,wkarr(isumd2f+65*ldsd),ntotpar2)          
c          call dblepr('sumD2f(,68)',-1,wkarr(isumd2f+67*ldsd),ntotpar2)          
c          call dblepr('sumEbb(,1)',-1,sumEbb(1,1),qb)          
c          call dblepr('sumEbb(,4)',-1,sumEbb(1,4),qb)          
c          call dblepr('sumMct(,1)',-1,summct(1,1),pvpq)          
c          call dblepr('sumMct(,3)',-1,summct(1,3),pvpq)          
c          return
c          end if
          
          offset1 = offset1 + tidp1
          offset2 = offset2 + tid
          offset3 = offset3 + ni
c          call ifprintf(0,0,0)
          if (mod(i,5).eq.0) call ifprintf(0,0,0)
 430    continue
        
        call updatep(sumebb,ldsb3,qb,m,wkarr(isumdf),ne1,px,qe,pv,
     $       wkarr(isumd2f),ldsd,summct,ldsmc,h0t,betas,
     $       sige,ldse,wkarr(iP2),ldp2,wkarr(inewh0t),newbetas,
     $       newsigb,ldsb4,newsige,ldse4,mcbetat,ldmcbet,mcvar,
     $       wkarr(iJ2),ldj2,wkarr(idJ2),lddj2,wkarr(inatrpar),
     $       wkarr(innatpar),wkarr(inatgrad),wkarr(inathess),ldnh,
     $       wkarr(inathinv),ldni,wkarr(idJwk),lddjwk)
c        if (iter.eq.2) then
        call dblepr('newh0t',-1,wkarr(inewh0t),ne1)
        call dblepr('newbetas',-1,newbetas,nfixpar)      
        do 423 k=1,qb
          call int2str(k,kstr,klen)
          labstr = 'newsigb(,'//kstr(1:klen)//')'
          call dblepr(labstr,10+klen,newsigb(1,k),qb)
 423    continue
        do 425 k=1,qe
          call int2str(k,kstr,klen)
          labstr = 'newsige(,'//kstr(1:klen)//')'
          call dblepr(labstr,10+klen,newsige(1,k),qe)
 425    continue        
c        call dblepr('mcbetat(,1)',-1,mcbetat(1,1),pvpq)
c        call dblepr('mcbetat(,3)',-1,mcbetat(1,3),pvpq)
c        return
c        end if
        
        call pctchg(h0t,ne1,betas,nfixpar,px,qe,sigb,ldsb,qb,
     $       sige,ldse,wkarr(inewh0t),newbetas,newsigb,ldsb4,
     $       newsige,ldse4,eps2,mcvar,hconv,sconv,pconv,dconv,mmcdel,
     $       dmcdel)
c        if (iter.eq.2) then        
        call dblepr('sconv',-1,sconv,ne1)
        call dblepr('hconv',-1,hconv,1)
        call dblepr('dconv',-1,dconv,nfixpar+qb*(qb+1)/2+qeu)
        call dblepr('pconv',-1,pconv,1)
        call dblepr('dmcdel',-1,dmcdel,pvpq)
        call dblepr('mmcdel',-1,mmcdel,1)
c        if (iter.eq.2) return
c        end if
        
        if (pconv.lt.eps1) go to 450
        
        if (mmcdel.lt.1.15d0) then
          nsize = min(nsize*2,maxnsize)
          call intpr('nsize increased to',-1,nsize,1)
        end if
c        call intpr('nsize',-1,nsize,1)
        
        call pcopy(wkarr(inewh0t),ne1,newbetas,nfixpar,newsigb,
     $       ldsb4,qb,newsige,ldse4,qe,h0t,betas,sigb,ldsb,sige,
     $       ldse)
c        call dblepr('h0t',-1,h0t,ne1)
c        call dblepr('betas',-1,betas,nfixpar)
c        call dblepr('sigb(,1)',-1,sigb(1,1),qb)
c        call dblepr('sigb(,4)',-1,sigb(1,4),qb)
c        call dblepr('sige(,1)',-1,sige(1,1),qe)
c        call dblepr('sige(,2)',-1,sige(1,2),qe)
c        return
        
        call genbtcom(betas,alphamat,ldam,px,qe,beta,pv,gamma)
c        call dblepr('alphamat(,1)',-1,alphamat(1,1),px)
c        call dblepr('alphamat(,2)',-1,alphamat(1,2),px)
c        call dblepr('beta',-1,beta,pv)
c        call dblepr('gamma',-1,gamma,qe)        
c        return
        
        call ifprintf(0,0,0)        
 440  continue
      iter = iter - 1
      
 450  continue
      
c     update the yic w.r.t. the updated alpha on the uncensored obs if exited due to max iter
      if (pconv.ge.eps1) then
      offset1 = 0
      do 480 i=1,m
        ni = iwkarr(ina+i-1)
        do 470 k=1,qe
          do 470 j=1,ni
            if (D(offset1+j,k).eq.0) then
              idx1 = offset1 + j
              eyc(idx1,k) = Q(idx1,k)
              do 460 l=1,px
                eyc(idx1,k)=eyc(idx1,k)-X(idx1,l)*alphamat(l,k)
 460          continue
            end if
 470    continue
        offset1 = offset1 + ni
 480  continue
c      call dblepr('eyc(,1)',-1,eyc(1,1),sumn)
c      call dblepr('eyc(,2)',-1,eyc(1,2),sumn)
      end if
      
      return
      end
