
C
C     Subroutines for survival analysis with left censoring covariate(s). 
C     This part of codes draws inference on the parameters.
C
      
C      
C     ChangeLog:
C
      
C     Metropolized Gibbs sampler
      subroutine mgibbs2(qic,ldqc,ni,q,dila,dil,di,lddi,sigbinv,ldsb,qb,
     $     deltai,vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,alphamat,ldam,
     $     gamma,zi,ldzi,sigeinv,ldse,gmat,ldgm,cstd,bimode,bicurv,
     $     ldbc,bhinv,ldbh,bim,yicm,ldycm,brarr,ldbra,nburnin,nsize,
     $     wkzb,ldwz)
      
c     Note: bim and yicm the initial bim0 and yic0 at the beginning and
c        are overwritten with the last bim and yic at the end of MCMC.
c        wkzb is a work array of the dimension q*maxni corresp. to t(Zi*bimm).
      parameter (MAXQ=4)
      integer ldqc,ni,q,dila,dil(*),lddi,ldsb,qb,deltai,tid,ldxt,px,
     $     ldzt,ldam,ldzi,ldse,ldgm,ldbc,ldbh,ldycm,ldbra,nburnin,nsize,
     $     ldwz
      integer di(lddi,*)
      double precision qic(ldqc,*),sigbinv(ldsb,*),vibeta,
     $     h0t(*),xit(ldxt,*),zit(ldzt,*),alphamat(ldam,*),gamma(*),
     $     zi(ldzi,*),sigeinv(ldse,*),gmat(ldgm,*),cstd(*),bimode(*),
     $     bicurv(ldbc,*),bhinv(ldbh,*),bim(*),yicm(ldycm,*),
     $     brarr(ldbra,*),wkzb(ldwz,*)
      
      integer j,k,l,s
      double precision r(2*MAXQ),bic(2*MAXQ),bimm(2,MAXQ),cmu,pp,ap,
     $     blogden,cpnorm,cqnorm
      
c      call dblepr('bimode',-1,bimode,qb)
c      call dblepr('bhinv(,1)',-1,bhinv(1,1),qb)
c      call dblepr('bhinv(,2)',-1,bhinv(1,2),qb)
c      call dblepr('bhinv(,3)',-1,bhinv(1,3),qb)
c      call dblepr('bhinv(,4)',-1,bhinv(1,4),qb)      
c      call dblepr('bim',-1,bim,qb)
c      call intpr('nburnin',-1,nburnin,1)
c      call intpr('nsize',-1,nsize,1)
c      return
      do 120 k=1,nburnin+nsize
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
          do 50 l=1,q
            do 50 s=1,2
              bimm(s,l) = bim((l-1)*2+s)
 50       continue
          
c        call intpr('ldwz',-1,ldwz,1)
c        return
          do 90 j=1,ni
            if ((k.gt.nburnin).or.(dil(j).eq.1)) then
              do 60 l=1,q
                wkzb(l,j) = 0.0d0
                do 60 s=1,2
                  wkzb(l,j) = wkzb(l,j) + zi(j,s)*bimm(s,l)
 60           continue
          
              if (dil(j).eq.1) then            
                do 80 l=1,q
                  if (di(j,l).eq.1) then
                    cmu = wkzb(l,j)
                    do 70 s=1,q
                      if (s.ne.l) cmu = cmu + 
     $                      gmat(s,l)*(yicm(j,s)-wkzb(s,j))
 70                 continue
                    call random(1,r)
                    pp = cpnorm(qic(j,l),cmu,cstd(l))*r(1)
                    yicm(j,l)=cqnorm(pp,cmu,cstd(l))
                  end if
 80             continue
              end if
            end if
 90       continue
          
          if (k.gt.nburnin) then
            do 100 l=1,qb
              brarr(l,k-nburnin)=bim(l)
 100        continue
            
            do 110 j=1,ni
              do 110 l=1,q
                brarr(qb+(j-1)*q+l,k-nburnin) = yicm(j,l)-wkzb(l,j)
 110        continue
          end if         
c         call dblepr('brarr(,1)',-1,brarr(1,1),qb)
c         call dblepr('biarr(,nsize)',-1,brarr(1,nsize),qb)
c         return
        end if
 120  continue
      
      return
      end
      
      
C     Evaluate Monte-Carlo estimates (active version)
      subroutine evalwpb2(brarr,ldbra,intvec,zi,ldzi,deltai,
     $     zit,ldzt,gamma,Eb,Ebb,ldeb,Eb4,ldeb4,Er,lder,Er2c,lder2c,
     $     Er2,lder2,Er3,lder3,Er4,lder4,Erbb,lderbb,Eb2r2,ldeb2r2,
     $     Erzb,lderzb,Ebzzb,ldebzzb,Ebzbb,ldebzbb,Ebzrr,ldebzrr,
     $     Eexpz,Eexpzb,ldexpzb,Eexz3btb,ldz3btb,Eexpzz,Eexpzr,ldexpzr,
     $     Eexpzzb,ldexpzzb,Eexpzrzb,ldzrzb,Eexpz4b2,ldz4b2,
     $     Eexpzbb,ldexpzbb,Eexpz2b3,ldz2b3,Eexpzrr,ldexpzrr,
     $     Eexz2br2,ldz2br2,Eexpz3b2,ldz3b2,tmpvec1,tmpmat1,ldtm1,
     $     tmpmat2,ldtm2)
c     total of 65 arguments (the maximum arguments allowed in R Fortran code)
      
      parameter (MAXQ=4)
      parameter (MAXQB=2*MAXQ)
      parameter (MAXQBA=MAXQB*MAXQB)
      integer ldbra,intvec(4),ldzi,deltai,ldzt,ldeb,ldeb4,lder,
     $     lder2c,lder2,lder3,lder4,lderbb,ldeb2r2,lderzb,ldebzzb,
     $     ldebzbb,ldebzrr,ldexpzb,ldz3btb,ldexpzr,ldexpzzb,ldzrzb,
     $     ldz4b2,ldexpzbb,ldz2b3,ldexpzrr,ldz2br2,ldz3b2,ldtm1,ldtm2
      double precision brarr(ldbra,*),zi(ldzi,*),zit(ldzt,*),gamma(*),
     $     Eb(*),Ebb(ldeb,*),Eb4(ldeb4,*),Er(lder,*),
     $     Er2c(lder2c,*),Er2(lder2,*),Er3(lder3,*),Er4(lder4,*),
     $     Erbb(lderbb,*),Eb2r2(ldeb2r2,*),Erzb(lderzb,*),
     $     Ebzzb(ldebzzb,*),Ebzbb(ldebzbb,*),Ebzrr(ldebzrr,*),
     $     Eexpz(*),Eexpzb(ldexpzb,*),Eexz3btb(ldz3btb,*),
     $     Eexpzz(*),Eexpzr(ldexpzr,*),Eexpzzb(ldexpzzb,*),
     $     Eexpzrzb(ldzrzb,*),Eexpz4b2(ldz4b2,*),Eexpzbb(ldexpzbb,*),
     $     Eexpz2b3(ldz2b3,*),Eexpzrr(ldexpzrr,*),Eexz2br2(ldz2br2,*),
     $     Eexpz3b2(ldz3b2,*),tmpvec1(*),tmpmat1(ldtm1,*),
     $     tmpmat2(ldtm2,*)
	  
	  integer qb,nsize,ni,tid
	  
c      return
c      end
      
      
C     Evaluate Monte-Carlo estimates (debug version)
c      subroutine evalwpbd(brarr,ldbra,qb,nsize,zi,ldzi,ni,deltai,tid,
c     $     zit,ldzt,gamma,Eb,Ebb,ldeb,Eb4,ldeb4,Er,lder,Er2c,lder2c,
c     $     Er2,lder2,Er3,lder3,Er4,lder4,Erbb,lderbb,Eb2r2,ldeb2r2,
c     $     Erzb,lderzb,Ebzzb,ldebzzb,Ebzbb,ldebzbb,Ebzrr,ldebzrr,
c     $     Eexpz,Eexpzb,ldexpzb,Eexz3btb,ldz3btb,Eexpzz,Eexpzr,ldexpzr,
c     $     Eexpzzb,ldexpzzb,Eexpzrzb,ldzrzb,Eexpz4b2,ldz4b2,
c     $     Eexpzbb,ldexpzbb,Eexpz2b3,ldz2b3,Eexpzrr,ldexpzrr,
c     $     Eexz2br2,ldz2br2,Eexpz3b2,ldz3b2,tmpvec1,tmpmat1,ldtm1,
c     $     tmpmat2,ldtm2)
c     total of 68 arguments

      
c      parameter (MAXQ=4)
c      parameter (MAXQB=2*MAXQ)
c      parameter (MAXQBA=MAXQB*MAXQB)
c      integer ldbra,qb,nsize,ldzi,ni,deltai,tid,ldzt,ldeb,ldeb4,lder,
c     $     lder2c,lder2,lder3,lder4,lderbb,ldeb2r2,lderzb,ldebzzb,
c     $     ldebzbb,ldebzrr,ldexpzb,ldz3btb,ldexpzr,ldexpzzb,ldzrzb,
c     $     ldz4b2,ldexpzbb,ldz2b3,ldexpzrr,ldz2br2,ldz3b2,ldtm1,ldtm2
c      double precision brarr(ldbra,*),zi(ldzi,*),zit(ldzt,*),gamma(*),
c     $     Eb(*),Ebb(ldeb,*),Eb4(ldeb4,*),Er(lder,*),
c     $     Er2c(lder2c,*),Er2(lder2,*),Er3(lder3,*),Er4(lder4,*),
c     $     Erbb(lderbb,*),Eb2r2(ldeb2r2,*),Erzb(lderzb,*),
c     $     Ebzzb(ldebzzb,*),Ebzbb(ldebzbb,*),Ebzrr(ldebzrr,*),
c     $     Eexpz(*),Eexpzb(ldexpzb,*),Eexz3btb(ldz3btb,*),
c     $     Eexpzz(*),Eexpzr(ldexpzr,*),Eexpzzb(ldexpzzb,*),
c     $     Eexpzrzb(ldzrzb,*),Eexpz4b2(ldz4b2,*),Eexpzbb(ldexpzbb,*),
c     $     Eexpz2b3(ldz2b3,*),Eexpzrr(ldexpzrr,*),Eexz2br2(ldz2br2,*),
c     $     Eexpz3b2(ldz3b2,*),tmpvec1(*),tmpmat1(ldtm1,*),
c     $     tmpmat2(ldtm2,*)
	 
c      return
      
      integer i,j,k,l,s,t,q,qba,qea,qecb,qni,qeani,nisq,tidp1,tid2u,
     $     idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8
      double precision tmpsum1,tmpsum2,tmpval1,tmpvec3(MAXQBA),
     $     bim(2,MAXQ),dexp
      
c      return
      qb = intvec(1)
      nsize = intvec(2)
      ni = intvec(3)
      tid = intvec(4)
      
      q = qb/2
      qba = qb**2
      qea = q**2
      qecb = q**3
      qni = q*ni
      qeani = qea*ni
      nisq = ni**2
      tidp1 = tid + 1
      tid2u = tid*(tid+1)/2
      
c      call intpr('qb',-1,qb,1)
c      call intpr('nsize',-1,nsize,1)
c      call intpr('ni',-1,ni,1)
c      call intpr('tid',-1,tid,1)

c      call intpr('q',-1,q,1)
c      call intpr('qba',-1,qba,1)
c      call intpr('qea',-1,qea,1)
c      call intpr('qecb',-1,qecb,1)
c      call intpr('qni',-1,qni,1)
c      call intpr('qeani',-1,qeani,1)
c      call intpr('nisq',-1,nisq,1)
c      call intpr('tidp1',-1,tidp1,1)
c      call intpr('tid2u',-1,tid2u,1)	  
c      return
      
c
c     Sum over all mc samples
c     
      
c     tid indepedent quantities
      
      call veczero(qb,Eb)
      call matzero(ldeb,qb,qb,Ebb)
      call matzero(ldeb4,qba,qba,Eb4)
      call matzero(lder,q,ni,Er)
      call matzero(lder2c,qea,nisq,Er2c)
      call matzero(lder3,qecb,nisq,Er3)
      call matzero(lder4,qea**2,nisq,Er4)
      call matzero(lderbb,q*qba,ni,Erbb)
      call matzero(ldeb2r2,qba*qea,ni,Eb2r2)
c      return
	  
      do 60 k=1,nsize
        do 10 l=1,qb
          Eb(l) = Eb(l) + brarr(l,k)
          do 10 s=1,l
            idx1 = (l-1)*qb + s
            tmpvec3(idx1)=brarr(s,k)*brarr(l,k)
            Ebb(s,l)=Ebb(s,l)+tmpvec3(idx1)
            if (s.ne.l) tmpvec3((s-1)*qb+l)=tmpvec3(idx1)
 10     continue
        
        do 20 l=1,qba
          do 20 s=1,l
            Eb4(s,l)=Eb4(s,l)+tmpvec3(s)*tmpvec3(l)
 20     continue
c        return
		
        do 50 j=1,ni
          do 50 l=1,q
            idx1 = qb + (j-1)*q + l
            Er(l,j)=Er(l,j)+brarr(idx1,k)
            
            do 30 s=1,qb
              do 30 i=1,s
                idx2 = (s-1)*qb + i
                idx3 = (idx2-1)*q + l
                Erbb(idx3,j)=Erbb(idx3,j)+brarr(idx1,k)*tmpvec3(idx2)
                do 30 t=1,l
                  idx4 = ((l-1)*q+t-1)*qba + idx2
                  Eb2r2(idx4,j)=Eb2r2(idx4,j)+
     $                 tmpvec3(idx2)*brarr(qb+(j-1)*q+t,k)*brarr(idx1,k)
 30         continue
            
            do 48 s=1,ni
              idx2 = (j-1)*ni + s
              if (s.eq.j) then
                do 42 t=1,l
                  idx3 = qb + (s-1)*q + t
                  idx4 = (l-1)*q + t
                  Er2c(idx4,idx2)=Er2c(idx4,idx2)+
     $                 brarr(idx3,k)*brarr(idx1,k)                
 42             continue
              else
                do 45 t=1,q
                  idx3 = qb + (s-1)*q + t
                  idx4 = (l-1)*q + t
                  Er2c(idx4,idx2)=Er2c(idx4,idx2)+
     $                 brarr(idx3,k)*brarr(idx1,k)                
 45             continue
              end if  
 48         continue
            
            do 50 i=1,l
              idx2 = (l-1)*q + i
              tmpmat1(idx2,j) = brarr(idx1,k)*brarr(qb+(j-1)*q+i,k)
              do 50 s=1,ni
                do 50 t=1,q
                  idx3 = (l-1)*qea + (i-1)*q + t
                  idx4 = (j-1)*ni + s
                  Er3(idx3,idx4)=Er3(idx3,idx4)+
     $                 brarr(qb+(s-1)*q+t,k)*tmpmat1(idx2,j)                  
 50     continue                  
c        return
		
        do 60 idx1=1,ni
          do 60 l=1,q
            do 60 s=1,l
              idx3 = (l-1)*q + s
              
              do 60 idx2=1,ni
                if (idx2.eq.idx1) then
                  do 52 t=1,l
                    if (t.lt.l) then
                      j = t
                    else
                      j = s
                    end if
                    
                    do 52 i=1,j
                      idx4 = (t-1)*q + i
                      
                      idx5 = (idx1-1)*ni + idx2
                      idx6 = (idx3-1)*qea + idx4
                      Er4(idx6,idx5)=Er4(idx6,idx5)+
     $                     tmpmat1(idx4,idx2)*tmpmat1(idx3,idx1)
 52               continue				
                else
                  do 55 t=1,q
                    do 55 i=1,t
                      idx4 = (t-1)*q + i
                      
                      idx5 = (idx1-1)*ni + idx2
                      idx6 = (idx3-1)*qea + idx4
                      Er4(idx6,idx5)=Er4(idx6,idx5)+
     $                     tmpmat1(idx4,idx2)*tmpmat1(idx3,idx1)
 55               continue                  
                end if				
 60   continue    
c      call dblepr('Eb',-1,Eb,qb)	  
c      call dblepr('Ebb(,1)',-1,Ebb(1,1),qb)	  
c      call dblepr('Ebb(,4)',-1,Ebb(1,4),qb)	  
c      return
      
	  
c     tid dependent quantities   
      
      call matzero(lderzb,qea,ni,Erzb)
      call matzero(ldebzzb,q,q,Ebzzb)
      call matzero(ldebzbb,q,qba,Ebzbb)
      call matzero(ldebzrr,qecb,ni,Ebzrr)
      
      if (tid.gt.0) then         
        call veczero(tid,Eexpz)
        call matzero(ldexpzb,qb,tid,Eexpzb)
        call matzero(ldz3btb,qea,tid,Eexz3btb)
        call veczero(tid2u,Eexpzz)
        call matzero(ldexpzr,qni,tid,Eexpzr)
        call matzero(ldexpzzb,qb,tid2u,Eexpzzb)
        call matzero(ldzrzb,qeani,tid,Eexpzrzb)
        call matzero(ldz4b2,qea,tid2u,Eexpz4b2)
        call matzero(ldexpzbb,qba,tid,Eexpzbb)
        call matzero(ldz2b3,q*qba,tid,Eexpz2b3)
        call matzero(ldexpzrr,qeani,tid,Eexpzrr)
        call matzero(ldz2br2,qecb*ni,tid,Eexz2br2)
        
        call matzero(ldz3b2,qea,tid,Eexpz3b2)
        
        do 240 k=1,nsize
          do 80 l=1,q
            do 80 s=1,2
              bim(s,l) = brarr((l-1)*2+s,k)
 80       continue
          
          if (deltai.eq.1) then
            do 100 l=1,q
              tmpsum1 = 0.0d0
              do 90 s=1,2
                tmpsum1 = tmpsum1 + zit(tidp1,s)*bim(s,l)
 90           continue
              tmpmat2(l,tidp1) = tmpsum1
 100        continue            
          end if
          
c          idx1 = 0
          do 200 l=1,tid
            tmpsum1 = 0.0d0 
            do 120 s=1,q
              tmpsum2 = 0.0d0
              do 110 t=1,2
                tmpsum2 = tmpsum2 + zit(l,t)*bim(t,s)
 110          continue
              tmpmat2(s,l) = tmpsum2
              tmpsum1 = tmpsum1 + tmpsum2*gamma(s)
 120        continue
            tmpvec1(l) = dexp(tmpsum1)
            Eexpz(l)=Eexpz(l)+tmpvec1(l)
            
            do 130 s=1,qb
              Eexpzb(s,l)=Eexpzb(s,l)+tmpvec1(l)*brarr(s,k)
              do 130 t=1,s
                idx2 = (s-1)*qb + t
                tmpval1 = brarr(t,k)*brarr(s,k)
                Eexpzbb(idx2,l)=Eexpzbb(idx2,l)+tmpval1*tmpvec1(l)              
                do 130 i=1,q
                  idx3 = (idx2-1)*q + i
                  Eexpz2b3(idx3,l)=Eexpz2b3(idx3,l)+
     $                 tmpval1*tmpvec1(l)*tmpmat2(i,l)
 130        continue
            
            do 150 s=1,q
              do 140 t=1,s
                idx2 = (s-1)*q + t
                Eexz3btb(idx2,l)=Eexz3btb(idx2,l)+
     $               tmpvec1(l)*tmpmat2(t,l)*tmpmat2(s,l)
 140          continue
              
              do 150 j=1,ni
                do 150 t=1,q
                  idx2 = (j-1)*qea + (s-1)*q + t
                  idx3 = qb + (j-1)*q + t
                  Eexpzrzb(idx2,l)=Eexpzrzb(idx2,l)+
     $                 tmpvec1(l)*brarr(idx3,k)*tmpmat2(s,l)
                  
                  tmpval1=brarr(idx3,k)*brarr(qb+(j-1)*q+s,k)
                  Eexpzrr(idx2,l)=Eexpzrr(idx2,l)+tmpval1*tmpvec1(l)
                  
                  if (t.le.s) then
                    do 145 i=1,q
                      idx4 = (j-1)*qecb+(s-1)*qea+(t-1)*q+i
                      Eexz2br2(idx4,l)=Eexz2br2(idx4,l)+
     $                     tmpval1*tmpvec1(l)*tmpmat2(i,l)
 145                continue
                  end if
 150        continue
            
            do 170 s=1,l
c              idx1 = idx1 + 1
              idx1 = l*(l-1)/2 + s
              tmpval1 = tmpvec1(s)*tmpvec1(l)
              Eexpzz(idx1)=Eexpzz(idx1)+tmpval1
              
              do 160 t=1,qb
                Eexpzzb(t,idx1)=Eexpzzb(t,idx1)+tmpval1*brarr(t,k)
 160          continue
              
              if (s.eq.l) then
                do 162 t=1,q
                  do 162 i=1,t
                    idx2 = (t-1)*q + i
                    Eexpz4b2(idx2,idx1)=Eexpz4b2(idx2,idx1)+
     $                   tmpval1*tmpmat2(i,s)*tmpmat2(t,l)	 
 162            continue			  
              else
                do 165 t=1,q
                  do 165 i=1,q
                    idx2 = (t-1)*q + i
                    Eexpz4b2(idx2,idx1)=Eexpz4b2(idx2,idx1)+
     $                   tmpval1*tmpmat2(i,s)*tmpmat2(t,l)	 
 165            continue
              end if	 
 170        continue
            
            do 180 s=1,qni
              Eexpzr(s,l)=Eexpzr(s,l)+brarr(qb+s,k)*tmpvec1(l)
 180        continue
            
            
            if (deltai.eq.1) then
              do 190 s=1,q
                do 190 t=1,q
                  idx2 = (s-1)*q + t
                  Eexpz3b2(idx2,l)=Eexpz3b2(idx2,l)+
     $                 tmpvec1(l)*tmpmat2(t,tidp1)*tmpmat2(s,l)
 190          continue               
            end if
 200      continue           
          
          
          if (deltai.eq.1) then
            do 210 j=1,ni
              do 210 l=1,q
                do 210 s=1,q
                  idx2 = (l-1)*q + s
                  idx3 = qb + (j-1)*q + s
                  Erzb(idx2,j)=Erzb(idx2,j)+
     $                 brarr(idx3,k)*tmpmat2(l,tidp1)
                  do 210 t=1,s
                    idx4 = (s-1)*qea + (t-1)*q + l
                    Ebzrr(idx4,j)=Ebzrr(idx4,j)+tmpmat2(l,tidp1)*
     $                   brarr(qb+(j-1)*q+t,k)*brarr(idx3,k)
 210        continue
            
            do 230 l=1,q
              do 220 s=1,l
                Ebzzb(s,l)=Ebzzb(s,l)+tmpmat2(s,tidp1)*tmpmat2(l,tidp1)
 220          continue
              
              do 230 s=1,qb
                do 230 t=1,s
                  idx2 = (s-1)*qb + t
                  Ebzbb(l,idx2)=Ebzbb(l,idx2)+
     $                 tmpmat2(l,tidp1)*brarr(t,k)*brarr(s,k)
 230        continue
          end if      
 240    continue      
      end if
c      return
      
	  
c      
c     divide the sums by the number of mc samples
c
      
c     tid independent quantities
      
      do 250 l=1,qb
        Eb(l) = Eb(l)/nsize
        do 250 s=1,l
          idx1 = (l-1)*qb + s
          Ebb(s,l)=Ebb(s,l)/nsize
          if (s.ne.l) Ebb(l,s)=Ebb(s,l)
 250  continue
      
      do 260 l=1,qba
        do 260 s=1,l
          Eb4(s,l)=Eb4(s,l)/nsize
          if (s.ne.l) Eb4(l,s)=Eb4(s,l)
 260  continue
      
      do 290 j=1,ni
        do 290 l=1,q
          Er(l,j)=Er(l,j)/nsize
          
          do 270 s=1,qb
            do 270 i=1,s
              idx1 = (s-1)*qb + i
              idx2 = (idx1-1)*q + l
              Erbb(idx2,j)=Erbb(idx2,j)/nsize
			  
              idx3 = (i-1)*qb + s
              if (i.ne.s) Erbb((idx3-1)*q+l,j)=Erbb(idx2,j)
              
              do 270 t=1,l
                idx4 = (l-1)*q + t - 1
                idx5 = idx4*qba + idx1
                Eb2r2(idx5,j)=Eb2r2(idx5,j)/nsize
				
                idx6 = (t-1)*q + l - 1
                if (i.ne.s) Eb2r2(idx4*qba+idx3,j)=Eb2r2(idx5,j)
                if (t.ne.l) Eb2r2(idx6*qba+idx1,j)=Eb2r2(idx5,j)
                if ((i.ne.s).and.(t.ne.l)) Eb2r2(idx6*qba+idx3,j)=
     $               Eb2r2(idx5,j)				        
 270      continue
          
          
          do 280 s=1,ni
            idx1 = (j-1)*ni + s
			
            if (s.eq.j) then
              do 272 t=1,l
                idx2 = (l-1)*q + t
                Er2c(idx2,idx1)=Er2c(idx2,idx1)/nsize
                Er2(idx2,j)=Er2c(idx2,idx1)
				
                if (t.ne.l) then
                  idx3 = (t-1)*q + l
                  Er2c(idx3,idx1)=Er2c(idx2,idx1)
                  Er2(idx3,j)=Er2c(idx3,idx1)
                end if
 272          continue
            else
              do 275 t=1,q
                idx2 = (l-1)*q + t
                Er2c(idx2,idx1)=Er2c(idx2,idx1)/nsize
 275          continue				
            end if
 280      continue
          
          do 290 i=1,l
            do 290 s=1,ni
              do 290 t=1,q
                idx1 = (l-1)*qea + (i-1)*q + t
                idx2 = (j-1)*ni + s
                Er3(idx1,idx2)=Er3(idx1,idx2)/nsize
                if (i.ne.l) Er3((i-1)*qea+(l-1)*q+t,idx2)=Er3(idx1,idx2)
 290  continue
      
	  
      do 300 idx1=1,ni
        do 300 l=1,q
          do 300 s=1,l
            idx3 = (l-1)*q + s
            idx7 = (s-1)*q + l
            
            do 300 idx2=1,ni
              if (idx2.eq.idx1) then
                do 292 t=1,l
                  if (t.lt.l) then
                    j = t
                  else
c				    t.eq.l
                    j = s
                  end if
			      
                  do 292 i=1,j
                    idx4 = (t-1)*q + i
                    
                    idx5 = (idx1-1)*ni + idx2
                    idx6 = (idx3-1)*qea + idx4
                    Er4(idx6,idx5)=Er4(idx6,idx5)/nsize
					
                    if ((t.ne.l).or.(i.ne.s)) 
     $                Er4((idx4-1)*qea+idx3,idx5)=Er4(idx6,idx5)
 292            continue 
              else
                do 295 t=1,q
                  do 295 i=1,t
                    idx4 = (t-1)*q + i
                    
                    idx5 = (idx1-1)*ni + idx2
                    idx6 = (idx3-1)*qea + idx4
                    Er4(idx6,idx5)=Er4(idx6,idx5)/nsize
                    
                    idx8 = (i-1)*q + t
                    if (s.ne.l) Er4((idx7-1)*qea+idx4,idx5)=
     $                   Er4(idx6,idx5)
                    if (i.ne.t) Er4((idx3-1)*qea+idx8,idx5)=
     $                   Er4(idx6,idx5)
                    if ((s.ne.l).and.(i.ne.t)) Er4((idx7-1)*qea+idx8,
     $                   idx5)=Er4(idx6,idx5)
 295            continue                  
              end if				
 300  continue

      do 305 idx1=1,ni
        do 305 l=1,q
          do 305 s=1,l
            idx3 = (l-1)*q + s
            idx7 = (s-1)*q + l
            
            do 305 t=1,q
              do 305 i=1,t
                idx4 = (t-1)*q + i
                
                idx5 = (idx1-1)*ni + idx1
                idx6 = (idx3-1)*qea + idx4
                
                idx8 = (i-1)*q + t
                if (s.ne.l) Er4((idx7-1)*qea+idx4,idx5)=
     $                   Er4(idx6,idx5)
                if (i.ne.t) Er4((idx3-1)*qea+idx8,idx5)=
     $                   Er4(idx6,idx5)
                if ((s.ne.l).and.(i.ne.t)) Er4((idx7-1)*qea+idx8,
     $                   idx5)=Er4(idx6,idx5)
 305   continue                  

c      call dblepr('Eb',-1,Eb,qb)	  
c      call dblepr('Ebb(,1)',-1,Ebb(1,1),qb)	  
c      call dblepr('Ebb(,4)',-1,Ebb(1,4),qb)	  
c      call dblepr('Eb4(,1)',-1,Eb4(1,1),16)	  
c      call dblepr('Eb4(,16)',-1,Eb4(1,16),16)	  
c      call dblepr('Er(,1)',-1,Er(1,1),q)	  
c      call dblepr('Er(,4)',-1,Er(1,4),q)	  
c      call dblepr('Erbb(,1)',-1,Erbb(1,1),32)	  
c      call dblepr('Erbb(,4)',-1,Erbb(1,4),32)	  
c      call dblepr('Eb2r2(,1)',-1,Eb2r2(1,1),64)	  
c      call dblepr('Eb2r2(,4)',-1,Eb2r2(1,4),64)	  
c      call dblepr('Er2c(,1)',-1,Er2c(1,1),4) 
c      call dblepr('Er2c(,16)',-1,Er2c(1,16),4)
c      call dblepr('Er2(,1)',-1,Er2(1,1),4)
c      call dblepr('Er2(,4)',-1,Er2(1,4),4)
c      call dblepr('Er3(,1)',-1,Er3(1,1),8)
c      call dblepr('Er3(,16)',-1,Er3(1,16),8)
c      call dblepr('Er4(,1)',-1,Er4(1,1),16)
c      call dblepr('Er4(,2)',-1,Er4(1,2),16)
c      call dblepr('Er4(,3)',-1,Er4(1,3),16)
c      call dblepr('Er4(,4)',-1,Er4(1,4),16)
c      call dblepr('Er4(,5)',-1,Er4(1,5),16)
c      call dblepr('Er4(,6)',-1,Er4(1,6),16)
c      call dblepr('Er4(,7)',-1,Er4(1,7),16)
c      call dblepr('Er4(,8)',-1,Er4(1,8),16)
c      call dblepr('Er4(,9)',-1,Er4(1,9),16)
c      call dblepr('Er4(,10)',-1,Er4(1,10),16)
c      call dblepr('Er4(,11)',-1,Er4(1,11),16)
c      call dblepr('Er4(,12)',-1,Er4(1,12),16)
c      call dblepr('Er4(,13)',-1,Er4(1,13),16)
c      call dblepr('Er4(,14)',-1,Er4(1,14),16)
c      call dblepr('Er4(,15)',-1,Er4(1,15),16)
c      call dblepr('Er4(,16)',-1,Er4(1,16),16)
c      return
      
      
c     tid dependent quantities
      
      if (tid.gt.0) then
c        idx1 = 0
        do 380 l=1,tid
          Eexpz(l)=Eexpz(l)/nsize
          
          do 310 s=1,qb
            Eexpzb(s,l)=Eexpzb(s,l)/nsize
            do 310 t=1,s
              idx2 = (s-1)*qb + t
              Eexpzbb(idx2,l)=Eexpzbb(idx2,l)/nsize
              idx3 = (t-1)*qb + s
              if (t.ne.s) Eexpzbb(idx3,l)=Eexpzbb(idx2,l)
              do 310 i=1,q
                idx4 = (idx2-1)*q + i
                Eexpz2b3(idx4,l)=Eexpz2b3(idx4,l)/nsize                
                if (t.ne.s) Eexpz2b3((idx3-1)*q+i,l)=Eexpz2b3(idx4,l)
 310      continue
          
          do 330 s=1,q
            do 320 t=1,s
              idx2 = (s-1)*q + t
              Eexz3btb(idx2,l)=Eexz3btb(idx2,l)/nsize
              if (t.ne.s) Eexz3btb((t-1)*q+s,l)=Eexz3btb(idx2,l)
 320        continue
            
            do 330 j=1,ni
              do 330 t=1,q
                idx2 = (j-1)*qea + (s-1)*q + t
                Eexpzrzb(idx2,l)=Eexpzrzb(idx2,l)/nsize
                
                Eexpzrr(idx2,l)=Eexpzrr(idx2,l)/nsize
                
                if (t.le.s) then
                  do 325 i=1,q
                    idx3 = (j-1)*qecb+(s-1)*qea+(t-1)*q+i
                    Eexz2br2(idx3,l)=Eexz2br2(idx3,l)/nsize
					
                    idx4 = (j-1)*qecb+(t-1)*qea+(s-1)*q+i
                    if (t.ne.s) Eexz2br2(idx4,l)=Eexz2br2(idx3,l)
 325              continue
                end if 
 330      continue
          
          do 350 s=1,l
c            idx1 = idx1 + 1
            idx1 = l*(l-1)/2 + s
            Eexpzz(idx1)=Eexpzz(idx1)/nsize
            
            do 340 t=1,qb
              Eexpzzb(t,idx1)=Eexpzzb(t,idx1)/nsize
 340        continue
            
            if (s.eq.l) then
              do 342 t=1,q
                do 342 i=1,t
                  idx2 = (t-1)*q + i
                  Eexpz4b2(idx2,idx1)=Eexpz4b2(idx2,idx1)/nsize
                  if (i.ne.t) Eexpz4b2((i-1)*q+t,idx1)=
     $                 Eexpz4b2(idx2,idx1)
 342          continue			  
            else
              do 345 t=1,q
                do 345 i=1,q
                  idx2 = (t-1)*q + i
                  Eexpz4b2(idx2,idx1)=Eexpz4b2(idx2,idx1)/nsize
 345          continue
            end if
 350      continue
          
          
          do 360 s=1,qni
            Eexpzr(s,l)=Eexpzr(s,l)/nsize
 360      continue
          
          
          if (deltai.eq.1) then
            do 370 s=1,q
              do 370 t=1,q
                idx2 = (s-1)*q + t
                Eexpz3b2(idx2,l)=Eexpz3b2(idx2,l)/nsize
 370        continue               
          end if
 380    continue           
c        call dblepr('Eexpz',-1,Eexpz,tid)
c        return
        
		
        if (deltai.eq.1) then
          do 390 j=1,ni
            do 390 l=1,q
              do 390 s=1,q
                idx2 = (l-1)*q + s
                Erzb(idx2,j)=Erzb(idx2,j)/nsize
                do 390 t=1,l
                  idx3 = (l-1)*qea + (t-1)*q + s
                  Ebzrr(idx3,j)=Ebzrr(idx3,j)/nsize
                  if (t.ne.l) Ebzrr((t-1)*qea+idx2,j)=
     $                 Ebzrr(idx3,j)
 390      continue
          
          do 410 l=1,q
            do 400 s=1,l
              Ebzzb(s,l)=Ebzzb(s,l)/nsize
              if (s.ne.l) Ebzzb(l,s)=Ebzzb(s,l)
 400        continue
            
            do 410 s=1,qb
              do 410 t=1,s
                idx2 = (s-1)*qb + t
                Ebzbb(l,idx2)=Ebzbb(l,idx2)/nsize
                if (t.ne.s) Ebzbb(l,(t-1)*qb+s)=Ebzbb(l,idx2)
  410     continue
        end if      
      end if
c      call dblepr('Erbb(,1)',-1,Erbb(1,1),32)
c      call dblepr('Erbb(,4)',-1,Erbb(1,4),32)
c      call dblepr('Eb2r2(,1)',-1,Eb2r2(1,1),64)
c      call dblepr('Eb2r2(,4)',-1,Eb2r2(1,4),64)
c      call dblepr('Erzb(,1)',-1,Erzb(1,1),4)
c      call dblepr('Erzb(,4)',-1,Erzb(1,4),4)
c      call dblepr('Ebzzb(,1)',-1,Ebzzb(1,1),2)
c      call dblepr('Ebzzb(,2)',-1,Ebzzb(1,2),2)
c      call dblepr('Ebzbb(,1)',-1,Ebzbb(1,1),2)
c      call dblepr('Ebzbb(,16)',-1,Ebzbb(1,16),2)
c      call dblepr('Ebzrr(,1)',-1,Ebzrr(1,1),8)
c      call dblepr('Ebzrr(,4)',-1,Ebzrr(1,4),8)
c      call dblepr('Eexpzb(,1)',-1,Eexpzb(1,1),4)
c      call dblepr('Eexpzb(,23)',-1,Eexpzb(1,23),4)
c      call dblepr('Eexz3btb(,1)',-1,Eexz3btb(1,1),4)
c      call dblepr('Eexz3btb(,23)',-1,Eexz3btb(1,23),4)
c      call dblepr('Eexpzz',-1,Eexpzz,276)
c      call dblepr('Eexpzr(,1)',-1,Eexpzr(1,1),8)
c      call dblepr('Eexpzr(,23)',-1,Eexpzr(1,23),8)
c      call dblepr('Eexpzzb(,1)',-1,Eexpzzb(1,1),4)
c      call dblepr('Eexpzzb(,276)',-1,Eexpzzb(1,276),4)
c      call dblepr('Eexpzrzb(,1)',-1,Eexpzrzb(1,1),16)
c      call dblepr('Eexpzrzb(,23)',-1,Eexpzrzb(1,23),16)
c      call dblepr('Eexpz4b2(,1)',-1,Eexpz4b2(1,1),4)
c      call dblepr('Eexpz4b2(,276)',-1,Eexpz4b2(1,276),4)
c      call dblepr('Eexpzbb(,1)',-1,Eexpzbb(1,1),16)
c      call dblepr('Eexpzbb(,23)',-1,Eexpzbb(1,23),16)
c      call dblepr('Eexpz2b3(,1)',-1,Eexpz2b3(1,1),32)
c      call dblepr('Eexpz2b3(,23)',-1,Eexpz2b3(1,23),32)
c      call dblepr('Eexpzrr(,1)',-1,Eexpzrr(1,1),16)
c      call dblepr('Eexpzrr(,23)',-1,Eexpzrr(1,23),16)
c      call dblepr('Eexz2br2(,1)',-1,Eexz2br2(1,1),32)
c      call dblepr('Eexz2br2(,23)',-1,Eexz2br2(1,23),32)
c      call dblepr('Eexpz3b2(,1)',-1,Eexpz3b2(1,1),4)
c      call dblepr('Eexpz3b2(,23)',-1,Eexpz3b2(1,23),4)
      
      return
      end
      
      
C     Get gradient and hessian contribution from subject i (active version)
      subroutine getghcpd(deltai,vibeta,h0t,tid,xit,ldvec,px,zit,
     $     alphamat,q,gamma,vi,pv,zitz,xi,ni,
     $     sigbinv,sigeinv,Eb,Ebb,Eb4,Er,
     $     Er2c,Er2,Er3,Er4,Erbb,
     $     Eb2r2,Erzb,Ebzzb,Ebzbb,
     $     Ebzrr,Eexpz,Eexpzb,Eexz3btb,
     $     Eexpzz,Eexpzr,Eexpzzb,Eexpzrzb,
     $     Eexpz4b2,Eexpzbb,Eexpz2b3,
     $     Eexpzrr,Eexz2br2,Eexpz3b2,
     $     vecwb,krwb2,vecwe,krwe2,dqb,dqe,
     $     tqbqb,tqeqe,retGrad,ne1,retInfo,retCpd,
     $     tmpvec1,tmpvec2,xitalp)
c     total of 60 arguments (5 fewer than the maximum arguments allowed in R Fortran code)
      
      parameter (MAXQ=4,MAXPX=50)
      parameter (MAXQB=2*MAXQ,MAXPXTQ=MAXPX*MAXQ,MAXQEA=MAXQ**2)
      parameter (MAXQBA=MAXQB*MAXQB)
      
c      integer deltai,tid,ldxt,px,ldzt,ldam,q,pv,ldzz,ldxi,ni,ldsb,ldse,
c     $     ldeb,ldeb4,lder,lder2c,lder2,lder3,lder4,lderbb,ldeb2r2,
c     $     lderzb,ldebzzb,ldebzbb,ldebzrr,ldexpzb,ldz3btb,ldexpzr,
c     $     ldexpzzb,ldzrzb,ldz4b2,ldexpzbb,ldz2b3,ldexpzrr,ldz2br2,
c     $     ldz3b2,ldkb,ldke,lddqb,lddqe,ldtqb,ldtqe,ne1,ldif,ldcp,
c     $     ldxa
c      double precision vibeta,h0t(*),xit(ldxt,*),zit(ldzt,*),
c     $     alphamat(ldam,*),gamma(*),vi(*),zitz(ldzz,*),xi(ldxi,*),
c     $     sigbinv(ldsb,*),sigeinv(ldse,*),Eb(*),Ebb(ldeb,*),
c     $     Eb4(ldeb4,*),Er(lder,*),Er2c(lder2c,*),Er2(lder2,*),
c     $     Er3(lder3,*),Er4(lder4,*),Erbb(lderbb,*),Eb2r2(ldeb2r2,*),
c     $     Erzb(lderzb,*),Ebzzb(ldebzzb,*),Ebzbb(ldebzbb,*),
c     $     Ebzrr(ldebzrr,*),Eexpz(*),Eexpzb(ldexpzb,*),
c     $     Eexz3btb(ldz3btb,*),Eexpzz(*),Eexpzr(ldexpzr,*),
c     $     Eexpzzb(ldexpzzb,*),Eexpzrzb(ldzrzb,*),Eexpz4b2(ldz4b2,*),
c     $     Eexpzbb(ldexpzbb,*),Eexpz2b3(ldz2b3,*),Eexpzrr(ldexpzrr,*),
c     $     Eexz2br2(ldz2br2,*),Eexpz3b2(ldz3b2,*),vecwb(*),
c     $     krwb2(ldkb,*),vecwe(*),krwe2(ldke,*),dqb(lddqb,*),
c     $     dqe(lddqe,*),tqbqb(ldtqb,*),tqeqe(ldtqe,*),retGrad(*),
c     $     retInfo(ldif,*),retCpd(ldcp,*),tmpvec1(*),tmpvec2(*),
c     $     xitalp(ldxa,*)
      
      integer deltai,tid,ldvec(40),px,q,pv,ni,ne1
      double precision vibeta,h0t(*),xit(ldvec(1),*),zit(ldvec(2),*),
     $     alphamat(ldvec(3),*),gamma(*),vi(*),zitz(ldvec(4),*),
     $     xi(ldvec(5),*),sigbinv(ldvec(6),*),sigeinv(ldvec(7),*),
     $     Eb(*),Ebb(ldvec(8),*),Eb4(ldvec(9),*),Er(ldvec(10),*),
     $     Er2c(ldvec(11),*),Er2(ldvec(12),*),Er3(ldvec(13),*),
     $     Er4(ldvec(14),*),Erbb(ldvec(15),*),Eb2r2(ldvec(16),*),
     $     Erzb(ldvec(17),*),Ebzzb(ldvec(18),*),Ebzbb(ldvec(19),*),
     $     Ebzrr(ldvec(20),*),Eexpz(*),Eexpzb(ldvec(21),*),
     $     Eexz3btb(ldvec(22),*),Eexpzz(*),Eexpzr(ldvec(23),*),
     $     Eexpzzb(ldvec(24),*),Eexpzrzb(ldvec(25),*),
     $     Eexpz4b2(ldvec(26),*),Eexpzbb(ldvec(27),*),
     $     Eexpz2b3(ldvec(28),*),Eexpzrr(ldvec(29),*),
     $     Eexz2br2(ldvec(30),*),Eexpz3b2(ldvec(31),*),vecwb(*),
     $     krwb2(ldvec(32),*),vecwe(*),krwe2(ldvec(33),*),
     $     dqb(ldvec(34),*),dqe(ldvec(35),*),tqbqb(ldvec(36),*),
     $     tqeqe(ldvec(37),*),retGrad(*),retInfo(ldvec(38),*),
     $     retCpd(ldvec(39),*),tmpvec1(*),tmpvec2(*),
     $     xitalp(ldvec(40),*)
      
c      call intpr('ldvec',-1,ldvec(1),40)
c      call dblepr('xit(,1)',-1,xit(1,1),24)
c      call dblepr('xit(,3)',-1,xit(1,3),24)
c      call dblepr('zit(,1)',-1,zit(1,1),24)
c      call dblepr('zit(,2)',-1,zit(1,2),24)
      
c      return
c      end
      
      
C     Get gradient and hessian contribution from subject i (debug version)
c      subroutine getghcpm(deltai,vibeta,h0t,tid,xit,ldxt,px,zit,ldzt,
c     $     alphamat,ldam,q,gamma,vi,pv,zitz,ldzz,xi,ldxi,ni,
c     $     sigbinv,ldsb,sigeinv,ldse,Eb,Ebb,ldeb,Eb4,ldeb4,Er,lder,
c     $     Er2c,lder2c,Er2,lder2,Er3,lder3,Er4,lder4,Erbb,lderbb,
c     $     Eb2r2,ldeb2r2,Erzb,lderzb,Ebzzb,ldebzzb,Ebzbb,ldebzbb,
c     $     Ebzrr,ldebzrr,Eexpz,Eexpzb,ldexpzb,Eexz3btb,ldz3btb,
c     $     Eexpzz,Eexpzr,ldexpzr,Eexpzzb,ldexpzzb,Eexpzrzb,ldzrzb,
c     $     Eexpz4b2,ldz4b2,Eexpzbb,ldexpzbb,Eexpz2b3,ldz2b3,
c     $     Eexpzrr,ldexpzrr,Eexz2br2,ldz2br2,Eexpz3b2,ldz3b2,
c     $     vecwb,krwb2,ldkb,vecwe,krwe2,ldke,dqb,lddqb,dqe,lddqe,
c     $     tqbqb,ldtqb,tqeqe,ldtqe,retGrad,ne1,retInfo,ldif,retCpd,ldcp,
c     $     tmpvec1,tmpvec2,xitalp,ldxa)
      
c      parameter (MAXQ=4,MAXPX=50)
c      parameter (MAXQB=2*MAXQ,MAXPXTQ=MAXPX*MAXQ,MAXQEA=MAXQ**2)
c      parameter (MAXQBA=MAXQB*MAXQB)
      
c      integer deltai,tid,ldxt,px,ldzt,ldam,q,pv,ldzz,ldxi,ni,ldsb,ldse,
c     $     ldeb,ldeb4,lder,lder2c,lder2,lder3,lder4,lderbb,ldeb2r2,
c     $     lderzb,ldebzzb,ldebzbb,ldebzrr,ldexpzb,ldz3btb,ldexpzr,
c     $     ldexpzzb,ldzrzb,ldz4b2,ldexpzbb,ldz2b3,ldexpzrr,ldz2br2,
c     $     ldz3b2,ldkb,ldke,lddqb,lddqe,ldtqb,ldtqe,ne1,ldif,ldcp,
c     $     ldxa
c      double precision vibeta,h0t(*),xit(ldxt,*),zit(ldzt,*),
c     $     alphamat(ldam,*),gamma(*),vi(*),zitz(ldzz,*),xi(ldxi,*),
c     $     sigbinv(ldsb,*),sigeinv(ldse,*),Eb(*),Ebb(ldeb,*),
c     $     Eb4(ldeb4,*),Er(lder,*),Er2c(lder2c,*),Er2(lder2,*),
c     $     Er3(lder3,*),Er4(lder4,*),Erbb(lderbb,*),Eb2r2(ldeb2r2,*),
c     $     Erzb(lderzb,*),Ebzzb(ldebzzb,*),Ebzbb(ldebzbb,*),
c     $     Ebzrr(ldebzrr,*),Eexpz(*),Eexpzb(ldexpzb,*),
c     $     Eexz3btb(ldz3btb,*),Eexpzz(*),Eexpzr(ldexpzr,*),
c     $     Eexpzzb(ldexpzzb,*),Eexpzrzb(ldzrzb,*),Eexpz4b2(ldz4b2,*),
c     $     Eexpzbb(ldexpzbb,*),Eexpz2b3(ldz2b3,*),Eexpzrr(ldexpzrr,*),
c     $     Eexz2br2(ldz2br2,*),Eexpz3b2(ldz3b2,*),vecwb(*),
c     $     krwb2(ldkb,*),vecwe(*),krwe2(ldke,*),dqb(lddqb,*),
c     $     dqe(lddqe,*),tqbqb(ldtqb,*),tqeqe(ldtqe,*),retGrad(*),
c     $     retInfo(ldif,*),retCpd(ldcp,*),tmpvec1(*),tmpvec2(*),
c     $     xitalp(ldxa,*)
      
c      integer ldxt,ldzt,ldam,ldzz,ldxi,ldsb,ldse,
c     $     ldeb,ldeb4,lder,lder2c,lder2,lder3,lder4,lderbb,ldeb2r2,
c     $     lderzb,ldebzzb,ldebzbb,ldebzrr,ldexpzb,ldz3btb,ldexpzr,
c     $     ldexpzzb,ldzrzb,ldz4b2,ldexpzbb,ldz2b3,ldexpzrr,ldz2br2,
c     $     ldz3b2,ldkb,ldke,lddqb,lddqe,ldtqb,ldtqe,ldif,ldcp,ldxa
      integer ldif,ldcp
      
      integer pxtq,nfixpar,qb,qbu,qeu,ntotpar,qba,qea,nisq,qbcb,qecb,
     $     qe4p,tidp1,i,j,k,l,s,t,ldtm2,ldtm3,ldtm4,ldtm5,offset1,
     $     offset2,offset3,idx1,idx2,idx3,idx4,k1,k2,kmin,kmax,idxrow,
     $     min,max
      double precision tmpsum1,tmpsum2,tmpval1,tmpval2,tmpval3,
     $     tmpvec3(MAXQBA+MAXPX),tmpvec4(MAXQBA),tmpvec5(MAXQBA),
     $     tmpvec6(MAXQ),tmpvec7(MAXPXTQ),tmpvec8(MAXQ),
     $     tmpmat1(MAXQBA,MAXQBA),tmpmat2(MAXQBA,MAXQBA),
     $     tmpmat3(MAXQ,MAXPXTQ),tmpmat4(MAXQBA,MAXPXTQ),
     $     tmpmat5(MAXQ,MAXQ),tmpmat6(MAXQBA,MAXPXTQ),dexp
      
	  ldif = ldvec(38)
	  ldcp = ldvec(39)
c      call intpr('ldif',-1,ldif,1)
c      call intpr('ldcp',-1,ldcp,1)	  
c      return
	  
      pxtq = px*q
      nfixpar = pxtq + pv + q
      qb = 2*q
      qbu = qb*(qb+1)/2
      qeu = q*(q+1)/2
      ntotpar = ne1 + nfixpar + qbu + qeu
      
      qba = qb**2
      qea = q**2
      nisq = ni**2
      qbcb = qb**3
      qecb = q**3
      qe4p = q**4
      
      tidp1 = tid + 1
      
      ldtm2 = MAXQBA
      ldtm3 = MAXQ
      ldtm4 = MAXQBA
      ldtm5 = MAXQ
c      call intpr('pxtq',-1,pxtq,1)
c      call intpr('nfixpar',-1,nfixpar,1)	  
c      call intpr('qb',-1,qb,1)
c      call intpr('qbu',-1,qbu,1)
c      call intpr('qeu',-1,qeu,1)
c      call intpr('ntotpar',-1,ntotpar,1)	  
c      call intpr('qba',-1,qba,1)
c      call intpr('qea',-1,qea,1)
c      call intpr('nisq',-1,nisq,1)
c      call intpr('qbcb',-1,qbcb,1)
c      call intpr('qecb',-1,qecb,1)
c      call intpr('qe4p',-1,qe4p,1)
c      call intpr('tidp1',-1,tidp1,1)
c      return
      
      
c     generate the constant vectors tmpvec1 and tmpvec2      
      if (tid.gt.0) then
        do 30 k=1,tid
          tmpsum1 = vibeta
          do 20 l=1,q
            tmpsum2=0.0d0
            do 10 s=1,px
              tmpsum2=tmpsum2+xit(k,s)*alphamat(s,l)
 10         continue
            xitalp(k,l)=tmpsum2
            tmpsum1 = tmpsum1 + tmpsum2*gamma(l)
 20       continue
          tmpvec1(k)=dexp(tmpsum1)
          tmpvec2(k)=h0t(k)*tmpvec1(k)
 30     continue
      end if
c      call dblepr('tmpvec1',-1,tmpvec1,23)
c      call dblepr('tmpvec2',-1,tmpvec2,23)
c      call dblepr('xit(,1)',-1,xit(1,1),24)
c      call dblepr('xit(,3)',-1,xit(1,3),24)
c      call dblepr('zit(,1)',-1,zit(1,1),24)
c      call dblepr('zit(,2)',-1,zit(1,2),24)
c      call dblepr('xitalp',-1,xitalp(1,1),24)
c      call dblepr('xitalp',-1,xitalp(1,2),24)
c      return
      
      do 35 l=1,q
        tmpsum2=0.0d0
        do 32 s=1,px
          tmpsum2=tmpsum2+xit(tidp1,s)*alphamat(s,l)
 32     continue
        xitalp(tidp1,l)=tmpsum2
 35   continue
c      call dblepr('xitalp',-1,xitalp(1,1),24)
c      call dblepr('xitalp',-1,xitalp(1,2),24)
c      return
      
      
      call veczero(ntotpar,retGrad)
      call matzero(ldif,ntotpar,ntotpar,retInfo)
      call matzero(ldcp,ntotpar,ntotpar,retCpd)
      
c     
c     First derivative
c
      
c     tid independent sums
      
c     palpha      
      do 50 j=1,ni
        do 50 k=1,q
          tmpsum1 = 0.0d0
          do 40 s=1,q
            tmpsum1=tmpsum1+Er(s,j)*sigeinv(s,k)
 40       continue
          
          do 50 l=1,px
            idx1 = ne1 + (k-1)*px + l
            retGrad(idx1)=retGrad(idx1) + tmpsum1*xi(j,l)
 50   continue
      
      
c     psigb      
      do 60 k=1,qb
        do 60 l=1,qb
          idx1 = (k-1)*qb + l
          tmpvec3(idx1)=vecwb(idx1)
          do 60 s=1,qb
            do 60 t=1,qb
              tmpvec3(idx1)=tmpvec3(idx1) -
     $              Ebb(t,s)*krwb2((s-1)*qb+t,idx1)
 60   continue              
      
      offset1 = ne1 + nfixpar
      do 80 k=1,qbu
        tmpsum1 = 0.0d0
        do 70 l=1,qba
          tmpsum1 = tmpsum1 + tmpvec3(l)*Dqb(l,k)
 70     continue
        idx1 = offset1 + k
        retGrad(idx1)=retGrad(idx1)-tmpsum1/2.0d0
 80   continue
      
      
c     psige      
      do 100 k=1,q
        do 100 l=1,k
          idx1 = (k-1)*q + l
          tmpvec4(idx1)=0.0d0
          do 90 j=1,ni
            tmpvec4(idx1)=tmpvec4(idx1)+Er2(idx1,j)
 90       continue
          if (l.ne.k) tmpvec4((l-1)*q+k)=tmpvec4(idx1)
 100  continue
c      call dblepr('Er2(,+)',-1,tmpvec4,4)
c      call dblepr('vecwe',-1,vecwe,4)
c      return
      
      do 110 k=1,q
        do 110 l=1,q
          idx1 = (k-1)*q + l
          tmpvec3(idx1)=ni*vecwe(idx1)
          do 110 s=1,q
            do 110 t=1,q
              idx2 = (s-1)*q + t
              tmpvec3(idx1)=tmpvec3(idx1)-tmpvec4(idx2)*krwe2(idx2,idx1)
 110  continue              
      
      offset1 = ne1 + nfixpar + qbu
      do 130 k=1,qeu
        tmpsum1 = 0.0d0
        do 120 l=1,qea
          tmpsum1 = tmpsum1 + tmpvec3(l)*Dqe(l,k)
 120    continue
        idx1 = offset1 + k
        retGrad(idx1)=retGrad(idx1)-tmpsum1/2.0d0
 130  continue
c      return
	  
      
c     tid dependent sums
      
      call veczero(px,tmpvec3(maxqba+1))
      call veczero(q,tmpvec5)
      call veczero(q,tmpvec6)
      
      if (tid.gt.0) then
        
        tmpsum1 = 0.0d0
        do 150 k=1,tid
c         ph0t		
          retGrad(k)=retGrad(k)-tmpvec1(k)*Eexpz(k)

          tmpval1 = tmpvec2(k)*Eexpz(k)
          do 140 l=1,px
            idx1 = maxqba + l
            tmpvec3(idx1)=tmpvec3(idx1)+tmpval1*xit(k,l)
 140      continue
          
          tmpsum1 = tmpsum1 + tmpval1
          
          do 150 l=1,q
            tmpvec5(l)=tmpvec5(l)+tmpval1*xitalp(k,l)
            
            do 150 s=1,2
              tmpvec6(l)=tmpvec6(l)+
     $              tmpvec2(k)*zit(k,s)*Eexpzb((l-1)*2+s,k)
 150    continue
        
        offset1 = ne1 + pxtq + pv
        do 170 k=1,q
c          tmpsum2 = 0.0d0
           
c         palpha		  
          do 160 l=1,px		  
            idx1 = ne1 + (k-1)*px + l
            idx2 = maxqba + l
            retGrad(idx1)=retGrad(idx1)-gamma(k)*tmpvec3(idx2)
c            tmpsum2 = tmpsum2 + tmpvec3(idx2)*alphamat(l,k)
 160      continue
          
c         pgamma
          idx2 = offset1 + k
          retGrad(idx2)=retGrad(idx2) - tmpvec5(k) - tmpvec6(k)
c          tmpvec6(k) = tmpsum2
 170    continue
        
c       pbeta 
        offset1 = ne1 + pxtq
        do 180 k=1,pv
          idx1 = offset1 + k
          retGrad(idx1)=retGrad(idx1) - tmpsum1*vi(k)
 180    continue
        tmpsum2 = tmpsum1
        
        
        if (deltai.eq.1) then
c         ph0t		
          retGrad(tid) = retGrad(tid) + 1.0d0/h0t(tid)
          
          offset1 = ne1 + pxtq + pv
          do 210 k=1,q
             
c           palpha			
            do 190 l=1,px
              idx1 = ne1 + (k-1)*px + l
              retGrad(idx1)=retGrad(idx1) + gamma(k)*xit(tidp1,l)
c              tmpsum1 = tmpsum1 + xit(tidp1,l)*alphamat(l,k)
 190        continue
            
c           pgamma 
            tmpsum1 = xitalp(tidp1,k)
            do 200 l=1,2
              tmpsum1 = tmpsum1 + zit(tidp1,l)*Eb((k-1)*2+l)
 200        continue
            
            idx2 = offset1 + k
            retGrad(idx2) = retGrad(idx2) + tmpsum1
 210      continue
          
c         pbeta 
          offset1 = ne1 + pxtq
          do 220 k=1,pv
            idx1 = offset1 + k
            retGrad(idx1) = retGrad(idx1) + vi(k)
 220      continue          
        end if
        
      end if
c      call dblepr('retGrad',-1,retGrad,ntotpar)
c      return
      
	  
c      
c     Second derivative (information matrix)
c
      
c     tid independent sums
      
c     palptalp      
      do 230 k=1,q
        do 230 l=1,px
          idx1 = ne1 + (k-1)*px + l
          do 230 s=1,q
            do 230 t=1,px
              idx2 = ne1 + (s-1)*px + t
              do 230 j=1,ni
                retInfo(idx2,idx1)=retInfo(idx2,idx1)+
     $                xi(j,t)*sigeinv(s,k)*xi(j,l)
 230  continue
      
      
c     psigbsigb      
      do 250 k=1,qb
        do 250 l=1,k
          idx1 = (k-1)*qb + l          
          tmpvec3(idx1)=0.0d0
          do 240 s=1,qb
            do 240 t=1,qb
              tmpvec3(idx1)=tmpvec3(idx1)+
     $              sigbinv(l,t)*Ebb(t,s)*sigbinv(s,k)
 240      continue
          if (l.ne.k) tmpvec3((l-1)*qb+k)=tmpvec3(idx1)
 250  continue
      
      do 260 k=1,qba
        idx1 = (k-1)/qb + 1
        idx2 = k - (idx1-1)*qb
        do 260 l=1,k
          idx3 = (l-1)/qb + 1
          idx4 = l - (idx3-1)*qb
          tmpmat1(l,k)=tmpvec3((idx1-1)*qb+idx3)*sigbinv(idx4,idx2)+
     $         sigbinv(idx3,idx1)*tmpvec3((idx2-1)*qb+idx4)
        if (l.ne.k) tmpmat1(k,l)=tmpmat1(l,k)
 260  continue
      
      do 270 k=1,qba
        do 270 l=1,qba
          tmpmat2(l,k)=0.0d0
          do 270 s=1,qba
            tmpmat2(l,k)=tmpmat2(l,k)+Tqbqb(l,s)*krwb2(s,k)-
     $            tmpmat1(l,s)*Tqbqb(s,k)
 270  continue
      
      offset1 = ne1 + nfixpar
      do 290 k=1,qbu
        idx1 = offset1 + k
        do 290 l=1,k
          idx2 = offset1 + l
          tmpsum1 = 0.0d0
          do 280 s=1,qba
            do 280 t=1,qba
              tmpsum1=tmpsum1+Dqb(t,l)*tmpmat2(t,s)*Dqb(s,k)
 280      continue
          retInfo(idx2,idx1)=retInfo(idx2,idx1)-tmpsum1/2.0d0
 290  continue
      
      
c     psigesige      
      do 310 k=1,q
        do 310 l=1,k
          idx1 = (k-1)*q + l          
          tmpvec3(idx1)=0.0d0
          do 300 s=1,q
            do 300 t=1,q
              tmpvec3(idx1)=tmpvec3(idx1)+
     $              sigeinv(l,t)*tmpvec4((s-1)*q+t)*sigeinv(s,k)
 300      continue
          if (l.ne.k) tmpvec3((l-1)*q+k)=tmpvec3(idx1)
 310  continue
      
      do 320 k=1,qea
        idx1 = (k-1)/q + 1
        idx2 = k - (idx1-1)*q
        do 320 l=1,k
          idx3 = (l-1)/q + 1
          idx4 = l - (idx3-1)*q
          tmpmat1(l,k)=tmpvec3((idx1-1)*q+idx3)*sigeinv(idx4,idx2)+
     $         sigeinv(idx3,idx1)*tmpvec3((idx2-1)*q+idx4)
        if (l.ne.k) tmpmat1(k,l)=tmpmat1(l,k)
 320  continue
      
      do 330 k=1,qea
        do 330 l=1,qea
          tmpmat2(l,k)=0.0d0
          do 330 s=1,qea
            tmpmat2(l,k)=tmpmat2(l,k)+ni*Tqeqe(l,s)*krwe2(s,k)-
     $            tmpmat1(l,s)*Tqeqe(s,k)
 330  continue
      
      offset1 = ne1 + nfixpar + qbu
      do 350 k=1,qeu
        idx1 = offset1 + k
        do 350 l=1,k
          idx2 = offset1 + l
          tmpsum1 = 0.0d0
          do 340 s=1,qea
            do 340 t=1,qea
              tmpsum1=tmpsum1+Dqe(t,l)*tmpmat2(t,s)*Dqe(s,k)
 340      continue
          retInfo(idx2,idx1)=retInfo(idx2,idx1)-tmpsum1/2.0d0
 350  continue
      
      
c     palpsige      
      call matzero(ldtm4,qea,pxtq,tmpmat4)
      do 380 j=1,ni
        do 370 k=1,q
          tmpsum1 = 0.0d0
          do 360 l=1,q
            tmpsum1=tmpsum1+sigeinv(k,l)*Er(l,j)
            do 360 s=1,px
              tmpmat3(l,(k-1)*px+s)=sigeinv(l,k)*xi(j,s)
 360      continue
          tmpvec5(k) = tmpsum1          
 370    continue
c        call dblepr('tmpvec5',-1,tmpvec5,q)
c        call dblepr('tmpmat3(,1)',-1,tmpmat3(1,1),q)
c        call dblepr('tmpmat3(,3)',-1,tmpmat3(1,3),q)
c        call dblepr('tmpmat3(,4)',-1,tmpmat3(1,4),q)
c        call dblepr('tmpmat3(,6)',-1,tmpmat3(1,6),q)
c        return
        		
        do 380 k=1,q
          do 380 l=1,q
            idx1 = (k-1)*q + l
            do 380 s=1,pxtq
              tmpmat4(idx1,s)=tmpmat4(idx1,s)+tmpvec5(k)*tmpmat3(l,s)+
     $              tmpmat3(k,s)*tmpvec5(l)
 380  continue
c 385  continue
c        call dblepr('tmpmat4(,1)',-1,tmpmat4(1,1),4)
c        call dblepr('tmpmat4(,2)',-1,tmpmat4(1,2),4)
c        call dblepr('tmpmat4(,3)',-1,tmpmat4(1,3),4)
c        call dblepr('tmpmat4(,4)',-1,tmpmat4(1,4),4)
c        call dblepr('tmpmat4(,5)',-1,tmpmat4(1,5),4)
c        call dblepr('tmpmat4(,6)',-1,tmpmat4(1,6),4)
c        return 
      
	  
      offset1 = ne1 + nfixpar + qbu
      do 400 k=1,qeu
        idx1 = offset1 + k
        do 400 l=1,pxtq
          idx2 = ne1 + l
          tmpsum1 = 0.0d0
          do 390 s=1,qea
c            do 390 t=1,qea
              tmpsum1=tmpsum1+tmpmat4(s,l)*Dqe(s,k)
 390      continue
          retInfo(idx2,idx1)=retInfo(idx2,idx1)+tmpsum1/2.0d0
 400  continue
c      return
	  
      
c     tid dependent sums
      
      if (tid.gt.0) then
        
c       summation over tid for some hessian parts
         
        tmpsum1 = 0.0d0
        call veczero(px,tmpvec3)
        call veczero(q,tmpvec6)
        
        offset1 = ne1 + pxtq
        offset2 = ne1 + pxtq + pv
        do 470 k=1,tid
           
          tmpval1 = tmpvec1(k)*Eexpz(k)
          tmpval2 = tmpvec2(k)*Eexpz(k)
		  
          do 410 l=1,px
            tmpvec3(l)=tmpvec3(l)+tmpval2*xit(k,l)		  
 410      continue      
		  
          tmpsum1=tmpsum1 + tmpval2
          
          do 430 l=1,q
            idx1 = offset2 + l
            
            tmpsum2 = 0.0d0
            do 420 s=1,2
              tmpsum2=tmpsum2+zit(k,s)*Eexpzb((l-1)*2+s,k)
 420        continue
            tmpvec8(l)=tmpsum2*tmpvec2(k)
			
            tmpvec6(l)=tmpvec6(l) + tmpval2*xitalp(k,l) + tmpvec8(l)
			
			
c           ph0tgam			
            retInfo(k,idx1)=retInfo(k,idx1)+tmpval1*xitalp(k,l)+
     $           tmpsum2*tmpvec1(k)
            
            
            do 430 s=1,px
              idx2 = (l-1)*px + s
              tmpvec7(idx2)=gamma(l)*xit(k,s)
              
              idx3 = ne1 + idx2
c             ph0talp			  
              retInfo(k,idx3)=retInfo(k,idx3)+tmpval1*tmpvec7(idx2)
              
 430      continue
          
          
          do 440 l=1,pv
            idx1 = offset1 + l
c           ph0tbet
            retInfo(k,idx1)=retInfo(k,idx1)+tmpval1*vi(l)
            
 440      continue      
          
          
c         palpalp			
          do 450 l=1,pxtq
            idx1 = ne1 + l
            do 450 s=1,l
              idx2 = ne1 + s
              retInfo(idx2,idx1)=retInfo(idx2,idx1)+
     $             tmpval2*tmpvec7(s)*tmpvec7(l)
 450      continue      
          
          
          do 470 l=1,q
            idx1 = offset2 + l
			
c           palpgam (more summation on palpgam at the end (outside the tid loop))
            do 460 s=1,pxtq
              idx2=ne1+s
              retInfo(idx2,idx1)=retInfo(idx2,idx1)+
     $             tmpvec7(s)*(tmpval2*xitalp(k,l)+tmpvec8(l))
 460        continue      
            
            
c           pgamgam			
            do 470 s=1,l
              idx2 = offset2 + s
              retInfo(idx2,idx1)=retInfo(idx2,idx1)+
     $             tmpval2*xitalp(k,s)*xitalp(k,l)+
     $             tmpvec2(k)*Eexz3btb((l-1)*q+s,k)+
     $             xitalp(k,s)*tmpvec8(l)+
     $             tmpvec8(s)*xitalp(k,l)
 
 470    continue         
        
        
c       use the summation from the previous part to calculate the rest of hessian parts 

        do 490 k=1,pv
          idx1 = offset1 + k
		  
          do 480 l=1,q
c           palpbet
            do 480 s=1,px
              idx2=ne1+(l-1)*px + s
              retInfo(idx2,idx1)=retInfo(idx2,idx1)+
     $             gamma(l)*tmpvec3(s)*vi(k)
 480      continue
          
		  
c         pbetbet
          do 490 l=1,k
            idx2 = offset1 + l
            retInfo(idx2,idx1)=retInfo(idx2,idx1)+tmpsum1*vi(l)*vi(k)
 490    continue
        
        
        do 510 k=1,q
          idx1 = offset2 + k
		  
c         palpgam (second part of summation for palpgam)
          do 500 s=1,px
            idx2=ne1+(k-1)*px + s
            retInfo(idx2,idx1)=retInfo(idx2,idx1)+tmpvec3(s)
 500      continue
          
          
c         pbetgam 		  
          do 510 l=1,pv
            idx2 = offset1 + l
            retInfo(idx2,idx1)=retInfo(idx2,idx1)+vi(l)*tmpvec6(k)
 510    continue
        
        
c       hessian parts for deltai==1
        
        if (deltai.eq.1) then
c         ph0th0t           
          retInfo(tid,tid)=retInfo(tid,tid)+1.0d0/h0t(tid)**2
          
          
c         palpgam
          offset1 = ne1 + pxtq + pv
          do 520 k=1,q
            idx1 = offset1 + k
            do 520 l=1,px
              idx2 = ne1 + (k-1)*px + l
              retInfo(idx2,idx1)=retInfo(idx2,idx1)-xit(tidp1,l)
 520      continue
        end if
        
      end if
c      call dblepr('retInfo(,1)',-1,retInfo(1,1),ntotpar)
c      call dblepr('retInfo(,tid)',-1,retInfo(1,tid),ntotpar)
c      call dblepr('retInfo(,57)',-1,retInfo(1,57),ntotpar)
c      call dblepr('retInfo(,62)',-1,retInfo(1,62),ntotpar)
c      call dblepr('retInfo(,63)',-1,retInfo(1,63),ntotpar)
c      call dblepr('retInfo(,64)',-1,retInfo(1,64),ntotpar)
c      call dblepr('retInfo(,65)',-1,retInfo(1,65),ntotpar)
c      call dblepr('retInfo(,66)',-1,retInfo(1,66),ntotpar)
c      call dblepr('retInfo(,75)',-1,retInfo(1,75),ntotpar)
c      call dblepr('retInfo(,76)',-1,retInfo(1,76),ntotpar)
c      call dblepr('retInfo(,77)',-1,retInfo(1,77),ntotpar)
c      call dblepr('retInfo(,78)',-1,retInfo(1,78),ntotpar)
c      return
      
      
	  
c      
c     Cross product of first derivative
c
      
c     tid independent sums

c     palptalp
      do 550 j=1,ni
        do 530 k=1,q
          do 530 l=1,px
            do 530 s=1,q
              tmpmat3(s,(k-1)*px+l)=sigeinv(s,k)*xi(j,l)
 530    continue
      
        do 550 i=1,ni
          do 540 k=1,q
            do 540 l=1,px
              do 540 s=1,q
                tmpmat4(s,(k-1)*px+l)=sigeinv(s,k)*xi(i,l)
 540      continue
          
          idx1 = (j-1)*ni + i
          
          do 550 k=1,pxtq
            idx2 = ne1 + k
            do 550 l=1,k
              idx3 = ne1 + l
              do 550 s=1,q
                do 550 t=1,q
                  retCpd(idx3,idx2)=retCpd(idx3,idx2)+
     $                  tmpmat4(t,l)*Er2c((s-1)*q+t,idx1)*tmpmat3(s,k)
 550  continue
      
      
c     palptsigb
      call veczero(pxtq,tmpvec7)
      call matzero(ldtm4,qba,pxtq,tmpmat4)
      do 580 j=1,ni
        do 560 k=1,q
          do 560 l=1,px
            idx1 = (k-1)*px + l
            do 560 s=1,q
              tmpmat3(s,idx1)=sigeinv(s,k)*xi(j,l)
 560    continue
        
        do 570 k=1,pxtq
          do 570 l=1,q
            tmpvec7(k)=tmpvec7(k)+Er(l,j)*tmpmat3(l,k)
 570    continue
        
        do 580 k=1,qba
          do 580 l=1,pxtq
            do 580 s=1,q
              tmpmat4(k,l)=tmpmat4(k,l)+tmpmat3(s,l)*Erbb((k-1)*q+s,j)
 580  continue
c      call dblepr('tmpvec7',-1,tmpvec7,pxtq)
c      call dblepr('tmpmat4(,1)',-1,tmpmat4(1,1),qba)
c      call dblepr('tmpmat4(,6)',-1,tmpmat4(1,6),qba)
c	  return      
	  
      do 590 k=1,qba
        do 590 l=1,pxtq
          tmpmat6(k,l) = tmpvec7(l)*vecwb(k)
          do 590 s=1,qba
            tmpmat6(k,l)=tmpmat6(k,l)-tmpmat4(s,l)*krwb2(s,k)
 590  continue
      
      offset1 = ne1 + nfixpar
      do 610 k=1,qbu
        do 610 l=1,pxtq
          tmpsum1 = 0.0d0
          do 600 s=1,qba
            tmpsum1=tmpsum1+tmpmat6(s,l)*Dqb(s,k)
 600      continue
          idx1 = offset1 + k
          idx2 = ne1 + l
          retCpd(idx2,idx1)=retCpd(idx2,idx1)-tmpsum1/2.0d0
 610  continue
c      call dblepr('retCpd(,66)',-1,retCpd(1,66),ntotpar)
c      call dblepr('retCpd(,75)',-1,retCpd(1,75),ntotpar)
c	  return
      
      
c     psigbtsigb
      do 630 k=1,qba
        do 630 l=1,k          
          tmpmat1(l,k)=vecwb(l)*vecwb(k)
          do 620 s=1,qba
             idx1 = (s-1)/qb + 1
             idx2 = s - (idx1-1)*qb
             tmpmat1(l,k)=tmpmat1(l,k) -
     $            vecwb(l)*Ebb(idx2,idx1)*krwb2(s,k) -
     $            krwb2(l,s)*Ebb(idx2,idx1)*vecwb(k)
            
            do 620 t=1,qba
              tmpmat1(l,k)=tmpmat1(l,k)+krwb2(l,t)*Eb4(t,s)*krwb2(s,k)
 620      continue
      
          if (l.ne.k) tmpmat1(k,l)=tmpmat1(l,k)
 630  continue
      
      do 650 k=1,qbu
        idx1 = offset1 + k
        do 650 l=1,k
          idx2 = offset1 + l
          
          tmpsum1 = 0.0d0
          do 640 s=1,qba
            do 640 t=1,qba
              tmpsum1=tmpsum1+Dqb(t,l)*tmpmat1(t,s)*Dqb(s,k)
 640      continue
      
          retCpd(idx2,idx1)=retCpd(idx2,idx1)+tmpsum1/4.0d0
 650  continue
      
      
c     palptsige
      call matzero(ldtm4,qea,pxtq,tmpmat4)
      do 690 j=1,ni
        do 680 k=1,q
          do 660 l=1,px
            idx1 = (k-1)*px + l
            do 660 s=1,q
              tmpmat3(s,idx1)=sigeinv(s,k)*xi(j,l)
 660      continue
      
          do 680 l=1,qea
            tmpsum1 = 0.0d0
            idx1 = (l-1)*q + k
            do 670 i=1,ni
              tmpsum1=tmpsum1+Er3(idx1,(i-1)*ni+j)
 670        continue
            tmpmat1(k,l)=tmpsum1
 680    continue

        do 690 k=1,pxtq
          do 690 l=1,qea
            do 690 s=1,q
              tmpmat4(l,k)=tmpmat4(l,k)+tmpmat3(s,k)*tmpmat1(s,l)
 690  continue
      
      do 700 k=1,qea
        do 700 l=1,pxtq
          tmpmat6(l,k)=ni*tmpvec7(l)*vecwe(k)
          do 700 s=1,qea
            tmpmat6(l,k)=tmpmat6(l,k)-tmpmat4(s,l)*krwe2(s,k)
 700  continue
      
      offset1 = ne1 + nfixpar + qbu
      do 720 k=1,qeu
        idx1 = offset1 + k
        do 720 l=1,pxtq
          idx2 = ne1 + l
          tmpsum1=0.0d0
          do 710 s=1,qea
            tmpsum1=tmpsum1+tmpmat6(l,s)*Dqe(s,k)
 710      continue
          retCpd(idx2,idx1)=retCpd(idx2,idx1)-tmpsum1/2.0d0
 720  continue
      
      
c     psigbtsige      
      do 730 k=1,qea
        tmpvec4(k)=0.0d0
        do 730 j=1,ni
          tmpvec4(k)=tmpvec4(k)+Er2(k,j)
 730  continue
      
      do 750 k=1,qea
        do 750 l=1,qba
          idx1 = (k-1)*qba + l
          tmpsum1 = 0.0d0
          do 740 j=1,ni
            tmpsum1=tmpsum1+Eb2r2(idx1,j)
 740      continue
          tmpmat1(l,k)=tmpsum1
 750  continue
      
      do 780 k=1,qea
        do 780 l=1,qba
          tmpmat2(l,k)=ni*vecwb(l)*vecwe(k)
          do 760 s=1,qea
            tmpmat2(l,k)=tmpmat2(l,k)-vecwb(l)*tmpvec4(s)*krwe2(s,k)
 760      continue
          
          do 770 s=1,qba
            idx1 = (s-1)/qb + 1
            idx2 = s - (idx1-1)*qb
            tmpmat2(l,k)=tmpmat2(l,k)-
     $            ni*krwb2(l,s)*Ebb(idx2,idx1)*vecwe(k)
 770      continue
          
          do 780 s=1,qea
            do 780 t=1,qba
              tmpmat2(l,k)=tmpmat2(l,k)+
     $              krwb2(l,t)*tmpmat1(t,s)*krwe2(s,k)
 780  continue          
      
      offset1 = ne1 + nfixpar + qbu
      offset2 = ne1 + nfixpar
      do 800 k=1,qeu
        idx1 = offset1 + k
        do 800 l=1,qbu
          idx2 = offset2 + l
          tmpsum1=0.0d0
          do 790 s=1,qea
            do 790 t=1,qba
              tmpsum1=tmpsum1+Dqb(t,l)*tmpmat2(t,s)*Dqe(s,k)
 790      continue
          retCpd(idx2,idx1)=retCpd(idx2,idx1) + tmpsum1/4.0d0
 800  continue
c      call dblepr('retCpd(,76)',-1,retCpd(1,76),ntotpar)
c      call dblepr('retCpd(,77)',-1,retCpd(1,77),ntotpar)
c      call dblepr('retCpd(,78)',-1,retCpd(1,78),ntotpar)
c      return      
      
      
c     psigetsige
      do 810 k=1,qea
        do 810 l=1,qea
          idx1 = (k-1)*qea + l
          tmpmat1(l,k) = 0.0d0
          do 810 s=1,nisq
            tmpmat1(l,k)=tmpmat1(l,k)+Er4(idx1,s)
 810  continue
      
      do 820 k=1,qea
        do 820 l=1,qea
          tmpmat2(l,k)=ni**2*vecwe(l)*vecwe(k)

          do 820 s=1,qea
            tmpmat2(l,k)=tmpmat2(l,k)-ni*vecwe(l)*tmpvec4(s)*krwe2(s,k)-
     $            ni*krwe2(l,s)*tmpvec4(s)*vecwe(k)

            do 820 t=1,qea
              tmpmat2(l,k)=tmpmat2(l,k)+
     $              krwe2(l,t)*tmpmat1(t,s)*krwe2(s,k)
 820  continue          
      
      do 840 k=1,qeu
        idx1 = offset1 + k
        do 840 l=1,qeu
          idx2 = offset1 + l
          tmpsum1=0.0d0
          do 830 s=1,qea
            do 830 t=1,qea
              tmpsum1=tmpsum1+Dqe(t,l)*tmpmat2(t,s)*Dqe(s,k)
 830      continue
          retCpd(idx2,idx1)=retCpd(idx2,idx1) + tmpsum1/4.0d0
 840  continue
c      return
	  
      
c     tid dependent sums

c      call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c      return
      
      if (tid.gt.0) then         
        
c
c       deltai dependent quantities
c
        if (deltai.eq.1) then
          
          offset1 = ne1 + pxtq
          offset2 = ne1 + pxtq + pv
          do 900 k=1,tid
            tmpval1 = tmpvec1(k)*Eexpz(k)
            
c           ph0tth0t
            if (k.lt.tid) then
              retCpd(k,tid)=retCpd(k,tid)-tmpval1/h0t(tid)
            else
              retCpd(k,tid)=retCpd(k,tid)-2.0d0*tmpval1/h0t(tid)
			end if
			
c           ph0ttalp
            do 850 l=1,q
              do 850 s=1,px
                idx1 = ne1 + (l-1)*px + s             
                retCpd(k,idx1)=retCpd(k,idx1)-
     $               tmpval1*gamma(l)*xit(tidp1,s)
 850        continue   

c           ph0ttbet
            do 860 l=1,pv
              idx1 = offset1 + l
              retCpd(k,idx1)=retCpd(k,idx1)-tmpval1*vi(l)
 860        continue
      
c           ph0ttgam
            do 880 l=1,q
              idx1 = offset2 + l
              
              tmpsum1=0.0d0
              do 870 s=1,2
                tmpsum1=tmpsum1+zit(tidp1,s)*Eexpzb((l-1)*2+s,k)
 870          continue               
      
              retCpd(k,idx1)=retCpd(k,idx1)-tmpval1*xitalp(tidp1,l)-
     $                  tmpsum1*tmpvec1(k)
 880        continue
            
            
c     
c           deltai==1 and (k equal to tid)
c     
            
            if (k.eq.tid) then
c             ph0tth0t
              retCpd(tid,tid)=retCpd(tid,tid)+1/h0t(tid)**2
              
c             ph0ttalp, ph0tbet, ph0ttgam, ph0ttsigb, ph0ttsige
              do 890 l=ne1+1,ntotpar
                retCpd(tid,l)=retCpd(tid,l)+retGrad(l)/h0t(tid)
 890          continue
            end if
            
 900      continue           
c         end of the loop over k=1,tid
c          call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c          return
		  
          
ccccccccccccccccccccccccccccccccccccc
      
c         components for calculating cross-products for tid > 0 and deltai=1
          
          call veczero(q,tmpvec4)
          call veczero(q,tmpvec6)
          call matzero(ldtm3,q,pxtq,tmpmat3)
          call matzero(ldtm5,q,q,tmpmat5)
          tmpsum2=0.0d0
          do 950 k=1,tid
            
            tmpval1 = tmpvec2(k)*Eexpz(k)
            tmpsum2 = tmpsum2 + tmpval1
            
            do 940 l=1,q              
              tmpsum1=0.0d0
              tmpval2=0.0d0
              do 910 s=1,2
                idx1=(l-1)*2+s
                tmpsum1=tmpsum1+zit(tidp1,s)*Eexpzb(idx1,k)
                tmpval2=tmpval2+zit(k,s)*Eexpzb(idx1,k)
 910          continue
              tmpsum1=tmpsum1*tmpvec2(k)
              tmpvec6(l)=tmpvec6(l)+tmpsum1
              
              do 930 s=1,q
                do 920 t=1,px
                  idx1 = (s-1)*px + t
                  tmpmat3(l,idx1)=tmpmat3(l,idx1)+
     $                 tmpsum1*gamma(s)*xit(k,t)
 920            continue
      
                tmpmat5(l,s)=tmpmat5(l,s)-tmpsum1*xitalp(k,s)-
     $             tmpvec2(k)*Eexpz3b2((s-1)*q+l,k)
 930          continue
      
              tmpvec4(l)=tmpvec4(l)+tmpval1*xitalp(k,l)+
     $             tmpval2*tmpvec2(k)
 940        continue
            
 950      continue
c          call dblepr('tmpvec6',-1,tmpvec6,q)
c          call dblepr('tmpvec4',-1,tmpvec4,q)
c          call dblepr('zit(tidp1,1)',-1,zit(tidp1,1),1)
c          call dblepr('zit(tidp1,2)',-1,zit(tidp1,2),1)
c          call dblepr('Eb(,1)',-1,Eb(1),2)
c          call dblepr('Eb(,2)',-1,Eb(3),2)
c          return
          
          
          do 970 k=1,q
            tmpsum1=0.0d0
            do 960 l=1,2
              tmpsum1=tmpsum1+zit(tidp1,l)*Eb((k-1)*2+l)
 960        continue
c            call dblepr('tmpsum1',-1,tmpsum1,1)
            tmpsum1=tmpsum1-tmpvec4(k)
c            call dblepr('tmpsum1',-1,tmpsum1,1)
            
            do 970 l=1,q
              tmpmat5(l,k)=tmpmat5(l,k)+tmpsum1*xitalp(tidp1,l)
 970      continue
c          call dblepr('tmpsum1',-1,tmpsum1,1)
c          call dblepr('tmpmat5(,1)',-1,tmpmat5(1,1),q)
c          call dblepr('tmpmat5(,2)',-1,tmpmat5(1,2),q) 
c          return
          
          
c         calculate the cross-products for tid > 0 and deltai=1
          
          offset1 = ne1 + pxtq
          offset2 = ne1 + pxtq + pv
          do 1020 k=1,q
            do 1010 l=1,px
              idx1 = ne1 + (k-1)*px + l
              tmpval1 = gamma(k)*xit(tidp1,l)
              
c             palptalp
              do 980 s=1,q
                do 980 t=1,px
                  idx2 = ne1 + (s-1)*px + t
                  retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $                 gamma(s)*xit(tidp1,t)*retGrad(idx1)+
     $                 tmpval1*(retGrad(idx2)-gamma(s)*xit(tidp1,t))
 980          continue
      
c             palptbet      
              do 990 s=1,pv
                idx2 = offset1 + s
                retCpd(idx1,idx2)=retCpd(idx1,idx2)+
     $               tmpval1*(retGrad(idx2)-vi(s))+
     $               retGrad(idx1)*vi(s)
 990          continue
              
c             palptgam
              idx3 = (k-1)*px + l
              do 1010 s=1,q
                idx2 = offset2 + s
                
                tmpsum1=0.0d0
                do 1000 j=1,ni
                  do 1000 t=1,q
                    tmpsum1=tmpsum1+
     $                    sigeinv(t,k)*xi(j,l)*Erzb((s-1)*q+t,j)
 1000           continue
      
                retCpd(idx1,idx2)=retCpd(idx1,idx2)+
     $             tmpval1*(retGrad(idx2)-xitalp(tidp1,s))+
     $             retGrad(idx1)*xitalp(tidp1,s)+tmpsum1-tmpmat3(s,idx3)
 1010       continue
            
c           pbettgam
            idx1 = offset2 + k
            do 1020 l=1,pv
              idx2 = offset1 + l
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $             vi(l)*(retGrad(idx1)-tmpvec6(k))-
     $             tmpsum2*vi(l)*xitalp(tidp1,k)
 1020     continue          
c          call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c          call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)
c		  return
          
          
c         pbettbet
          do 1030 k=1,pv
            idx1 = offset1 + k
            
            do 1030 l=1,pv
              idx2 = offset1 + l
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+vi(l)*retGrad(idx1)+
     $             retGrad(idx2)*vi(k)-vi(l)*vi(k)
 1030     continue
          
          
c         pgamtgam
          do 1040 k=1,q
            idx1 = offset2 + k
            do 1040 l=1,q
              idx2 = offset2 + l
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $             xitalp(tidp1,l)*xitalp(tidp1,k)+Ebzzb(l,k)+
     $             tmpmat5(l,k)+tmpmat5(k,l)
 1040     continue
c          call dblepr('xitalp(tidp1,1)',-1,xitalp(tidp1,1),1)
c          call dblepr('xitalp(tidp1,2)',-1,xitalp(tidp1,2),1)
c          call dblepr('Ebzzb(,1)',-1,Ebzzb(1,1),q)
c          call dblepr('Ebzzb(,2)',-1,Ebzzb(1,2),q)
c          call dblepr('tmpmat5(,1)',-1,tmpmat5(1,1),q)
c          call dblepr('tmpmat5(,2)',-1,tmpmat5(1,2),q) 
c          call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c          call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)
c		  return
          
          
c         more component quantities for (palp, pbet, pgam) x (psigb, psige)
      
          do 1070 k=1,q
            tmpsum1=0.0d0
            do 1050 l=1,2
              tmpsum1=tmpsum1+zit(tidp1,l)*Eb((k-1)*2+l)
 1050       continue
            
            do 1060 l=1,qba
              tmpmat1(k,l)=tmpsum1*vecwb(l)
              do 1060 s=1,qba
                tmpmat1(k,l)=tmpmat1(k,l)-Ebzbb(k,s)*krwb2(s,l)
 1060       continue
            
            do 1070 l=1,qea
              tmpmat2(k,l)=ni*tmpsum1*vecwe(l)
 1070     continue
          
          
          do 1090 k=1,q
            do 1090 l=1,qea
              tmpsum1=0.0d0
              do 1080 j=1,ni
                tmpsum1=tmpsum1+Ebzrr((l-1)*q+k,j)
 1080         continue
      
              do 1090 s=1,qea
                tmpmat2(k,s)=tmpmat2(k,s)-tmpsum1*krwe2(l,s)              
 1090     continue
          
          
c         calculate the cross prod (palp, pbet, pgam) x (psigb, psige)
          
          offset1 = ne1 + nfixpar
          offset2 = ne1 + pxtq          
          offset3 = ne1 + pxtq + pv
          
          do 1130 k=1,qbu
            idx1 = offset1 + k
            
            do 1120 l=1,q
c             palptsigb
              do 1100 s=1,px
                idx2 = ne1 + (l-1)*px + s
                retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $               gamma(l)*xit(tidp1,s)*retGrad(idx1)
 1100         continue
              
c             pgamtsigb
              tmpsum1=0.0d0
              do 1110 s=1,qba
                tmpsum1=tmpsum1+tmpmat1(l,s)*Dqb(s,k)
 1110         continue
              idx2 = offset3 + l            
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $             xitalp(tidp1,l)*retGrad(idx1)-tmpsum1/2.0d0
 1120       continue
            
            
c           pbettsigb
            do 1130 l=1,pv
              idx2 = offset2 + l
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+vi(l)*retGrad(idx1)
 1130     continue
c		  return
          
          
          offset1 = ne1 + nfixpar + qbu
          do 1170 k=1,qeu
            idx1 = offset1 + k
            
            do 1160 l=1,q
c             palptsige
              do 1140 s=1,px
                idx2 = ne1 + (l-1)*px + s
                retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $               gamma(l)*xit(tidp1,s)*retGrad(idx1)
 1140         continue
              
c             pgamtsige
              tmpsum1=0.0d0
              do 1150 s=1,qea
                tmpsum1=tmpsum1+tmpmat2(l,s)*Dqe(s,k)
 1150         continue
              idx2 = offset3 + l            
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $             xitalp(tidp1,l)*retGrad(idx1)-tmpsum1/2.0d0
 1160       continue
            
            
c           pbettsige
            do 1170 l=1,pv
              idx2 = offset2 + l
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+vi(l)*retGrad(idx1)
 1170     continue
          
        end if
c       end of if (deltai==1)
c        call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c        call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c        call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)
c        call dblepr('retCpd(,66)',-1,retCpd(1,66),ntotpar)
c        call dblepr('retCpd(,75)',-1,retCpd(1,75),ntotpar)
c        call dblepr('retCpd(,76)',-1,retCpd(1,76),ntotpar)
c        call dblepr('retCpd(,77)',-1,retCpd(1,77),ntotpar)
c        call dblepr('retCpd(,78)',-1,retCpd(1,78),ntotpar)
c        return
		
        
ccccccccccccccccccccccccccccccccccccccccccccc
        
c         
c       deltai independent quantities
c        
     
c       deltai indep otherwise tid > 0 - 1 dimensional

        do 1470 k=1,tid

c         ph0ttalp (and some component quantities)
          call veczero(pxtq,tmpvec7)
          call matzero(ldtm4,q,pxtq,tmpmat4)
          do 1210 j=1,ni
            do 1180 l=1,q
              do 1180 s=1,px
                do 1180 t=1,q
                  tmpmat3(t,(l-1)*px+s)=sigeinv(t,l)*xi(j,s)
 1180       continue
            
            do 1210 l=1,q
              do 1210 s=1,px
                idx1 = (l-1)*px + s
				idx2 = ne1 + idx1
                tmpsum1 = 0.0d0
                do 1200 t=1,q
                  tmpsum1=tmpsum1+Eexpzr((j-1)*q+t,k)*tmpmat3(t,idx1)
                  
                  tmpsum2=0.0d0
                  do 1190 i=1,q
                    tmpsum2=tmpsum2+
     $                   tmpmat3(i,idx1)*Eexpzrzb((j-1)*qea+(t-1)*q+i,k)
 1190             continue
                  tmpmat4(t,idx1)=tmpmat4(t,idx1)+tmpsum2*tmpvec2(k)
 1200           continue
                retCpd(k,idx2)=retCpd(k,idx2)-tmpsum1*tmpvec1(k)
                tmpvec7(idx1)=tmpvec7(idx1)+tmpsum1*tmpvec2(k)                     
 1210     continue
c          call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c          return
          
          
          offset1 = ne1 + pxtq
          offset2 = ne1 + pxtq + pv
          do 1240 l=1,q
            do 1240 s=1,px
              idx1 = (l-1)*px + s
              idx2 = ne1 + idx1
              
              do 1230 t=1,q
c               palptalp
                do 1220 i=1,px
                  idx3 = (t-1)*px + i
                  idx4 = ne1 + idx3
                  retCpd(idx4,idx2)=retCpd(idx4,idx2)-
     $                 tmpvec7(idx3)*gamma(l)*xit(k,s)-
     $                 tmpvec7(idx1)*gamma(t)*xit(k,i)
 1220           continue
                
c               palptgam
                idx3 = offset2 + t
                retCpd(idx2,idx3)=retCpd(idx2,idx3)-
     $               tmpvec7(idx1)*xitalp(k,t)-tmpmat4(t,idx1)
 1230         continue
              
c             palptbet
              do 1240 t=1,pv
                idx3 = offset1 + t
                retCpd(idx2,idx3)=retCpd(idx2,idx3)-tmpvec7(idx1)*vi(t)
 1240     continue      
c          call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c          return
          
          
c         quantities for cross products of (ph0t, palp, pbet, pgam) by psigb
          do 1260 l=1,q
            tmpsum1=0.0d0
            do 1250 s=1,2
              tmpsum1=tmpsum1+zit(k,s)*Eexpzb((l-1)*2+s,k)
 1250       continue
            tmpvec6(l)=tmpsum1
 1260     continue
          
          do 1300 l=1,qba
            tmpsum1=Eexpz(k)*vecwb(l)
            do 1270 s=1,qba
              tmpsum1=tmpsum1-Eexpzbb(s,k)*krwb2(s,l)
 1270       continue
            tmpvec4(l)=tmpsum1
            
            do 1300 s=1,q
              do 1280 t=1,px
                tmpmat4(l,(s-1)*px+t)=
     $                gamma(s)*xit(k,t)*tmpvec4(l)
 1280         continue
      
              tmpsum2=xitalp(k,s)*tmpvec4(l)+tmpvec6(s)*vecwb(l)
              do 1290 t=1,qba
                tmpsum2=tmpsum2-Eexpz2b3((t-1)*q+s,k)*krwb2(t,l)
 1290         continue
              tmpmat1(s,l)=tmpsum2
 1300     continue
          
          offset1 = ne1 + nfixpar
          offset2 = ne1 + pxtq
          offset3 = ne1 + pxtq + pv
          do 1360 l=1,qbu
            idx1=offset1+l
            
c           ph0ttsigb
            tmpsum1=0.0d0
            do 1310 s=1,qba
              tmpsum1=tmpsum1+tmpvec4(s)*Dqb(s,l)
 1310       continue
            retCpd(k,idx1)=retCpd(k,idx1)+tmpsum1*tmpvec1(k)/2.0d0
            
c           pbettsigb            
            tmpval1=tmpsum1*tmpvec2(k)/2.0d0
            do 1320 s=1,pv
              idx2 = offset2+s
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+tmpval1*vi(s)
 1320       continue
            
c           palptsigb
            do 1340 s=1,pxtq
              idx2=ne1+s
              tmpsum1=0.0d0
              do 1330 t=1,qba
                tmpsum1=tmpsum1+tmpmat4(t,s)*Dqb(t,l)
 1330         continue
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $           tmpsum1*tmpvec2(k)/2.0d0
 1340       continue
            
c           pgamtsigb
            do 1360 s=1,q
              idx2=offset3 + s
              tmpsum1=0.0d0              
              do 1350 t=1,qba
                tmpsum1=tmpsum1+tmpmat1(s,t)*Dqb(t,l)
 1350         continue
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $             tmpsum1*tmpvec2(k)/2.0d0      
 1360     continue
c          call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c          call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c          call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)
c          return
          
          
c         quantities for cross product of (ph0t, palp, pbet, pgam) by psige
          
          call veczero(qea,tmpvec3)
          call matzero(ldtm2,q,qea,tmpmat2)
          do 1370 j=1,ni
            do 1370 l=1,qea
              tmpvec3(l)=tmpvec3(l)+Eexpzrr((j-1)*qea+l,k)
              do 1370 s=1,q
                 tmpmat2(s,l)=tmpmat2(s,l)+
     $                Eexz2br2((j-1)*qecb+(l-1)*q+s,k)
 1370     continue
          
c          do ? l=1,q
c            tmpsum1=0.0d0
c            do ? s=1,2
c              tmpsum1=tmpsum1+zit(k,s)*Eexpzb((l-1)*2+s,k)
c ?          continue
c            tmpvec6(l)=tmpsum1
c ?        continue
          
          do 1410 l=1,qea
            tmpsum1=ni*Eexpz(k)*vecwe(l)
            do 1380 s=1,qea
              tmpsum1=tmpsum1-tmpvec3(s)*krwe2(s,l)
 1380       continue
            tmpvec4(l)=tmpsum1
            
            do 1410 s=1,q
              do 1390 t=1,px
                tmpmat4(l,(s-1)*px+t)=
     $                gamma(s)*xit(k,t)*tmpvec4(l)
 1390         continue
              
              tmpsum2=xitalp(k,s)*tmpvec4(l)+ni*tmpvec6(s)*vecwe(l)
              do 1400 t=1,qea
                tmpsum2=tmpsum2-tmpmat2(s,t)*krwe2(t,l)
 1400         continue
              tmpmat1(s,l)=tmpsum2
 1410     continue
          
          
          offset1 = ne1 + nfixpar + qbu
          offset2 = ne1 + pxtq
          offset3 = ne1 + pxtq + pv
          do 1470 l=1,qeu
            idx1=offset1+l
            
c           ph0ttsige
            tmpsum1=0.0d0
            do 1420 s=1,qea
              tmpsum1=tmpsum1+tmpvec4(s)*Dqe(s,l)
 1420       continue
            retCpd(k,idx1)=retCpd(k,idx1)+tmpsum1*tmpvec1(k)/2.0d0
            
c           pbettsige            
            tmpval1=tmpsum1*tmpvec2(k)/2.0d0
            do 1430 s=1,pv
              idx2 = offset2+s
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+tmpval1*vi(s)
 1430       continue
            
c           palptsige
            do 1450 s=1,pxtq
              idx2=ne1+s
              tmpsum1=0.0d0
              do 1440 t=1,qea
                tmpsum1=tmpsum1+tmpmat4(t,s)*Dqe(t,l)
 1440         continue
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $           tmpsum1*tmpvec2(k)/2.0d0
 1450       continue
            
c           pgamtsige
            do 1470 s=1,q
              idx2=offset3 + s
              tmpsum1=0.0d0              
              do 1460 t=1,qea
                tmpsum1=tmpsum1+tmpmat1(s,t)*Dqe(t,l)
 1460         continue
              retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $             tmpsum1*tmpvec2(k)/2.0d0
        
 1470   continue
c        call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c        call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)
c        call dblepr('retCpd(,66)',-1,retCpd(1,66),ntotpar)
c        call dblepr('retCpd(,75)',-1,retCpd(1,75),ntotpar)
c        call dblepr('retCpd(,76)',-1,retCpd(1,76),ntotpar)
c        call dblepr('retCpd(,77)',-1,retCpd(1,77),ntotpar)
c        call dblepr('retCpd(,78)',-1,retCpd(1,78),ntotpar)
c        return
		
        
ccccccccccccccccccccccccccccccccc
      
c       deltai indep otherwise tid > 0 - 2 dimensional
        
c        call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c        call dblepr('tmpvec1',-1,tmpvec1,ne1)
c        call dblepr('Eexpzz',-1,Eexpzz,tid*(tid+1)/2)
c        call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c        call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)
c		return
		
        offset1 = ne1 + pxtq
        offset2 = ne1 + pxtq + pv
        do 1580 k1=1,tid
          do 1580 k2=1,tid
            kmin = min(k1,k2)
            kmax = max(k1,k2)
            idxrow = kmax*(kmax-1)/2 + kmin
            
            tmpval1 = tmpvec1(k1)*tmpvec2(k2)*Eexpzz(idxrow)
            tmpval2 = tmpvec2(k1)*tmpvec2(k2)*Eexpzz(idxrow)

c           ph0tth0t            
            retCpd(k1,k2)=retCpd(k1,k2)+
     $           tmpvec1(k1)*tmpvec1(k2)*Eexpzz(idxrow)
            
            do 1490 l=1,q
              do 1480 s=1,px
                idx1 = ne1 + (l-1)*px + s
                
c               ph0ttalp
                retCpd(k1,idx1)=retCpd(k1,idx1)+
     $               tmpval1*gamma(l)*xit(k2,s)
                
c               palptalp
                do 1480 t=1,q
                  do 1480 i=1,px
                    idx2 = ne1 + (t-1)*px + i
                    retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $                   tmpval2*gamma(t)*xit(k1,i)*gamma(l)*xit(k2,s)
 1480         continue
              
c             ph0ttgam (more quantities are added later)
              idx1 = offset2 + l
              retCpd(k1,idx1)=retCpd(k1,idx1)+tmpval1*xitalp(k2,l)
 1490       continue
            
            do 1510 l=1,pv
              idx1 = offset1 + l
c             ph0ttbet
              retCpd(k1,idx1)=retCpd(k1,idx1)+tmpval1*vi(l)
              
c             palptbet
              do 1500 s=1,q
                do 1500 t=1,px
                  idx2 = ne1 + (s-1)*px + t
                  retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $                 tmpval2*gamma(s)*xit(k1,t)*vi(l)
 1500         continue
              
c             pbettbet
              do 1510 s=1,pv
                idx2 = offset1 + s
                retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $                 tmpval2*vi(s)*vi(l)
 1510       continue
            
            
c           quantities for pgam related cross products      
            do 1530 l=1,q
              tmpsum1=0.0d0
              tmpsum2=0.0d0
              do 1520 s=1,2
                tmpsum1=tmpsum1+zit(k2,s)*Eexpzzb((l-1)*2+s,idxrow)
                tmpsum2=tmpsum2+zit(k1,s)*Eexpzzb((l-1)*2+s,idxrow)
 1520         continue
              tmpvec6(l) = tmpsum1
              tmpvec4(l) = tmpsum2
 1530       continue
            
            if (k1.le.k2) then
              do 1540 l=1,q
                do 1540 s=1,q
                  tmpmat5(s,l)=Eexpz4b2((l-1)*q+s,idxrow)
 1540         continue
            else
              do 1550 l=1,q
                do 1550 s=1,q
                  tmpmat5(s,l)=Eexpz4b2((s-1)*q+l,idxrow)
 1550         continue
            end if
            
            
            tmpval1 = tmpvec1(k1)*tmpvec2(k2)
            tmpval3 = tmpvec2(k1)*tmpvec2(k2)
            do 1580 l=1,q
              idx1 = offset2 + l

c             ph0ttgam
              retCpd(k1,idx1)=retCpd(k1,idx1)+tmpval1*tmpvec6(l)
              
c             palptgam
              do 1570 s=1,q
                do 1560 t=1,px
                  idx2 = ne1 + (s-1)*px + t
                  retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $                 tmpval2*gamma(s)*xit(k1,t)*xitalp(k2,l)+
     $                 tmpval3*gamma(s)*xit(k1,t)*tmpvec6(l)
 1560           continue
                
c               pgamtgam
                idx3 = offset2 + s
                retCpd(idx3,idx1)=retCpd(idx3,idx1)+
     $               tmpval2*xitalp(k1,s)*xitalp(k2,l)+
     $               tmpval3*(xitalp(k1,s)*tmpvec6(l)+
     $               tmpvec4(s)*xitalp(k2,l)+
     $               tmpmat5(s,l))
 1570         continue
              
c             pbettgam
              do 1580 s=1,pv
                idx2 = offset1 + s
                retCpd(idx2,idx1)=retCpd(idx2,idx1)+
     $               tmpval2*vi(s)*xitalp(k2,l)+
     $               tmpval3*vi(s)*tmpvec6(l)
                
 1580   continue
        
      end if
c     end of if (tid > 0)
c      call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c      call dblepr('retCpd(,tid)',-1,retCpd(1,tid),ntotpar)
c      call dblepr('retCpd(,57)',-1,retCpd(1,57),ntotpar)
c      call dblepr('retCpd(,62)',-1,retCpd(1,62),ntotpar)
c      call dblepr('retCpd(,63)',-1,retCpd(1,63),ntotpar)
c      call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c      call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)
c      call dblepr('retCpd(,66)',-1,retCpd(1,66),ntotpar)
c      call dblepr('retCpd(,75)',-1,retCpd(1,75),ntotpar)
c      call dblepr('retCpd(,76)',-1,retCpd(1,76),ntotpar)
c      call dblepr('retCpd(,77)',-1,retCpd(1,77),ntotpar)
c      call dblepr('retCpd(,78)',-1,retCpd(1,78),ntotpar)
c      return      
      
	  
c     flip the upper triangles to lower triangles on retInfo and retCpd
      do 1590 k=1,ntotpar
        do 1590 l=1,k
          if (l.ne.k) then
            retInfo(k,l) = retInfo(l,k)
            retCpd(k,l) = retCpd(l,k)
          end if
 1590 continue			
c      call dblepr('retCpd(,1)',-1,retCpd(1,1),ntotpar)
c      call dblepr('retCpd(,tid)',-1,retCpd(1,tid),ntotpar)  	  
c      call dblepr('retCpd(,57)',-1,retCpd(1,57),ntotpar)
c      call dblepr('retCpd(,62)',-1,retCpd(1,62),ntotpar)	  
c      call dblepr('retCpd(,63)',-1,retCpd(1,63),ntotpar)	  
c      call dblepr('retCpd(,64)',-1,retCpd(1,64),ntotpar)
c      call dblepr('retCpd(,65)',-1,retCpd(1,65),ntotpar)	  	  
c      call dblepr('retCpd(,66)',-1,retCpd(1,66),ntotpar)
c      call dblepr('retCpd(,75)',-1,retCpd(1,75),ntotpar)	  
c      call dblepr('retCpd(,76)',-1,retCpd(1,76),ntotpar)
c      call dblepr('retCpd(,77)',-1,retCpd(1,77),ntotpar)
c      call dblepr('retCpd(,78)',-1,retCpd(1,78),ntotpar)      
      
      return
      end
      
      
C     Sum the contributions from each subject to get the gradient/information
      
      subroutine sumgrhs(retGrad,ntotpar,retInfo,ldif,retCpd,ldcp,
     $     sumGrad,sumInfo,ldsi)
      
      integer ntotpar,ldif,ldcp,ldsi
      double precision retGrad(*),retInfo(ldif,*),retCpd(ldcp,*),
     $     sumGrad(*),sumInfo(ldsi,*)
      
      integer k,l
      
      do 10 k=1,ntotpar
        sumGrad(k)=sumGrad(k)+retGrad(k)
        do 10 l=1,k
          sumInfo(l,k)=sumInfo(l,k)+retInfo(l,k)-retCpd(l,k)+
     $          retGrad(l)*retGrad(k)
          if (l.ne.k) sumInfo(k,l)=sumInfo(l,k)
 10   continue
      
      return
      end

      
c     Get the Jacobian and derivative of Jacobian on log(h0t)
      subroutine getjdlh2(h0t,ne1,lh0t,Jc,ldj,ntotpar,dJ,lddj)
      
      integer ne1,ldj,lddj,ntotpar
      double precision h0t(*),lh0t(*),Jc(ldj,*),dJ(lddj,*)
      
      integer j,dlog
      
c     save the matrices J and dJ for the final deriv Jacobian data (ntotpar by ntotpar)
      call matzero(ldj,ntotpar,ntotpar,Jc)
      call matzero(lddj,ntotpar**2,ntotpar,dJ)
      do 10 j=1,ne1
        lh0t(j) = dlog(h0t(j))
        Jc(j,j) = h0t(j)
        dJ((j-1)*ntotpar+j,j) = h0t(j)
 10   continue
      
      return
      end
      
      
c     Get the Jacobian and derivative of Jacobian on log(chol) transformed
c     parameters
      
      subroutine getjdlc2(vmat,ldvm,qbe,vmat2,ldvm2,Jc,ldj,dJ,lddj)
      
      parameter (MAXQ=4)
      parameter (MAXQB=2*MAXQ)
      integer ldvm,qbe,ldvm2,ldj,lddj
      double precision vmat(ldvm,*),vmat2(ldvm2,*),Jc(ldj,*),dJ(lddj,*)
      
      integer qbeu,jpvt(MAXQB),ierr,i,j,k,idx1,dlog,dexp
      double precision tmp1(MAXQB)
      
      qbeu = qbe*(qbe+1)/2
      
c     log(chol) transformation on the parameters
      call mtxcopy(vmat,ldvm,qbe,qbe,vmat2,ldvm2)
      call dchdc(vmat2,ldvm2,qbe,tmp1,jpvt,0,ierr)
      
      do 10 i=1,qbe
        vmat2(i,i)=dlog(vmat2(i,i))
 10   continue
      
c     get the first derivative (Jacobian)
      call matzero(ldj,qbeu,qbeu,Jc)
      
      do 40 j=1,qbe
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
      call matzero(lddj,qbeu**2,qbeu,dJ)
      
      do 70 j=1,qbe
        do 70 i=1,j
          idx1=j*(j-1)/2+i
          if (i.lt.j) then
            do 50 k=1,i
              if (k.lt.i) then
                dJ((i*(i-1)/2+k-1)*qbeu+idx1,j*(j-1)/2+k)=1.0d0
                dJ((j*(j-1)/2+k-1)*qbeu+idx1,i*(i-1)/2+k)=1.0d0
              else
                dJ((i*(i-1)/2+i-1)*qbeu+idx1,i*(i-1)/2+i)=
     $                dexp(vmat2(i,i))*vmat2(i,j)
                dJ((i*(i-1)/2+i-1)*qbeu+idx1,j*(j-1)/2+i)=
     $                dexp(vmat2(i,i))
                dJ((j*(j-1)/2+i-1)*qbeu+idx1,i*(i-1)/2+i)=
     $                dexp(vmat2(i,i))
              end if
 50         continue
          else
            do 60 k=1,i
              if (k.lt.i) then
                dJ((i*(i-1)/2+k-1)*qbeu+idx1,i*(i-1)/2+k)=2.0d0
              else
                dJ((i*(i-1)/2+i-1)*qbeu+idx1,i*(i-1)/2+i)=
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
      
      
c     Get the gradient and information on the full parameter vector with h0t log-transformed
c     and variance parameters log(chol) transformed
      
      subroutine getgdiff(h0t,ne1,coefs,nfixpar,sigb,ldsb,qb,lcholb,
     $     ldlb,sige,ldse,lchole,ldle,Jb,ldjb,dJb,lddjb,Je,ldje,dJe,
     $     lddje,origParm,origGrad,origInfo,ldoi,P2,ldp2,natrParm,
     $     natrGrad,natrInfo,ldni,J2,ldj2,dJ2,lddj2,wkinfo,ldwi)
      
      integer ne1,nfixpar,ldsb,qb,ldlb,ldse,ldle,ldjb,lddjb,
     $     ldje,lddje,ldoi,ldp2,ldni,ldj2,lddj2,ldwi
      double precision h0t(*),coefs(*),sigb(ldsb,*),lcholb(ldlb,*),
     $     sige(ldse,*),lchole(ldle,*),
     $     Jb(ldjb,*),dJb(lddjb,*),Je(ldje,*),dJe(lddje,*),origParm(*),
     $     origGrad(*),origInfo(ldoi,*),P2(ldp2,*),natrParm(*),
     $     natrGrad(*),natrInfo(ldni,*),J2(ldj2,*),dJ2(lddj2,*),
     $     wkinfo(ldwi,*)
      
      integer q,qbu,qeu,ntotpar,i,j,k,l,offset1,idx1,idx2
      double precision tmpsum1,dlog
      
      q = qb / 2
      qbu = qb*(qb+1)/2
      qeu  = q*(q+1)/2
      ntotpar = ne1 + nfixpar + qbu + qeu
      
c     assign the original and natural parameters
      do 10 i=1,ne1+nfixpar
        if (i.le.ne1) then
          origParm(i)=h0t(i)
          natrParm(i)=dlog(h0t(i))
        else  
          origParm(i)=coefs(i-ne1)
          natrParm(i)=coefs(i-ne1)
        end if
 10   continue
      
      offset1 = ne1 + nfixpar
      do 20 j=1,qb 
        do 20 k=1,j
          idx1=offset1+k*(2*qb+1-k)/2-qb+j
c         count on lower triangular matrix position (j,k) from bottom to top and right to left
          origParm(idx1)=sigb(k,j)
          natrParm(idx1)=lcholb(k,j)
 20   continue
      
      offset1 = ne1 + nfixpar + qbu
      do 30 j=1,q 
        do 30 k=1,j
          idx1=offset1+k*(2*q+1-k)/2-q+j
c         count on lower triangular matrix position (j,k) from bottom to top and right to left
          origParm(idx1)=sige(k,j)
          natrParm(idx1)=lchole(k,j)
 30   continue
      
c      call dblepr('natrParm',-1,natrParm,24)
      
      
c     assign the full Jacobian
c     J2 is pre-filled with Jacobian for lh0t
      
      do 40 i=ne1+1,ne1+nfixpar
        J2(i,i)=1.0d0
 40   continue
      
      offset1 = ne1 + nfixpar
      do 50 j=1,qbu
        do 50 i=1,qbu
          J2(offset1+i,offset1+j)=Jb(i,j)
 50   continue
      
      offset1 = ne1 + nfixpar + qbu
      do 60 j=1,qeu
        do 60 i=1,qeu
          J2(offset1+i,offset1+j)=Je(i,j)
 60   continue      
c      call dblepr('J2(,1)',-1,J2(1,1),ntotpar)
c      call dblepr('J2(,56)',-1,J2(1,56),ntotpar)
c      call dblepr('J2(,57)',-1,J2(1,57),ntotpar)
c      call dblepr('J2(,65)',-1,J2(1,65),ntotpar)
c      call dblepr('J2(,66)',-1,J2(1,66),ntotpar)
c      call dblepr('J2(,75)',-1,J2(1,75),ntotpar)
c      call dblepr('J2(,76)',-1,J2(1,76),ntotpar)
c      call dblepr('J2(,78)',-1,J2(1,78),ntotpar)
c      return
      
c     assign the full derivative of Jacobian
c     dJ2 is pre-filled with dJ for lh0t
      
      offset1 = ne1 + nfixpar            
      do 70 j=1,qbu
        do 70 i=1,qbu**2
          idx1 = (i-1)/qbu + 1
          idx2 = i - (idx1-1)*qbu
          dJ2((offset1+idx1-1)*ntotpar+offset1+idx2,offset1+j)=dJb(i,j)
 70   continue
      
      offset1 = ne1 + nfixpar + qbu            
      do 80 j=1,qeu
        do 80 i=1,qeu**2
          idx1 = (i-1)/qeu + 1
          idx2 = i - (idx1-1)*qeu
          dJ2((offset1+idx1-1)*ntotpar+offset1+idx2,offset1+j)=dJe(i,j)
 80   continue      
c      call dblepr('dJ2(,1)',-1,dJ2(1,1),ntotpar**2)
c      call dblepr('dJ2(,56)',-1,dJ2(1,56),ntotpar**2)
c      call dblepr('dJ2(,57)',-1,dJ2(1,57),ntotpar**2)
c      call dblepr('dJ2(,65)',-1,dJ2(1,65),ntotpar**2)
c      call dblepr('dJ2(,66)',-1,dJ2(1,66),ntotpar**2)
c      call dblepr('dJ2(,75)',-1,dJ2(1,75),ntotpar**2)
c      call dblepr('dJ2(,76)',-1,dJ2(1,76),ntotpar**2)
c      call dblepr('dJ2(,78)',-1,dJ2(1,78),ntotpar**2)
c      return
      
      
c     get the gradient and information on natural scale
      call matzero(ldwi,ntotpar,ntotpar,wkInfo)
      
      do 130 j=1,ntotpar
        tmpsum1=0.0d0
        do 90 k=1,ntotpar
          do 90 l=1,ntotpar
            tmpsum1=tmpsum1+origGrad(l)*J2(l,k)*P2(k,j)
 90     continue
        natrGrad(j)=tmpsum1
        
        do 130 i=1,j
          tmpsum1=0.0d0
          do 120 k=1,ntotpar
            if (k.eq.i) then
              do 100 l=1,ntotpar
                tmpsum1=tmpsum1+J2(l,i)*origInfo(l,k)*J2(k,j)-
     $                origGrad(l)*dJ2((i-1)*ntotpar+l,j)
 100          continue
            else
              do 110 l=1,ntotpar
                tmpsum1=tmpsum1+J2(l,i)*origInfo(l,k)*J2(k,j)
 110          continue
            end if
 120      continue
          wkinfo(i,j)=tmpsum1
          if (i.ne.j) wkinfo(j,i)=tmpsum1
 130  continue
c      call dblepr('wkinfo(,1)',-1,wkinfo(1,1),ntotpar)
c      call dblepr('wkinfo(,78)',-1,wkinfo(1,78),ntotpar)
c      return
      
      do 150 j=1,ntotpar
        do 150 i=1,j
          tmpsum1=0.0d0
          do 140 k=1,ntotpar
            do 140 l=1,ntotpar
              tmpsum1=tmpsum1+P2(l,i)*wkinfo(l,k)*P2(k,j)
 140      continue      
          natrInfo(i,j)=tmpsum1
          if (i.ne.j) natrInfo(j,i)=tmpsum1
 150  continue
c      call dblepr('natrInfo(,1)',-1,natrInfo(1,1),ntotpar)
c      call dblepr('natrInfo(,78)',-1,natrInfo(1,78),ntotpar)


c     permute the origGrad and origInfo to lower triangular "major" format
      do 170 j=1,ntotpar
        tmpsum1=0.0d0
        do 160 k=1,ntotpar
          tmpsum1=tmpsum1+origGrad(k)*P2(k,j)
 160    continue
        wkinfo(j,1)=tmpsum1
 170  continue
      call veccopy(wkinfo(1,1),ntotpar,origGrad)
	  
      do 190 j=1,ntotpar
        do 190 i=1,j
          tmpsum1=0.0d0
          do 180 k=1,ntotpar
            do 180 l=1,ntotpar
              tmpsum1=tmpsum1+P2(l,i)*origInfo(l,k)*P2(k,j)
 180      continue      
          wkinfo(i,j)=tmpsum1
          if (i.ne.j) wkinfo(j,i)=tmpsum1
 190  continue
      call mtxcopy(wkinfo,ldwi,ntotpar,ntotpar,origInfo,ldoi)
	  
      return
      end
      
      
C     Main subroutine for getting information matrix on survival outcome with left-censored 
c     covariates
      
      subroutine survlinf(time,event,id1,m,Q,ldq,sumn,qe,D,ldd,id2,repl,
     $     V,ldv,pv,X,ldx,px,Z,ldz,samp_pd,eb,ldeb,eyc,ldeyc,h0t,
     $     alphamat,ldam,beta,gamma,sigb,ldsb,sige,ldse,
     $     maxiter2,nburnin,nsize,origParm,origGrad,origInfo,ldoi,
     $     natrParm,natrGrad,natrInfo,ldni,wkarr,iwkarr)
      
      parameter (MAXQ=4,MAXPX=50,MAXPV=20)
      parameter (MAXQEA=MAXQ**2,MAXQEU=MAXQ*(MAXQ+1)/2,MAXQB=2*MAXQ,
     $     MAXPXTQ=MAXPX*MAXQ,MAXPVPQ=MAXPV+MAXQ)
      parameter (MAXQBA=MAXQB*MAXQB,MAXQBU=MAXQB*(MAXQB+1)/2)
      integer event(*),id1(*),m,ldq,sumn,qe,ldd,id2(*),repl(*),ldv,pv,
     $     ldx,px,ldz,ldeb,ldeyc,ldam,ldsb,ldse,maxiter2,
     $     nburnin,nsize,ldoi,ldni,iwkarr(*)
      integer D(ldd,*)
      double precision time(*),Q(ldq,*),V(ldv,*),X(ldx,*),Z(ldz,*),
     $     samp_pd,eb(ldeb,*),eyc(ldeyc,*),h0t(*),alphamat(ldam,*),
     $     beta(*),gamma(*),sigb(ldsb,*),sige(ldse,*),origParm(*),
     $     origGrad(*),origInfo(ldoi,*),natrParm(*),natrGrad(*),
     $     natrInfo(ldni,*),wkarr(*)
      
      integer i,j,k,l,qb,pxtq,pvpq,nfixpar,qba,qbu,qea,qeu,qecb,ne0,ne1,
     $     ne1p1,ne12u,ntotpar,maxiter3,cid,idx1,idx2,offset1,offset2,
     $     offset3,maxrepl,maxni,maxnisq,sumtd,sumtdpm,tid,tidp1,ni,
     $     dila,dilj,replidx,replidx2,sseqb(MAXQB),sseqe(MAXQ),ierr,
     $     klen,ks,kt,int,min,max,mod
      integer ina,ietimes0,ietimes1,itd,idla,idl,it2imap,ixt,izt,iztz,
     $     ip2,idi,ivbeta,
     $     iebp,ieycp,ixi,izi,iqic,ieyic,iyic0,ixit,izit,izitz,
     $     ibrarr,iwkzb,ier,ier2c,ier2,ier3,ier4,ierbb,ieb2r2,
     $     ierzb,iebzrr,iexpz,iexpzb,iexz3btb,iexpzz,iexpzr,iexpzzb,
     $     iexpzrzb,iexpz4b2,iexpzbb,iexpz2b3,iexpzrr,iexz2br2,iexpz3b2,
     $     itmpvec1,itmpvec2,itmpmat1,itmpmat2,iretGrad,iretInfo,
     $     iretCpd,ixitalp,iJ2,idJ2,ilh0t,iwkinfo
      integer ldxt,ldzt,ldztz,ldp2,lddi,ldebp,ldeycp,ldxi,
     $     ldzi,ldqic,ldeyic,ldyic0,ldxit,ldzit,ldzitz,ldbrarr,
     $     ldwkzb,lder,lder2c,lder2,lder3,lder4,lderbb,ldeb2r2,lderzb,
     $     ldebzrr,ldexpzb,ldz3btb,ldexpzr,ldexpzzb,ldzrzb,ldz4b2,
     $     ldexpzbb,ldz2b3,ldexpzrr,ldz2br2,ldz3b2,ldtm1,ldtm2,ldif,
     $     ldcp,ldxa,ldj2,lddj2,ldwi
      integer ldsb2,ldse2,ldtqb,ldtqe,lddqb,lddqe,ldpqb,ldpqe,ldgm,
     $     ldbc,ldbh,ldeb2,ldeb4,ldebzzb,ldebzbb,ldlb,ldjb,lddjb,ldle,
     $     ldje,lddje,ldkb,ldke,intvec(4),ldvec(40)
      double precision eps1b,tmp1,etimes1j,tqbqb(MAXQBA,MAXQBA),
     $     tqeqe(MAXQEA,MAXQEA),dqb(MAXQBA,MAXQBU),dqe(MAXQEA,MAXQEU),
     $     pqb(MAXQBU,MAXQBU),pqe(MAXQEU,MAXQEU),sigbinv(MAXQB,MAXQB),
     $     vecwb(MAXQBA),sigeinv(MAXQ,MAXQ),vecwe(MAXQEA),
     $     krwb2(MAXQBA,MAXQBA),krwe2(MAXQEA,MAXQEA),gmat(MAXQ,MAXQ),
     $     cstd(MAXQ),ebi0(MAXQB),bimode(MAXQB),bicurv(MAXQB,MAXQB),
     $     bhinv(MAXQB,MAXQB),bim0(MAXQB),vi(MAXPV),ebi(MAXQB),
     $     ebbi(MAXQB,MAXQB),eb4(MAXQBA,MAXQBA),ebzzb(MAXQ,MAXQ),
     $     ebzbb(MAXQ,MAXQBA),lcholb(MAXQB,MAXQB),Jb(MAXQBU,MAXQBU),
     $     dJb(MAXQBU*MAXQBU,MAXQBU),lchole(MAXQ,MAXQ),
     $     Je(MAXQEU,MAXQEU),dJe(MAXQEU*MAXQEU,MAXQEU),
     $     betas(MAXPXTQ+MAXPVPQ),det(2)
      character kstr*2,labstr*20
      
c      return
	  
      qb = 2*qe
      pxtq = px*qe
      pvpq = pv + qe
      nfixpar = pxtq + pvpq
      qba = qb**2
      qbu = qb*(qb+1)/2
      qea = qe**2
      qeu = qe*(qe+1)/2
      qecb = qe**3
      
      ldsb2 = MAXQB
      ldse2 = MAXQ
      
      ldtqb = MAXQBA
      ldtqe = MAXQEA
      lddqb = MAXQBA
	  ldkb = MAXQBA
      lddqe = MAXQEA
      ldpqb = MAXQBU	  
      ldpqe = MAXQEU
      ldke = MAXQEA	  
      ldgm = MAXQ
      ldbc = MAXQB
      ldbh = MAXQB
      ldeb2 = MAXQB
      ldeb4 = MAXQBA
      ldebzzb = MAXQ
      ldebzbb = MAXQ
      ldlb = MAXQB
      ldjb = MAXQBU
      lddjb = MAXQBU*MAXQBU
      ldle = MAXQ
      ldje = MAXQEU
      lddje = MAXQEU*MAXQEU
      
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
c      call intpr('maxiter2',-1,maxiter2,1)
c      call intpr('nburnin',-1,nburnin,1)
c      call intpr('nsize',-1,nsize,1)
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
c      return
      
      ietimes1 = ietimes0 + ne0
      ne1 = 1
      wkarr(ietimes1+ne1-1) = wkarr(ietimes0)
      do 50 k=2,ne0
        if (wkarr(ietimes0+k-1).ne.wkarr(ietimes1+ne1-1)) then
          ne1 = ne1 + 1
          wkarr(ietimes1+ne1-1) = wkarr(ietimes0+k-1)
        end if
 50   continue
      ne1p1 = ne1 + 1
      ne12u = ne1*(ne1+1)/2
c      call intpr('ne1',-1,ne1,1)
c      call dblepr('etimes1',-1,wkarr(ietimes1),ne1)
c      return
      
      ntotpar = ne1 + nfixpar + qbu + qeu
c      call intpr('ntotpar',-1,ntotpar,1)
c      return
      
      
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
c      return
      
      
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
      maxnisq = maxni**2
c      call intpr('maxni',-1,maxni,1)
c      call intpr('dla',-1,iwkarr(idla),m)
c      call intpr('dl',-1,iwkarr(idl),sumn)
c      return
      
      
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
c      return
      
      
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
c      return
      
      
c     generate the communication matrix for sigb
      call matzero(ldtqb,qba,qba,tqbqb)
      do 210 i=1,qb
        do 210 j=1,qb
          tqbqb((i-1)*qb+j,(j-1)*qb+i) = 1.0d0
 210  continue
c      call dblepr('tqbqb(,1)',-1,tqbqb(1,1),qba)
c      call dblepr('tqbqb(,16)',-1,tqbqb(1,16),qba)
c      return
      
c     generate the communication matrix for sige
      call matzero(ldtqe,qea,qea,tqeqe)
      do 215 i=1,qe
        do 215 j=1,qe
          tqeqe((i-1)*qe+j,(j-1)*qe+i) = 1.0d0
 215  continue
c      call dblepr('tqeqe(,1)',-1,tqeqe(1,1),qea)
c      call dblepr('tqeqe(,2)',-1,tqeqe(1,2),qea)
c      call dblepr('tqeqe(,3)',-1,tqeqe(1,3),qea)
c      call dblepr('tqeqe(,4)',-1,tqeqe(1,4),qea)
c      return
      
      
c     generate the duplication matrix for sigb
      call matzero(lddqb,qba,qbu,dqb)
      sseqb(1)=0
      do 220 j=1,qb-1
        sseqb(j+1)=sseqb(j)+j
 220  continue
      
      do 230 j=1,qb
        do 230 i=1,qb
          idx1 = (j-1)*qb + i
          if (i.le.j) then
            dqb(idx1,sseqb(j)+i) = 1.0d0
          else
            dqb(idx1,sseqb(i)+j) = 1.0d0
          end if
 230  continue
c      call intpr('sseqb',-1,sseqb,qb)
c      call dblepr('dqb(,1)',-1,dqb(1,1),qba)
c      call dblepr('dqb(,2)',-1,dqb(1,2),qba)
c      call dblepr('dqb(,10)',-1,dqb(1,10),qba)
c      return
      
c     generate the duplication matrix for sige
      call matzero(lddqe,qea,qeu,dqe)
      sseqe(1)=0
      do 232 j=1,qe-1
        sseqe(j+1)=sseqe(j)+j
 232  continue
      
      do 235 j=1,qe
        do 235 i=1,qe
          idx1 = (j-1)*qe + i
          if (i.le.j) then
            dqe(idx1,sseqe(j)+i) = 1.0d0
          else
            dqe(idx1,sseqe(i)+j) = 1.0d0
          end if
 235  continue
c      call intpr('sseqe',-1,sseqe,qe)
c      call dblepr('dqe(,1)',-1,dqe(1,1),qea)
c      call dblepr('dqe(,2)',-1,dqe(1,2),qea)
c      call dblepr('dqe(,3)',-1,dqe(1,3),qea)
c      return
      
      
c     generate the permutation matrix for sigb
      call matzero(ldpqb,qbu,qbu,pqb)
      do 240 j=1,qb
        do 240 i=1,j
          kt = (j-1)*j/2 + i
          ks = qb*(i-1) - (i-1)*(i-2)/2 + (j-i+1)
          pqb(kt,ks) = 1.0d0
 240  continue
c      call dblepr('pqb(,1)',-1,pqb(1,1),qbu)
c      call dblepr('pqb(,2)',-1,pqb(1,2),qbu)
c      call dblepr('pqb(,10)',-1,pqb(1,10),qbu)
c      return
      
c     generate the permutation matrix for sige
      call matzero(ldpqe,qeu,qeu,pqe)
      do 245 j=1,qe
        do 245 i=1,j
          kt = (j-1)*j/2 + i
          ks = qe*(i-1) - (i-1)*(i-2)/2 + (j-i+1)
          pqe(kt,ks) = 1.0d0
 245  continue
c      call dblepr('pqe(,1)',-1,pqe(1,1),qeu)
c      call dblepr('pqe(,2)',-1,pqe(1,2),qeu)
c      call dblepr('pqe(,3)',-1,pqe(1,3),qeu)
c      return
      
      
c     populate the combined permutation matrix
      ip2 = iztz + ldztz*4
      ldp2 = ntotpar
      call matzero(ldp2,ntotpar,ntotpar,wkarr(ip2))
      do 250 j=1,ne1+nfixpar
        wkarr(ip2+(j-1)*ntotpar+j-1) = 1.0d0
 250  continue
      
      offset1 = ne1+nfixpar
      do 260 j=1,qbu
        do 260 i=1,qbu
          wkarr(ip2+(offset1+j-1)*ntotpar+offset1+i-1) = pqb(i,j)
 260  continue
      
      offset1 = ne1+nfixpar+qbu
      do 265 j=1,qeu
        do 265 i=1,qeu
          wkarr(ip2+(offset1+j-1)*ntotpar+offset1+i-1) = pqe(i,j)
 265  continue
c      call dblepr('p2(,1)',-1,wkarr(ip2),ntotpar)
c      call dblepr('p2(,56)',-1,wkarr(ip2+55*ntotpar),ntotpar)
c      call dblepr('p2(,57)',-1,wkarr(ip2+56*ntotpar),ntotpar)
c      call dblepr('p2(,65)',-1,wkarr(ip2+64*ntotpar),ntotpar)
c      call dblepr('p2(,66)',-1,wkarr(ip2+65*ntotpar),ntotpar)
c      call dblepr('p2(,75)',-1,wkarr(ip2+74*ntotpar),ntotpar)
c      call dblepr('p2(,76)',-1,wkarr(ip2+75*ntotpar),ntotpar)
c      call dblepr('p2(,78)',-1,wkarr(ip2+77*ntotpar),ntotpar)
c      return
      
      
c     Start the iteration      
      idi = idl + sumn
      lddi = maxni
c     dimension of di is maxni by qe (a matrix)
      
      ivbeta = ip2 + ldp2*ntotpar
      iebp = ivbeta + m
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
      ldxit = ne1p1
      izit = ixit + ldxit*px
      ldzit = ne1p1
      izitz = izit + ldzit*2
      ldzitz = ne1
      ibrarr = izitz + ldzitz*4
      ldbrarr = qb + qe*maxni
c      ieycbi = ibiarr + ldbrarr*maxnsize
c      ldeycbi = maxni
      iwkzb = ibrarr + ldbrarr*nsize
      ldwkzb = qe	  
      ier = iwkzb + ldwkzb*maxni
      lder = qe
      ier2c = ier + lder*maxni
      lder2c = qea
      ier2 = ier2c + lder2c*maxnisq
      lder2 = qea
      ier3 = ier2 + lder2*maxni
      lder3 = qecb
      ier4 = ier3 + lder3*maxnisq
      lder4 = qe**4
      ierbb = ier4 + lder4*maxnisq
      lderbb = qe*qba
      ieb2r2 = ierbb + lderbb*maxni
      ldeb2r2 = qba*qea
      ierzb = ieb2r2 + ldeb2r2*maxni
      lderzb = qea
      iebzrr = ierzb + lderzb*maxni
      ldebzrr = qecb
      iexpz = iebzrr + ldebzrr*maxni
      iexpzb = iexpz + ne1
      ldexpzb = qb
      iexz3btb = iexpzb + ldexpzb*ne1
      ldz3btb = qea
      iexpzz = iexz3btb + ldz3btb*ne1
      iexpzr = iexpzz + ne12u
      ldexpzr = qe*maxni
      iexpzzb = iexpzr + ldexpzr*ne1
      ldexpzzb = qb
      iexpzrzb = iexpzzb + ldexpzzb*ne12u
      ldzrzb = qea*maxni
      iexpz4b2 = iexpzrzb + ldzrzb*ne1
      ldz4b2 = qea
      iexpzbb = iexpz4b2 + ldz4b2*ne12u
      ldexpzbb = qba
      iexpz2b3 = iexpzbb + ldexpzbb*ne1
      ldz2b3 = qe*qba
      iexpzrr = iexpz2b3 + ldz2b3*ne1
      ldexpzrr = qea*maxni
      iexz2br2 = iexpzrr + ldexpzrr*ne1
      ldz2br2 = qecb*maxni
      iexpz3b2 = iexz2br2 + ldz2br2*ne1
      ldz3b2 = qea
      itmpvec1 = iexpz3b2 + ldz3b2*ne1
      itmpvec2 = itmpvec1 + ne1
      itmpmat1 = itmpvec2 + ne1
      ldtm1 = qea
      itmpmat2 = itmpmat1 + ldtm1*maxni
      ldtm2 = qe
      iretGrad = itmpmat2 + ldtm2*ne1p1
      iretInfo = iretGrad + ntotpar
      ldif = ntotpar
      iretCpd = iretInfo + ldif*ntotpar
      ldcp = ntotpar
      ixitalp = iretCpd + ldcp*ntotpar
      ldxa = ne1p1
      iJ2 = ixitalp + ldxa*qe
      ldj2 = ntotpar
      idJ2 = iJ2 + ldj2*ntotpar
      lddj2 = ntotpar**2
      ilh0t = idJ2 + lddj2*ntotpar
      iwkinfo = ilh0t + ne1
      ldwi = ntotpar
c     dimension of wkinfo is ntotpar by ntotpar (a matrix)
      
c      do 440 iter=1,maxiter
c        call ifprintf(iter,0,1)
      call intpr('Integrating the information matrix...',-1,0,0)
      
      offset1 = pxtq + pv
      do 270 k=1,qe
        betas(offset1+k)=gamma(k)
        do 270 l=1,px
          idx1 = (k-1)*px + l
          betas(idx1)=alphamat(l,k)
 270  continue
      
      do 280 k=1,pv
        betas(pxtq+k)=beta(k)
 280  continue
c      call dblepr('betas',-1,betas,nfixpar)
c      return
      
      do 290 i=1,m
        idx1 = ivbeta + i - 1 
        wkarr(idx1)=0.0d0
        do 290 k=1,pv
          wkarr(idx1)=wkarr(idx1)+V(i,k)*beta(k)
 290  continue
c      call dblepr('Vbeta',-1,wkarr(ivbeta),m)
c      return
      
      call mtxcopy(sigb,ldsb,qb,qb,sigbinv,ldsb2)
      call invert(ldsb2,qb,sigbinv,det,ierr)
c      call dblepr('sigbinv(,1)',-1,sigbinv(1,1),qb)
c      call dblepr('sigbinv(,2)',-1,sigbinv(1,2),qb)
c      call dblepr('sigbinv(,3)',-1,sigbinv(1,3),qb)
c      call dblepr('sigbinv(,4)',-1,sigbinv(1,4),qb)
c      return
      
      call mtxcopy(sige,ldse,qe,qe,sigeinv,ldse2)
      call invert(ldse2,qe,sigeinv,det,ierr)
c      call dblepr('sigeinv(,1)',-1,sigeinv(1,1),qe)
c      call dblepr('sigeinv(,2)',-1,sigeinv(1,2),qe)
c      return
      
c      call matzero(ldkb,qba,qba,krwb2)
      do 300 j=1,qb
        do 300 i=1,qb
          vecwb((j-1)*qb+i)=sigbinv(i,j)
          do 300 l=1,qb
            do 300 k=1,qb
              krwb2((i-1)*qb+k,(j-1)*qb+l) = sigbinv(k,l)*sigbinv(i,j)
 300  continue
c      call dblepr('vecwb',-1,vecwb,qba)
c      call dblepr('krwb2(,1)',-1,krwb2(1,1),qba)
c      call dblepr('krwb2(,16)',-1,krwb2(1,16),qba)
c      return
      
c      call matzero(ldke,qea,qea,krwe2)
      do 305 j=1,qe
        do 305 i=1,qe
          vecwe((j-1)*qe+i)=sigeinv(i,j)
          do 305 l=1,qe
            do 305 k=1,qe
              krwe2((i-1)*qe+k,(j-1)*qe+l) = sigeinv(k,l)*sigeinv(i,j)
 305  continue
c      call dblepr('vecwe',-1,vecwe,qea)
c      call dblepr('krwe2(,1)',-1,krwe2(1,1),qea)
c      call dblepr('krwe2(,2)',-1,krwe2(1,2),qea)
c      call dblepr('krwe2(,3)',-1,krwe2(1,3),qea)
c      call dblepr('krwe2(,4)',-1,krwe2(1,4),qea)
c      return
      
      call setgibbs(sige,ldse,qe,gmat,ldgm,cstd)
c      call dblepr('gmat(,1)',-1,gmat(1,1),qe)
c      call dblepr('gmat(,2)',-1,gmat(1,2),qe)
c      call dblepr('cstd',-1,cstd,qe)
c      return
      
      call veczero(ntotpar,origGrad)
      call matzero(ldoi,ntotpar,ntotpar,origInfo)
      
      call mtxcopy(eb,ldeb,m,2,wkarr(iebp),ldebp)
      call mtxcopy(eyc,ldeyc,sumn,qe,wkarr(ieycp),ldeycp)
      
c     assign the dimension vectors	  
      intvec(1) = qb
      intvec(2) = nsize
      
      ldvec(1) = ldxit
      ldvec(2) = ldzit
      ldvec(3) = ldam
      ldvec(4) = ldzitz
      ldvec(5) = ldxi
      ldvec(6) = ldsb2
      ldvec(7) = ldse2
      ldvec(8) = ldeb2
      ldvec(9) = ldeb4
      ldvec(10) = lder
      ldvec(11) = lder2c
      ldvec(12) = lder2
      ldvec(13) = lder3
      ldvec(14) = lder4
      ldvec(15) = lderbb
      ldvec(16) = ldeb2r2
      ldvec(17) = lderzb
      ldvec(18) = ldebzzb
      ldvec(19) = ldebzbb
      ldvec(20) = ldebzrr
      ldvec(21) = ldexpzb
      ldvec(22) = ldz3btb
      ldvec(23) = ldexpzr
      ldvec(24) = ldexpzzb
      ldvec(25) = ldzrzb
      ldvec(26) = ldz4b2
      ldvec(27) = ldexpzbb
      ldvec(28) = ldz2b3
      ldvec(29) = ldexpzrr
      ldvec(30) = ldz2br2
      ldvec(31) = ldz3b2
      ldvec(32) = ldkb
      ldvec(33) = ldke
      ldvec(34) = lddqb
      ldvec(35) = lddqe
      ldvec(36) = ldtqb
      ldvec(37) = ldtqe
      ldvec(38) = ldif
      ldvec(39) = ldcp
      ldvec(40) = ldxa
      
      offset1 = 0
      offset2 = 0
      offset3 = 0
c      do 430 i=1,2
      do 430 i=1,m
        call ifprintf(i,1,0)
c        if (mod(i,5).eq.0) call ifprintf(i,1,0)
        
        ni = iwkarr(ina+i-1)  
        
        do 310 k=1,qb
          ebi0(k) = eb(i,k)
 310    continue
c        call dblepr('ebi0',-1,ebi0,qb)
c        return
        
        do 350 j=1,ni
          idx1 = offset3 + j
            
          do 320 k=1,px
             wkarr(ixi+(k-1)*ldxi+j-1) = X(idx1,k)
 320      continue
          
          do 330 k=1,2
            wkarr(izi+(k-1)*ldzi+j-1) = Z(idx1,k)
 330      continue
          
          do 350 k=1,qe
            idx2 = iqic + (k-1)*ldqic + j - 1
            wkarr(idx2) = Q(idx1,k)
            do 340 l=1,px
              wkarr(idx2)=wkarr(idx2)-X(idx1,l)*alphamat(l,k)
 340        continue
            
            wkarr(ieyic+(k-1)*ldeyic+j-1) = eyc(idx1,k)
              
            iwkarr(idi+(k-1)*lddi+j-1) = D(idx1,k)
 350    continue
c        call dblepr('xi(,1)',-1,wkarr(ixi),ldxi)
c        call dblepr('xi(,2)',-1,wkarr(ixi+ldxi),ldxi)
c        call dblepr('xi(,3)',-1,wkarr(ixi+2*ldxi),ldxi)
c        call dblepr('zi(,1)',-1,wkarr(izi),ldzi)
c        call dblepr('zi(,2)',-1,wkarr(izi+ldzi),ldzi)
c        call dblepr('qic(,1)',-1,wkarr(iqic),ldqic)
c        call dblepr('qic(,2)',-1,wkarr(iqic+ldqic),ldqic)
c        call dblepr('eyic(,1)',-1,wkarr(ieyic),ldeyic)
c        call dblepr('eyic(,2)',-1,wkarr(ieyic+ldeyic),ldeyic)
c        call intpr('di(,1)',-1,iwkarr(idi),lddi)
c        call intpr('di(,2)',-1,iwkarr(idi+lddi),lddi)
c        return
        
        tid = iwkarr(itd+i-1)
        tidp1 = tid + 1
        do 380 j=1,tidp1
          idx1 = offset1 + j
          
          do 360 k=1,px
            wkarr(ixit+(k-1)*ldxit+j-1)=wkarr(ixt+(k-1)*ldxt+idx1-1)
 360      continue
          
          do 370 k=1,2
            wkarr(izit+(k-1)*ldzit+j-1)=wkarr(izt+(k-1)*ldzt+idx1-1)
 370      continue          
 380    continue
        
        if (tid.gt.0) then
          do 390 j=1,tid
            idx2 = offset2 + j
            do 390 k=1,4
              wkarr(izitz+(k-1)*ldzitz+j-1)=
     $              wkarr(iztz+(k-1)*ldztz+idx2-1)
 390      continue
        end if
        
        do 400 k=1,pv
          vi(k) = V(i,k)
 400    continue
c        call dblepr('vi',-1,vi,pv)
c        return
        
        
        call getycbm0(wkarr(iqic),ldqic,ni,qe,iwkarr(idla+i-1),
     $       iwkarr(idl+offset3),iwkarr(idi),lddi,ebi0,qb,sigbinv,
     $       ldsb2,event(i),wkarr(ivbeta+i-1),h0t,tid,wkarr(ixit),
     $       ldxit,px,wkarr(izit),ldzit,alphamat,ldam,gamma,
     $       wkarr(ieyic),ldeyic,wkarr(izi),ldzi,sigeinv,ldse2,
     $       maxiter2,eps1b,maxiter3,gmat,ldgm,cstd,bimode,bicurv,
     $       ldbc,bhinv,ldbh,bim0,wkarr(iyic0),ldyic0)
c        if (i.eq.2) then
c          call ifprintf(0,0,0)
c          call dblepr('bimode',-1,bimode,qb)
c          call dblepr('bicurv(,1)',-1,bicurv(1,1),qb)
c          call dblepr('bicurv(,4)',-1,bicurv(1,4),qb)
c          call dblepr('bim0',-1,bim0,qb)
c          call dblepr('yic0(,1)',-1,wkarr(iyic0),ldyic0)
c          call dblepr('yic0(,2)',-1,wkarr(iyic0+ldyic0),ldyic0)          
c          return
c        end if
        
        
        call mgibbs2(wkarr(iqic),ldqic,ni,qe,iwkarr(idla+i-1),
     $       iwkarr(idl+offset3),iwkarr(idi),lddi,sigbinv,ldsb2,qb,
     $       event(i),wkarr(ivbeta+i-1),h0t,tid,wkarr(ixit),ldxit,px,
     $       wkarr(izit),ldzit,alphamat,ldam,gamma,wkarr(izi),ldzi,
     $       sigeinv,ldse2,gmat,ldgm,cstd,bimode,bicurv,ldbc,bhinv,ldbh,
     $       bim0,wkarr(iyic0),ldyic0,wkarr(ibrarr),ldbrarr,nburnin,
     $       nsize,wkarr(iwkzb),ldwkzb)
c        if (i.eq.2) then
c          call ifprintf(0,0,0)
c          call dblepr('brarr(,1)',-1,wkarr(ibrarr),qb+8)
c          call dblepr('brarr(,1000)',-1,wkarr(ibrarr+(1000-1)*ldbrarr),
c     $         qb+8)
c	      return
c        end if
        
		
c          call intpr('ldbrarr',-1,ldbrarr,1)
c          call intpr('qb',-1,qb,1)
c          call intpr('nsize',-1,nsize,1)
c          call dblepr('brarr(,1)',-1,wkarr(ibrarr),ldbrarr)
c          call dblepr('brarr(,1000)',-1,wkarr(ibrarr+(1000-1)*ldbrarr),
c     $         ldbrarr)
c          call intpr('ldzi',-1,ldzi,1)
c          call intpr('ni',-1,ni,1)
c          call dblepr('zi(,1)',-1,wkarr(izi),ldzi)
c          call dblepr('zi(,2)',-1,wkarr(izi+ldzi),ldzi)
c          call intpr('deltai',-1,event(i),1)
c          call intpr('tid',-1,tid,1)
c          call intpr('ldzit',-1,ldzit,1)
c          call dblepr('zit(,1)',-1,wkarr(izit),ldzit)
c          call dblepr('zit(,2)',-1,wkarr(izit+ldzit),ldzit)
c          call dblepr('gamma',-1,gamma,2)
c          call dblepr('ebi',-1,Ebi,qb)
c          call intpr('ldeb2',-1,ldeb2,1)
c          call dblepr('ebbi(,1)',-1,Ebbi(1,1),4)
c          call dblepr('ebbi(,4)',-1,Ebbi(1,4),4)
c          call intpr('ldeb4',-1,ldeb4,1)
c          call dblepr('Eb4(,1)',-1,Eb4(1,1),16)
c          call dblepr('Eb4(,16)',-1,Eb4(1,16),16)
c          call intpr('lder',-1,lder,1)
c          call dblepr('Er(,1)',-1,wkarr(ier),qe)
c          call dblepr('Er(,4)',-1,wkarr(ier+3*lder),qe)
c          call intpr('lder2c',-1,lder2c,1)
c          call dblepr('Er2c(,1)',-1,wkarr(ier2c),qea)
c          call dblepr('Er2c(,16)',-1,wkarr(ier2c+15*lder2c),qea)
c          call intpr('lder2',-1,lder2,1)
c          call dblepr('Er2(,1)',-1,wkarr(ier2),qea)
c          call dblepr('Er2(,4)',-1,wkarr(ier2+3*lder2),qea)
c          call intpr('lder3',-1,lder3,1)
c          call dblepr('Er3(,1)',-1,wkarr(ier3),lder3)
c          call dblepr('Er3(,16)',-1,wkarr(ier3+15*lder3),lder3)
c          call intpr('lder4',-1,lder4,1)
c          call dblepr('Er4(,1)',-1,wkarr(ier4),lder4)
c          call dblepr('Er4(,16)',-1,wkarr(ier4+15*lder4),lder4)
c          call intpr('lderbb',-1,lderbb,1)
c          call dblepr('erbb(,1)',-1,wkarr(ierbb),lderbb)
c          call dblepr('erbb(,4)',-1,wkarr(ierbb+3*lderbb),lderbb)
c          call intpr('ldeb2r2',-1,ldeb2r2,1)
c          call dblepr('eb2r2(,1)',-1,wkarr(ieb2r2),ldeb2r2)
c          call dblepr('eb2r2(,4)',-1,wkarr(ieb2r2+3*ldeb2r2),ldeb2r2)
c          call intpr('lderzb',-1,lderzb,1)
c          call dblepr('erzb(,1)',-1,wkarr(ierzb),lderzb)
c          call dblepr('erzb(,4)',-1,wkarr(ierzb+3*lderzb),lderzb)
c          call intpr('ldebzzb',-1,ldebzzb,1)
c          call dblepr('Ebzzb(,1)',-1,Ebzzb(1,1),qe)
c          call dblepr('Ebzzb(,2)',-1,Ebzzb(1,2),qe)
c          call intpr('ldebzbb',-1,ldebzbb,1)
c          call dblepr('Ebzbb(,1)',-1,Ebzbb(1,1),qe)
c          call dblepr('Ebzbb(,16)',-1,Ebzbb(1,16),qe)
c          call intpr('ldebzrr',-1,ldebzrr,1)
c          call dblepr('ebzrr(,1)',-1,wkarr(iebzrr),ldebzrr)
c          call dblepr('ebzrr(,4)',-1,wkarr(iebzrr+3*ldebzrr),ldebzrr)
c          call dblepr('expz',-1,wkarr(iexpz),ne1)
c          call intpr('ldexpzb',-1,ldexpzb,1)
c          call dblepr('expzb(,1)',-1,wkarr(iexpzb),ldexpzb)
c          call dblepr('expzb(,23)',-1,wkarr(iexpzb+22*ldexpzb),ldexpzb)
c          call intpr('ldz3btb',-1,ldz3btb,1)
c          call dblepr('exz3btb(,1)',-1,wkarr(iexz3btb),ldz3btb)
c          call dblepr('exz3btb(,23)',-1,wkarr(iexz3btb+22*ldz3btb),
c     $		  ldz3btb)
c          call dblepr('expzz',-1,wkarr(iexpzz),ne12u)
c          call intpr('ldexpzr',-1,ldexpzr,1)
c          call dblepr('expzr(,1)',-1,wkarr(iexpzr),ldexpzr)
c          call dblepr('expzr(,23)',-1,wkarr(iexpzr+22*ldexpzr),
c     $		  ldexpzr)
c          call intpr('ldexpzzb',-1,ldexpzzb,1)
c          call dblepr('expzzb(,1)',-1,wkarr(iexpzzb),ldexpzzb)
c          call dblepr('expzzb(,1596)',-1,wkarr(iexpzzb+1595*ldexpzzb),
c     $		  ldexpzzb)
c          call intpr('ldzrzb',-1,ldzrzb,1)
c          call dblepr('expzrzb(,1)',-1,wkarr(iexpzrzb),ldzrzb)
c          call dblepr('expzrzb(,23)',-1,wkarr(iexpzrzb+22*ldzrzb),
c     $		  ldzrzb)
c          call intpr('ldz4b2',-1,ldz4b2,1)
c          call dblepr('expz4b2(,1)',-1,wkarr(iexpz4b2),ldz4b2)
c          call dblepr('expz4b2(,1596)',-1,wkarr(iexpz4b2+1595*ldz4b2),
c     $		  ldz4b2)
c          call intpr('ldexpzbb',-1,ldexpzbb,1)
c          call dblepr('expzbb(,1)',-1,wkarr(iexpzbb),ldexpzbb)
c          call dblepr('expzbb(,23)',-1,wkarr(iexpzbb+22*ldexpzbb),
c     $		  ldexpzbb)
c          call intpr('ldz2b3',-1,ldz2b3,1)
c          call dblepr('expz2b3(,1)',-1,wkarr(iexpz2b3),ldz2b3)
c          call dblepr('expz2b3(,23)',-1,wkarr(iexpz2b3+22*ldz2b3),
c     $		  ldz2b3)
c          call intpr('ldexpzrr',-1,ldexpzrr,1)
c          call dblepr('expzrr(,1)',-1,wkarr(iexpzrr),ldexpzrr)
c          call dblepr('expzrr(,23)',-1,wkarr(iexpzrr+22*ldexpzrr),
c     $		  ldexpzrr)
c          call intpr('ldexpzrr',-1,ldexpzrr,1)
c          call dblepr('expzrr(,1)',-1,wkarr(iexpzrr),ldexpzrr)
c          call dblepr('expzrr(,23)',-1,wkarr(iexpzrr+22*ldexpzrr),
c     $		  ldexpzrr)
c          call intpr('ldz2br2',-1,ldz2br2,1)
c          call dblepr('exz2br2(,1)',-1,wkarr(iexz2br2),ldz2br2)
c          call dblepr('exz2br2(,23)',-1,wkarr(iexz2br2+22*ldz2br2),
c     $		  ldz2br2)
c          call intpr('ldz3b2',-1,ldz3b2,1)
c          call dblepr('expz3b2(,1)',-1,wkarr(iexpz3b2),ldz3b2)
c          call dblepr('expz3b2(,23)',-1,wkarr(iexpz3b2+22*ldz3b2),
c     $		  ldz3b2)
c          call dblepr('tmpvec1',-1,wkarr(itmpvec1),ne1)
c          call intpr('ldtm1',-1,ldtm1,1)
c          call dblepr('tmpmat1(,1)',-1,wkarr(itmpmat1),ldtm1)
c          call dblepr('tmpmat1(,4)',-1,wkarr(itmpmat1+3*ldtm1),ldtm1)
c          call intpr('ldtm2',-1,ldtm2,1)
c          call dblepr('tmpmat2(,1)',-1,wkarr(itmpmat2),ldtm2)
c          call dblepr('tmpmat2(,57)',-1,wkarr(itmpmat2+56*ldtm2),ldtm2)		  
c        return
        
        intvec(3) = ni
        intvec(4) = tid
        call evalwpb2(wkarr(ibrarr),ldbrarr,intvec,wkarr(izi),ldzi,
     $       event(i),wkarr(izit),ldzit,gamma,
     $       Ebi,Ebbi,ldeb2,Eb4,ldeb4,wkarr(ier),lder,wkarr(ier2c),
     $       lder2c,wkarr(ier2),lder2,wkarr(ier3),lder3,wkarr(ier4),
     $       lder4,wkarr(ierbb),lderbb,wkarr(ieb2r2),ldeb2r2,
     $       wkarr(ierzb),lderzb,Ebzzb,ldebzzb,Ebzbb,ldebzbb,
     $       wkarr(iebzrr),ldebzrr,wkarr(iexpz),wkarr(iexpzb),ldexpzb,
     $       wkarr(iexz3btb),ldz3btb,wkarr(iexpzz),wkarr(iexpzr),
     $       ldexpzr,wkarr(iexpzzb),ldexpzzb,wkarr(iexpzrzb),ldzrzb,
     $       wkarr(iexpz4b2),ldz4b2,wkarr(iexpzbb),ldexpzbb,
     $       wkarr(iexpz2b3),ldz2b3,wkarr(iexpzrr),ldexpzrr,
     $       wkarr(iexz2br2),ldz2br2,wkarr(iexpz3b2),ldz3b2,
     $       wkarr(itmpvec1),wkarr(itmpmat1),ldtm1,
     $       wkarr(itmpmat2),ldtm2)
c        total of 65 arguments (the maximum arguments allowed in R Fortran code)
        
c        call evalwpb2(wkarr(ibrarr),ldbrarr,qb,nsize,wkarr(izi),ldzi,ni,
c     $       event(i),tid,wkarr(izit),ldzit,gamma,
c     $       Ebi,Ebbi,ldeb2,Eb4,ldeb4,wkarr(ier),lder,wkarr(ier2c),
c     $       lder2c,wkarr(ier2),lder2,wkarr(ier3),lder3,wkarr(ier4),
c     $       lder4,wkarr(ierbb),lderbb,wkarr(ieb2r2),ldeb2r2,
c     $       wkarr(ierzb),lderzb,Ebzzb,ldebzzb,Ebzbb,ldebzbb,
c     $       wkarr(iebzrr),ldebzrr,wkarr(iexpz),wkarr(iexpzb),ldexpzb,
c     $       wkarr(iexz3btb),ldz3btb,wkarr(iexpzz),wkarr(iexpzr),
c     $       ldexpzr,wkarr(iexpzzb),ldexpzzb,wkarr(iexpzrzb),ldzrzb,
c     $       wkarr(iexpz4b2),ldz4b2,wkarr(iexpzbb),ldexpzbb,
c     $       wkarr(iexpz2b3),ldz2b3,wkarr(iexpzrr),ldexpzrr,
c     $       wkarr(iexz2br2),ldz2br2,wkarr(iexpz3b2),ldz3b2,
c     $       wkarr(itmpvec1),wkarr(itmpmat1),ldtm1,
c     $       wkarr(itmpmat2),ldtm2)

c        if (i.eq.2) then
c          call ifprintf(0,0,0)
c          call dblepr('Ebi',-1,Ebi,qb)
c          call dblepr('Ebbi(,1)',-1,Ebbi(1,1),qb)
c          call dblepr('Ebbi(,4)',-1,Ebbi(1,4),qb)
c          call dblepr('Eb4(,1)',-1,Eb4(1,1),16)	  
c          call dblepr('Eb4(,16)',-1,Eb4(1,16),16)	  
c          call dblepr('Er(,1)',-1,wkarr(ier),qe)	  
c          call dblepr('Er(,4)',-1,wkarr(ier+3*lder),qe)	  
c          call dblepr('Er2c(,1)',-1,wkarr(ier2c),4)
c          call dblepr('Er2c(,16)',-1,wkarr(ier2c+15*lder2c),4)
c          call dblepr('Er2(,1)',-1,wkarr(ier2),4)
c          call dblepr('Er2(,4)',-1,wkarr(ier2+3*lder2),4)
c          call dblepr('Er3(,1)',-1,wkarr(ier3),8)
c          call dblepr('Er3(,16)',-1,wkarr(ier3+15*lder3),8)
c          call dblepr('Er4(,1)',-1,wkarr(ier4),16)
c          call dblepr('Er4(,16)',-1,wkarr(ier4+15*lder4),16)
c          call dblepr('Erbb(,1)',-1,wkarr(ierbb),32)	  
c          call dblepr('Erbb(,4)',-1,wkarr(ierbb+3*lderbb),32)	  
c          call dblepr('Eb2r2(,1)',-1,wkarr(ieb2r2),64)	  
c          call dblepr('Eb2r2(,4)',-1,wkarr(ieb2r2+3*ldeb2r2),64)
c          call dblepr('Erzb(,1)',-1,wkarr(ierzb),4)
c          call dblepr('Erzb(,ni)',-1,wkarr(ierzb+(ni-1)*lderzb),4)
c          call dblepr('Ebzzb(,1)',-1,Ebzzb(1,1),2)
c          call dblepr('Ebzzb(,2)',-1,Ebzzb(1,2),2)
c          call dblepr('Ebzbb(,1)',-1,Ebzbb(1,1),2)
c          call dblepr('Ebzbb(,16)',-1,Ebzbb(1,16),2)
c          call dblepr('Ebzrr(,1)',-1,wkarr(iebzrr),8)
c          call dblepr('Ebzrr(,ni)',-1,wkarr(iebzrr+(ni-1)*ldebzrr),8)
c          call dblepr('Eexpz',-1,wkarr(iexpz),tid)
c          call dblepr('Eexpzb(,1)',-1,wkarr(iexpzb),4)
c          call dblepr('Eexpzb(,tid)',-1,wkarr(iexpzb+55*ldexpzb),4)
c          call dblepr('Eexz3btb(,1)',-1,wkarr(iexz3btb),4)
c          call dblepr('Eexz3btb(,tid)',-1,wkarr(iexz3btb+55*ldz3btb),4)
c          call dblepr('Eexpzz',-1,wkarr(iexpzz),tid*(tid+1)/2)
c          call dblepr('Eexpzr(,1)',-1,wkarr(iexpzr),8)
c          call dblepr('Eexpzr(,tid)',-1,wkarr(iexpzr+55*ldexpzr),8)
c          call dblepr('Eexpzzb(,1)',-1,wkarr(iexpzzb),4)
c          call dblepr('Eexpzzb(,tid2u)',-1,wkarr(iexpzzb+1595*ldexpzzb)
c     $       ,4)
c          call dblepr('Eexpzrzb(,1)',-1,wkarr(iexpzrzb),16)
c          call dblepr('Eexpzrzb(,tid)',-1,wkarr(iexpzrzb+55*ldzrzb),16)
c          call dblepr('Eexpz4b2(,1)',-1,wkarr(iexpz4b2),4)
c          call dblepr('Eexpz4b2(,tid2u)',-1,wkarr(iexpz4b2+1595*ldz4b2)
c     $       ,4)
c          call dblepr('Eexpzbb(,1)',-1,wkarr(iexpzbb),16)
c          call dblepr('Eexpzbb(,tid)',-1,wkarr(iexpzbb+55*ldexpzbb),16)
c          call dblepr('Eexpz2b3(,1)',-1,wkarr(iexpz2b3),32)
c          call dblepr('Eexpz2b3(,tid)',-1,wkarr(iexpz2b3+55*ldz2b3),32)
c          call dblepr('Eexpzrr(,1)',-1,wkarr(iexpzrr),16)
c          call dblepr('Eexpzrr(,tid)',-1,wkarr(iexpzrr+55*ldexpzrr),16)
c          call dblepr('Eexz2br2(,1)',-1,wkarr(iexz2br2),32)
c          call dblepr('Eexz2br2(,tid)',-1,wkarr(iexz2br2+55*ldz2br2)
c     $       ,32)
c          call dblepr('Eexpz3b2(,1)',-1,wkarr(iexpz3b2),4)
c          call dblepr('Eexpz3b2(,tid)',-1,wkarr(iexpz3b2+55*ldz3b2),4)
c	      return
c        end if
        
        
        call getghcpd(event(i),wkarr(ivbeta+i-1),h0t,tid,wkarr(ixit),
     $       ldvec,px,wkarr(izit),alphamat,qe,gamma,vi,pv,
     $       wkarr(izitz),wkarr(ixi),ni,sigbinv,sigeinv,
     $       Ebi,Ebbi,Eb4,wkarr(ier),wkarr(ier2c),
     $       wkarr(ier2),wkarr(ier3),wkarr(ier4),
     $       wkarr(ierbb),wkarr(ieb2r2),
     $       wkarr(ierzb),Ebzzb,Ebzbb,
     $       wkarr(iebzrr),wkarr(iexpz),wkarr(iexpzb),
     $       wkarr(iexz3btb),wkarr(iexpzz),wkarr(iexpzr),
     $       wkarr(iexpzzb),wkarr(iexpzrzb),
     $       wkarr(iexpz4b2),wkarr(iexpzbb),
     $       wkarr(iexpz2b3),wkarr(iexpzrr),
     $       wkarr(iexz2br2),wkarr(iexpz3b2),
     $       vecwb,krwb2,vecwe,krwe2,dqb,dqe,
     $       tqbqb,tqeqe,wkarr(iretGrad),ne1,
     $       wkarr(iretInfo),wkarr(iretCpd),wkarr(itmpvec1),
     $       wkarr(itmpvec2),wkarr(ixitalp))
c        total of 60 arguments (5 fewer than the maximum arguments allowed in R Fortran code)
        
c        call getghcpd(event(i),wkarr(ivbeta+i-1),h0t,tid,wkarr(ixit),
c     $       ldxit,px,wkarr(izit),ldzit,alphamat,ldam,qe,gamma,vi,pv,
c     $       wkarr(izitz),ldzitz,wkarr(ixi),ldxi,ni,sigbinv,ldsb2,
c     $       sigeinv,ldse2,
c     $       Ebi,Ebbi,ldeb2,Eb4,ldeb4,wkarr(ier),lder,wkarr(ier2c),
c     $       lder2c,wkarr(ier2),lder2,wkarr(ier3),lder3,wkarr(ier4),
c     $       lder4,wkarr(ierbb),lderbb,wkarr(ieb2r2),ldeb2r2,
c     $       wkarr(ierzb),lderzb,Ebzzb,ldebzzb,Ebzbb,ldebzbb,
c     $       wkarr(iebzrr),ldebzrr,wkarr(iexpz),wkarr(iexpzb),ldexpzb,
c     $       wkarr(iexz3btb),ldz3btb,wkarr(iexpzz),wkarr(iexpzr),
c     $       ldexpzr,wkarr(iexpzzb),ldexpzzb,wkarr(iexpzrzb),ldzrzb,
c     $       wkarr(iexpz4b2),ldz4b2,wkarr(iexpzbb),ldexpzbb,
c     $       wkarr(iexpz2b3),ldz2b3,wkarr(iexpzrr),ldexpzrr,
c     $       wkarr(iexz2br2),ldz2br2,wkarr(iexpz3b2),ldz3b2,
c     $       vecwb,krwb2,ldkb,vecwe,krwe2,ldke,dqb,lddqb,dqe,lddqe,
c     $       tqbqb,ldtqb,tqeqe,ldtqe,wkarr(iretGrad),ne1,
c     $       wkarr(iretInfo),ldif,wkarr(iretCpd),ldcp,wkarr(itmpvec1),
c     $       wkarr(itmpvec2),wkarr(ixitalp),ldxa)
        
c        if (i.eq.2) then		  
c          call ifprintf(0,0,0)
c          call dblepr('retGrad',-1,wkarr(iretGrad),ntotpar)		  
c          call dblepr('retInfo(,1)',-1,wkarr(iretInfo),ntotpar)
c          call dblepr('retInfo(,tid)',-1,wkarr(iretInfo+55*ldif),78)
c          call dblepr('retInfo(,57)',-1,wkarr(iretInfo+56*ldif),78)
c          call dblepr('retInfo(,62)',-1,wkarr(iretInfo+61*ldif),78)	  
c          call dblepr('retInfo(,63)',-1,wkarr(iretInfo+62*ldif),78)	  
c          call dblepr('retInfo(,64)',-1,wkarr(iretInfo+63*ldif),78)
c          call dblepr('retInfo(,65)',-1,wkarr(iretInfo+64*ldif),78)	  	  
c          call dblepr('retInfo(,66)',-1,wkarr(iretInfo+65*ldif),78)
c          call dblepr('retInfo(,75)',-1,wkarr(iretInfo+74*ldif),78)	  
c          call dblepr('retInfo(,76)',-1,wkarr(iretInfo+75*ldif),78)
c          call dblepr('retInfo(,77)',-1,wkarr(iretInfo+76*ldif),78)
c          call dblepr('retInfo(,78)',-1,wkarr(iretInfo+77*ldif),78)
c          call dblepr('retCpd(,1)',-1,wkarr(iretCpd),ntotpar)
c          call dblepr('retCpd(,tid)',-1,wkarr(iretCpd+55*ldcp),78)
c          call dblepr('retCpd(,57)',-1,wkarr(iretCpd+56*ldcp),78)
c          call dblepr('retCpd(,62)',-1,wkarr(iretCpd+61*ldcp),78)	  
c          call dblepr('retCpd(,63)',-1,wkarr(iretCpd+62*ldcp),78)	  
c          call dblepr('retCpd(,64)',-1,wkarr(iretCpd+63*ldcp),78)
c          call dblepr('retCpd(,65)',-1,wkarr(iretCpd+64*ldcp),78)	  	  
c          call dblepr('retCpd(,66)',-1,wkarr(iretCpd+65*ldcp),78)
c          call dblepr('retCpd(,75)',-1,wkarr(iretCpd+74*ldcp),78)	  
c          call dblepr('retCpd(,76)',-1,wkarr(iretCpd+75*ldcp),78)
c          call dblepr('retCpd(,77)',-1,wkarr(iretCpd+76*ldcp),78)
c          call dblepr('retCpd(,78)',-1,wkarr(iretCpd+77*ldcp),78)      
c          return
c        end if
        
        call sumgrhs(wkarr(iretGrad),ntotpar,wkarr(iretInfo),ldif,
     $       wkarr(iretCpd),ldcp,origGrad,origInfo,ldoi)
c        if (i.eq.2) then
c          call ifprintf(0,0,0)
c          call dblepr('origGrad',-1,origGrad,ntotpar)
c          call dblepr('origInfo(,1)',-1,origInfo(1,1),ntotpar)
c          call dblepr('origInfo(,78)',-1,origInfo(1,78),ntotpar)
c          return
c        end if
        
        offset1 = offset1 + tidp1
        offset2 = offset2 + tid
        offset3 = offset3 + ni
        call ifprintf(0,0,0)
c        if (mod(i,5).eq.0) call ifprintf(0,0,0)
 430  continue
      
      
c      call matzero(ldj2,ntotpar,ntotpar,wkarr(iJ2))
c      call matzero(lddj2,ntotpar**2,ntotpar,wkarr(idJ2))
      
      call getjdlh2(h0t,ne1,wkarr(ilh0t),wkarr(iJ2),ldj2,ntotpar,
     $     wkarr(idJ2),lddj2)
c      call dblepr('lh0t',-1,wkarr(ilh0t),ne1)
c      call dblepr('J2(,1)',-1,wkarr(iJ2),ne1)
c      call dblepr('J2(,ne1)',-1,wkarr(iJ2+(ne1-1)*ldj2),ne1)
c      call dblepr('dJ2(,1)',-1,wkarr(idJ2),ne1**2)
c      call dblepr('dJ2(,ne1)',-1,wkarr(idJ2+(ne1-1)*lddj2),ne1**2)
c      return
      
      call getjdlc2(sigb,ldsb,qb,lcholb,ldlb,Jb,ldjb,dJb,lddjb)
c      call dblepr('lcholb(,1)',-1,lcholb(1,1),qb)
c      call dblepr('lcholb(,4)',-1,lcholb(1,4),qb)
c      call dblepr('Jb(,1)',-1,Jb(1,1),qbu)
c      call dblepr('Jb(,10)',-1,Jb(1,10),qbu)
c      call dblepr('dJb(,1)',-1,dJb(1,1),qbu**2)
c      call dblepr('dJb(,10)',-1,dJb(1,10),qbu**2)
c      return
      
      call getjdlc2(sige,ldse,qe,lchole,ldle,Je,ldje,dJe,lddje)
c      call dblepr('lchole(,1)',-1,lchole(1,1),qe)
c      call dblepr('lchole(,2)',-1,lchole(1,2),qe)
c      call dblepr('Je(,1)',-1,Je(1,1),qeu)
c      call dblepr('Je(,3)',-1,Je(1,3),qeu)
c      call dblepr('dJe(,1)',-1,dJe(1,1),qeu**2)
c      call dblepr('dJe(,3)',-1,dJe(1,3),qeu**2)
c      return
      
      call getgdiff(h0t,ne1,betas,nfixpar,sigb,ldsb,qb,lcholb,
     $     ldlb,sige,ldse,lchole,ldle,Jb,ldjb,dJb,lddjb,Je,ldje,
     $     dJe,lddje,origParm,origGrad,origInfo,ldoi,wkarr(iP2),
     $     ldp2,natrParm,natrGrad,natrInfo,ldni,wkarr(iJ2),ldj2,
     $     wkarr(idJ2),lddj2,wkarr(iwkinfo),ldwi)
	  
c      call dblepr('origParm',-1,origParm,ntotpar)
c      call dblepr('origGrad',-1,origGrad,ntotpar)
c      call dblepr('origInfo(,1)',-1,origInfo(1,1),ntotpar)
c      call dblepr('origInfo(,78)',-1,origInfo(1,78),ntotpar)
c      call dblepr('natrParm',-1,natrParm,ntotpar)
c      call dblepr('natrGrad',-1,natrGrad,ntotpar)
c      call dblepr('natrInfo(,1)',-1,natrInfo(1,1),ntotpar)
c      call dblepr('natrInfo(,78)',-1,natrInfo(1,78),ntotpar)
      
      
      call ifprintf(0,0,0)        
      
      return
      end
