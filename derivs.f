	subroutine derivs(p,q,F,alpha,beta,lambda,cnu,df)
        implicit none

        integer p,q,alpha,beta,mini
        integer N
        integer nmax

        parameter(nmax=100)

        double precision F(0:nmax,0:nmax,3,3)
        double precision cnu,df,retshift,ds,nu,lambda
        double precision sum,sumret,sumcon,sumccr,sumrep
        double precision pi,taue, extdot,dssqinv,ds2inv
        double precision pp,pm,qp,qm,epsilon,eps2inv
        double precision fpq,fpp1q,fpm1q,fpqp1,fpqm1
        double precision fpp1qp1,fpm1qm1
        double precision feqpq,feqpp1q,feqpm1q,feqpqp1,feqpqm1
        double precision trpp,trqq
        double precision trppp1,trppm1,trqqp1,trqqm1
        double precision trppsr,trqqsr
        double precision trppp1sr,trppm1sr,trqqp1sr,trqqm1sr
        double precision flampp,flamqq
        double precision flamppp1,flamppm1,flamqqp1,flamqqm1
        double precision flambylpp,flambylqq
        double precision flambylppp1,flambylppm1
        double precision flambylqqp1,flambylqqm1
        double precision trmsr,trmp1sr,trmm1sr,dclfpq

        double precision TrF,Feq,Dclf,flambda,flambdabyl
	integer findMini
        external TrF,Feq,Dclf,flambda,flambdabyl, findMini

        common/retractionshift/RetShift
        common/segment/ds,dssqinv,ds2inv
        common/epsilon/epsilon,eps2inv
        common/taue/taue       
        common/Pi/pi 
        common/extension/extdot
        common/points/N

        nu=Cnu*lambda
        sum=0.0
        sumret = 0.0
        sumcon = 0.0
        sumccr = 0.0
        sumrep = 0.0

        fpq = f(p,q,alpha,beta)
        fpp1q = f(p+1,q,alpha,beta)
        fpm1q = f(p-1,q,alpha,beta)
        fpqp1 = f(p,q+1,alpha,beta)
        fpqm1 = f(p,q-1,alpha,beta)

        feqpq = feq(p,q,alpha,beta)
        feqpp1q = feq(p+1,q,alpha,beta)
        feqpm1q = feq(p-1,q,alpha,beta)
        feqpqp1 = feq(p,q+1,alpha,beta)
        feqpqm1 = feq(p,q-1,alpha,beta)

        trpp = TrF(p,p,F)
        trqq = TrF(q,q,F)
        trppp1 = TrF(p+1,p+1,F)
        trppm1 = TrF(p-1,p-1,F)
        trqqp1 = TrF(q+1,q+1,F)
        trqqm1 = Trf(q-1,q-1,F)
        
        trppsr = dsqrt(TrF(p,p,F))
        trqqsr = dsqrt(TrF(q,q,F))
        trppp1sr = dsqrt(TrF(p+1,p+1,F))
        trppm1sr = dsqrt(TrF(p-1,p-1,F))
        trqqp1sr = dsqrt(TrF(q+1,q+1,F))
        trqqm1sr = dsqrt(Trf(q-1,q-1,F))

        flampp = flambda(p,p,F)
        flamqq = flambda(q,q,F)
        flamppp1 = flambda(p+1,p+1,F)
        flamppm1 = flambda(p-1,p-1,F)
        flamqqp1 = flambda(q+1,q+1,F)
        flamqqm1 = flambda(q-1,q-1,F)

        flambylpp = flambdabyl(p,p,F)
        flambylqq = flambdabyl(q,q,F)
        flambylppp1 = flambdabyl(p+1,p+1,F)
        flambylppm1 = flambdabyl(p-1,p-1,F)
        flambylqqp1 = flambdabyl(q+1,q+1,F)
        flambylqqm1 = flambdabyl(q-1,q-1,F)


C     =======Retraction term======================
   
C       ###  (df/dp)*(1/lambda_p)*(dflambda/dp) term ###

        sumret=sumret+
     &     ds2inv*(fpp1q-fpm1q)
     &    *1.0/trppsr
     &    *ds2inv*(flamppp1-flamppm1) 
									    
        sumret=sumret+
     &     ds2inv*(fpqp1-fpqm1)
     &    *1.0/trqqsr
     &    *ds2inv*(flamqqp1-flamqqm1)  

									
C       ###  f*d(1/lambda_p)/dp*(dflambda/dp) term ###

        sumret=sumret+fpq
     &    *ds2inv*(1.0/trppp1sr-1.0/trppm1sr)
     &    *ds2inv*(flamppp1-flamppm1)

        sumret=sumret+fpq
     &   *ds2inv*(1.0/trqqp1sr-1.0/trqqm1sr)
     &   *ds2inv*(flamqqp1-flamqqm1)  


C       ###  f*(1/lambda_p)*(d2/dp^2)flambda term  ###

        sumret=sumret+fpq
     &   *1.0/trppsr
     &   *dssqinv*(flamppp1+flamppm1-2.0*flampp)     

        sumret=sumret+fpq
     &   *1.0/trqqsr
     &   *dssqinv*(flamqqp1+flamqqm1-2.0*flamqq)

        
        sumret=sumret/(Pi**2*taue)*RetShift
      

C     =======Convection term for extension======================

C       ###  Fxx  ###

        !if (alpha.EQ.1) then
	!   if (beta.EQ.1) then
        !  sumcon=sumcon+2.0*extdot*F(p,q,1,1)
	!endif
	!endif

C       ###  Fyy  ###

        !if (alpha.EQ.2) then
        ! if (beta.EQ.2) then
        !  sumcon=sumcon-extdot*F(p,q,2,2)
        ! endif
        !endif 
        
	 !==shear==
      !Fxx
      if (alpha.EQ.1) then
         if (beta.EQ.1) then
            sumcon=sumcon+2.0*extdot*F(p,q,1,2)
         endif
      endif
      

      !Fxy,yx term
      if (alpha.EQ.1) then
         if (beta.EQ.2) then
            sumcon=sumcon+extdot*F(p,q,2,2)
         endif
      endif


C     =======CCr term======================
 
C       ###   f(s,s') terms   ###

C       ###   (Flambda/(lambdap2))*(d^2/ds^2)F term   ###

        sumccr=sumccr+
     &     dssqinv*(fpp1q+fpm1q-2.0*fpq)
     &     *flampp/Trpp
        
        sumccr=sumccr+
     &     dssqinv*(fpqp1+fpqm1-2.0*fpq)
     &     *flamqq/Trqq


C       ###   flambda*(d(1/lambdap2)/ds)*d/ds(F)   ###

        sumccr=sumccr+flampp
     &     *ds2inv*(1.0/Trppp1-1.0/Trppm1) 
     &     *ds2inv*(fpp1q-fpm1q)
    
        sumccr=sumccr+flamqq
     &     *ds2inv*(1.0/Trqqp1-1.0/Trqqm1) 
     &     *ds2inv*(fpqp1-fpqm1)


C       ###   d(flambda)/ds*(1/lambdap2)*d/ds(F)   ###

        sumccr=sumccr+
     &      ds2inv*(flamppp1-flamppm1)
     &     *1.0/Trpp 
     &     *ds2inv*(fpp1q-fpm1q)

        sumccr=sumccr+
     &      ds2inv*(flamqqp1-flamqqm1)
     &     *1.0/Trqq
     &     *ds2inv*(fpqp1-fpqm1)


C       ###   Extra terms  ###

C       ###   f/lambdap*d2(flambdabyl)/dp2

        sumccr=sumccr+fpq/trppsr
     &     *dssqinv*(flambylppp1+flambylppm1-2.0*flambylpp)

        sumccr=sumccr+fpq/trqqsr
     &     *dssqinv*(flambylqqp1+flambylqqm1-2.0*flambylqq)
                

C       ###  (df/ds)*(1/lambdap)*(dflambdabyl/dp)

        sumccr=sumccr+1.0/trppsr
     &     *ds2inv*(fpp1q-fpm1q)
     &     *ds2inv*(flambylppp1-flambylppm1)

        sumccr=sumccr+1.0/trqqsr
     &     *ds2inv*(fpqp1-fpqm1)
     &     *ds2inv*(flambylqqp1-flambylqqm1)


C       ###  f*(d(1/lambdap)/dp)*(dflambdabyl/dp)

        sumccr=sumccr+fpq
     &   *ds2inv*(1.0/trppp1sr-1.0/trppm1sr)
     &   *ds2inv*(flambylppp1-flambylppm1)

        sumccr=sumccr+fpq
     &   *ds2inv*(1.0/trqqp1sr-1.0/trqqm1sr)
     &   *ds2inv*(flambylqqp1-flambylqqm1)


C       ###   feq terms   ###

C       ###   1/sqrt(TrFpp)(d^2/ds^2)Feq term   ###

        sumccr=sumccr-
     &     dssqinv*(Feqpp1q+Feqpm1q-2.0*Feqpq)
     &     *1.0/trppsr
     
        sumccr=sumccr-
     &     dssqinv*(Feqpqp1+Feqpqm1-2.0*Feqpq)
     &     *1.0/trqqsr


C       ###   d/ds(1/sqrt(TrFpp)*d/ds(Feq)   ###

        sumccr=sumccr-
     &     ds2inv*(1.0/trppp1sr-1.0/trppm1sr)
     &     *ds2inv*(Feqpp1q-Feqpm1q)
     
        sumccr=sumccr-
     &     ds2inv*(1.0/trqqp1sr-1.0/trqqm1sr) 
     &     *ds2inv*(Feqpqp1-Feqpqm1)


        sumccr=sumccr*1.5*nu


C     =======Reptation + CLF term======================
 
C       find point closest to chain end
        mini=findMini(p,q)

        pp=ds*p+epsilon
        pm=ds*p-epsilon 
        qp=ds*q+epsilon
        qm=ds*q-epsilon
 
        fpp1qp1 = f(p+1,q+1,alpha,beta)
        fpm1qm1 = f(p-1,q-1,alpha,beta)

        trmsr = dsqrt(trF(mini,mini,F))
        trmp1sr = dsqrt(TrF(mini+1,mini+1,F))
        trmm1sr = dsqrt(TrF(mini-1,mini-1,F))

        dclfpq = Dclf(ds*p,ds*q)

C       ###   (d/dp+d/dq)Dclf*(1/sqrt(Trfmin))*(d/dp+d/dq)F   ###

        sumrep=sumrep+
     &     eps2inv*(Dclf(pp,qp)-Dclf(pm,qm))
     &     *1.0/trmsr
     &     *ds2inv*(fpp1qp1-fpm1qm1)


C       ###  Dclf*(d/dp+d/dq)(1/sqrt(Trfmin))*(d/dp+d/dq)F

        sumrep=sumrep+
     &     Dclfpq
     &     *ds2inv*(1.0/trmp1sr-1.0/trmm1sr) 
     &     *ds2inv*(fpp1qp1-fpm1qm1)


C       ###  Dclf*(1/sqrt(Trfmin))*(d/dq+d/dq)**2F  ###

        sumrep=sumrep+
     &     Dclfpq/trmsr
     &     *dssqinv*(fpp1qp1+fpm1qm1-2.0*fpq)

        sumrep=sumrep*1.0/(3.0*Pi**2*taue)*(1.0/trmsr)


        df = sumret + sumcon + sumccr + sumrep

        return
        end


        Double Precision function TrF(p,q,F)
        implicit none

        integer p,q
        integer nmax

        parameter(nmax=100)

        double precision  F(0:nmax,0:nmax,3,3)

        TrF=F(p,q,1,1)+2.0*F(p,q,2,2)

        end function

        
        Double precision function flambda(p,q,F)
        implicit none
  
        integer p,q
        integer nmax
        double precision term,lambdam2,const,trace

        double precision TrF
        external TrF

        common/lambdamaxsq/lambdam2
        common/constant/const
        
        parameter(nmax=100)

        double precision F(0:nmax,0:nmax,3,3)
 
        trace = TrF(p,q,F)
        if (trace .gt. lambdam2) then
         write(*,*) "Lambda exceeding lambda_max !!!!!!"
         write(*,*) lambdam2,trace
         stop
        endif
        term=(lambdam2-trace/3.0)/(lambdam2-trace)
        term=dsqrt(trace)*term
        flambda=term*const
        end function

        
        Double precision function flambdabyl(p,q,F)
        implicit none
  
        integer p,q
        integer nmax
        double precision term,lambdam2,const,trace

        double precision TrF
        external TrF

        common/lambdamaxsq/lambdam2
        common/constant/const
        
        parameter(nmax=100)

        double precision F(0:nmax,0:nmax,3,3)
 




        trace = TrF(p,q,F)
        term=(lambdam2-trace/3.0)/(lambdam2-trace)
        flambdabyl=term*const
        end function        


        Double Precision function Feq(p,q,alpha,beta)
        implicit none

        integer p,q,alpha,beta
        double precision dum,ds,preal,qreal,dssqinv,ds2inv

        common/segment/ds,dssqinv,ds2inv

        dum=0.0
        preal = 1.0*p
        qreal = 1.0*q
        if(dabs(preal-qreal).lt.ds2inv) then
         if(alpha.EQ.beta) then
            dum=1.0/3.0
         endif
        endif
        Feq=dum
    
        end function

       
        Double Precision Function Dclf(p,q)
        implicit none

        Double Precision  alphad,s,p,q,zreal
        integer z
        
        common/entanglements/Z

        zreal = 1.0*z
        alphad=1.15
        s=Min(1.0*p,1.0*q,1.0*(Z-p),1.0*(Z-q))
 
        if (s.lt.alphaD*dsqrt(zreal)) then
         Dclf=(alphaD/s)**2
        else
         if(s.gt.Zreal-alphaD*dsqrt(zreal)) then
          Dclf=(alphaD/(s-Zreal))**2
         else
          Dclf=1.0/Zreal
         endif
        endif

        end function


	Integer function findMini(p,q)

	implicit none
	integer :: p,q, answ,N

	common/points/N
	
	
	
	
	
	if(p==q) then
	   answ=p
	endif
	
	if( (1.0*p< 1.0*N/2.0) .AND. (1.0*q<1.0*N/2.0)) then
	   answ = min( p,q)
	endif
      
	
	
	
	if( (1.0*p > 1.0*N/2.0) .AND. (1.0*q > 1.0*N/2.0)) then
	   answ = max(p,q)
	endif
	
	
	
	if( (1.0*p < 1.0*N/2.0) .AND. (1.0*q > 1.0*N/2.0)) then
	   
	   if( p< N-q) then
	      answ = p
	   else
	      answ = q
	   endif
	   
	endif
	
	
	if( (1.0*p > 1.0*N/2.0) .AND. (1.0*q < 1.0*N/2.0)) then
	   
	   if( q< N-p) then
	      answ = q
	   else
	      answ = p
	   endif
	   
	endif
	
	
	
	
	if(p==q) then
	   answ=p
	endif
	
	
	if(p==N/2) then
	   answ = q
	endif
	
	if(q==N/2) then
	   answ = p
	endif
	
	!print*,"mini",p,q,N, N/2, answ
	
	findMini = answ
	
	
	end function findMini
	
	
