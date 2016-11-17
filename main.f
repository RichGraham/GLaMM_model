	Program main
        implicit none

        integer Z,N
        integer nstep, nsave, nprint
        integer p,q,alpha,beta,k
        integer nmax

        parameter(nmax=100)

        double precision F(0:nmax,0:nmax,3,3)
        double precision k11(nmax,nmax),k22(nmax,nmax),k12(nmax,nmax)
        double precision lambdahistory, lambda, ax, dtrfpp
        double precision zeff, zreal, betarcr, retshift
        double precision start, finish, h, t, sum, fractaue
        double precision taue, taur, taud, ds, ge, cnu, df
        double precision pi, extdot, epsilon,dssqinv,ds2inv,eps2inv
        double precision N1, shearStress
        double precision lambdam,lambdam2,const
        double precision trace,term
	double precision angle, length, Rx, Ry, a,b,c

	character(2) :: ZString
	character(9) :: rateString


        double precision Feq,rethistory,TrF,flambdabyl
	integer findMini
        external Feq, rethistory, TrF,flambdabyl, findMini

        common/retractionshift/RetShift
        common/segment/ds,dssqinv,ds2inv
        common/taue/taue       
        common/Pi/pi 
        common/epsilon/epsilon,eps2inv
        common/extension/extdot
        common/entanglements/Z
        common/points/N
        common/lambdamaxsq/lambdam2
        common/constant/const

C       ### N has to be an odd multiple of Z ###
C       ### N = (2*m+1)Z ###

        open(unit=1,file='input',status='old')
        read(1,*) taue
        read(1,*) ge
        read(1,*) Z
        read(1,*) N
        read(1,*) cnu
        read(1,*) extdot
	read(1,*) rateString
        read(1,*) lambdam
        read(1,*) start
        read(1,*) finish
        read(1,*) fractaue
        read(1,*) nsave
        close(unit=1)


	ZString=char(Z/10+48)//Char(Z-10*(Z/10)+48)


        pi = 3.14159265359

C       ### Rs ###
        retshift = 2.0






        ds=(1.0*Z)/(1.0*N)
        dssqinv = 1.0/(ds**2)
        ds2inv = 0.5/ds
        zreal = 1.0*z
        BetaRCR = (2.13-8.91/(dsqrt(Zreal))+12.29/Zreal)
     &           *(1.0+0.46*dlog10(Cnu))

        taur = z**2*taue
        taud = 3*z*taur
        h = fractaue*taue
        nprint = nsave*10
        epsilon = 0.001
        eps2inv = 0.5/epsilon
        lambdam2 = lambdam*lambdam
        const = (lambdam2 - 1.0)/(lambdam2 - 1.0/3.0)

        write(*,*) " "
        write(*,*) " ## All parameters in real units ## "
        write(*,*) " "
        write(*,*) "Taue          =   ",taue
        write(*,*) "Taur          =   ",taur
        write(*,*) "Taud          =   ",taud
        write(*,*) "Extensiondot  =   ",extdot
        write(*,*) "Lambda-max    =   ",lambdam
        write(*,*) "Ge            =   ",ge 
        write(*,*) "Cnu           =   ",cnu
        write(*,*) "Entanglements =   ",Z
        write(*,*) "No. of points =   ",N
        write(*,*) "Timestep      =   ",h
        write(*,*) "Final time    =   ",finish
        write(*,*) "Total steps   =   ",dnint((finish-start)/h)
        write(*,*) "#########################################"
        write(*,*) " "
        write(*,*) "Computing stress every  ",nsave,"  steps"
        write(*,*) " "


        nstep = 0
        t = start

        do p = 1,nmax
         do q = 1,nmax
          k11(p,q) = 0.0
          k22(p,q) = 0.0
	  k12(p,q) = 0.0
         enddo
        enddo

        do p = 0,nmax
         do q = 0,nmax
          do alpha = 1,3
           do beta = 1,3
            F(p,q,alpha,beta) = 0.0
           enddo
          enddo
         enddo
        enddo

        do p = 0,N
         do q = 0,N
          do alpha = 1,3
           F(p,q,alpha,alpha) = Feq(p,q,alpha,alpha)
          enddo
         enddo
        enddo
         

C	open(unit=1,file='StdyFpq50/StdyRs23-e-5.dat',status='old')
C	do p = 0,N
C
C	   read(1,*) k, a,b,c
C	   F(p,p, 1,2) =a
C	   F(p,p, 1,1) =b
C	   F(p,p, 2,2) =c
C	   F(p,p, 3,3) =c
C
C	   print*,p,a,b,c
C        enddo

C        close(unit=1)
C	sum=0.0
C         do p=1,N-1
C	    sum=sum+ds*flambdabyl(p,p,F)*F(p,p,1,2)
C         enddo
C        
C	 shearStress=3.0/Z*(4.0/5.0)*Ge*sum
C	 print*,shearStress/3e-5
C
C	stop





        open(unit=1,file='Trans'//ZString//
     &  '/Rs2'//rateString//'.dat',status='unknown')
C        open(unit=2,file='Trans'//ZString//
C     &  '/lam'//rateString//'.dat',status='unknown')
C        open(unit=3,file='Trans'//ZString//
C     & '/taue'//rateString//'.dat',status='unknown')
        open(unit=4,file=
     &   'Fpq'//ZString//'/Length'//rateString//'.dat',status='unknown')


	write(*,*) "Full theory calculation begins.... "
        write(*,*) "Step #      ","Time       ","Z_effective"




C       ### Main loop begins ###
        do while (t.le.finish)

C	   if(t.ge.100) then
C	      extdot =0.0
C	   endif

C       ### Zeff (effective number of entanglements)  ###

C       include two end lengths and then others
        sum=(ds/2.0)*2.0          
        do p=1,N-1
         sum=sum+ds*dsqrt(TrF(p,p,F)) 
        enddo
        Zeff=sum

C       ###  Reptation CCR ###

        ax=(1.0)/(3.0*Z*Z*taue*BetaRCR*Zeff)

C       ### Retraction CCR (lambda) ###
    
        sum=0.0
        do p=1,N-1 
         dtrFpp=RetHistory(p,p,F,1,1,Zeff)
     &         +2.0*RetHistory(p,p,F,2,2,Zeff)
         sum=sum+ds*0.5*dtrFpp/(Zeff*(dsqrt(TrF(p,p,F))))
        enddo
        LambdaHistory=-1.0*sum
							
C       ### add both contributions ###
        
        lambda = lambdahistory+ax


C       ### compute all components ###
     
        do p=1,N-1
         do q=1,N-1
          call derivs(p,q,F,1,1,lambda,cnu,df)
          k11(p,q) = df
          call derivs(p,q,F,2,2,lambda,cnu,df)
          k22(p,q) = df
          call derivs(p,q,F,1,2,lambda,cnu,df)
          k12(p,q) = df
	  !write(*,*) dF, k12(p,q)
         enddo
        enddo		

C       ### Euler advancement in time ###

        do p=1,N-1
         do q=1,N-1
          f(p,q,1,1) = f(p,q,1,1) + h*k11(p,q)
          f(p,q,2,2) = f(p,q,2,2) + h*k22(p,q)
	  f(p,q,1,2) = f(p,q,1,2) + h*k12(p,q)
         enddo
        enddo

        t=t+h
        nstep=nstep+1

C       ### compute stress every few time steps ###

        if (mod(nstep,nsave).eq.0) then

C       ### N1 ###
         sum=0.0
         do p=1,N-1
          sum=sum+ds*flambdabyl(p,p,F)*(F(p,p,1,1)-F(p,p,2,2))
         enddo
         N1=3.0/Z*(4.0/5.0)*Ge*sum

	 sum=0.0
         do p=1,N-1
	    sum=sum+ds*flambdabyl(p,p,F)*F(p,p,1,2)
         enddo
         shearStress=3.0/Z*(4.0/5.0)*Ge*sum

         trace = Trf(25,25,F)
C         write(1,*) t,Zeff,shearStress/extdot, n1
 20	     format(i4,3f12.8)
 55	     format(4e20.12)
         write(1,55) t,shearStress/extdot,Zeff,n1





	 p=5
	 print*,p,F(p,p,1,2), F(p,p,1,1), 
     &   F(p,p,2,2), F(p,p,1,1)+2.0*F(p,p,2,2)



C         do p=0,N
C          trace = Trf(p,p,F)
C          write(2,*)t,p,dsqrt(trace),flambdabyl(p,p,F)
C          term= const*(lambdam2*lambdam2 + trace*trace/3.0)
C          term = term/(lambdam2 - trace)**2
C          write(3,*) t,p,taue/term
C         enddo

C       ###save "steady" data#####
C	  open(unit=7,file='StdyFpq'//ZString//'/StdyRs2'
C     &   //rateString//'.dat',status='unknown')

C	   do p=0,N
C	      write(7,20) (p+1),
C     &	      F(p,p,1,2), F(p,p,1,1), 
C     &   F(p,p,2,2)
C	   enddo
C	  close(unit=7)


        endif
 
C       ### Display on screen ###
        if (mod(nstep,nprint).eq.0) then
C        if (mod(nstep,1).eq.0) then
         write(*,*) nstep,"   ",t,"   ",zeff,"   ",
     &       dsqrt(Trf(25,25,f)),"   ",shearStress,flambdabyl(25,25,F)
        endif

        enddo
 
C       ###  Loop ends ###

        close(unit=1)
C        close(unit=2)
	close(unit=4)

        stop
        end


        Double Precision Function RetHistory(p,q,F,alpha,beta,Zeff)
        implicit none

        Integer p,q,i,j,alpha,beta,N,mini 
        integer nmax

        parameter(nmax=100)

        Double Precision F(0:nmax,0:nmax,3,3)
        double precision sumret,sumrep,zeff, ds,dssqinv,ds2inv
        double precision retshift, taue, pi
        double precision pp,pm,qp,qm,epsilon,eps2inv
        double precision fpq,fpp1q,fpm1q,fpqp1,fpqm1
        double precision fpp1qp1,fpm1qm1
        double precision trppsr,trqqsr
        double precision trppp1sr,trppm1sr,trqqp1sr,trqqm1sr
        double precision flampp,flamqq
        double precision flamppp1,flamppm1,flamqqp1,flamqqm1
        double precision trmsr,trmp1sr,trmm1sr,dclfpq
        
        double precision TrF,Dclf,flambda
	integer findMini
        external TrF,Dclf,flambda, findMini

        common/retractionshift/RetShift
        common/segment/ds,dssqinv,ds2inv
        common/taue/taue       
        common/Pi/pi 
        common/epsilon/epsilon,eps2inv
        common/points/N

        sumret=0.0
        sumrep=0.0

        fpq = f(p,q,alpha,beta)
        fpp1q = f(p+1,q,alpha,beta)
        fpm1q = f(p-1,q,alpha,beta)
        fpqp1 = f(p,q+1,alpha,beta)
        fpqm1 = f(p,q-1,alpha,beta)
        
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

        Rethistory = sumret+sumrep

	end function 
	
