	Program k_E
	implicit none
    integer :: date_time(8)
    character*10 :: b(3),t
    double precision :: redmass,Energy_MECP,EMECP,J_cm,EPS,Energy_min
    double precision :: QSIMP,Psh,P_sh
    double precision,dimension(1000)::FREQ,FREQ_MECP,ENERGY
    double precision,dimension(500) :: GradS,GradT
 	DOUBLE PRECISION RHO_SADDLE_POINT,RHO_REACTANT,rate,NMECP
	DOUBLE PRECISION E_MIN  
	EXTERNAL RHO_SADDLE_POINT
    integer i,NumGrad,Ncalc,Natoms,Nfreq,Nfreq_MECP
    double precision :: Planck,NA,H12,deltaF,Bohr,coef,sqr
    
	 	Planck = 6.6260693D-34  ! Planck's constant (J * s)
		NA = 6.0221415E23       ! Avogadro constant (mol-1)
      	BOHR = 0.5291772108E-10 ! Bohr radius (m)
		J_cm = 1.98644561E-23 	! J to cm-1 convertion factor

8   format(30x,'The output file is printed at ',I2,':',I2,1x,A2,' on ',I2,'/',I2,'/',I4,' in Taiwan')
!9   format(5x,F8.2,20x,E8.5,10x,E8.5,10x,E8.5,10x,E8.5)
10  format(40x,'THE RESULTS OF k(E) ARE LISTED BELOW')
11  format(5x,'E(Internal),kcal/mol',15x,'k(E),s-1',15x,'P(sh)',15x,'NMECP',15x,'RHO(Reactant)')
12  format(5x,'-------------------------------------------------------------------------------------------------------------')

    open(5,file='input.dat')
	open(6,file='output.log')
	open(7,file='energy.out')
	open(8,file='Ncr(E).out')
!	open(9,file='Psh.out')
    
	read(5,*)
    read(5,*) Natoms	!Number of atoms
    read(5,*) NumGrad	!Read the number of gradients    
    read(5,*) H12		!Electronic matrix element (cm-1) for coupling between two surfaces
    read(5,*) redmass	!Read the reduced mass of the system as it moves along the hopping coordinate
	READ(5,*) EPS 		! Convergence criterion for numerical integration
	READ(5,*)    
 	READ(5,*) E_Min  	!Relative energy of local minima (kcal/mol)
    READ(5,*)    
    read(5,*) EMECP		!Relative energy of MECP (kcal/mol)
    read(5,*)
    read(5,*) Ncalc		!Read number of energies to calculate
    read(5,*)
    
    Do i=1,Ncalc
      	read(5,*) Energy(i)	!Read internal energy (kcal/mol) above initial reactant
    End do
    read(5,*)
    Do i=1,NumGrad
      	read(5,*) GradS(i)	!Read singlet gradients from input file
    End do
    read(5,*)
    Do i=1,NumGrad
      	read(5,*) GradT(i)	!Read triplet gradients from input file
    End do
      NFREQ = 3*Natoms - 6
	  NFREQ_MECP = 3*Natoms - 7
	read(5,*)      
    DO i=1,NFREQ
      READ(5,*) FREQ(i) 
	End do
      READ(5,*) 
	DO i=1,NFREQ_MECP
      READ(5,*) FREQ_MECP(i) 
	End do
      ! =============== Calculation block ============    
	H12 = H12*11.9626/NA	!convert from cm-1 to J; 1cm-1 = 11.9626 J/mol; H12 = 63*11.9626/6.0221415E+23=1.25145E-21(J)
    redmass = redmass*1.66054E-27 !convert from a.u to kg; redmass =13.5411*1.66054E-27=2.24855E-26 (kg)

	! --- Calculate norm of gradients on two surfaces
    	deltaF = 0
    do i=1,NumGrad
      	deltaF = deltaF + (GradS(i)-GradT(i))*(GradS(i)-GradT(i))	!GradE = dE/dx, dE/dy, dE/dz (J/m)
    end do
    	coef = 62750.94717E-2*1000*4.184/NA/Bohr !convert from a.u to J/m; coefficient =62750.94717E-2*1000*4.184/6.0221415E+23/0.5291772108E-10=8.23872E-08
        sqr = sqrt(deltaF) !1 au = 627.5094717 kcal/mol.
    	deltaF = sqr*coef	!=0.1211*8.23872E-08=9.97709E-09 (J/m) 
	
	write(6,10)
    write(6,*)
    write(6,11)
	write(6,12)
    
    DO i=1,Ncalc
    	Energy_MECP = abs(Energy(i)-EMECP)	!Available engergy above MECP   
   		Energy_Min  = abs(Energy(i)-E_Min)  !Available Energy above local minima 
       write(7,*)
	   write(7,'(A18,F6.2,A10)') 'Internal energy = ',Energy(i),' kcal/mol' 
	   write(7,'(A29,F6.2,A10)') 'Available energy above MECP = ',Energy_MECP,' kcal/mol' 
	   write(7,'(A31,F6.2,A10)') 'Available energy above minima = ',Energy_Min,' kcal/mol'       

       write(8,*)
	   write(8,'(A58)') 'Calculating the integrated density of states Ncr(E) ...'
	   write(7,*)
       
         NMECP = qsimp(Energy_MECP,H12,RedMass,DeltaF,FREQ_MECP,NFREQ_MECP,EPS)
         !NMECP:= Ncr(E), the effective integrated density of states in the crossing seam between two surfaces.
	     RHO_REACTANT = RHO_SADDLE_POINT(Energy_MIN,FREQ,NFREQ)
         !RHO_Reactant:= Rho(E),the density of rovibrational states of the reactant.
         rate = NMECP/(RHO_REACTANT*(Planck/J_CM))
         !rate:= k(E) = Ncr(E)/[h*Rho(E)]
         P_sh = Psh(Energy_MECP,H12,Redmass,deltaF)
         
	   write(7,*)
	   write(7,'(A12,E10.5)') 'NMECP(E) = ', NMECP
	   write(7,'(A12,E10.5)') 'RHO_R    = ', RHO_REACTANT  
       write(7,'(A12,E10.5,A4)') 'k(E)     = ', rate, ' s-1' 
	   write(7,*) 
	   write(6,'(F15.2,E35.5,E22.5,E20.5,E20.5)') ENERGY(i),rate,P_sh,NMECP,RHO_REACTANT 	 
    END DO
	   write(6,*)
       write(6,12)
       call date_and_time(b(1),b(2),b(3),date_time)
       if (date_time(5)<12) then
         t = 'AM'
       else
         t = 'PM'
       end if
       write(6,'(50x,A21)')'Calculation finished!'
       write(6,8)date_time(5),date_time(6),t,date_time(2),date_time(3),date_time(1)
       write(7,*) 'Calculation finished! Open output file to get results.'
    
    End program k_E
!*******************************************************************************************    
      DOUBLE PRECISION FUNCTION qsimp(ENERGY,H12,RedMass,DeltaF,FREQ,NFREQ,EPS)
      IMPLICIT NONE
      INTEGER JMAX, i, j,it, NFREQ, maxfreq
      DOUBLE PRECISION ENERGY,RedMass
	  DOUBLE PRECISION PSH,RHO_SADDLE_POINT,EPS,b,a,a1,b1,c,d
      DOUBLE PRECISION os,ost,st,s1
	  DOUBLE PRECISION del,sum,tnm,x, y0,y1
      PARAMETER (maxfreq=1000) 
      DOUBLE PRECISION FREQ(maxfreq)
!	  EXTERNAL PSH, RHO_SADDLE_POINT
      PARAMETER (JMAX=200)
      double precision :: Cal_cm,s,H12,DeltaF
      	Cal_cm = 349.755!1153   ! kcal to cm-1 convertion factor
		a = 0.0    !  
		b = 0.0
		st=-1.e30
      	os= -1.e30
        ost= 0.0 !I'm not sure this value because origin version didn't mention
      DO i=1,JMAX
         ! ===== Trapezoid rule block ===
	   if (i.eq.1) then
       		c = 0.5*(ENERGY-a)*Psh(a,H12,RedMass,DeltaF)
            d = RHO_SADDLE_POINT((ENERGY-a),FREQ,NFREQ)
       		a1 = c*d
            b1 = Psh(ENERGY,H12,RedMass,DeltaF)*RHO_SADDLE_POINT(b,FREQ,NFREQ)
            st = a1+b1
            st=st*Cal_cm
	   else
            it=2**(i-2)! i = 20 --> it = 1048576
            tnm=it
            del=(ENERGY-a)/tnm !Energy=35.9 --> del = 35.9/1048576=3.42369E-05
            x=a+0.5*del!=1.71185E-05
            sum=0.
            DO j=1,it
      		 	y0 = Psh(x,H12,RedMass,DeltaF) 
				!if (y0==1.0) then
				!   write(9,'(5x,A6,F5.2)')'Psh = ',y0
				!end if
	         	y1 = RHO_SADDLE_POINT((ENERGY-x),FREQ,NFREQ)
			 	sum= sum + y0*y1*Cal_cm
			 	x  = x + del
            ENDDO
            	st=0.5*(st+(ENERGY-a)*sum/tnm)
       endif
         ! ==== End of Trapezoid block ====

         	s1 = 4*st-ost !actually there is no value of ost
	   		s=s1/3.
         	qsimp = s
	   write(8,13)'Iter. ', i, '   NMECP(E)= ', s,'   Grad = ', abs(s-os)
       if (abs(s-os).lt.EPS) return 
       if (s==0. .and. os==0. .and. i>6) return
       !if (i>=18) return
         
         os=s
         ost=st
      ENDDO
  13  format(A8,I3,A14,E10.5,A10,E10.3)      
!      write(*,*) 'too many steps in QSIMP'
!      pause
	  END FUNCTION QSIMP

!*******************************************************************************************  
    double precision function Psh(Eh,H12,redmass,deltaF)
    double precision :: h,NA,Pi
    double precision :: P,Eh,redmass,deltaE
    double precision :: H12,deltaF,A,B,C    
    	h = 6.6260693D-34 !Planck's constant (J.s)
   	 	NA = 6.0221415E23 !Avogadro constant (mol-1)
    	Pi = 4*atan(1.)
    if (Eh<=0) then	!when Eh <= Ec ---> Psh(Eh-Ec) = 0, Harvey pp.
      Psh = 0.0
	else
    deltaE = Eh*1000*4.184/NA	!Convert from kcal/mol to J; deltaE=(E(i)-EMECP)*10^3*4.184/6.0221415E23=(0.0+25.9)*10^3*4.184/6.0221415E23=1.79945E-19(J)
    A = -2*Pi*H12*H12 !=-2*3.14*1.25145E-21^2 = -9.83528E-42 (J^2)
    B = sqrt(redmass/(2*deltaE))!=sqrt(2.24855E-26(kg)/(2*1.79945E-19 J))=2.49958E-4 (s/m);1 J = 1 kg m2/s2.
    C = h*deltaF  !=6.6260693D-34(J.s)*9.97709E-09 (J/m)=6.61089E-42 (J^2.s/m) 
    P = exp(A*B/C) !=exp(-9.83528E-42*2.49958E-4/6.61089E-42)=0.999628197 (no units)
    Psh = (1-P)*(1+P)	!The probability of hopping, Psh = (1-0.999628197)*(1+0.999628197)=0.000743468
    end if
    end function Psh  
!*******************************************************************************************
      ! === Calculate dencity of states using Saddle point method	===
      ! 
	  DOUBLE PRECISION FUNCTION RHO_SADDLE_POINT(ENERGY, FREQ, NFREQ)
      IMPLICIT NONE
	  EXTERNAL zriddr, expsafe
	  INTEGER maxfreq 
      PARAMETER (maxfreq=1000)  
      DOUBLE PRECISION eng,ENERGY,FREQ(maxfreq)
	  INTEGER i,NFREQ
	  DOUBLE PRECISION sum,prod,temp1,pi,expsafe
	  DOUBLE PRECISION xacc,zriddr 
      double precision :: x1,x2,b2

			PI = acos(-1.) 
      		eng = ENERGY * 349.755!1153  ! Convert energy from kcal to cm-1
		if(eng.le.1.) then
	                 RHO_SADDLE_POINT = 0.0
		else
 	   		x1=1.E-20
	   		x2=1.E+3
	   		xacc=1.E-12
	   		b2 = zriddr(x1,x2,xacc,eng,FREQ,NFREQ)
	   		sum=0.
	   		prod=1.
	   
	   DO i=1,NFREQ
	  	  	temp1=FREQ(i)/(expsafe(b2*FREQ(i))-1.)
		  	temp1=temp1 * temp1
		  	sum=sum+expsafe(b2*FREQ(i))*temp1
		  	prod = prod / (1.0 - expsafe(-b2*FREQ(i)))
	   ENDDO
	   		RHO_SADDLE_POINT=expsafe(b2*eng) * prod / sqrt(2.0 * pi * sum)
	  endif
	  END FUNCTION RHO_SADDLE_POINT 

	! ============================================================================================
      double precision FUNCTION ZRIDDR(x1,x2,xacc,eng, FREQ, NFREQ)
      IMPLICIT NONE
      INTEGER MAXIT, maxfreq
	  PARAMETER (maxfreq=1000) 
      DOUBLE PRECISION eng,FREQ(maxfreq)
	  INTEGER NFREQ
      DOUBLE PRECISION xacc,funcrzo,UNUSED
      PARAMETER (MAXIT=60,UNUSED=-1.11E30)
      EXTERNAL funcrzo
      INTEGER j
      DOUBLE PRECISION fh,fl,fm,fnew,s,xh,xl,xnew,a,b
      double precision:: xn,xm,x1,x2
      
			fl=funcrzo(x1,eng,FREQ,NFREQ)
      		fh=funcrzo(x2,eng,FREQ,NFREQ)
      	if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        	xl=x1
        	xh=x2
        	zriddr=UNUSED
        DO j=1,MAXIT
          xn = xl+xh
          xm=0.5*xn
          fm=funcrzo(xm,eng,FREQ,NFREQ)
          s=sqrt(fm * fm - fl * fh)
          if(s==0.) return
          a = 1.0
          b = fl-fh
          xnew=xm+(xm-xl)*(sign(a,b)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          	zriddr=xnew
          	fnew=funcrzo(zriddr,eng,FREQ,NFREQ)
          if (fnew==0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr
            fl=fnew
          else
            print*, 'never get here in zriddr'	!old version replace print = pause
          end if
          if(abs(xh-xl).le.xacc) return
        END DO
        	print*, 'zriddr exceed maximum iterations'
        else if (fl==0.) then
        	zriddr=x1
        else if (fh==0.) then
        	zriddr=x2
        else
        	print*, 'root must be bracketed in zriddr'
        end if
        END FUNCTION ZRIDDR
	! ============================================================================================
	DOUBLE PRECISION FUNCTION funcrzo(x, eng, FREQ, NFREQ)
      IMPLICIT NONE
      INTEGER maxfreq 
	  PARAMETER (maxfreq=1000)
	  DOUBLE PRECISION eng, FREQ(maxfreq)
      DOUBLE PRECISION sum,expsafe,divsafe
      INTEGER i,NFREQ
      EXTERNAL expsafe, divsafe 
      double precision :: x,y,z
	
      sum=0.
      DO i=1, NFREQ
        	y = x*FREQ(i)
        	z = expsafe(y)-1.0
         	sum = sum+FREQ(i)/divsafe(z)
      ENDDO
      		funcrzo = eng-sum
      END FUNCTION funcrzo

	! ============================================================================================     
      double precision FUNCTION expsafe(arg)
      IMPLICIT NONE
      double precision:: arg
		if (arg > 700.) then
	                  arg = 700.
		endif 
      		expsafe = exp(arg)
	  END FUNCTION expsafe 
	! ============================================================================================     

	  DOUBLE PRECISION FUNCTION divsafe(arg)
      IMPLICIT NONE
      double precision:: arg
	  
      if (abs(arg) < 1.d-100) then
  	  	arg = 1.d-100
	  endif 
      	divsafe = arg
	  END FUNCTION divsafe 
	! ============================================================================================