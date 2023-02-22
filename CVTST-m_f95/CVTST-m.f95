      PROGRAM VTST

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 TITLE
      CHARACTER*20 SPECIES
      CHARACTER*3 ZPE, RES, ASK, TUN
      DIMENSION TEMP(100),GEO(100),DGT(100,100),W(100),E(100),GDAT(100)
      DIMENSION RCV(100), TV(100), EZ(100), EZT(100), GEOM(100),M(10)
      DIMENSION AA(500), BB(500), X(100), WW(100),GG(100),dGdR(100)


!        character*32 file1, file2
!        write(*,*)'input file name'
!        read(*,*) file1
!      OPEN (UNIT=10,FILE=file1,STATUS='unknown')

      OPEN (UNIT=10,FILE='bimo.inp')
      OPEN (UNIT=11,FILE='bimo.out',STATUS='UNKNOWN')
      OPEN (UNIT=12,FILE='bimo.paste',STATUS='UNKNOWN')
      REWIND 10
      REWIND 11
      REWIND 12
  101 format (2x,'Geometrical Parameters:')
  102 format (2x,'--------- dG at' ,F10.2,' Degrees Kelvin:')
  103 format (4x,'R [Angstrom]    G [kcals/mol]    dG/dR')  
  104 format (F15.4,F15.6,5x,'   NAN')     
  105 format (2x,'Greatest Gibbs Free-Energy Barriers [kcals/mol]')
  106 format (2x,'at the temperatures specified:')  
  107 format (2x,'T [Kelvin]  R [Angstrom]  G [kcals/mol]')
  108 format (2x,'Same Barriers from quadratic fit:')  
  109 format (2x,'T [Kelvin]  G [kcals/mol]')
  110 format (2x,'dG(FIT) = ',E12.6,' T^2 ',E12.6,' T + ',E12.6,'  [Calories])')
  111 format (2x,'Rate Constant:')
  112 format (2x,'k = ',E12.6,' T^',I1,' exp(-dG/RT)  [cc/molecule]^',I1,'/second')
  113 format (20x,'dG = Raw Data       dG = Quad. Fit')  
  114 format (2x,'at T =' ,F8.2,'    k = ',E12.6,'    k = ',E12.6)  
  

      READ(10,"(A80)")TITLE
      WRITE(11,"(A80)") TITLE
      READ(10,"(28X,I3)") NP !Total number of structures
      READ(10,"(21X,I3)") NR !Number of reactants
      READ(10,"(A80)")TITLE
      READ(10,"(10I5)") M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8)
      NN = M(1)+M(2)+M(3)+M(4)+M(5)+M(6)+M(7)+M(8)
      READ(10,"(24X,I3)") NT !Number of temperatures
      READ(10,"(26X,A3)") ZPE
      READ(10,"(30X,A3)") TUN
      READ(10,"(A80)")TITLE
      RES = "bar"
      READ(10,"(23X,A3)") ASK
      IF(ASK.EQ."yes") RES="atm"
      READ(10,"(16X,A3)") ASK
      IF(ASK.EQ."yes") RES="bar"
      READ(10,"(A80)")TITLE
!-------------------------- CONSTANTS ----------------------------------C
!     Physical constants are from Rev. Mod. Phys., 59, 1121 (1987).
      FAC =  0.46499834D0
!         = (g.cm.cm/au)*10^40   [conversion factor]
!         = ((0.5291772E-8)^2/6.02214E23)*10^40

      R=1.9872D0
!      = Nk = "gas constant"  [cal/(mole.kelvin)]

      RCC1= 2.08367541593D10
!         = k/h  = 1.38066E-16/6.62608E-27
      IF(RES.EQ."atm") RCC2= 1.3626055D-22
!                          = k/pressure factor, = 1.38066E-16/1.01325E6
      IF(RES.EQ."bar") RCC2= 1.38066D-22
!                          = k/pressure factor, = 1.38066E-16/1.0E6
      RCC = RCC1 * (RCC2)**(NN-1)

      C1=1.438767553D0
!       = hc/k = 6.62608E-27 *2.99792458E+10 /1.38066E-16 [cm.molecule.kelvin]
!--------------------------- Read Temperatures -------------------------C

      DO 10 I=1,NT
         READ(10,"(F10.3)") TEMP(I)
   10 CONTINUE
      DO 20 I=1,10
         READ(10,"(A80)")TITLE
   20 CONTINUE
!
!------------------------ Loop over structures -------------------------C
      DO 30 L=1,NP
         READ(10,*) TITLE
         READ(10,*) GEO(L), E(L), SPECIES
         E(L) = E(L)*627.50595D0
           WRITE(11,"(2F15.8,5X,A10)") GEO(L), E(L)/627.50595, SPECIES
           BE = 0.0
           AE = 0.0
           WE = 0.0
           WEX = 0.0

           IF(SPECIES.EQ."diatomic") Then 
             READ(10,*) WT,G0,G1,EE
             READ(10,*) BE,AE,WE,WEX
            ELSE
             READ(10,*) WT,G0,G1,EE
           ENDIF
 
!        A,B,C and F must be read as atomic units:
         READ(10,*) A,B,C,SN,SF,F,NF,NV
         A = A * FAC
         B = B * FAC
         C = C * FAC
         F = F * FAC
         IF(NV.NE.0) READ(10,*)(W(j),j=1,NV)
         IF(TUN.EQ."yes") READ(10,*) RCV(L)
         CALL GIBBS (NT,TEMP,L,SPECIES,WT,G0,G1,EE,BE,AE,WE,WEX,SN,SF,A,B,C,F,NF,W,NV,DGT,E(L),RES,ZPE,EZ)
   30 CONTINUE

!------------- List "Delta G" over range of temperatures ---------------C
      IF(NR.GT.1) THEN
        DO 32 J=1,NT
           DGT(NR,J) = M(NR)*DGT(NR,J)
   32   CONTINUE
        EZ(NR) = M(NR)*EZ(NR)
        DO 40 I=1,NR-1
           EZ(NR) = EZ(NR) + M(I)*EZ(I)
           DO 35 J=1,NT
              DGT(NR,J) = DGT(NR,J) + M(I)*DGT(I,J)
   35      CONTINUE
   40   CONTINUE
      ENDIF
      WRITE(11,101)
      DO 50 I=NR,NP
         WRITE(11,"(F15.8)") GEO(I)
   50 CONTINUE

!     Save Geometric Parameters of dG/dR
 
      DO 80 I=1,NT
         WRITE(11,102) TEMP(I) 
         GDAT(I) = -10.0D99
         DO 70 J=NR,NP
            DD = DGT(J,I) - DGT(NR,I)
!            WRITE(11,"(F15.6)") DD
!            WRITE(11,"(2F15.6)") DGT(J,I) , DGT(NR,I)
            GG(J)=DD   
            IF(DD.GT.GDAT(I)) THEN
               GDAT(I) = DD
               IF(TUN.EQ."yes") THEN
               	 TV(I)   = RCV(J)
               ELSE
                 TV(I)   = 0
               ENDIF
               EZT(I)  = EZ(J) - EZ(NR)
	         GEOM(I)=GEO(J)
            ENDIF
   70    CONTINUE

      WRITE(11,103) 
!     Do dG/dR (Central Difference)
         DO 71 J=NR,NP  
            IF(J==NR .OR. J==NP) THEN
             WRITE(11,104) GEO(J), GG(J)
            ELSE 
		   dGdRA=(GG(J+1)-GG(J))/(GEO(J+1)-GEO(J))
	       dGdRB=(GG(J)-GG(J-1))/(GEO(J)-GEO(J-1))
		   dGdR(J)=(dGdRA+dGdRB)/2
		   WRITE(11,"(F15.4, F15.6, F15.6)") GEO(J), GG(J), dGdR(J)
		  ENDIF
   71    CONTINUE 

   80 CONTINUE

      WRITE(11,"(1X)") 
      WRITE(11,105) 
      WRITE(11,106) 
      WRITE(11,107) 
    	DO 90 I=1,NT 
		  WRITE(11,"(F10.4, F15.4,  F15.6)") TEMP(I), GEOM(I), GDAT(I)
!            WRITE(11,"(F15.6)") GDAT(I)
   90 CONTINUE

!-----------------------------------------------------------------------C
!- Do quadratic least squares fit to GDAT as function of temperature ---C

      DO 100 I=1,NT
         AA(I) = TEMP(I)*TEMP(I)
         AA(I+NT) = TEMP(I)
         AA(I+(2*NT)) = 1.0D0
         BB(I) = GDAT(I)
         !WRITE(6,*) AA(I), AA(I + NT), AA(I + 2*NT)
  100 CONTINUE

      CALL HOUSEFIT (AA,NT,3,BB,NT,1,X,WW,ISING)

      WRITE(11,"(1X)") 
	  WRITE(11,108) 
      WRITE(11,109) 
	DO 130 I=1,NT
         PRED= X(1)*TEMP(I)*TEMP(I) + X(2)*TEMP(I) + X(3)
         WRITE(11,"(F10.4, F15.6)") TEMP(I), PRED
  130 CONTINUE
!     Convert to Calories:
      X(1) = X(1) * 1000.0D0
      X(2) = X(2) * 1000.0D0
      X(3) = X(3) * 1000.0D0

      WRITE(11,"(1X)") 
      WRITE(11,110)  X(1), X(2), X(3)
      WRITE(11,"(1X)") 
!------------- Print out Rate Constant equation and k @ temps ----------C
      DO 135 I=1,NT
         WRITE(12,"(F10.3)") TEMP(I)
  135 CONTINUE
      WRITE(12,"(1X)") 

      WRITE(11,"(1X)") 
      WRITE(11,111) 
      WRITE(11,"(1X)") 
      WRITE(11,112) RCC, NN, NN-1
      WRITE(11,"(1X)") 
      WRITE(11,113)
      WRITE(11,"(1X)") 
      DO 140 I=1,NT
         PRED= X(1)*TEMP(I)*TEMP(I) + X(2)*TEMP(I) + X(3)
         FRC = RCC * (TEMP(I)**NN) * DEXP(-PRED/(R*TEMP(I)))
         RRC = RCC * (TEMP(I)**NN) * DEXP(-GDAT(I)*1000.0D0/(R*TEMP(I)))
         WRITE(11,114) TEMP(I), RRC, FRC
         WRITE(12,"(E12.6)") RRC
  140 CONTINUE
      IF(TUN.NE."yes") GOTO 900
      WRITE(11,"(1X)") 
      WRITE(12,"(1X)") 
      WRITE(11,"(47x)") 'Tunneling factor = [1 +(hcw/kT)^2(1+RT/Eo)/24])'
      WRITE(11,"(49x)") 'where w(cm^-1) = negative frequency (along path))'
      WRITE(11,"(48x)") 'Note: This factor only valid for small barriers)'
      WRITE(11,"(47x)") 'with a width of less than about 1.8 Angstroms.)'
      WRITE(11,"(1X)") 
      WRITE(11,"(33x)") 'Rates with Tunneling Correction:)'
      WRITE(11,"(1X)") 
      DO 150 I=1,NT
         PRED= X(1)*TEMP(I)*TEMP(I) + X(2)*TEMP(I) + X(3)
         FRC = RCC * (TEMP(I)**NN) * DEXP(-PRED/(R*TEMP(I)))
         RRC = RCC * (TEMP(I)**NN) * DEXP(-GDAT(I)*1000.0D0/(R*TEMP(I)))
         TF  = (C1*TV(I)/TEMP(I))**2.0D0
         TF  = TF * (1.0D0 + R*TEMP(I)/(EZT(I)*1000.0D0))/24.0D0
         TF  = TF + 1.0D0
         RRC = TF * RRC
         FRC = TF * FRC
         WRITE(11,114) TEMP(I), RRC, FRC
         WRITE(12,"(E12.6)") RRC
  150 CONTINUE
!-----------------------------------------------------------------------C
  900 CLOSE (UNIT=10)
      CLOSE (UNIT=11)

      END

!///////////////////////////////////////////////////////////////////////C

      SUBROUTINE GIBBS (NT,TEMP,NN,SPECIES,WT,G0,G1,EE,BE,AE,WE,WEX,SN,SF,A,B,C,F,NF,W,NV,DGT,E,RES,ZPE,EZ)
!=======================================================================C
!                                                                       C
!     Calculates Gibbs Free Energy using standard Statistical Thermo.   C
!                                                                       C
!                                                                       C
!     NT      = Number of Temperatures.                                 C
!     SPECIES = "monatomic" or "diatomic" or "linear" or "nonlinear".   C
!     WT      = MOLECULAR WEIGHT in amu                                 C
!     SN      = Rotational Symmetry Number.                             C
!     SF      = SYMMERTRY NUMBER OF FREE ROTATION, IF NO FR, SET TO 0.  C
!                                                                       C
!         ELECTRONIC PARTITION FUNCTION = G0 + G1*EXP(-EE/RT)           C
!     G0, G1  =  DEGENERACIES FOR GROUND STATE & FIRST EXCITED STATE    C
!         EE in cm-1, WHEN EE IS LARGE ENOUGH, PUT ZERO TO IGNORE.      C
!                                                                       C
!         FOR DIATOMIC MOLECULES, SET THE FOLLOWING PARAMETERS :        C
!     BE      = ROTATIONAL CONSTANT IN cm-1                             C
!     AE      = ROTATION-VIBRATION INTERACTION CONSTANT IN cm-1         C
!     WE      = FUNDAMENTAL VIBRATONAL CONSTANT IN cm-1                 C
!     WEX     = VIBRATIONAL ANHARMONICITY CONSTANT IN cm-1              C
!                                                                       C
!        FOR POLYATOMIC MOLECULES, SET THE FOLLOWING PARAMETERS  :      C
!     A       = MOMENT OF INERTIA FOR LINEAR CASE  (au,ie g.bohrs^2)    C
!     A,B,C   = MOMENT OF INERTIA FOR NONLINEAR CASE  (au)              C
!     F       = MOMENT OF INERTIA FOR FREE ROTATION  (au)               C
!     NF      = DEGREES OF FREEDOM FOR FREE ROTATION                    C
!     W(I)    = VIBRATIONAL FREQUENCIES IN cm-1                         C
!     NV      = NUMBER OF VIBRATIONAL FREQUENCIES                       C
!                                                                       C
!=======================================================================C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10 SPECIES
      CHARACTER*3 RES,ZPE
      DIMENSION W(100),V(100)
      DIMENSION TEMP(100), DGT(100,100), EZ(100)


!-------------------------- CONSTANTS ----------------------------------C
!     Physical constants are from Rev. Mod. Phys., 59, 1121 (1987).
      C1=1.438767553D0
!       = hc/k = 6.62608E-27 *2.99792458E+10 /1.38066E-16 [cm.molecule.kelvin]
      C2=40.27512295D0
!       = (h*h/k.8pi*pi) * 10^40   (erg.sec*sec.kelvin)
      R=1.9872D0
!      = Nk = "gas constant"  [cal/(mole.kelvin)]
      PI=DACOS(-1.0D0)
      IF(RES.EQ."atm") TC= -3.664854875D0
!       = translational constant = ln [ (k/1.01325E6)(2PI*k/N*h*h)^1.5 ]
!            i.e. in atmospheres.
      IF(RES.EQ."bar") TC= -3.651691889D0
!       = translational constant = ln [ (k/1.0E6)(2PI*k/N*h*h)^1.5 ]
!            i.e. in Bars.
!-----------------------------------------------------------------------C
!----------------------- Loop over Temperatures ------------------------C
      DO 400 II=1,NT
      T =TEMP(II)

!     Electronic Contribution to (Gibbs Free Energy)/-RT:
      IF(EE==0)  CELECT = DLOG(G0)
      IF(EE.NE.0)  CELECT = DLOG(G0 + G1*DEXP(C1*EE/T))

!     Translational Contribution to (Gibbs Free Energy)/-RT:
      CTRANS = 1.5D0 *DLOG(WT) +2.5D0 *DLOG(T) +TC

!----------------------- Monatomic Species -----------------------------C
      IF(SPECIES.EQ."monatomic ") G = (CELECT + CTRANS)*R*T
      IF(SPECIES.EQ."monatomic ") E0= 0.0D0

!----------------------- Diatomic Species ------------------------------C
      IF(SPECIES.EQ."diatomic  ") THEN
        Y=C1*(BE-AE/2)/T
        U=C1*(WE-2*WEX)/T
        XE=WEX/WE
        D=AE/BE
        GG=BE/WE
        IF(U.GT.80)THEN
          EU=EXP(80.)
          U1=0.
         ELSE
          EU=DEXP(U)-1
          U1= U/EU/EU
        ENDIF      
!   20 CONTINUE
       IF(ZPE.NE."yes") E0 =  0.0D0
       IF(ZPE.EQ."yes") E0 =  - U/2.0D0
!                          =  Zero-Point Energy correction + Corrections for anharmonic vibrations:
       G= CTRANS  - DLOG(SN*Y) + Y/3 + Y**2/90 - DLOG(1-DEXP(-U)) + CELECT + E0 + (2*XE*U1) + D/EU + 8*GG/U
       G = G*R*T
       E0= E0*R*T
      ENDIF

!----------------------------- Polyatomic Species ---------------------C
!     Vibrational Contribution to (Gibbs Free Energy)/-RT:
!     NOTE WELL: The zero-point energy correction MAY BE INCLUDED;
      E0   = 0.0D0
      CVIB = 0.0D0
      DO 30 I=1,NV
         V(I)=C1*W(I)/T
!            = hv/kT  i.e. hcw/kT , w in reciprocal centimeters.
         IF(ZPE.EQ."yes") E0 = E0 - V(I)/2.0D0
!                              =      +  Zero-Point Energy correction.
         CVIB = CVIB - DLOG(1.0D0 -DEXP(-V(I)))
   30 CONTINUE
      CVIB = CVIB + E0

!     Linear Polyatomic Species:
      IF(SPECIES.EQ."linear    ")THEN
        Y=C2/A/T
      IF(F.LE.0)THEN
        FR=0.
        GO TO 32
      ENDIF 
        YY=C2/F/T
      IF(NF-2)35,36,37
   35   FR=.5*DLOG(SF**2*YY/PI)
        GO TO 32
   36   FR=DLOG(SF*YY)
        GO TO 32
   37   FR=.5*DLOG(SF**2*YY**3/PI)     
   32   G =CTRANS - FR - DLOG(SN*Y) + Y/3 + Y**2/90 + CVIB + CELECT
      G = G*R*T
      E0= E0*R*T
      ENDIF

!     Nonlinear Polyatomic Species:
      IF(SPECIES.EQ."nonlinear ")THEN
        Y=(C2/A/T)*(C2/B/T)*(C2/C/T)
      IF(F.LE.0)THEN
        FR=0.
        GO TO 33
      ENDIF
        YY=C2/F/T
      IF(NF-2)38,39,40  
   38   FR=.5*DLOG(SF**2*YY/PI)
         GO TO 33
   39   FR=DLOG(SF*YY)
         GO TO 33
   40   FR=.5*DLOG(SF**2*YY**3/PI)
   33   G= CTRANS - FR - .5*DLOG(SN**2*Y/PI) + CVIB + CELECT
        G = G*R*T
        E0= E0*R*T
      ENDIF
!         Get G in kcals:
         DGT(NN,II) = E  -G/1000.0D0
         EZ(NN) = E - E0/1000.0D0
  400 CONTINUE
!-----------------------------------------------------------------------C
      RETURN
      END

!///////////////////////////////////////////////////////////////////////

      SUBROUTINE HOUSEFIT(A,M,N,B,K,L,X,W,ISING)
!======================================================================C
!     Solve a system of simultaneous equations, AX=B. If B is a matrix C
!   then X, the  matrix of multiple solutions,  must  have  the same   C
!   dimensions as B.                                                   C
!     The system is allowed to be  overdetermined, in  which case you  C
!   will get the least-squares solution(s).                            C
!======================================================================C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), B(K,L), X(K,L), W(M)
      CALL HOUSETRANS(A,M,N,B,K,L,W,ISING)
      CALL BACKSUB(A,M,N,B,K,L,X)
      RETURN
      END

!///////////////////////////////////////////////////////////////////////

      SUBROUTINE HOUSETRANS(A,M,N,B,KB,LB,W,ISING)
!======================================================================C
!                 Householder transformation.                          C
!======================================================================C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(M,N) , B(KB,LB) , W(M)
      ISING=0
      DO 300 K=1,N
         MX=IDAMAX(M-K+1,A(K,K),1) + K-1
         write(6,*) MX
         IF(A(MX,K).EQ.0) ISING=K
!                   i.e. A is a singular matrix; Kth column became zero.
!         IF(A(MX,K).EQ.0) THEN
!           ISING=K
!           GOTO 400
!         ENDIF
         RMS=0.0
         DO 50 I=K,M
            W(I)=A(I,K)/DABS(A(MX,K))
            RMS=RMS+W(I)*W(I)
   50    CONTINUE
         RMS=DSQRT(RMS)
         BK=1/(RMS*(RMS + DABS(W(K))))
         ZZ=W(K)
         W(K)=W(K) + DSIGN(RMS,ZZ)
         DO 100 J=1,N
            S=DDOT(M-K+1,W(K),1,A(K,J),1)
            S=BK*S
            CALL DAXPY(M-K+1,-S,W(K),1,A(K,J),1)
  100    CONTINUE
         DO 200 J=1,LB
            S=DDOT(M-K+1,W(K),1,B(K,J),1)
            S=BK*S
            CALL DAXPY(M-K+1,-S,W(K),1,B(K,J),1)
  200    CONTINUE
  300 CONTINUE
!  400 RETURN
      RETURN
      END

!///////////////////////////////////////////////////////////////////////

      SUBROUTINE BACKSUB(A,M,N,B,KB,LB,X)
!======================================================================C
!                       Execute backsubstitution.                      C
!======================================================================C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N) , B(KB,LB) , X(KB,LB)
      CALL DCOPY(KB*LB,B,1,X,1)
      DO 300 L=1,LB
         DO 200 J=N,2,-1
            X(J,L)=X(J,L)/A(J,J)                        
            CALL DAXPY(J-1,-X(J,L),A(1,J),1,X(1,L),1)  
  200    CONTINUE                                   
         X(1,L)=X(1,L)/A(1,1)                 
  300 CONTINUE                                  
      RETURN    
      END

!/////////////////////////////////////////////////////////////////////C

!         Basic Linear Algebra Subroutines

      subroutine dswap(n,dx,incx,dy,incy)
      implicit double precision (a-h,o-z)
      dimension dx(*),dy(*)
!     swaps the vectors dx and dy
      if (incx.eq.1.and.incy.eq.1) then
         do 100 i=1,n
           s=dx(i)
           dx(i)=dy(i)
           dy(i)=s
 100     continue
      else
         ix=1
         iy=1
         do 200 i=1,n
           s=dx(ix)
           dx(ix)=dy(iy)
            dy(iy)=s
           ix=ix+incx
           iy=iy+incy
 200     continue
      end if
      end

      double precision function ddot(n,dx,incx,dy,incy)
      implicit double precision (a-h,o-z)
      dimension dx(*),dy(*)
!     scalar (dot) productt of two vectors

      s=0.0d0
      if(incx.eq.1.and.incy.eq.1) then
         do 100 i=1,n
           s=s+dx(i)*dy(i)
 100     continue
      else
         ix=1
         iy=1
         do 200 i=1,n
           s=s+dx(ix)*dy(iy)
           ix=ix+incx
           iy=iy+incy
 200     continue
      end if
      ddot=s    
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
      implicit double precision (a-h,o-z)
      dimension dx(*),dy(*)
!     this is the famous daxpy operation  (outer product)
!     - add to the vector dy da times the vector dx
      if (incx.eq.1.and.incy.eq.1) then
         do 100 i=1,n
            dy(i)=dy(i)+da*dx(i)
 100     continue
      else
         ix=1
         iy=1
         do 200 i=1,n
            dy(iy)=dy(iy)+da*dx(ix)
            iy=iy+incy
            ix=ix+incx
 200     continue
      end if
      end

      subroutine dscal(n,da,dx,incx)
      implicit double precision (a-h,o-z)
      dimension dx(*)
!     multiply a vector by a scalar
      if(incx.eq.1) then
         do 100 i=1,n
            dx(i)=da*dx(i)
 100     continue
      else
         ix=1
         do 200 i=1,n
            dx(ix)=da*dx(ix)
            ix=ix+incx
 200     continue
      end if
      end

      integer function idamax(n,dx,incx)
      implicit double precision (a-h,o-z)
      dimension dx(*)
!     the value of the function is the index of the largest
!     element (in absolute value) of the vector dx
         ii=1
         xmax=abs(dx(1)) 
      if(incx.eq.1) then
         do 100 i=2
           WRITE(6,*) abs(dx(i)),n
           if (abs(dx(i)).gt.xmax) then
               xmax=abs(dx(i))
               ii=i
           end if
 100     continue
      else
         ix=1
         do 200 i=2,n
            ix=ix+incx
            if(abs(dx(ix)).gt.xmax) then
               xmax=abs(dx(ix))
               ii=i
            end if
 200     continue
      end if
      idamax=ii
      end
      subroutine dcopy(n, dx, incx, dy, incy)
      implicit double precision (a-h,o-z)
      dimension dx(*),dy(*)
!     performs the elementary vector operation dy=dx
      if (incx.eq.1.and.incy.eq.1) then
         do 100 i=1,n
           dy(i)=dx(i)
 100     continue
      else
         ix=1
         iy=1
         do 200 i=1,n
           dy(iy)=dx(ix)
           iy=iy+incy
           ix=ix+incx
 200     continue
      end if
      end

      subroutine  drot (n,dx,incx,dy,incy,c,s)
      implicit double precision(a-h,o-z)
      dimension dx(*),dy(*)

      if (abs(c).gt.abs(s)) then
        c1=c
        s1=s
        s2=-s1/c1
        c2=c1-s2*s1
      else
        c1=s
        s1=c
        call dswap (n,dx,incx,dy,incy)
        s2=s1/c1
        c2=-s2*s1-c1
      endif
      call dscal (n,c1,dx,incx)
      call daxpy (n,s1,dy,incy,dx,incx)
      call dscal (n,c2,dy,incy)
      call daxpy (n,s2,dx,incx,dy,incy)
      end
