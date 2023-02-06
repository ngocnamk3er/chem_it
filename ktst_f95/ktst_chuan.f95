	  Program main_program
	  implicit real*8(a-h,o-z)
      Dimension FVA(50),FVB(50),FVX(99),air(50),dira(50)
      Dimension bir(50),dirb(50),XIR(50),DIRx(50),FVXIR(50),VIRX(50)
      Common /FVS/ FVX,FVA,FVB,NA,NB,NX,NIRx,FVXIR,DIRx,VIRX
      Common /HRA/ H,R,AN,RT,WA,WB,EA
      Common /ABX/ AI1,AI2,AI3,BI1,BI2,BI3,XI1,XI2,XI3,XIR
	  real*16::T,Eaa,Vimag,Qtun
      
      open (unit=5, file='input.dat')
      open (unit=21, file='output.log')

  100 format (20A4)
  101 format (3I10)
  103 format (8F10.4)
  105 format (//,5X,20A4,/,10X,'INPUT DATA:',/)
  106 format (5x,'Temperature',10x,'k(T),cm^3/mole.s',10x,'k(T),cm^3/molecule.s',10x,'Eckart coefficients')
  107 format ('-----------------------------------------------------------------------------------------------------')
  108 format (6x,F8.2,7x,'*',7x,E10.4,8x,'*',9x,E10.4,10x,'*',9x,F10.6)
  109 format (6F10.4)
  112 format (33x,'THE RATE CONSTANTS ARE LISTED BELOW')
  113 format (3I10,5F10.4)
  114 format(/,14x,a)
  115 format(/,'VHR=',F10.4,1x,'kcal/mol',2x,'EA=',F10.4,1x,'kcal/mol'/)
  116 format(14x,a)
      AN=6.0225E+23
      H=2.8592E-03
      R=1.9872E-03

    1 read (5,100) Title
      print *,"------------------------"
      print *, Title
      print *,"------------------------"
      read (5,*)
      read (5,*) na,nb,nx
      read (5,*)
      read (5,*) DGN
      read (5,*)
      read (5,*) VIMAG
      read (5,*)
      read (5,*) FHR,DIHR,RSHR
      read (5,*)
      IF(NA.EQ.0) GO TO 99
      read (5,*) WA,WB
      read (5,*)
      read (5,*) AI1,AI2,AI3
      read (5,*)
      read (5,*) BI1,BI2,BI3
      read (5,*)
      read (5,*) XI1,XI2,XI3
      read (5,*)      
      read (5,*) DL
      read (5,*)
      read (5,*) EA,EAA
      read (5,*)
      read (5,*) Nirx,Nira,Nirb
	  read (5,*)
      IF(NIRx.EQ.0) GO TO 460
      read (5,*) (XIR(I), FVXIR(I),I=1,NIRx)
      READ (5,*) (DIRx(I),I=1,NIRx) 
  460 IF(NIRa.EQ.0) GO TO 470
      read (5,*) (aIR(I),I=1,NIRa)
      READ (5,*) (DIRa(I),I=1,NIRa) 
  470 IF(NIRb.EQ.0) GO TO 500
      read (5,*) (bIR(I),I=1,NIRb)
      READ (5,*) (DIRb(I),I=1,NIRb) 
  500 read (5,*) (FVX(i),i=1,NX)
  	  read (5,*)
      read (5,*) (FVA(i),i=1,NA)
      IF(NB.EQ.0) GO TO 17
      read (5,*)
      read (5,*) (FVB(i),i=1,NB)
   17 write (6,105) TITLE
      write (6,113) NA,NB,NX,DGN,VIMAG,FHR,DIHR,RSHR
      write (6,103) WA,WB,AI1,AI2,AI3,BI1,BI2,BI3
      write (6,109) DL,XI1,XI2,XI3,EA,EAA
      write (6,101) NIRx,nira,nirb
      IF (NIRx.EQ.0) GO TO 510
      WRITE (6,103) (XIR(I),FVXIR(I),I=1,NIRx)
      WRITE (6,103) (DIRx(I),I=1,NIRx)
  510 IF (NIRa.EQ.0) GO TO 520
      WRITE (6,103) (aIR(I),I=1,NIRa)
      WRITE (6,103) (DIRa(I),I=1,NIRa)
  520 IF (NIRb.EQ.0) GO TO 530
      WRITE (6,103) (bIR(I),I=1,NIRb)
      WRITE (6,103) (DIRb(I),I=1,NIRb)
  530   continue
      WRITE(6,103) (FVX(I),I=1,NX)
      WRITE(6,103) (FVA(I),I=1,NA)
      WRITE(6,103) (FVB(I),I=1,NB)

	  VHR = 0.0d0
      write(6,115) VHR,EA
      write(21,*)
	  write(21,112)
      write(21,*)
      write(21,107)
      write(21,106)
      write(21,107)
      read (5,*)
      read(5,*) k
      read (5,*)
      Do j=1,k
      Read(5,*)T
      IF(T)1,1,2
2     RT=R*T
      QHR=1.0
	  if(VIMAG)199,199,200 
199	     Qtun = 1.0
200   Call ECKART(T,EA,EAA,VIMAG,Qtun)
      T1 = T
      Call PARTF(T1,Q,nira,dira,air,nirb,dirb,bir,QX)
      rkcal=DL*(2.083E10)*T*Q*DEXP(-EA/RT)*Qtun*QHR
      rkcal2=rkcal/6.022E+23
      write(21,108) T,  RKCAL,  rkcal2,  QTUN 
  	  End do
     
      write(21,*)
      write(21,107)
      write(21,114)'Thank you for KTST code. First writing by Prof. M.C. Lin - 1965'
      write(21,116)'Then it was modified by J. Park on August 20,1998 and I. Tokmakov on September 2,1998'
      write(21,116)'Now it is modified by P.V. Tien on November 10,2018'

   99 Stop
      End program main_program

!************************************************************************************************************
      SUBROUTINE ECKART(T,E0,E1,VIMAG,QTUN)
      IMPLICIT real*16(a-h,o-z)
	  common /eint/ Emin,gamma,V1,V22,delta,RT      
	  real*8::E0,Emax,RT,Pi
      
      pi = 3.14
      H=9.5377077e-14
      R=1.9872e-3
      CL=2.99792e10

      EPSIL=1.e-6
      freq=CL*VIMAG
      V1=E0-E1
      V22=sqrt(E1)+sqrt(E0)

      gamma=2*sqrt(E0*E1)/(H*freq)
      IF (gamma<0.5) THEN
      	write(6,*) '    E1*E0 IS TOO SMALL. THE PROGRAM WILL STOP'
        stop
      END IF
     	delta=sqrt(gamma**2-.25)
      	RT=R*T
      	EMAX=E0+14*RT
        N=20
        if (E0.LE.E1) then 
        Emin=0
      else 
        Emin=E0-E1
      end if
      ST=(Emax-Emin)/N
      alfa=gamma*sqrt(EMAX)/V22
      beta=gamma*sqrt(EMAX-V1)/V22
      x = 2*Pi*(alfa+beta)
      y = 2*pi*(alfa-beta)
      z = 2*Pi*delta
        c = (exp(x)+exp(-x))/2
        d = (exp(y)+exp(-y))/2
        e = (exp(z)+exp(-z))/2
      
      EDG=(c-d)/2
      EDG=EDG/(c+e)
      EDG=EDG*dexp((E0-EMAX)/RT)/RT

4     call total(N-1,ST,ST,E0,SUM1)
      SUM1=SUM1+EDG
      X0=ST/2

      call total(N,X0,ST,E0,SUM2)
      SUM2=SUM2+SUM1
      N=N*2
      ST=(EMAX-EMIN)/N
      X0=ST/2
      if (abs((SUM2*ST-SUM1*2*ST)/3).LE.EPSIL) go to 7
      SUM1=SUM2
      GO TO 4
7     Qtun=SUM2*ST
      return
	  end subroutine eckart
!************************************************************************************************************
      Subroutine Total(n,x0,ST,E0,Sum)
      Implicit real*16(a-h,o-z)
	  common /eint/ Emin,gamma,V1,V22,delta,RT
      real*8::E0,E,RT,Pi
      pi = 3.14      
      SUM=0
      E=Emin+x0
      Do i=1,n
      alfa=gamma*sqrt(E)/V22
      beta=gamma*sqrt(E-V1)/V22
      x1 = 2*Pi*(alfa+beta)
      y1 = 2*pi*(alfa-beta)
      z1 = 2*Pi*delta
      
      ! x1 = x1/10
      ! y1 = y1/10
      ! z1 = z1/10

        c1 = (exp(x1)+exp(-x1))/2
        d1 = (exp(y1)+exp(-y1))/2
        f1 = (exp(z1)+exp(-z1))/2
      
      S = (c1 - d1)/(c1+f1)
      S=S*dexp((E0-E)/RT)/RT
      SUM=SUM+S
      E=E+ST
      End do
      Return
      end subroutine Total
!************************************************************************************************************
      SUBROUTINE PARTF(T,Q,nira,dira,air,nirb,dirb,bir,QX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FVA(50),FVB(50),FVX(99),air(50),dira(50),dirb(50),bir(50),XIR(50),FVXIR(50),DIRx(50),VIRX(50)
      COMMON /FVS/ FVX,FVA,FVB,NA,NB,NX,NIRx,FVXIR,DIRx,VIRX
      COMMON /HRA/ H,R,AN,RT,WA,WB,EA
      COMMON /ABX/ AI1,AI2,AI3,BI1,BI2,BI3,XI1,XI2,XI3,XIR
      real*4 :: CIR,CF3      

	  !open(10,file='qv.out',status='unknown')
      CT=3.11985E-04
      CR=2.48320E-02
      CRP=6.93580E-03
      QTA=CT*(WA*T)**1.5
      QTB=CT*(WB*T)**1.5
      if (qtb) 888,888,889
  888 qtb=1.
  889 QTX=CT*((WA+WB)*T)**1.5
      IF(XI3) 1,1,2
    1 QRX=CR*DSQRT(XI1*XI2)*T
      GO TO 3
    2 QRX=CRP*DSQRT(XI1*XI2*XI3)*T**1.5
    3 QVA=1.
      QVB=1.
      QRB=1.
      QVX=1.
      DO 4 I=1,NA
    4 QVA=QVA/(1.-DEXP(-H*FVA(I)/RT))
	!write(10,*) T, QVA
      IF(AI1) 5,5,6
    5 QRA=CR*DSQRT(AI3*AI2)*T
      GO TO 7
    6 QRA=CRP*DSQRT(AI1*AI2*AI3)*T**1.5
    7 IF(NB.EQ.0) GO TO 11
      DO 8 I=1,NB
    8 QVB=QVB/(1.-DEXP(-H*FVB(I)/RT))
      IF(BI3) 9,9,10
    9 QRB=CR*DSQRT(BI1*BI2)*T
      GO TO 11
   10 QRB=CRP*DSQRT(BI1*BI2*BI3)*T**1.5
   11 DO 12 I=1,NX
   12 QVX=QVX/(1.-DEXP(-H*FVX(I)/RT))
      Q=QTX/(QTA*QTB)*QRX/(QRA*QRB)*QVX/(QVA*QVB)

      CIR= 0.279321

      IF(NIRx.EQ.0) GO TO 20

      DO 15 I=1,NIRx
      QX=CIR*(DSQRT(XIR(I)*T))/DIRx(I)

      If(FVXIR(I)) 16,16,17

   16 VIRx(I)=0.
      CF=1.
      GO TO 18

   17 VIRx(I)=1.0215E-04*XIR(I)*FVXIR(I)**2/DIRx(I)**2


      CF1=(1.-exp(-RT/VIRx(I)))**1.2
      CF2=exp(-1.2*RT/VIRx(I))
      CF3=(Qx*(1.-exp(-DIRx(I)/QX*SQRT(3.14159*VIRx(I)/RT))))
      CF=CF1+(CF2/CF3)

   18 Q=Q*QX*CF
    
   15 CONTINUE

   20 IF(NIRa.EQ.0) GO TO 25
      DO 22 I=1,NIRa 
      QA=(CIR*(DSQRT(aIR(I)*T))/DIRa(I))
   22 Q=Q/QA

   25 IF(NIRb.EQ.0) GO TO 30
      DO 27 I=1,NIRb 
      QB=(CIR*(DSQRT(bIR(I)*T))/DIRb(I))
   27 Q=Q/QB

   30 RETURN
      END
!************************************************************************************************************




