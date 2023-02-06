import math
import venv
import numpy
# funtion drot
def drot(n,dx,incx,dy,incy,c,s):
    if(abs(c) > abs(s)):
        c1 = c
        s1 = s
        s2 = -s1/c1
        c2 = c1 - s2*s1
    else:
        c1 = s
        s1 = c
        dswap(n,dx,incx,dy,incy)
        s2 = s1/c1
        c2 = -s2*s1 - c1
    dscal(n,c1,dx,incx)
    daxpy(n,s1,dy,incy,dx,incx)
    dscal(n,c2,dy,incy)
    daxpy(n,s2,dx,incx,dy,incy)
#end funtion drot

#funtion dcopy
def dcopy(n, dx, incx, dy, incy):
    if incx == 1 and incy == 1:
        i = 1
        while i <= n:
            dy[i] = dx[i]
            i += 1
    else:
        ix = 1
        iy = 1
        i = 1
        while i <= n:
            dy[iy] = dx[ix]
            iy = iy + incy
            ix = ix + incx 
            i += 1
#end funtion dcopy    

#function idamax
def idamax(n,dx,incx):
    ii = 1
    xmax = abs(dx[1])
    if (incx == 1):
        i = 1
        while(i <= n):
            if (abs(dx[i]) > xmax):
                xmax = abs(dx[i])
                ii = i
            i += 1
    else:
        ix = 1
        i = 1
        while (i <= n):
            ix = ix + incx
            if(abs(dx[ix]) > xmax):
                xmax = abs(dx[ix])
                ii = i
            i += 1 
    idamax = ii
    return idamax
# end function idamax

# function dscal       
def dscal(n,da,dx,incx):
    if (incx == 1):
        i = 1
        while (i <= n):
            dx[i] = da*dx[i]
            i += 1
    else:
        ix = 1
        while (i <= n):
            dx[ix] = da*dx[ix]
            ix = ix + incx
            i += 1
#end function dscal

#function daxpy
def daxpy(n,da,dx,incx,dy,incy):
    if (incx == 1 and incy == 1):
        i = 1
        while (i <= n):
            dy[i] = dy[i] + da*dx[i]
            i += 1
    else:
        ix = 1
        iy = 1
        i = 1
        while (i <= n):
            dy[iy] = dy[iy] + da*dx[ix]
            iy = iy + incy
            ix = ix + incx
            i += 1
#end function daxpy

#function ddot
def ddot(n,dx,incx,dy,incy):
    s = 0
    if (incx == 1 and incy == 1):
        i = 1
        while (i <= n):
            s = s + dx[i]*dy[i]
            i += 1
    else:
        ix = 1
        iy = 1
        i = 1
        while (i <= n):
            s = s + dx[ix]*dy[iy]
            ix = ix + incx
            iy = iy + incy
            i += 1
    ddot = s        
    return ddot
#end function ddot

#funtion dswap
def dswap(n,dx,incx,dy,incy):
    if(incx == 1 and incy == 1):
        i = 1
        while(i <= n):
            s = dx[i]
            dx[i] = dy[i]
            dy[i] = s
            i += 1
    else:
        ix = 1
        iy = 1
        i = 1
        while (i <= n):
            s = dx[ix]
            dx[ix] = dy[iy]
            dy[iy] = s
            ix = ix + incx
            iy = iy + incy
            i += 1
#end funtion dswap

# function backsub
def BACKSUB(A,M,N,B,KB,LB,X):
    A[M,N],B[KB,LB],X[KB,LB]
    dcopy(KB*LB,B,1,X,1)
    L = 1
    while (L <= LB):
        J == N
        while (J >= 2):
            X[J,L] = X[J,L]/A[J,J]
            daxpy(J-1,-X[J,L],A[1,J],1,X[1,L],1)
            J -= 1
        X[1,L] = X[1,L]/A[1,1]
# end function backsub  

# function housetrans
def HOUSETRANS(A,M,N,B,KB,LB,W,ISING):
    A[M,N],B[KB,LB],W(M)
    ISING = 0
    K = 1 
    while (K <= N):
        MX = idamax(M-K+1,A[K,K],1) + K - 1
        if (A[MX,K] == 0): ISING = K
        RMS = 0.0
        I = K
        while (I <= M):
            W[I] = A[I,K]/abs(A[MX,K])
            RMS = RMS + W[I]*W[I]
            I += 1
        RMS = math.sqrt(RMS)
        BK = 1/(RMS*(RMS + abs(W[K])))
        ZZ = W[K]
        W[K] = W[K] + numpy.sign(RMS,ZZ)
        J = 1
        while (J <= N):
            S = ddot(M-K+1,W[K],1,A[K,J],1)
            S = BK*S
            daxpy(M-K+1,-S,W[K],1,A[K,J],1)
            J += 1
        J = 1
        while (J <= LB):
            S = ddot(M-K+1,W[K],1,B[K,J],1)
            S = BK*S
            daxpy(M-K+1,-S,W[K],1,B[K,J],1)
            J += 1
        K += 1
#end funtion housetrans

#housefit
def HOUSEFIT(A,M,N,B,K,L,X,W,ISING):
    A[M,N],B[K,L],X[K,L],W[M]
    HOUSETRANS(A,M,N,B,K,L,W,ISING)
    BACKSUB(A,M,N,B,K,L,X)
#end housefit

#gibbbs
def GIBBS (NT,TEMP,NN,SPECIES,WT,G0,G1,EE,BE,AE,WE,WEX,SN,SF,A,B,C,F,NF,W,NV,DGT,E,RES,ZPE,EZ,V):
    W[100],V[100],TEMP[100], DGT[100,100], EZ[100]
    C1=1.438767553
    C2=40.27512295
    R=1.9872
    PI=numpy.arccos(-1.0)
    if (RES == "atm"): TC = -3.664854875
    if (RES == "bar") : TC = -3.651691889
    #Loop over Temperatures
    II = 1
    while (II <= NT):
        T = TEMP(II)
        if (EE == 0): CELECT = numpy.log(G0)
        if (EE != 0): CELECT = numpy.log(G0 + G1*numpy.exp(C1*EE/T))
        
        CTRANS = 1.5*numpy.log(WT) + 2.5*numpy.log(T) + TC
        
        if(SPECIES == "monatomic"): G = (CELECT + CTRANS)*R*T
        if(SPECIES == "monatomic"): E0 = 0
        #Diatomic Species    
        if(SPECIES == "diatomic"):
            Y=C1*(BE-AE/2)/T
            U=C1*(WE-2*WEX)/T
            XE=WEX/WE
            D=AE/BE
            GG=BE/WE
            if(U > 80):
                EU = numpy.exp(80)
                U1 = 0
            else: 
                EU = numpy.exp(U) - 1
                U1 = U/EU/EU
        if(ZPE != "yes"): E0 = 0
        if(ZPE == "yes"): E0 = -U/2
        G= CTRANS  - numpy.log(SN*Y) + Y/3 + Y**2/90 - numpy.log(1-numpy.exp(-U)) + CELECT + E0 + (2*XE*U1) + D/EU + 8*GG/U
        G = G*R*T
        E0= E0*R*T
            
        #Polyatomic Species        
        E0   = 0
        CVIB = 0
        for I in range (1,NV+1):
            V[I] = C1*W[I]/T
            if (ZPE == "yes"):
                E0 = E0 -V[I]/2
                CVIB = CVIB - numpy.log(1-numpy.exp(-V[I]))
        CVIB = CVIB + E0
        if(SPECIES == "linear   "):
            Y= C2/A/T
            if(F<=0):
                FR = 0
            else: 
                YY=C2/F/T
                if(NF-2<0):
                    FR=0.5*math.log(math.pow(SF,2)*YY/PI)
                elif(NF-2 ==0):
                    FR=math.log(SF*YY)
                else:
                    FR= 0.5*math.log(math.pow(SF,2)*math.pow(YY,3)/PI)
            G = CTRANS - FR - math.log(SN*Y) + Y/3 + math.pow(Y,2)/90 + CVIB + CELECT
            G = G*R*T
            E0 = E0*R*T
        if(SPECIES == "nonlinear   "):
            Y=(C2/A/T)*(C2/B/T)*(C2/C/T)
            if(F<=0):
                FR=0
            else:
                YY=C2/F/T
            if(NF-2<0):
                    FR=0.5*math.log(math.pow(SF,2)*YY/PI)
            elif(NF-2==0):
                    FR= math.log(SF*YY)
            else:   
                    FR=0.5*math.log(math.pow(SF,2)*math.pow(YY,3)/PI)
            G = CTRANS -FR - 0.5*math.log(math.pow(SN,2)*Y/PI) + CVIB + CELECT
            G = G*R*T
            E0 = E0*R*T
        
        DGT[NN,II]= E - G/1000
        EZ[NN] = E - E0/1000
    II += 1
#end gibbs
def main():
    pass


if __name__ == "__main__":
    main()
    # inputArr = []
    # file1 = open('bimo.inp', 'r')
    # Lines = file1.readlines()
    # for line in Lines:
    #     inputArr.append(line.strip())
     
    # #Gán giá trị các biến 
    # Title = inputArr[0]
    # NP = int(inputArr[1][29:32])
    # NR = int(inputArr[2][22:25])
    # Title = inputArr[3]
    # NN = 0
    # for i in range (1,9):
    #     M[i]= list(map(int, inputArr[4].split('  ')))
    #     NN += M[i]
    # NT = int(inputArr[5][25:28])
    # ZPE = bool(inputArr[6][27:30])
    # TUN = bool(inputArr[7][31:33])
    # Title = inputArr[8]
    # RES = "bar"
    # ASK = bool(inputArr[9][24:27])
    # if(ASK == "yes"):
    #     RES = "atm"
    # ASK = bool(inputArr[10][17:19])
    # if(ASK == "yes"):
    #     RES = "bar"
    
    # #CONSTANTS
    # FAC = 0.46499834
    # R = 1.9872
    # RCC1 = 2.08367541593
    # if(RES == "atm"):
    #     RCC2= 1.3626055
    # if(RES=="bar"):
    #     RCC2=1.38066
    # RCC = RCC1 * numpy.pow(RCC2,(NN-1))
    # C1 = 1.438767553
    
    # #Read Temperatures
    
    # for i in range(1,NT+1):
    #     TEMP[i]= list(map(int, inputArr[i+11].split('.')))
        
    # #Loop over structures
    # for L in range(1,NP+1):
    #     GEO[L],E[L],SPECIES = list(map(float, inputArr[36].split('   ')))
    #     E[L] = E[L]*627.50595E0
    #     file1.writelines(GEO(L), E(L)/627.50595, SPECIES, '\n')
    #     BE = 0.0
    #     AE = 0.0
    #     WE = 0.0
    #     WEX = 0.0
    #     if (SPECIES == "diatomic"):
    #     # READ(10,*) WT,G0,G1,EE
    #     # READ(10,*) BE,AE,WE,WEX
    #         WT,G0,G1,EE,BE,AE,WE,WEX = list(map(float, inputArr[38].split(' ')))
    #     else:
    #     # READ(10,*) WT,G0,G1,EE
    #     # READ(10,*) A,B,C,SN,SF,F,NF,NV
    #         WT,G0,G1,EE = list(map(float, inputArr[37].split('    '))) 
    #     A,B,C,SN,SF,F,NF,NV = list(map(float, inputArr[38].split('')))
    #     A = A * FAC
    #     B = B * FAC
    #     C = C * FAC
    #     F = F * FAC
    #     #IF(NV.NE.0) READ(10,*)(W(j),j=1,NV)
    #     #IF(TUN.EQ."yes") READ(10,*) RCV(L)
    #     if(NV!=0):
    #         for j in range (1,NV+1):
    #             W[j]=list(map(float, inputArr[39].split('               ')))
    #     GIBBS (NT,TEMP,L,SPECIES,WT,G0,G1,EE,BE,AE,WE,WEX,SN,SF,A,B,C,F,NF,W,NV,DGT,E[L],RES,ZPE,EZ)
    # # Ghi ra file 
    # file1 = open('bimo.out', 'w')
    # file1.writelines(L)    
    
    # #List "Delta G"
    # if(NR > 1):
    #     for J in range(1,NT+1):
    #         DGT[NR,J]= M[NR]*DGT[NR,J]
    #     EZ[NR]=M[NR]*EZ[NR]
    #     for I in range(1,NR):
    #         EZ[NR]=EZ[NR]+M[I]*EZ[I]
    #         for J in range(1,NT+1):
    #             DGT[NR,J]=DGT[NR,J]+M[I]*DGT[I,J]
    # file1.writelines()
    # for I in range(NR,NP+1):
    #     file1.writelines(GEO[I])
    # #Save Geometric 
    # for I in range(1,NT+1):
    #     #WRITE(11,102) TEMP(I) 
    #     GDAT[I] = -10.0
    #     for J in range(NR,NP+1):
    #         DD=DGT[J,I]-DGT[NR,I]
    #         GG[J]=DD
    #         if(DD>GDAT[I]):
    #             GDAT[I]=DD
    #             if(TUN=="yes"):
    #                 TV[I]=RCV(J)
    #             else: TV[I]=0
    #             EZT[I]=EZ[J]-EZ[NR]
    #         GEOM[I]=GEO[J]
    #     #Do dG/dR (Central Difference)
    #     for J in range(NR,NP+1):
    #         if(J==NR or J==NP):
    #             file1.writelines(GEO[J],GG[J])
    #         else:
    #             dGdRA = (GG[J+1]-GG[J])/(GEO[J+1]-GEO[J])
    #             dGdRB = (GG[J]-GG[J-1])/(GEO[J]-GEO[J-1])
    #             dGdR[J]=(dGdRA+dGdRB)/2
    #             file1.writelines(GEO[J],GG[J],dGdR[J])
    #             # WRITE(11,"(F15.4, F15.6, F15.6)") GEO(J), GG(J), dGdR(J)
    # #WRITE(11,"(1X)") 
    # #WRITE(11,105) 
    # #WRITE(11,106) 
    # #WRITE(11,107)'
    # for I in range(1,NT+1):
    #     file1.writelines(TEMP[I], GEOM[I],GDAT[I])
    # # DO quadratic least
    # for I in range(1,NT+1):
    #     AA[I]= TEMP[I]*TEMP[I]
    #     AA[I+NT]=TEMP[I]
    #     AA[I+(2*NT)]=1.0
    #     BB[I] = GDAT[I]
    # HOUSEFIT(AA,NT,3,BB,NT,1,X,WW,ISING)
    # #WRITE(11,"(1X)") 
	# #WRITE(11,108) 
    # #WRITE(11,109)
    # for I in range(1,NT+1):
    #     PRED = X[1]*TEMP[I]*TEMP[I]+X[2]*TEMP[I]+X[3]
    #     file1.writelines(TEMP[I],PRED)
    # X[1]=X[1]*1000.0
    # X[2]=X[2]*1000.0
    # X[3]=X[3]*1000.0 
    # #WRITE(11,"(1X)") 
    # #WRITE(11,110)  X(1), X(2), X(3)
    # #WRITE(11,"(1X)") 
    # for I in range(1,NT):
    #     #WRITE(12,"(F10.3)") TEMP(I)
    # #WRITE(12,"(1X)") 
    # #WRITE(11,"(1X)") 
    # #WRITE(11,111) 
    # #WRITE(11,"(1X)") 
    # #WRITE(11,112) RCC, NN, NN-1
    # #WRITE(11,"(1X)") 
    # #WRITE(11,113)
    # #WRITE(11,"(1X)") 
    
    # for I in range(1,NT+1):
    #     PRED = X[1]*TEMP[I]*TEMP[I]+X[2]*TEMP[I]+X[3]
    #     FRC = RCC * numpy.pow(TEMP[I],NN) * numpy.exp(-PRED/(R*TEMP[I]))
    #     RRC = RCC * numpy.pow(TEMP[I],NN) * numpy.exp(-GDAT[I]*1000.0/(R*TEMP[I]))
    #     #WRITE(11,114) TEMP(I), RRC, FRC
    #     #WRITE(12,"(E12.6)") RRC
    # if(TUN != "yes"):
    #     #CLOSE (UNIT=10)
    #     #CLOSE (UNIT=11)
    # else:
    #    #WRITE(11,"(1X)") 
    #   #WRITE(12,"(1X)") 
    #   #WRITE(11,"(49x)") 'where w(cm^-1) = negative frequency (along path))'
    #   #WRITE(11,"(48x)") 'Note: This factor only valid for small barriers)'
    #   #WRITE(11,"(47x)") 'with a width of less than about 1.8 Angstroms.)'
    #   #WRITE(11,"(1X)") 
    #   #WRITE(11,"(33x)") 'Rates with Tunneling Correction:)'
    #   #WRITE(11,"(1X)") 
    # for I in range(1,NT+1):
    #         PRED = X[1]*TEMP[I]*TEMP[I]+X[2]*TEMP[I]+X[3]
    #         FRC = RCC * numpy.pow(TEMP[I],NN) * numpy.exp(-PRED/(R*TEMP[I]))
    #         RRC = RCC * numpy.pow(TEMP[I],NN) * numpy.exp(-GDAT[I]*1000.0/(R*TEMP[I]))
    #         TF = math.pow((C1*TV[I]/TEMP[I]),2)
    #         TF = TF * (1.0 + R*TEMP[I]/(EZT[I]*1000.0))/24.0
    #         TF = TF + 1.0
    #         RRC = TF * RRC
    #         FRC = TF * FRC
    #         #WRITE(11,114) TEMP(I), RRC, FRC
    #         #WRITE(12,"(E12.6)") RRC
          
        
    
                       