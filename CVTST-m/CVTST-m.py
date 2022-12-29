import math
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
def GIBBS (NT,TEMP,NN,SPECIES,WT,G0,G1,EE,BE,AE,WE,WEX,SN,SF,A,B,C,F,NF,W,NV,DGT,E,RES,ZPE,EZ):
    C1=1.438767553
    C2=40.27512295
    R=1.9872
    PI=numpy.arccos(-1.0)
    if (RES == "atm"): TC = -3.664854875
    if (RES == "bar") : TC = -3.651691889
    II = 1
    while (II <= NT):
        T = TEMP(II)
        if (EE == 0): CELECT = numpy.log(G0)
        if (EE != 0): CELECT = numpy.log(G0 + G1*numpy.exp(C1*EE/T))
        
        CTRANS = 1.5*numpy.log(WT) + 2.5*numpy.log(T) + TC
        
        if(SPECIES == "monatomic "): G = (CELECT + CTRANS)*R*T
        if(SPECIES == "monatomic "): E0 = 0
        
        if(SPECIES == "diatomic  "):
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
        E0   = 0
        CVIB = 0
        # dang code den dong 350 o fortran
#end gibbs
def main():
    print('main')

if __name__ == "__main__":
    main()