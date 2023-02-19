      SUBROUTINE HOUSETRANS(A,M,N,B,KB,LB,W,ISING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(M,N) , B(KB,LB) , W(M)
      ISING=0
      DO 300 K=1,N
         MX=IDAMAX(M-K+1,A(K,K),1) + K-1
  300 CONTINUE
      RETURN
      END
            integer function idamax(n,dx,incx)
      implicit double precision (a-h,o-z)
      dimension dx(*)
         ii=1
         xmax=abs(dx(1)) 
      if(incx.eq.1) then
         do 100 i=2,n
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