!
! Input  : N,A,B
! Output : X
!

!=========================================================================
!
! Title         : MatrixSolve_mod
! Application   : Solve linear equations
! Copyright     : Tata Consultancy Services
! Author        : Chetan Malhotra
! Creation Date : 
! Requirements  : None
! Description   : This file contains the module MatrixSolve which solves 
!				  a set of linear algebraic equations by LU 
!				  decomposition. The LU decomposition solver is taken from
!				  from Numerical Recipes Software.
! Limitations   : None
! Dependencies  : Fortran 90
! Modifications :       WHO           WHEN         WHY
!
!=========================================================================

module MatrixSolve_mod

  contains

    subroutine matrixsolver(N,A,B,X)

	implicit none

	integer N
	double precision A(N,N),B(N),X(N)
	integer dummy1(N)
	double precision A_bak(N,N),B_bak(N),dummy2
	
    A_bak=A
	B_bak=B
	call ludcmp(A,N,N,dummy1,dummy2)
	call lubksb(A,N,N,dummy1,B)
	call mprove(A_bak,A,N,N,dummy1,B_bak,B)
	X=B
	A=A_bak
	B=B_bak

	return

	end subroutine matrixsolver

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END SUBROUTINE ludcmp
!  (C) Copr. 1986-92 Numerical Recipes Software "0'!zQ-.
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END SUBROUTINE lubksb
!  (C) Copr. 1986-92 Numerical Recipes Software "0'!zQ-.
      SUBROUTINE mprove(a,alud,n,np,indx,b,x)
      INTEGER n,np,indx(n),NMAX
      REAL*8 a(np,np),alud(np,np),b(n),x(n)
      PARAMETER (NMAX=500)
!    USES lubksb
      INTEGER i,j
      REAL*8 r(NMAX)
      DOUBLE PRECISION sdp
      do 12 i=1,n
        sdp=-b(i)
        do 11 j=1,n
          sdp=sdp+dble(a(i,j))*dble(x(j))
11      continue
        r(i)=sdp
12    continue
      call lubksb(alud,n,np,indx,r)
      do 13 i=1,n
        x(i)=x(i)-r(i)
13    continue
      return
      END SUBROUTINE mprove
!  (C) Copr. 1986-92 Numerical Recipes Software "0'!zQ-.
end module MatrixSolve_mod