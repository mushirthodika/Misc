    Program main
!=============================================
! Integration of a function using Simpson rule 
!=============================================
implicit none
double precision f, a, b, integral
integer n, i
double precision, parameter:: pi = 6.28318530717958647692d0
external f

a = 0.0
b = pi

n = 101

write(*,100)

do i=1,16
   call simpson(f,a,b,integral,n)
   write (*,101) n, integral
   n = n*2
end do

100   format('     nint   Simpson')
101   format(i9,1pe15.6)
end
!----------------HCN-HNC ASYMMETRIC DOUBLE WELL POTENTIAL-----------------------!
 FUNCTION gn(x,t,alpha)

 Implicit None

!----------------------VARIABLE DECLARATION-------------------------------------!
real (kind=8)            :: gn,x,t,alpha

real (kind=8), parameter :: a1=0.076060d0,b1=-0.73950d0,c1=1.9920d0,a2=-0.02288d0

real (kind=8), parameter :: b2=2.90400d0,c2=1.31100d0,a3=-92.9500d0,b3=-1.78000d0

real (kind=8), parameter :: c3=106.400000d0,a4=0.021070000000d0,b4=-0.012320000d0

real (kind=8), parameter :: c4=1.054000000d0,a5=-0.34730000000d0,b5=5.850000000d0

real (kind=8), parameter :: c5=2.6490000d0
!-------------------------------------------------------------------------------!

 f  =  a1*dexp(-((x-(alpha*dcos(t))-b1)/c1)*((x-(alpha*dcos(t))-b1)/c1)) +&

    &  a2*dexp(-((x-(alpha*dcos(t))-b2)/c2)*((x-(alpha*dcos(t))-b2)/c2)) +&

    &  a3*dexp(-((x-(alpha*dcos(t))-b3)/c3)*((x-(alpha*dcos(t))-b3)/c3)) +&

    &  a4*dexp(-((x-(alpha*dcos(t))-b4)/c4)*((x-(alpha*dcos(t))-b4)/c4)) +&

    &  a5*dexp(-((x-(alpha*dcos(t))-b5)/c5)*((x-(alpha*dcos(t))-b5)/c5))

 END
!-------------------------------------------------------------------------------!


!  Function f(x)
!----------------------------------------
! Function for integration
!----------------------------------------
!implicit none
!double precision f, x
! f = sin(x)
! f = x*cos(10.0*x**2)/(x**2 + 1.0)
!return
!end

 Subroutine simpson(f,a,b,integral,n)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of integration
!==========================================================
implicit none
double precision f, a, b, integral,s
double precision h, x
integer nint
integer n, i

! if n is odd we add +1 to make it even
if((n/2)*2.ne.n) n=n+1

! loop over n (number of intervals)
s = 0.0
h = (b-a)/dfloat(n)
do i=2, n-2, 2
   x   = a+dfloat(i)*h
   s = s + 2.0*f(x) + 4.0*f(x+h)
end do
integral = (s + f(a) + f(b) + 4.0*f(a+h))*h/3.0
return
end subroutine simpson
