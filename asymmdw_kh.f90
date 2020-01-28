! This program calculates the time averaged KH potential for a given potential
! Author :: Mushir Ul Hasan
 
 PROGRAM timeavgKH

 IMPLICIT NONE

!------------------VARIABLE DECLARATION-----------------------------------------!
          integer (kind=8), parameter  :: n=1001, n1=500
          integer (kind=8)             :: i,j,k
          real (kind=8)                :: d,x,am,fn,gn,l,d0,t
          real (kind=8), parameter     :: alpha = 6.700000000000000d0
          real (kind=8)                :: w(n1),p(n1),v_mat  
          real (kind=8), parameter     :: a=-7.0000000000d0,b=7.00000000000d0
          real (kind=8), parameter     :: pi = 3.14159265358979323846260d0
!-------------------------------------------------------------------------------!
 
 d  = (b-a)/real(n - 1)
 d0 = 0.200000000d0
 open(6,file='tau_w.txt')
 open(7,file='tau_x.txt')

 do i=1,n1
           read(6,*) w(i)
           read(7,*) p(i)
 enddo

 close(6)
 close(7)

 x = 0.00000000d0
 v_mat = 0.00000000000d0

 open(5,file='zero_order.dat')
 do i = 1,n
           x = a + real(i - 1)*d
           fn=0.000000000d0
 do j = 1,n1
          fn = fn + (gn(x,p(j),alpha)*w(j))
 enddo
       ! write(5,*) x, fn/(2.00000000d0*pi)
       v_mat = fn/(2.0000000000d0*pi)
       write(5,*) x,v_mat
 enddo
 close(5)

 END

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

 gn =  a1*dexp(-((x-(alpha*dcos(t))-b1)/c1)*((x-(alpha*dcos(t))-b1)/c1)) +&

    &  a2*dexp(-((x-(alpha*dcos(t))-b2)/c2)*((x-(alpha*dcos(t))-b2)/c2)) +&

    &  a3*dexp(-((x-(alpha*dcos(t))-b3)/c3)*((x-(alpha*dcos(t))-b3)/c3)) +&

    &  a4*dexp(-((x-(alpha*dcos(t))-b4)/c4)*((x-(alpha*dcos(t))-b4)/c4)) +&

    &  a5*dexp(-((x-(alpha*dcos(t))-b5)/c5)*((x-(alpha*dcos(t))-b5)/c5))

 gn = (gn + 92.839470338523768) !* dcos(8.00d0*t)

 END

!-------------------AMMONIA DOUBLE WELL POTENTIAL-------------------------------!
 FUNCTION am(x,t,alpha)

 Implicit none

 real (kind=8) :: x,t,alpha,am
 real (kind=8), parameter :: A = 5.224852D-03, B = 1.39753572D-02

 am = ((x-(alpha*dcos(t)))**4) - (x -(alpha*dcos(t)))**2

 am = am*dcos(t)

 END
