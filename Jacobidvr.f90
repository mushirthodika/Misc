 ! Author :: Mushir Ul Hasan
 ! Program :: Discrete Variable Representation (2D)
 ! Registration No :: MS12047
 ! Version :: 0.1
 Program dvr
 implicit none

!--------------------VARIABLE DECLARATION---------------------------------------!
 integer (kind=4), parameter :: n=4
 real (kind=8),parameter     :: pi=3.141592653589793238462d0 !, eps=1.00000e-07
 real (kind=8), parameter    :: tmin=0.000000000d0,s1=0.250000000000000000d0 
 real (kind=8), parameter    :: a=-8.0000000d0, b=8.000000000000d0
 real (kind=8), parameter    :: ymin=-4.0000000d0, ymax=4.00000000000d0
 real (kind=8)               :: h(n*n,n*n),ht(n*n,n*n),t(n*n,n*n)
 real (kind=8)               :: ev(n*n,n*n),e(n*n),v(n*n,n*n)
 real (kind=8)               :: dx,s,x,w,k,c,xmin,tx,y,tdho,dy
 integer (kind=4)            :: i1,i2,j1,j2,k1,k2,l1,l2,i,j
!-------------------------------------------------------------------------------!

!write(*,*) "Enter the value of a & b"
!read(*,*) a,b

 dx = (b-a)/real(n)
 dy = (ymax-ymin/real(n))
 s  = (0.50000000000000000d0/((dx*dx)))
 c  = ((pi*pi)/3.0000000000d0)
 tx = (pi-tmin)/real(n)
 h  = 0.00000000000d0
 ht = 0.0000000000d0
 e  = 0.0000000000d0
 ev = 0.00000000000d0
!----------------------HAMILTONIAN CONSTRUCTION---------------------------------!
 do i1=0,n-1
 do j1=1,n
 
         h((n*i1)+j1,(n*i1)+j1) = (c - (0.5000000000d0/real(j1*j1)))*s
 
 do j2=j1+1,n
 
 if (mod(abs(j1-j2),2).eq.0) then
 
  h((n*i1)+j1,(n*i1)+j2) = s *&
 
 & (2.0000000d0*((1.0000000d0/real((j1-j2)*(j1-j2))) -&
 
 & (1.0000000d0/real((j1+j2)*(j1+j2)))))
 
 else
 
  h((n*i1)+j1,(n*i1)+j2) = s *&
 
 & (2.000000d0*((1.000000d0/real((j1-j2)*(j1-j2))) -&
 & (1.0000d0/real((j1+j2)*(j1+j2)))))*(-1)
         
 endif
 
            h((n*i1)+j2,(n*i1)+j1) = h((n*i1)+j1,(n*i1)+j2)
 enddo
 enddo
 enddo
 
 do i2=0,n-1
 
 do k1=1,n-1
 
 
 ht((n*i2)+k1,(n*i2)+k1) = (1.00000d0/(k1*dx))* s1 *&

 & ( (((2.0000d0*n*n) + 1.000d0)/3.00000d0) -&

 & (1.00000000d0/(dsin((pi*k1)/real(n))*(dsin((pi*k1)/real(n))))))

 do k2=k1+1,n-1
 
 
 if (mod((k1-k2),2)==0) then

 ht((n*i2)+k1,(n*i2)+k2) = (1.00000d0/(k1*dx)) *&

 & ((1.00d0/(dsin(0.5000d0*pi*(k1-k2)/real(n)))) *&
 
 & dsin(0.50000000d0*pi*(k1-k2)/(real(n)))) -&
 
 & (1.0000000d0/(dsin(0.5000000d0*pi*(k1-k2)/(real(n))) *&
 
 & dsin(0.500000000d0*pi*(k1-k2)/(real(n)))))*s1
 
 else 

 ht((n*i2)+k1,(n*i2)+k2) = (1.00000d0/(k1*dx))*s1 *&

 & ((1.000000d0/(dsin(0.5000000d0*pi*(k1-k2)/(real(n))) *&
 
 & dsin(0.5000000d0*pi*(k1-k2)/(real(n))))) -&
 
 & (1.0000000d0/(dsin(0.50000000d0*pi*(k1-k2)/(real(n))) *&
 
 & dsin(0.500000000d0*pi*(k1-k2)/(real(n))))))*(-1)
 
 endif
        ht((n*i2)+k2,(n*i2)+k1) = ht((n*i2)+k1,(n*i2)+k2)
 
 enddo
 enddo
 enddo
 open(12,file='thetamat.dat')
 do i=1,n*n
 write(12,*) (ht(i,j),j=1,n*n)
 enddo
 
 open(111,file='potval.dat')
 do i=1,n*n
 read(111,*) (v(i,j),j=1,n*n)
 enddo
 close(111)
 t = h + ht + v
!-------------------PRINTS THE HAMILTONIAN--------------------------------------
 open(13,file='hmat.dat')
 do i=1,n*n
 write(13,*) (t(i,j),j=1,n*n)
 enddo
!-------------------------------------------------------------------------------

 write(*,*) "Hamiltonian done"

!------------------MATRIX DIAGONALISATION---------------------------------------
 call CALL_DSYEV_WT_EIGVECS(h,n,e,ev)

 open(14,file='eigenvalues.dat')
 open(15,file='eigenvectors.dat')
 do i = 1,n*n
 write(14,*) e
 write(15,*) (ev(i,j),j=1,n)
 enddo
 close(14)
 close(15)
!-------------------------------------------------------------------------------

 END

!------------POTENTIAL FUNCTION-------------------------------------------------
 Function tdho(x,y)

 implicit none

 real (kind=8) :: x,y,tdho
!real (kind=8), parameter :: =1.0000000000d0

 tdho = 0.5000000000*((x*x) + (y*y))
 
 END
