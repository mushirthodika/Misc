! Program :: Time propagation (Kramers-Henneberger)
! Author  :: Mushir
! Regn No :: MS12047 

 PROGRAM TIME_PROPAGATION_KH

 IMPLICIT NONE

!------------------------PARAMETERS---------------------------------------------!

       integer (kind=4), parameter   :: n = 1001 , n1 = 1000
       real    (kind=8), parameter   :: zero = 0.00000000000000000000000000d0
       real    (kind=8), parameter   :: pi = 3.141592653589793238462d0
       real    (kind=8), parameter   :: E = 1.00000000d0
       real    (kind=8), parameter   :: xmin = -10.000000d0,xmax = 10.00000d0
       real    (kind=8), parameter   :: tmin = zero ,tmax = 10.00000000d0
       real    (kind=8), parameter   :: omega = 0.5000000000d0

!------------------------VARIABLES----------------------------------------------!

          integer                    :: i,j,i1
          real (kind=8)              :: dx,b,a,c2,x,w,k,c1,c,gn,dt,dw,alpha
          real (kind=8)              :: h(n,n),t(n,n),evec(n,n),eval(n),v(n,n)
          real (kind=8)              :: evec_t(n,n),evec_t_d(n,n),p(n,1),q(n,1)
          real (kind=8)              :: r(n,1),s(n,1),u(n),z(n,1),evec_1(n,1)
          real (kind=8)              :: eval_t(n),eval_n(n,n),vmat(n,n),tau
          complex                    :: iota

!-------------------------------------------------------------------------------!

 iota     = (0.000d0,1)
 dx       = (xmax - xmin)/real(n - 1)
 dt       = tmax/real(n1 - 1)
 c        = (0.5000000000d0/(dx*dx))
 c1       = ((pi*pi)/3.0000d0)*c
 c2       = ((n+1)*0.50000000000d0)
 t        = zero
 h        = zero
 v        = zero 
 p        = zero 
 q        = zero 
 r        = zero 
 s        = zero 
 evec     = zero 
 evec_t   = zero 
 evec_t_d = zero 
 eval     = zero 
 eval_n   = zero 
 eval_t   = zero 
 u        = zero 
 z        = zero 
 vmat     = zero 
 evec_1   = zero
 do i=1,n

          x      = xmin + real(i-1)*dx
          v(i,i) = gn(x,zero,zero)
          h(i,i) = c1 + v(i,i)
          w      = i - c2

 do j=i+1,n

          k = w + (j-i)

 if (mod((w-k),2.0000d0)==0) then

          h(i,j) = c*(2.00d0/((w - k)*(w - k)))
 else
          h(i,j) = c*(2.00d0/((w - k)*(w - k)))*(-1)

 endif
          h(j,i)=h(i,j)

 enddo

 enddo

 
 t = h

 do i = 1, n

           t(i,i) = h(i,i) - v(i,i) 
 enddo

!stop

!-------------------MATRIX DIAGONALISATION--------------------------------------!

 call CALL_DSYEV_WT_EIGVECS(h,n,eval,evec) ! Diagonalisation of Hamiltonian


 call CALL_DSYEV_WT_EIGVECS(t,n,eval_t,evec_t) ! Diagonalisation of KE Matrix

!-------------------------------------------------------------------------------!

 do i = 1,n                      ! Conversion of 1D array to a diagonal matrix
                                 ! form 
     eval_n(i,i) = eval_t(i)     
                                
 enddo                         


 evec_t_d = transpose(evec_t)    ! Transpose of the eigenvector matrix 
                                 ! (Obtained after diagonalising the KE MATRIX)

 open(12,file='tpkhwavefuncs.dat')
 open(13,file='tdkhpot.dat')

!----------------TIME PROPAGATION-----------------------------------------------!

 x = zero

 do i = 1,n1
            tau = dt*real(i)
            alpha = (E/(omega*omega))*sin(omega*tau)

      do j = 1,n

                         x = xmin + real(j-1)*dx

               evec_1(j,1) = evec(j,1)
     
                      u(j) = gn(x,tau,alpha)               ! Calculation of the  
                                                           ! exponential of a 
                 vmat(j,j) = exp(-iota*dt*u(j)*0.500000d0) ! matrix

               eval_n(j,j) = exp(-iota*dt*eval_n(j,j))    

               write(13,*) u(j)
 
      enddo

            write(13,*)
            write(13,*)
            write(13,*)

            p = matmul(vmat,evec_1)

            q = matmul(evec_t_d,p)

            r = matmul(eval_n,q)

            s = matmul(evec_t,r)

            z = matmul(vmat,s)

     do i1 = 1,n

               write(12,*) z(i1,1)

     enddo
          
           write(12,*)
           write(12,*)
           write(12,*)
!           write(*,*) i
         
              evec_1 = z
 enddo

 close (13)
 close (12) 

 END

!----------------HCN-HNC ASYMMETRIC DOUBLE WELL POTENTIAL-----------------------!
 FUNCTION gn(x,t,alpha)
!-------------------------------------------------------------------------------!

 IMPLICIT NONE

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

 END

!--------------DOUBLE WELL POTENTIAL (FEIT,1982)--------------------------------!
 FUNCTION dw(x)
!-------------------------------------------------------------------------------!

 IMPLICIT NONE

!-------------VARIABLE DECLARATION----------------------------------------------!

 real (kind=8), parameter :: k0 = -132.7074997d0, k2 = 7.0000000000d0
 real (kind=8), parameter :: k3 = 0.500000000d0, k4 = 1.0000000000d0
 real (kind=8)            :: x, dw

!-------------------------------------------------------------------------------!

 dw = k0 - (k2*x*x) + (k3*x*x*x) + (k4*x*x*x*x)

END

