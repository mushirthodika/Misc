! Program :: Kramer-Henneberger time propagation for symmetric double well potential
 
! Author  :: Mushir Ul hasan

! Regn No :: MS12047 

!-------------------------------------------------------------------------------!
 PROGRAM TIME_PROPAGATION_KH
!-------------------------------------------------------------------------------!

 IMPLICIT NONE

!------------------------VARIABLE DECLARATION-----------------------------------!
          integer, parameter         :: n = 201 , n1 = 1000
          integer                    :: i,j,i1
          real (kind=8)              :: dx,b,a,c2,x,w,k,c1,c,gn,dt,dw
          real (kind=8)              :: h(n,n),t(n,n),evec(n,n),eval(n),v0(n,n)
          real (kind=8)              :: evec_d(n,n),mat1(n,n),mat2(n,n)
          real (kind=8)              :: mat3(n,n),mat4(n,n),mat5(n,n),mat6(n,n)
          real (kind=8)              :: eval_(n),d_eval(n,n),tau
          real (kind=8)              :: eval_1(n),eval_2(n),eval_3(n)
          real (kind=8)              :: eval_a(n,n),eval_b(n,n),eval_c(n,n)
          real (kind=8)              :: evec_1(n,n),evec_1d(n,n),evec_2(n,n)
          real (kind=8)              :: evec_2d(n,n),evec_3(n,n),evec_3d(n,n)
          real (kind=8)              :: a_mat(n,n),b_mat(n,n),c_mat(n,n)
          real (kind=8)              :: d_mat(n,n),e_mat(n,n),f_mat(n,n)
          real (kind=8)              :: g_mat(n,n),h_mat(n,n),i_mat(n,n)
          real (kind=8)              :: j_mat(n,n),k_mat(n,n),l_mat(n,n)
          real (kind=8)              :: m_mat(n,n),n_mat(n,n),o_mat(n,n)
          real (kind=8)              :: p_mat(n,n),q_mat(n,n),r_mat(n,n)
          real (kind=8)              :: s_mat(n,n),t_mat(n,n),z_mat(n,n)
          real (kind=8)              :: v1(n,n),v2(n,n),v3(n,n)
          real (kind=8), parameter   :: pi=3.141592653589793238462d0
          real (kind=8), parameter   :: xmin=-10.000000d0,xmax=10.0000000d0
          complex                    :: iota
          real (kind=8), parameter   :: tmin=0.0000000d0,tmax=10.00000000d0
          real (kind=8), parameter   :: zero=0.000000000000000000000000000000d0
          real (kind=8), parameter   :: alpha=1.15600000000000000000000000000d0
          real (kind=8), parameter   :: omega=0.01000000000000000000000000d0
!-------------------------------------------------------------------------------!

 iota     = (0.000d0,1)
 dx       = (xmax - xmin)/real(n - 1)
 dt       = tmax/real(n1-1)
 c        = (0.5000000000d0/(dx*dx))
 c1       = ((pi*pi)/3.0000d0)*c
 c2       = ((n+1)*0.50000000000d0)
 t        = zero
 h        = zero
 v0       = zero 
 v1       = zero 
 v2       = zero 
 v3       = zero 
 evec     = zero 
 evec_d   = zero 
 eval     = zero 
 d_eval   = zero 
 eval_1=zero;eval_2=zero;eval_3=zero;eval_a=zero;eval_b=zero
 eval_c=zero;evec_1=zero;evec_2=zero;evec_3=zero;evec_1d=zero
 evec_2d=zero;evec_3d=zero
 mat1 = zero;mat2=zero;mat3=zero;mat4=zero;mat5=zero;mat6=zero
 a_mat=zero;b_mat=zero;c_mat=zero;d_mat=zero;e_mat=zero;f_mat=zero
 g_mat=zero;h_mat=zero;i_mat=zero;j_mat=zero;k_mat=zero;l_mat=zero
 m_mat=zero;n_mat=zero;o_mat=zero;p_mat=zero;q_mat=zero;r_mat=zero
 s_mat=zero;t_mat=zero;z_mat=zero

 do i=1,n

          x       = xmin+real(i-1)*dx

          v0(i,i) = ((x**4) + (((3.0000d0*alpha*alpha) - 1)*x*x) -& 
                   & (0.500000d0*alpha*alpha) + (0.375000000000d0*(alpha**4))) 

          v1(i,i) = ((-2.00000d0*x*x*x*alpha) - (3.0000d0*0.50000d0*&
                    & x*(alpha**3)) + (x*alpha))

          v2(i,i) = (((alpha**4)*0.250000d0) + (1.50000d0*x*x*alpha*alpha) -&
                    &  (alpha*alpha*0.5000000d0))

          v3(i,i) = (-x*(alpha**3)*0.500000d0)

          h(i,i)  = c1 + v0(i,i)

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
 
!-------------------MATRIX DIAGONALISATION--------------------------------------!

 call CALL_DSYEV_WT_EIGVECS(h,n,eval,evec) ! Diagonalisation of Hamiltonian

 evec_d = transpose(evec)                  ! Transpose of the eigenvector 
                                           ! matrix for the diagonalized
 mat1 = matmul(v1,evec)                    ! Hamiltonian
 mat2 = matmul(evec_d,mat1)

 call CALL_DSYEV_WT_EIGVECS(mat2,n,eval_1,evec_1)

 mat3 = matmul(v2,evec)
 mat4 = matmul(evec_d,mat3)

 call CALL_DSYEV_WT_EIGVECS(mat4,n,eval_2,evec_2)

 mat5 = matmul(v3,evec)
 mat6 = matmul(evec_d,mat5)

 call CALL_DSYEV_WT_EIGVECS(mat6,n,eval_3,evec_3)

 evec_1d = transpose(evec_1)

 evec_2d = transpose(evec_2)

 evec_3d = transpose(evec_3)

 do i = 1,n
 
          d_eval(i,i) = eval(i)  
          eval_a(i,i) = eval_1(i)
          eval_b(i,i) = eval_2(i)
          eval_c(i,i) = eval_3(i)
 
 enddo

 
 open(12,file='khtwavefuncs.dat')

!----------------TIME PROPAGATION-----------------------------------------------!

 do i = 1,n1

       tau = dt*real(i)

      do j = 1,n
               
       d_eval(j,j) = exp(-iota*dt*d_eval(j,j))       
       eval_a(j,j) = exp(-iota*dt*cos(omega*tau)*0.500000d0*eval_a(j,j)) 
       eval_b(j,j) = exp(-iota*dt*cos(2.000d0*omega*tau)*0.500000d0*eval_b(j,j)) 
       eval_c(j,j) = exp(-iota*dt*cos(3.000d0*omega*tau)*0.500000d0*eval_c(j,j)) 

      enddo

       a_mat = matmul(evec_3,evec) ; b_mat = matmul(eval_c,a_mat)

       c_mat = matmul(evec_3d,b_mat) ; d_mat = matmul(evec_2,c_mat)

       e_mat = matmul(eval_b,d_mat) ; f_mat = matmul(evec_2d,e_mat)

       g_mat = matmul(evec_1,f_mat) ; h_mat = matmul(eval_a,g_mat)
 
       i_mat = matmul(evec_1d,h_mat) ; j_mat = matmul(evec,i_mat)

       k_mat = matmul(d_eval,j_mat) ; l_mat = matmul(evec_d,k_mat)

       m_mat = matmul(evec_1,l_mat) ; n_mat = matmul(eval_a,m_mat)

       o_mat = matmul(evec_1d,n_mat) ; p_mat = matmul(evec_2,o_mat)      

       q_mat = matmul(eval_b,p_mat) ; r_mat = matmul(evec_2d,q_mat)      

       s_mat = matmul(evec_3,r_mat) ; t_mat = matmul(eval_c,s_mat)      

       z_mat = matmul(evec_3d,t_mat)      

     do i1 = 1,n

               write(12,*) (z_mat(i1,j),j=1,n)

     enddo
          
           write(12,*)
           write(12,*)
           write(12,*)

           write(*,*) i
         
              evec = z_mat
 enddo

 close (12) 

 END

