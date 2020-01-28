! This program calculates the ground state energy for the Helium atom using Hartree-Fock method
! Author :: Mushir Ul Hasan
! Date   :: 09/09/2016 (DD/MM/YYYY)

                          !****************************************************************!
                           !***Check for stop command and comment it before compilation***!
                          !****************************************************************!      


     !------------------------------------------!
      MODULE INPUT
     !------------------------------------------!

      IMPLICIT NONE

     !-------------------------------------------------------------------------!
      real (kind=8) , parameter   :: zero = 0.00000000000000000000000000000d0
      real (kind=8) , parameter   :: one = 1.00000000000000000000000000000d0
      real (kind=8) , parameter   :: pi = 3.1415926535897932384626433830000d0
      real (kind=8) , parameter   :: thby2 = 1.50000000000000000000000000000d0
      real (kind=8) , parameter   :: two = 2.000000000000000000000000000000d0
      real (kind=8) , parameter   :: fby2 = 2.50000000000000000000000000000d0
      real (kind=8) , parameter   :: three = 3.0000000000000000000000000000d0
      real (kind=8) , parameter   :: four = 4.00000000000000000000000000000d0
      real (kind=8) , parameter   :: thby4 = 0.750000000000000000000000000d0
     !-------------------------------------------------------------------------!

      END
     
   !-----------------------------------------!
    PROGRAM HELIUM_HF
   !-----------------------------------------!

       USE INPUT

      IMPLICIT NONE

     !-----------------PARAMETER DECLARATION------------------------------------!
      integer (kind=8) , parameter :: n=4,Niter=100
      real    (kind=8) , parameter :: eps=1.00e-08
     !--------------------------------------------------------------------------!

     !----------------ARRAY DECLARATAION----------------------------------------!
      real (kind=8), allocatable  :: TMAT(:,:), VMAT(:,:), SMAT(:,:), FMAT(:,:)
      real (kind=8), allocatable  :: QMAT(:,:,:,:),CMAT(:,:),ALPHA(:),EPRIME(:)
      real (kind=8), allocatable  :: CMAT_NEW(:,:), HMAT(:,:), XMAT(:,:)
      real (kind=8), allocatable  :: SDIAG(:), FPRIME(:,:),FPRIME1(:,:)
      real (kind=8), allocatable  :: FDIAG(:),XMAT_D(:,:),UMAT(:,:),UMAT_D(:,:)
      real (kind=8), allocatable  :: LDA(:,:), XMAT1(:,:)
      real (kind=8), allocatable  :: FMAT1(:,:),FMAT2(:,:),EPRIMEOLD(:)
     !--------------------------------------------------------------------------!
 
     !---------------VARIABLE DECLARATION---------------------------------------!
      integer (kind=4)  :: i,j,k,l,i1,i2,i3,j3
      real    (kind=8)  :: NORM,E,EOLD
     !--------------------------------------------------------------------------!

    !---------------------ARRAY ALLOCATION--------------------------------------!

      allocate(TMAT(n,n),VMAT(n,n),SMAT(n,n),FMAT(n,n),EPRIME(n),QMAT(n,n,n,n))
      allocate(ALPHA(n),HMAT(n,n),XMAT(n,n),FPRIME(n,n),FPRIME1(n,n),FDIAG(n))
      allocate(SDIAG(n),XMAT_D(n,n),XMAT1(n,n),LDA(n,n),UMAT(n,n),UMAT_D(n,n))
      allocate(CMAT_NEW(n,n),CMAT(n,n),EPRIMEOLD(n),FMAT1(n,n),FMAT2(n,n))

    !------------------INITIALIZATION OF ARRAYS---------------------------------!

      TMAT = zero; VMAT = zero; SMAT = zero; FMAT = zero; EPRIME = zero
     
      CMAT = zero; QMAT = zero; FPRIME = zero; FPRIME1 = zero; SDIAG = zero 

      HMAT = zero; XMAT = zero; FDIAG = zero; XMAT_D = zero; XMAT1 = zero

      LDA = zero; UMAT = zero; UMAT_D = zero; E = zero

      FMAT1 = zero; FMAT2 = zero
   
   !-------------------EXPONENT VALUES------------------------------------------!

      ALPHA(1) = 0.2980730000000000000d0
      ALPHA(2) = 1.2425670000000000000d0
      ALPHA(3) = 5.7829480000000000000d0
      ALPHA(4)  = 38.47497000000000000d0

   !***************************************************************************!
      ! Formation of the kinetic energy matrix, potential matrix
      ! and the overlap matrix
   !***************************************************************************!

      do 101 i = 1, n

             do 102 j = 1,n

                     TMAT(i,j) = ((three*ALPHA(i)*ALPHA(j)*((pi)**(thby2)))/&
                               & ((ALPHA(i) + ALPHA(j))**(fby2)))*&
                               & (((two*alpha(i))/pi)**thby4)*&
                               & (((two*alpha(j))/pi)**thby4)

                     VMAT(i,j) = -((four*pi)/(ALPHA(i) + ALPHA(j)))*&
                                  & (((two*alpha(i))/pi)**thby4)*&
                                  & (((two*alpha(j))/pi)**thby4)


                     SMAT(i,j) = (((pi)/(ALPHA(i) + ALPHA(j)))**thby2)*&
                                 & (((two*alpha(i))/pi)**thby4)*&
                                 & (((two*alpha(j))/pi)**thby4)

             102 enddo

      101 enddo


   !***************************************************************************!
    ! Formation of the Q matrix (Two - electron integrals)
   !***************************************************************************!
 
     do 103 i = 1,n

         do 104 j = 1,n

            do 105 k = 1,n

               do 106 l = 1,n


                QMAT(i,k,j,l) =  ((two*(pi**fby2))/((ALPHA(i) +&
                               & ALPHA(j))*(ALPHA(k) + ALPHA(l)) *&
                               & (sqrt(ALPHA(i)+ALPHA(j)+ALPHA(k)+ALPHA(l)))))*&
                               & (((two*alpha(i))/pi)**thby4)*&
                               & (((two*alpha(j))/pi)**thby4)*&
                               & (((two*alpha(k))/pi)**thby4)*&
                               & (((two*alpha(l))/pi)**thby4)

               106 enddo

            105 enddo

         104 enddo

     103 enddo

   !-----------ORTHOGONALISATION OF BASIS---------------------------------------!

     HMAT = TMAT + VMAT

     call CALL_DSYEV_WT_EIGVECS(SMAT,n,SDIAG,UMAT)

     UMAT_D = transpose(UMAT) 

     do i = 1, n

     LDA(i,i) = one/sqrt(SDIAG(i))

     enddo

     XMAT1 = matmul(LDA,UMAT_D)
  
     XMAT = matmul(UMAT,XMAT1)

     XMAT_D = transpose(XMAT)

     FPRIME1 = matmul(HMAT,XMAT)
 
     FPRIME = matmul(XMAT_D,FPRIME1)

     call CALL_DSYEV_WT_EIGVECS(FPRIME,n,FDIAG,CMAT)

     CMAT = matmul(XMAT,CMAT)
 
     write(*,*)

   !-------------SCF PROCEDURE--------------------------------------------------!
     
 222   do 107 i = 1,n

          do 108 j = 1,n

                         FMAT(i,j) = HMAT(i,j) 

               do 109 k = 1,n

                    do 110 l = 1,n

                         FMAT(i,j) = FMAT(i,j) + ((two*QMAT(i,k,j,l))-&
                                   QMAT(i,j,k,l))*CMAT(k,1)*CMAT(l,1)

                    110 enddo

               109 enddo

          108 enddo
     
    107 enddo

     FMAT1 = matmul(FMAT,XMAT)

     FMAT2 = matmul(XMAT_D,FMAT1)

     EPRIMEOLD = EPRIME
 
     EPRIME = zero
 
     call CALL_DSYEV_WT_EIGVECS(FMAT2,n,EPRIME,CMAT_NEW)

     CMAT = matmul(XMAT,CMAT_NEW)

   !-----------------CONVERGENCE CHECK------------------------------------------!
 
     if ((abs((EPRIME(1))-(EPRIMEOLD(1)))).lt.eps) then

      E = zero
     
     do i = 1,n

          do j = 1,n

               E  = E + two*CMAT(i,1)*CMAT(j,1)*HMAT(i,j)

               do k = 1,n

                    do l = 1,n

                         E = E + QMAT(i,k,j,l)*&
                           & CMAT(i,1)*CMAT(j,1) *&
                           & CMAT(k,1)*CMAT(l,1)

                    enddo

               enddo

          enddo

     enddo

     write(*,*) 
     write(*,*) "------The eigenvalues are as follows:------"
     write(*,*)
 
     do i2 = 1, n 

           write(*,'(F10.4)') EPRIME(i2)

     enddo

     write(*,*)
     write(*,*) "--------------------------------------------------------------------------------"
     write(*,'(A50,F15.10,A10)') "The ground state electronic energy of Helium is",E,"Hartrees"
     write(*,*) "--------------------------------------------------------------------------------"

     write(*,*) 
     write(*,*) "-------The Eigenvectors are as follows:------" 
     write(*,*)

     do i3 = 1, n

     write(*,'(4F10.6)') (CMAT(i3,j3),j3=1,n)

     enddo
     write(*,*)

     else

     goto 222
  
     endif

     END
   !----------------------------------------------------------------------------!

