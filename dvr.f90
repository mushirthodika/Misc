! DISCRETE VARIABLE REPRESENTATION IN ONE DIMENSION (DVR)
! AUTHOR :: MUSHIR UL HASAN


!***************************************************************************!
! COMMENT THE STOP COMMAND BEFORE COMPILATION
!***************************************************************************!


     PROGRAM DISCRETE_VAR_REP_1D

     IMPLICIT NONE

!-----------------VARIABLE DECLARATION--------------------------------------!
     INTEGER*4, PARAMETER        :: n=301
     REAL*8                      :: h(n,n),ev(n,n),e(n)
     INTEGER*4                   :: i,j
     REAL*8                      :: dx,s,x,v,w,k,c1,c,potential
     REAL*8, PARAMETER           :: pi=3.141592653589793238462d0
     REAL*8, PARAMETER           :: xmin=-20.000000d0,xmax=20.0000000000d0
!---------------------------------------------------------------------------!

     dx=(xmax-xmin)/n
     s=(0.50000000/(dx*dx))
     c=((pi*pi)/3.0000d0)*s
     c1=((n+1)*0.50000000000d0)

     OPEN(18,FILE='potential.dat')

     DO i = 1, n

          x=xmin+real(i-1)*dx

          v=potential(x)

          h(i,i) = c + v 

          w = i - c1

          DO j = i+1, n

              k = w + (j-i)

              IF (mod((w-k),2.0000d0)==0) THEN

                  h(i,j) = s*(2.00d0/((w - k)*(w - k)))

                  ELSE

                  h(i,j) = s*(2.00d0/((w - k)*(w - k)))*(-1)

                  ENDIF

              h(j,i)=h(i,j)

          ENDDO

          write(18,*) x,v

     ENDDO

     OPEN(13,FILE='KE.txt')

     DO i = 1, n

         WRITE(13,*) (h(i,j),j=1,n)

     ENDDO

     CLOSE(18)

!    STOP

!-------------------------------------------------------------------------------
     call CALL_DSYEV_WT_EIGVECS(h,n,e,ev) !***MATRIX DIAGONALIZATION***!
!-------------------------------------------------------------------------------

     OPEN(14,file='eval.dat')
     OPEN(15,file='evec.dat')
     OPEN(18,file='potential.dat')

     DO i = 1, n

        READ(18,*) x,v
        WRITE(14,*) e(i)
        WRITE(15,*) x,v,((3.000d0*ev(i,j)*ev(i,j))+e(j),j=1,n)

     ENDDO

     CLOSE(14)
     CLOSE(15)
     CLOSE(18)

     END

!--------------------------------------------------------------!
     FUNCTION potential(x)
!--------------------------------------------------------------!

     IMPLICIT NONE

     REAL*8          :: potential,x

     potential = ((x*x*0.50000d0)-0.800000d0)*dexp(-0.100000*x*x)
 
     !potential = 0.50000d0*x*x

     !potential = 10.00d0*(dexp(-4*(x-1))-(2*dexp(-2*(x-1))))

     END
