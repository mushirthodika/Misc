!--------------------------------------------------------------------------------
!      THIS PROGRAM CALCULATES THE ABSORPTION SPECTRA USING THE
!      ENERGETICS OF THE GROUND STATES AS WELL AS THE EXCITED STATES AND THE
!      TRANSITION DIPOLE MOMENTS.
!--------------------------------------------------------------------------------

!      AUTHOR   :: MUSHIR UL HASAN K T
!      DATE     :: 03/09/2018
!      VERSION  :: 0.01

!--------------------------------------------------------------------------------
       PROGRAM ABSORPTION_SPECTRA
!--------------------------------------------------------------------------------

       IMPLICIT NONE

!--------------------------VARIABLE DECLARATION----------------------------------
       INTEGER*4            :: I,J,I1,J1,I2
       INTEGER*4, PARAMETER :: NPOINTS=74, NGRID=1000

       REAL*8               :: E(NGRID),DELTA_E,MU_TOT,DE,C,F(NPOINTS,1)
       REAL*8               :: EN(NPOINTS,4),TDM_X(NPOINTS,3),TOTAL_SUM
       REAL*8               ::TDM_Y(NPOINTS,3),TDM_Z(NPOINTS,3),P(NGRID)
       REAL*8               :: G(NGRID,NPOINTS),S(NGRID,1),F1(NPOINTS,1)
       REAL*8, PARAMETER    :: TW_BY_THR=0.666666D0, ONE_BY_FOUR=0.25000000D0
       REAL*8, PARAMETER    :: EMIN=0.00000000D0, EMAX=10.00000000D0
       REAL*8, PARAMETER    :: ZERO=0.00000000D0, DELTA=0.05000000D0
       REAL*8, PARAMETER    :: PI=3.1415926535897D0
       REAL*8, PARAMETER    :: FACTOR=27.211600000000D0
!--------------------------------------------------------------------------------

      DE = (EMAX-EMIN)/REAL(NGRID - 1)

      OPEN(10,FILE="energetics.dat")

      DO I = 1, NPOINTS

                READ(10,*) EN(I,:)

      ENDDO

      CLOSE(10)

      OPEN(11,FILE="TDMs-x.dat")

      OPEN(12,FILE="TDMs-y.dat")

      OPEN(13,FILE="TDMs-z.dat")

      DO I = 1, NPOINTS

                READ(11,*) TDM_X(I,:)

                READ(12,*) TDM_Y(I,:)

                READ(13,*) TDM_Z(I,:)

      ENDDO

      CLOSE(11)

      CLOSE(12)

      CLOSE(13)

      OPEN(14,FILE='abs_spectra_unnorm.dat')




            TOTAL_SUM = ZERO

             J1 = 2!, 4    

                     S = ZERO

                  DO I1 = 1, NPOINTS 
 
                              DELTA_E = EN(I1,J1) - EN(I1,1)

                              F1(I1,1) = TW_BY_THR*(DELTA_E)

                             MU_TOT = (TDM_X(I1,J1-1)*TDM_X(I1,J1-1)) +&
                                   & (TDM_Y(I1,J1-1)*TDM_Y(I1,J1-1)) +&
                                   & (TDM_Z(I1,J1-1)*TDM_Z(I1,J1-1))

                             F(I1,1) = F1(I1,1)*MU_TOT

                        DO I2 = 1, NGRID

                              E(I2) = EMIN + REAL(I2-1)*DE


                G(I2,I1) = DELTA*(((E(I2) - (DELTA_E*FACTOR))*&
                   & (E(I2) - (DELTA_E*FACTOR))) +&
                   ((DELTA*DELTA)*(ONE_BY_FOUR)))**(-1)*F(I1,1)


                       ENDDO 

                  !  TOTAL_SUM = TOTAL_SUM + S/REAL(NPOINTS)

                  ENDDO

               !DO I = 1,NGRID

                    !S = MATMUL(G,F)

                    DO I = 1,NGRID

                      P(I) = SUM(G(I,:))

                    WRITE(14,'(2f20.10)') E(I),P(I)/NPOINTS

                    ENDDO

      !ENDDO

      CLOSE(14)

!--------------------------------------------------------------------------------
      END PROGRAM ABSORPTION_SPECTRA
!--------------------------------------------------------------------------------
