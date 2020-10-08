      PROGRAM WE_GONNA_DO_SOME_SHIT

              IMPLICIT NONE

              ! VARIABLE DECLARATION

              REAL*8               :: E_REF, RXCOORD, EN_HARTREE, FAC=27.2114
              INTEGER*4, PARAMETER :: NSTATES = 15, NPOINTS = 10
              INTEGER*4            :: I, J
              CHARACTER (LEN=80)   :: FLNOUT, FLNINP, B, FORMAT

              !-----------------------------------------

              FORMAT = '(2F12.6)'


                   OPEN(111,file="irc.dat")
              
                   OPEN(222,file="gs-hartree-1.dat")

                   READ(222,*) E_REF

                   CLOSE(222)

                   OPEN(222,file="gs-hartree-1.dat")
                   OPEN(333,file="gs-ev-1.dat")

                   DO J = 1, NPOINTS

                        READ(111,*) RXCOORD
                        READ(222,*) EN_HARTREE

                        WRITE(333,FORMAT) RXCOORD,(EN_HARTREE-E_REF)*FAC

                   ENDDO

              CLOSE(111)
              CLOSE(222)
              CLOSE(333)
              CLOSE(444)

      END
