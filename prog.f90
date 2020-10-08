      PROGRAM WE_GONNA_DO_SOME_SHIT

              IMPLICIT NONE

              ! VARIABLE DECLARATION

              REAL*8               :: E_REF, RXCOORD, EN_HARTREE, FAC=27.2114
              INTEGER*4, PARAMETER :: NSTATES = 15, NPOINTS = 10
              INTEGER*4            :: I, J
              CHARACTER (LEN=80)   :: FLNOUT, FLNINP, B, FORMAT

              !-----------------------------------------

              FORMAT = '(2F12.6)'


              DO I = 1, NSTATES

                   WRITE(B,'(i2)') I

                   B = ADJUSTL(B)

                   FLNOUT = 'state-ev-'//TRIM(B)//'.dat'
                   FLNINP = 'state-hartree-'//TRIM(B)//'.dat'

                   FLNOUT = TRIM(FLNOUT)
                   FLNINP = TRIM(FLNINP)

                   FLNOUT = ADJUSTL(FLNOUT)
                   FLNINP = ADJUSTL(FLNINP)

                   OPEN(111,file="irc.dat")
              
                   OPEN(222,file="gs-hartree-1.dat")

                   READ(222,*) E_REF

                   OPEN(333,file=FLNINP)
                   OPEN(444,file=FLNOUT)

                   DO J = 1, NPOINTS

                        READ(111,*) RXCOORD
                        READ(333,*) EN_HARTREE

                        WRITE(444,FORMAT) RXCOORD,(EN_HARTREE-E_REF)*FAC

                   ENDDO

              CLOSE(111)
              CLOSE(222)

              ENDDO

              CLOSE(333)
              CLOSE(444)

      END
