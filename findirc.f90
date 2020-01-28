! THIS PROGRAM PRINTS THE COORDINATES OF THE GEOMETRIES ALONG THE IRC
! AUTHOR :: MUSHIR UL HASAN
! DATE   :: 11/07/16

     !-------------------------------!
      PROGRAM IRC_GEOM_SEARCH
     !-------------------------------!

      IMPLICIT NONE

     !--------------------------------------------------------------------------!
       CHARACTER (LEN=80) :: TEXT, INPUT
       REAL (KIND=8)      :: X(3),Y(3),Z(3)
       INTEGER            :: SNO(3),ANO(3),ATYP(3),I,J
     !--------------------------------------------------------------------------!
    
      X = 0.000000000d0
      Y = 0.000000000d0
      Z = 0.000000000d0


      OPEN(111,FILE='irc.log')
      OPEN(222,FILE='coordinates.dat')
     
      DO

      READ(111,'(A80)') TEXT

      IF (TEXT(17:44) == 'Cartesian Coordinates (Ang):') THEN 
     
      DO 

      READ(111,'(A80)') TEXT

      IF ((TEXT(2:7)) == 'Number') THEN
    
      EXIT

      ENDIF

      ENDDO

      READ(111,'(A80)') TEXT 

      DO I = 1, 3

      READ(111,*) SNO(I),ANO(I),X(I),Y(I),Z(I)

      WRITE(222,11) SNO(I),ANO(I),X(I),Y(I),Z(I)
  11  format(i1,i1,f10.6,f10.6,f10.6)
      ENDDO

!      WRITE(222,*)

      ENDIF

      ENDDO

      CLOSE(111)
      CLOSE(222)

      END

