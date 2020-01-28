!This program generates the input files for GAMESS or GAUSSIAN calculation
! Author :: Mushir Ul Hasan
! Date   :: 09/07/2016

 PROGRAM gamess_or_gaussian_input

 IMPLICIT NONE

  character (len=80)  :: text(3),flninp,c
  character (len=10)  :: b,extension
  integer             :: i,j

 open(12,file='coordinates.dat')


 do i = 1,72


        do j = 1,3
 
               read(12,'(a80)') text(j)
 
        enddo

  write(b,'(i2)') i
  write(*,*) b

  b = adjustl(b)
  flninp = trim(b)//'.inp'
  flninp=trim(flninp)
  call gamessinp(text,flninp)
 
 enddo

 close(12)

 END 

!-------------------------------------------------------------------------------!
 SUBROUTINE gamessinp(text,flninp)
!-------------------------------------------------------------------------------!

 IMPLICIT NONE

 character (len=80)  :: text(3)
 character (len=10)  :: flninp

 
 flninp = adjustl(flninp) 

 open(111,file=flninp)

 write(111,*) '$CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG=0 MULT=1 $END'
! write(111,*) '$TDDFT NSTATE=10 IROOT=4 MULT=1 $END'
 write(111,*) '$SCF VVOS=.TRUE. $END'
 write(111,*) '$SYSTEM MWORDS=20 $END'
 write(111,*) '$DATA'
 write(111,*) 'Excited state calcualtion along IRC'
 write(111,*) 'C1'

 write(111,11) 'H' ,'', '1.0' ,'', text(1)
 11 format (a1,a2,a4,a2,a78)
 call basisset('h.txt')
 write(111,12) 'C' ,'', '6.0' ,'', text(2)
 12 format (a1,a2,a4,a2,a78)
 call basisset('c.txt')
 write(111,13) 'N' ,'', '7.0' ,'', text(3)
 13 format (a1,a2,a4,a2,a78)
 call basisset('n.txt')
 write(111,*) '$END'

 close(111)

END

!-------------------------------------------------------------------------------!
 SUBROUTINE gaussianinp(text,flninp)
!-------------------------------------------------------------------------------!

 IMPLICIT NONE

 character (len=80)  :: text(3)
 character (len=10)  :: flninp


 flninp = adjustl(flninp)

 open(111,file=flninp)

 write(111,10) '%mem=1gb'
 10 format(a8)
 write(111,15) '%nprocs=3'
 15 format(a9)
 write(111,16) '#p M06/gen td=(nstate=10,root=4,singlets)'
 16 format(a41)
 write(111,*) 
 write(111,18) 'Excited state calcualtion along IRC'
 18 format(a35)
 write(111,*) 
 write(111,20) '0 1'
 20 format(a3)

 write(111,11) 'H' ,'', text(1)
 11 format (a1,a2,a78)

 write(111,12) 'C' ,'', text(2)
 12 format (a1,a2,a78)

 write(111,13) 'N' ,'', text(3)
 13 format (a1,a2,a78)

write(111,*)

!call basisset('coemdref.dat')

close(111)

END

!-------------------------------------------------------------------------------!
 SUBROUTINE basisset(filename2)
!-------------------------------------------------------------------------------!

 IMPLICIT NONE

 character (len=80) :: text
 character (len=5) :: filename2
 integer            :: i

 open(14,file=filename2)

 do

  read(14,'(a80)',iostat=i) text

  if (i.lt.0) then
  exit
  endif

  write(111,*) text(1:80)
  !14 format (a78)

 enddo

 close(14)

 END
