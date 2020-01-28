! Discrete variable representation (DVR)
! Author :: Mushir Ul hasan
! Please remove the stop command before compilation
Program dvr

implicit none
integer (kind=8), parameter :: n=1001
real (kind=8)               :: h(n,n),ev(n,n),e(n)
integer (kind=8)            :: i,j
real (kind=8)               :: dx,s,x,v,w,k,c1,c
real (kind=8),parameter     :: pi=3.141592653589793238462d0!,k4=1.00000d0
!real (kind=8), parameter    :: k0 = -132.7074997d0,k2=7.0000d0,k3=0.500000d0
real (kind=8), parameter    :: xmin=-20.000000d0,xmax=20.0000000000d0
real (kind=8), parameter    ::  mu = 
dx=(xmax-xmin)/n
s=(0.50000000/(dx*dx*mu))
c=((pi*pi)/3.0000d0)*s
c1=((n+1)*0.50000000000d0)
open(18,file='symmdwpot.dat')

do i=1,n
        x=xmin+real(i-1)*dx
        !v = 0.5000000*x*x
        h(i,i) = c !+ v 
        w = i - c1

do j=i+1,n
        k = w + (j-i)
if (mod((w-k),2.0000d0)==0) then
        h(i,j) = s*(2.00d0/((w - k)*(w - k)))
else
        h(i,j) = s*(2.00d0/((w - k)*(w - k)))*(-1)
endif
        h(j,i)=h(i,j)
enddo

write(18,*) x,v

enddo
open(13,file='KE.txt')
do i=1,n
write(13,*) (h(i,j),j=1,n)
enddo
close(18)
stop
call CALL_DSYEV_WT_EIGVECS(h,n,e,ev)

open(14,file='eval_symmdw.dat')
open(15,file='evec_symmdw.dat')
open(18,file='symmdwpot.dat')
do i=1,n
        read(18,*) x,v
        write(14,*) e(i)
        write(15,*) x,v,(ev(i,j),j=1,n)
enddo
close(14)
close(15)
close(18)
END

