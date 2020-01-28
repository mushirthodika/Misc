      ! THIS PROGRAM CALCULATES THE ONE-ELECTRON CAP INTEGRALS IN
      ! THE ATOMIC ORBITAL BASIS

      ! AUTHOR            :: MUSHIR
      ! START DATE        :: 11/26/2019

      PROGRAM CAP_NUMERICAL_INTEGRATION

      IMPLICIT NONE

!------------------VARIABLE DECLARATION------------------------------------------ 
      INTEGER*4, PARAMETER :: npoints = 500
      INTEGER*4            :: i, i1, i2, i3

      REAL*8               :: alpha_mu, alpha_nu, A1, A2, A3, B1, B2, B3
      REAL*8               :: c1, c2, c3, ka1, ka2, ka3, kb1, kb2, kb3
      REAL*8               :: AB1, AB2, AB3, xi1, xi2, xi3, alpha, W
      REAL*8               :: cap_x(npoints), cap_w(npoints)
!      REAL*8               :: int1, int2_xtot
!      REAL*8               :: int2_x, main_int, W
!--------------------------------------------------------------------------------

      alpha_mu = 0.0500000000d0
 
      alpha_nu = 0.0500000000d0

      alpha = alpha_mu + alpha_nu

      A1 = 0.00000000d0
      A2 = 0.00000000d0
      A3 = 0.69688403d0

      B1 = 0.00000000d0
      B2 = 0.00000000d0
      B3 = -0.69688403d0

      c1 = 1.000000000d0
      c2 = 1.000000000d0
      c3 = 2.200000000d0

      ka1 = 0.00000000d0
      ka2 = 0.00000000d0
      ka3 = 0.00000000d0

      kb1 = 0.00000000d0
      kb2 = 0.00000000d0
      kb3 = 0.00000000d0

      open(111,file='cap_x.txt')
      open(222,file='cap_w.txt')

      do i = 1,500

           read(111,*) cap_x(i)
           read(222,*) cap_w(i)

      enddo

      close(111)
      close(222)

!------------------------------CALCULATION OF XI_1-------------------------------!

               xi1 = 0.000000000d0

                AB1 = ((alpha_mu*A1) + (alpha_nu*B1))/alpha

                do i1 = 1, npoints

                if (abs(cap_x(i1)) .le. c1) then

                     W = 0

                else

                     W = (abs(cap_x(i1)) - c1)*(abs(cap_x(i1)) - c1)

                endif

                xi1 = xi1 + (((cap_x(i1) - A1)**ka1)*&

                     & ((cap_x(i1) - B1)**kb1)*&

                     & dexp(-alpha*(cap_x(i1) - AB1)*&

                     & (cap_x(i1) - AB1))*W)*cap_w(i1)

                enddo
!--------------------------------------------------------------------------------!

!------------------------------CALCULATION OF XI_2-------------------------------!

               xi2 = 0.000000000d0

                AB2 = ((alpha_mu*A2) + (alpha_nu*B2))/alpha

                do i2 = 1, npoints

                if (abs(cap_x(i2)) .le. c2) then

                     W = 0

                else

                     W = (abs(cap_x(i2)) - c2)*(abs(cap_x(i2)) - c2)

                endif

                xi2 = xi2 + (((cap_x(i2) - A2)**ka2)*&

                     & ((cap_x(i2) - B2)**kb2)*&

                     & dexp(-alpha*(cap_x(i2) - AB2)*&

                     & (cap_x(i2) - AB2))*W)*cap_w(i2)

                enddo
!--------------------------------------------------------------------------------!

!------------------------------CALCULATION OF XI_3-------------------------------!

               xi3 = 0.000000000d0

                AB3 = ((alpha_mu*A3) + (alpha_nu*B3))/alpha

                do i3 = 1, npoints

                if (abs(cap_x(i3)) .le. c3) then

                     W = 0

                else

                     W = (abs(cap_x(i3)) - c3)*(abs(cap_x(i3)) - c3)

                endif

                xi3 = xi3 + (((cap_x(i3) - A3)**ka3)*&

                     & ((cap_x(i3) - B3)**kb3)*&

                     & dexp(-alpha*(cap_x(i3) - AB3)*&

                     & (cap_x(i3) - AB3))*W)*cap_w(i3)

                enddo
!--------------------------------------------------------------------------------!

          print*, xi1, xi2, xi3

!         int2_xtot = 1.0000000000d0

!                    do j2 = 1,2

!                        do i1 = 1, 500

!                      int2_x = int2_x + (((cap_x(i1) - A(j2))**ka(j2))*&
!                        & ((cap_x(i1) - B(j2))**kb(j2))*&
!                        & dexp(-alpha*(cap_x(i1) - AB(j2))*&
!                        & (cap_x(i1) - AB(j2))))*cap_w(i1)
!
!                        enddo
!
!                        int2_xtot = int2_xtot*int2_x
!
!                    enddo
!
!                    main_int = int1*int2_xtot
!
!                    print*, main_int
!
!                elseif (j == 2) then
!
!                    do j2 = 1,3
!
!                        do i1 = 1, 500
!
!                      int2_y = int2_y + (((cap_x(i1) - A(j2))**ka(j2))*&
!                        & ((cap_x(i1) - B(j2))**kb(2j))*&
!                        & dexp(-alpha*(cap_(i1) - AB(j2))*&
!                        & (cap_x(i1) - AB(j2))))*cap_w(i1)
!
!                        enddo
!
!                        int2_ytot = int2_ytot*int2_y
!
!                    enddo
!
!                elseif (j == 3) then
!
!                    do j2 = 1,2
!
!                        do i1 = 1, 500
!
!                      int2_z = int2_z + (((cap_x(i1) - A(j2))**ka(j2))*&
!                        & ((cap_x(i1) - B(j2))**kb(2j))*&
!                        & dexp(-alpha*(cap_(i1) - AB(j2))*&
!                        & (cap_x(i1) - AB(j2))))*cap_w(i1)
!
!                        enddo
!
!                        int2_ztot = int2_ztot*int2_z
!
!                    enddo
!
!                endif
!
!           enddo
!
!      enddo

      END 
