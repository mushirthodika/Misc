 PROGRAM tesselation

 implicit none

 ! declare variables
 integer (kind=8), parameter :: nx=101, ny=nx, n_intp=60, nr=nx*ny, max_iter=100000
 real    (kind=8), parameter :: xmin=-4.000d0, xmax=4.000d0, ymin=xmin, ymax=xmax
 real    (kind=8), parameter :: tol_conv=1e-3
 
 integer (kind=4)  :: i,j,ii,ipt,jj,k,iter
 integer (kind=4)  :: rseed,igrid, intp(n_intp)
 real    (kind=8)  :: dx,dy,bounds(2,2),newpts(2,n_intp),a(nx)
 real    (kind=8)  :: rho,oldpts(2,n_intp),dist,dist_tmp,diff
 real    (kind=8)  :: weightedpos(2),weights,mindist,minpts(2),vtmp(2)
 real    (kind=8), allocatable :: icenter(:),fullgrid(:,:)
 real    (kind=8), allocatable :: fullrho(:)

 allocate(icenter(nr),fullgrid(2,nr),fullrho(nr))

 open(111,file='cvt_debug.dat')

! generate the grid points

 dx=(xmax-xmin)/real(nx-1)
 dy=(ymax-ymin)/real(ny-1)

 write(111,*) dx,dy,xmin

 k = 1

 do i = 1,nx
    a(i) = xmin + real(i-1)*dx
    do j = 1,ny
       fullgrid(1,k) = a(i)
       fullgrid(2,k) = ymin + real(j-1)*dy
       k = k+1
    enddo
 enddo

! generate rho on the full grid

 do i = 1, nr
    fullrho(i) = rho(fullgrid(1,i),fullgrid(2,i))
 enddo

 ! find the lower and higher bounds of full-grid points
 do ii = 1,2
     bounds(1,ii) = minval(fullgrid(ii,:))
     bounds(2,ii) = maxval(fullgrid(ii,:))
 enddo

 rseed = 1518543090
 write(111,'(a,i20)') "# rseed for random number generator is: ", rseed

 call srand(rseed)

 write(111,*) rand(0)

! initialize the tesselation points and make sure they are within the bounds

 do ipt = 1, n_intp
      do ii = 1,2
          newpts(ii,ipt) = bounds(1,ii) + rand(0)*0.9*(bounds(2,ii)-bounds(1,ii)) 
      enddo
 enddo

! print the tesselation points in the debug file

 do ii = 1, 2
      write(111,*) (newpts(ii,jj),jj=1,n_intp)
 enddo

 ! perform centroidal voronoi tesselation (CVT) algorithm

 ! Grid points closest to each interpolation point are clustered below 

 oldpts = newpts
 do iter = 1, max_iter
    ! for each point in the full grid, find which interpolation points is 
    ! the closest it. Then put the index of intp point to icenter(:)
    do igrid = 1, nr
       ! set a very large initial value for dist
       dist = 100.0 * ( maxval(bounds(2,:)) - minval(bounds(1,:)) )**2  
       icenter(igrid) = -1
       ! find which intp point is closest to the current grid point
       do ipt = 1, n_intp
          vtmp = newpts(:,ipt) - fullgrid(:,igrid)
          dist_tmp = dot_product(vtmp, vtmp)
          if (dist_tmp < dist) then
             dist = dist_tmp
             icenter(igrid) = ipt
          endif
       enddo 
    enddo

    do igrid = 1, nr

       write(111,*) igrid, icenter(igrid)

    enddo

    ! Now update the interpolation pts
    diff = 0.0
    do ipt = 1, n_intp
       weightedpos(1:2) = 0.0
       weights = 0.0
       do igrid = 1, nr
          if (icenter(igrid) .ne. ipt) cycle
          weightedpos = weightedpos + fullgrid(:,igrid)*fullrho(igrid)
          weights = weights + fullrho(igrid)
       enddo
       ! update the new intp points with the centroid
       ! This is just a quick simple trick to avoid the case of weights == 0. It
       if(weights .lt. 1.0e-10) then
          newpts(:,ipt) = oldpts(:,ipt)
       else
          newpts(:,ipt) = weightedpos(:) / weights
       endif
       mindist = 1.e9 ! initialize a very large number
       ! Loop over all the grid points to find a grid point that is closest to newpts(:,ipt)
       do igrid = 1, nr
             vtmp = newpts(1:2,ipt) - fullgrid(1:2,igrid)
             dist = sqrt(dot_product(vtmp,vtmp)) 
             if (dist < mindist) then
                mindist = dist
                minpts(1:2) = fullgrid(1:2,igrid)
      !          minig = igrid
                intp(ipt) = igrid
             endif
       enddo ! igrid
       newpts(1:2,ipt) = minpts(1:2)
      ! select_grid(ipt) = minig
       vtmp = newpts(:,ipt) - oldpts(:,ipt)
       diff = diff + sqrt(dot_product(vtmp,vtmp))
       !write(outdbg,'(8f8.3)') newpts(1:3,ipt), vtmp(1:3), &
       !  sqrt(dot_product(vtmp,vtmp)),weights
    enddo ! loop ipt
   !write(*,*) select_grid(:)

    write(*,'(i8,a,f18.12)') iter, " diff (a.u.) ", diff/n_intp
    if (diff < tol_conv) then ! conv threshold is satisfied, break the loop!
       exit 
    endif
    oldpts = newpts
 enddo ! iter
 
 do ipt = 1, n_intp
    write(111,'(i8,2f12.6)') intp(ipt),fullgrid(1:2,intp(ipt)) !,fullgrid(1:2,intp(ipt))*0.529177+ 4.1011
 enddo

 close(111) ! close the dbg file

 END

!---------------------------------------------------------------------------------
 FUNCTION rho(x,y)
!---------------------------------------------------------------------------------

 implicit none

! declare variables
 
 real (kind=8), parameter :: c=0.500000d0
 real (kind=8)            :: x,y,rho


 rho = x*x*dexp(c*(-(x*x)-(y*y)))
 !rho = dexp(c*(-(x*x)-(y*y)))

 END
