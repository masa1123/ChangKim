!	----------------------------------------------------------------------
!	Module name: ConstructGrids.f90
!	----------------------------------------------------------------------

module GridsModule

	use Globals

	implicit none

	integer i
	real(8) step


contains


subroutine IndividualCapitalGrid()

	implicit none


	!	construct unequally spaced asset grids

	!	(0, a1)

	agrid(1) = amin
	do i = 2, na1
		agrid(i) = agrid(i-1) + s1
	end do

	!	(a1, a2)

	do i = na1 + 1, na2
		agrid(i) = agrid(i-1) + s2
	end do

	!	(a2, a3)

	do i = na2 + 1, na3
		agrid(i) = agrid(i-1) + s3
	end do

	!	(a3, a4)

	do i = na3 + 1, na4
		agrid(i) = agrid(i-1) + s4
	end do

	!	(a4, a5)

	do i = na4 + 1, na5
		agrid(i) = agrid(i-1) + s5
	end do

	!	(a5, a6)

	do i = na5 + 1, na6
		agrid(i) = agrid(i-1) + s6
	end do

	!	(a6, a7)

	do i = na6 + 1, na7
		agrid(i) = agrid(i-1) + s7
	end do


	! save grid

	open(1, file='Output\agrid.txt', status='unknown')
	write(1, '(f12.6)') agrid

end subroutine



subroutine AggregateCapitalGrid()

	implicit none

	real(8) step


	! step for kgrid construction

	step = (kmax - kmin)/(nk - 1)

	
	! aggregate capital grid

	kgrid(1) = kmin
	do i = 2, nk
		kgrid(i) = kgrid(i-1) + step
	end do


	! save grid

	open(1, file='Output\kgrid.txt', status='unknown')
	write(1, '(f12.6)') kgrid
	close(1)
	

	! step for k2grid construction

	step = (k2max - k2min)/(nk2 - 1)

	
	! aggregate capital grid

	k2grid(1) = k2min
	do i = 2, nk2
		k2grid(i) = k2grid(i-1) + step
	end do


	! save grid

	open(1, file='Output\k2grid.txt', status='unknown')
	write(1, '(f12.6)') k2grid
	close(1)
	

end subroutine



subroutine IdiosyncraticShockGrid()

	implicit none
	
	real(8) xmax, xmin
	

	! xgrid

	xmax = avex + 3.0*sigx
	xmin = avex - 3.0*sigx
	step = (xmax - xmin)/(nx - 1)

	xgrid(1) = xmin
	do i = 2, nx
		xgrid(i) = xgrid(i-1) + step
	end do

	exgrid = dexp(xgrid)


	! save grids

	open(1, file='Output\xgrid.txt', status='unknown')
	open(2, file='Output\exgrid.txt', status='unknown')

	write(1,'(f12.6)') xgrid
	write(2,'(f12.6)') exgrid


end subroutine



subroutine AggregateShockGrid()

	implicit none
	
	real(8) pmax, pmin
	

	! pgrid: more than 2 states

	pmax = avep + pbnd*sigp
	pmin = avep - pbnd*sigp
	step = (pmax - pmin)/(np - 1)

	pgrid(1) = pmin
	do i = 2, np
		pgrid(i) = pgrid(i-1) + step
	end do

	epgrid = dexp(pgrid)


	! save grids

	open(1, file='Output\pgrid.txt', status='unknown')
	open(2, file='Output\epgrid.txt', status='unknown')

	write(1,'(f12.6)') pgrid
	write(2,'(f12.6)') epgrid


end subroutine


end module
