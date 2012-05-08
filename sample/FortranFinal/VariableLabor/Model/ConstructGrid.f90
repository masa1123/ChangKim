!	----------------------------------------------------------------------
!	File name: CaptialGrid.cpp
!	----------------------------------------------------------------------


subroutine CapitalGrid()

	use Globals
	
	implicit none
	
	real(8) kmax, kmin, step
	integer i


	! max, min, step for grid construction

	kmax = kss*1.13
	kmin = kss*0.87
	step = (kmax - kmin)/(nk-1)

	kgrid(1) = kmin
	do i = 2, nk
		kgrid(i) = kgrid(i-1) + step
	end do

	! save capital grid

	open(1, file='Output\kgrid.txt', status='unknown')
	write(1,'(f12.6)') kgrid


end subroutine



subroutine ShockGrid()

	use Globals
	
	implicit none
	
	integer i
	real(8) pmax, pmin, step
	
	pmax = avep + 3.0D0*sigp
	pmin = avep - 3.0D0*sigp
	step = (pmax - pmin)/(np - 1)

	pgrid(1) = pmin
	do i = 2, np
		pgrid(i) = pgrid(i-1) + step
	end do

	epgrid = dexp(pgrid)

	open(1, file='Output\pgrid.txt', status='unknown')
	write(1,'(f12.6)') pgrid

	open(2, file='Output\epgrid.txt', status='unknown')
	write(2,'(f12.6)') epgrid


end subroutine
