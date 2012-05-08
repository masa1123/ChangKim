!	----------------------------------------------------------------------
!	Module name: GridsModule.f90
!	----------------------------------------------------------------------


module GridsModule

	use Globals

	implicit none

contains


subroutine CapitalGrid()

	implicit none

	integer i


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



subroutine ShockGrid()

	implicit none
	
	integer i
	real(8) xmax, xmin, step
	

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



end module
