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



subroutine CapitalGridMu()

	implicit none

	integer i


	!	construct unequally spaced asset grids

	!	(0, a1)

	agridmu(1) = amin
	do i = 2, namu1
		agridmu(i) = agridmu(i-1) + smu1
	end do

	!	(a1, a2)

	do i = namu1 + 1, namu2
		agridmu(i) = agridmu(i-1) + smu2
	end do

	!	(a2, a3)

	do i = namu2 + 1, namu3
		agridmu(i) = agridmu(i-1) + smu3
	end do

	!	(a3, a4)

	do i = namu3 + 1, namu4
		agridmu(i) = agridmu(i-1) + smu4
	end do

	!	(a4, a5)

	do i = namu4 + 1, namu5
		agridmu(i) = agridmu(i-1) + smu5
	end do

	!	(a5, a6)

	do i = namu5 + 1, namu6
		agridmu(i) = agridmu(i-1) + smu6
	end do

	!	(a6, a7)

	do i = namu6 + 1, namu7
		agridmu(i) = agridmu(i-1) + smu7
	end do


	! save grid

	open(1, file='Output\agridmu.txt', status='unknown')
	write(1, '(f12.6)') agridmu

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
