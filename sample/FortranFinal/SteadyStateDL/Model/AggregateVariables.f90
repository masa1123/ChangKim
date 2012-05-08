!	----------------------------------------------------------------------
!	File name: AggregateCapital.cpp
!
!	Calculates the aggregate capital by numerically integrating 
!	the measure, and then interest rate.
!	----------------------------------------------------------------------


subroutine AggregateVariables()

	use Globals
	use dfport
	use Numerical_Libraries

	implicit none

	real(8) AK, AH, AE, AKE, irateb, wageb, Prod
	integer i, k, kk
	real(8):: break(na,nx), cscoef(4,na,nx)


	!	interpolate HR using cubic spline and calculate efficiency units

	do ix = 1, nx
		
		call dcsint(na, agrid, HR(:,ix), break(:,ix), cscoef(:,:,ix))

		do ia = 1, namu
			HRmu(ia,ix) = dcsval(agridmu(ia), namu-1, break(:,ix), cscoef(:,:,ix))
		end do

		EFmu(:,ix) = HRmu(:,ix)*exgrid(ix)

	end do


	!	aggregate capital, hours

	AK = dot_product(agridmu, sum(mu, dim=2))
	AH = sum(HRmu*mu)
	AE = sum(EFmu*mu)
	AKE = AK/AE

	
	!	real interest rate, wage

	irateb = alpha*AKE**(alpha-1.0D0) - delta
	wageb  = (1.0D0-alpha)*AKE**alpha


	!	output

	Prod = AK**alpha*AE**(1.0D0-alpha)


	!	print results on screen

	write(*, '(A,F13.10)') "beta     = ", bta
	write(*, '(A,F10.4)') "B        = ", B
	write(*, '(A,F10.4)') "ave_x    = ", avex
	write(*, '(A,F10.4)') "rho_x    = ", rhox
	write(*, '(A,F10.4)') "sig_x    = ", sigex
	write(*, *) " "
	print*, "Aggregate Capital                     = ", AK
	print*, "Aggregate Hours Worked                = ", AH
	print*, "Aggregate Effeciency Units of Hours   = ", AE
	print*, "Aggregate Capital per Effeciency Unit = ", AKE
	print*, "Aggregate Output                      = ", Prod
	print*, ""
	print*, "Discount Factor                       = ", bta
	print*, "Real Interest Rate                    = ", irateb
	print*, "Wage Rate per Effeciency Unit         = ", wageb
	
	
	!	save results on file

	open(1, file='Output\Aggregates.txt', status='unknown')

	write(1, '(A/)') "Parameter Values:"
	write(1, '(A,F13.10)') "beta     = ", bta
	write(1, '(A,F10.4)') "B        = ", B
	write(1, '(A,F10.4)') "gamma    = ", gama
	write(1, '(A,F10.4)') "ave_x    = ", avex
	write(1, '(A,F10.4)') "rho_x    = ", rhox
	write(1, '(A,F10.4)') "sig_x    = ", sigex
	write(1, '(/)') 
	write(1, '(A/)') "Simulated Aggregate Variables:"
	write(1, '(A,F15.10)') "Aggregate Capital                     = ", AK
	write(1, '(A,F15.10)') "Aggregate Hours Worked                = ", AH
	write(1, '(A,F15.10)') "Aggregate Effeciency Units of Hours   = ", AE
	write(1, '(A,F15.10)') "Aggregate Capital per Effeciency Unit = ", AKE
	write(1, '(A,F15.10)') "Aggregate Output                      = ", Prod
	write(1, '(/)') 
	write(1, '(A,F15.10)') "Real Interest Rate                    = ", irateb
	write(1, '(A,F15.10)') "Wage Rate per Effeciency Unit         = ", wageb


	!	set the finishing time

	call date(finishdate)
	finishtime = clock()


	!	print the starting and finishing times

	write(*, '(/)') 
	write(*, '(A,A,A)') "The program started at  ", startdate, starttime
	write(*, '(A,A,A)') "The program finished at ", finishdate, finishtime

	write(1, '(/)') 
	write(1, '(A,A,A)') "The program started at  ", startdate, starttime
	write(1, '(A,A,A)') "The program finished at ", finishdate, finishtime

end subroutine
