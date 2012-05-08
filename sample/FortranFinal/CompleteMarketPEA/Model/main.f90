!	----------------------------------------------------------------------
!	Main Program: CompleteMarketPEA
!
!	Model:	   Complete market with worker heterogeneity 
!			   (benchmark specification: rhox=0.95, sgix=0.225)
!	Algorithm: Parameterized Expectation Algorithm a la Marcet (1999)
!	----------------------------------------------------------------------


program CompleteMarketPEA

    include 'link_f90_dll.h'

	use Globals
	use InitializationModule
	use SimulationModule
	use dfport
	
	implicit none

	integer iter
	real(8) err
	character(8) starttime, finishtime
	character(10) startdate, finishdate


	!	set the varaibles for numerical integration

	errabs = 1.0D-6
	errrel = 0.0D0
	irule  = 2


	!	set the starting time

	call date(startdate)
	starttime = clock()


	!	initialize the model parameters

	call InitializeModelParameters()
	call InitializeExpectationParameters()

	
	!	solve for the steady state equilibrium

	call SteadyState()


	!	generate shocks

	call GenerateShocks()


	!	iteration starts

	final = .false.
	iter = 0
	err = 10000.0D0

	do while (err > tol_Iter .and. iter < Max_Iter)

		iter = iter + 1

		call Simulation()

		call Regression()

		err = maxval(dabs(coefNew - coef))
		coef = sfac*coefNew + (1.0D0-sfac)*coef

		if (mod(iter,10) == 0) then
			print*, iter, err
		end if

	end do


	!	final simulation

	final = .true.
	call Simulation()


	!	save results

	print*, iter
	print*, coef
	call SaveTimeSeries()


	!	set the finishing time

	call date(finishdate)
	finishtime = clock()


	!	print the starting and finishing times

	write(*, '(/)') 
	write(*, '(A,A,A)') "The program started at  ", startdate, starttime
	write(*, '(A,A,A)') "The program finished at ", finishdate, finishtime

end program CompleteMarketPEA