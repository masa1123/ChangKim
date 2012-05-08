!	----------------------------------------------------------------------
!	File name : Main.f90
!
!	Main program that solves the model with aggregate fluctuations.
!	----------------------------------------------------------------------


program Fluctuations

    include 'link_f90_dll.h'

	use Globals
	use InitializationModule
	use GridsModule
	use TransitionModule
	use ValueModule
	use SimulationModule
	use RegressionModule
	use dfport

	implicit none

	integer iterLOM
	real(8) errK, errR, errW, err
	character(8) starttime, finishtime, solvetime, simultime
	character(10) startdate, finishdate, today


	!	open a log file

	open(0, file='Output\Logfile.txt', status = 'unknown')


	!	set the starting time

	call date(startdate)
	starttime = clock()


	!	initialize parameters and global variables

	call InitializeParameters()
	call InitializeCoefficients()


	!	construct grids

	call IndividualCapitalGrid()
	call IdiosyncraticShockGrid()
	call AggregateCapitalGrid()
	call AggregateShockGrid()


	!	initialize transition probability matrix

	call IdiosyncraticTransition()
	call AggregateTransition()


	!	initialize aggregate shock for simulation

	call AggregateShockSeries()


	!	initialize value function and decision rule for asset

	call InitializeValue()


	!	iteration on laws of motion

	do iterLOM = 1, maxiterLOM


		!	solve value function

		call date(today)
		solvetime = clock()

		write(0, '(//A,I2)') "Coefficients of Equations in Iteration ", iterLOM
		write(0, '(/A, 3F10.5)') "Law of Motion for Kapital:", Kcoef
		write(0, '(/A, 3F10.5)') "Interest Rate Equation   :", Rcoef
		write(0, '(/A, 3F10.5)') "Wage Rate Equation       :", Wcoef
		write(0, '(//A,I2,A,A,A)') "SolveValueFunction in Iteration ",		&
								  iterLOM, " started at ", today, solvetime
	
		call SolveValueFunction()


		!	generate artificial data through simulation

		call date(today)
		simultime = clock()
		write(0, '(/A,I2,A,A,A)') "SimulateData in Iteration ",			&
								 iterLOM, " started at ", today, simultime

		call SimulateData()


		!	updata laws of motion through regressions of simulated data

		call RegressLOM()

		errK = maxval(dabs(NewKcoef - Kcoef))
		errR = maxval(dabs(NewRcoef - Rcoef))
		errW = maxval(dabs(NewWcoef - Wcoef))
		err  = max(errK, errR, errW)

		Kcoef = sfac*NewKcoef + (1.0D0 - sfac)*Kcoef
		Rcoef = sfac*NewRcoef + (1.0D0 - sfac)*Rcoef
		Wcoef = sfac*NewWcoef + (1.0D0 - sfac)*Wcoef


		!	convergence test

		if (errK < tol_LOM) then
!		if (err < tol_LOM) then
			write(0, '(//A,I2)') "Errors in LOM Iteration ", iterLOM
			write(0, '(/A,F10.6,A,F10.6,A,F10.6)')			&
				"errK = ", errK, "    errR = ", errR, "    errW = ", errW
			write(0, '(//A,I2)') "Law of Motion converged at iteration ", iterLOM
			write(0, '(//A)') "Converged Coefficients of Equations"
			write(0, '(/A, 3F10.5)') "Law of Motion for Kapital:", Kcoef
			write(0, '(/A, 3F10.5)') "Interest Rate Equation   :", Rcoef
			write(0, '(/A, 3F10.5)') "Wage Rate Equation       :", Wcoef
			exit 
		else
			write(0, '(//A,I2)') "Errors in LOM Iteration ", iterLOM
			write(0, '(/A,F10.6,A,F10.6,A,F10.6)')			&
				"errK = ", errK, "    errR = ", errR, "    errW = ", errW
		end if

	end do


	!	whether to construct panel data in final simulation?

	panel = .true.


	!	final simulation

	call date(today)
	simultime = clock()
	write(0, '(//A,A,A)') "FinalSimulation started at ", today, simultime
	call FinalSimulation()
	call FinalRegressLOM()
	call FinalRegressLOM2()


	!	save results

	call SaveValueFunction2()
	call SaveCoefficients()
	call SaveSimulatedTimeSeriesData()


	!	set the finishing time

	call date(finishdate)
	finishtime = clock()

	
	!	print the starting and finishing times

	write(0, '(/)') 
	write(0, '(A,A,A)') "The program started at ", startdate, starttime
	write(0, '(A,A,A)') "The program finished at ", finishdate, finishtime


end program Fluctuations
