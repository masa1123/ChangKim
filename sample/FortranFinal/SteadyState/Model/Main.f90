!	----------------------------------------------------------------------
!	File name : SteadyStateFixedDiscrete.f90
!
!	Solves heterogeneous agent model in which worker's productivity (x,z)
!	is changing over time.
!	----------------------------------------------------------------------


program SteadyState

    include 'link_f90_dll.h'

	use Globals
	use dfport
	use ReservationWageModule
	use InitializationModule
	use GridsModule
	use ValueModule
	use MeasureModule
	use SaveResultsModule

	implicit none


	!	set the starting time

	call date(startdate)
	starttime = clock()


	!	initialize parameters and global variables

	call InitializeParameters()
	call GlobalVariables()


	!	construct grids

	call CapitalGrid()
	call ShockGrid()
	
	
	!	initialize transition probability matrix

	call InitializeTransition()


	!	initialize value function and decision rule for asset

	call InitializeValue()


	!	solve value function

	call SolveValueFunction()


	!	calculate invariant measure
	
	call InitializeMeasure()
	call InvariantMeasure()


	!	calculate aggregate variables

	call AggregateVariables()


	!	distribution
	
	call AssetEarningsDistribution()


	!	find reservation wages

	call FindReservationWages()


	!	find labor supply curve

	call LaborSupply()


	!	save results

	call SaveEquilibrium()
	call SaveDistribution()
	call SaveReservationWages()
	Call SaveLaborSupply()


end program SteadyState
