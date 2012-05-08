!	----------------------------------------------------------------------
!	Module name: SimulationModule.f90
!
!	collection of sub-programs to simulate artificial data
!	----------------------------------------------------------------------



Module SimulationModule

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use Numerical_Libraries
	
	implicit none

	integer, parameter:: xseed = 2
	real(8), dimension(:,:), allocatable:: ASkp, HRkp
	type Individual
		real(8) :: A, AP, H, W, C
		integer :: X
	end type Individual
	type(Individual), dimension(:), allocatable :: Worker

contains



!	----------------------------------------------------------------------
!	File name: SimulateData.f90
!
!	generate artificial data for aggregate capital, interest rate, and 
!	wage rate through simulating individuals' decision rules. These data
!	will be used in RegressLOM().
!	----------------------------------------------------------------------


subroutine SimulateData()

	implicit none


	!	initialize time series data

	Kdata = 0.0D0
	Ldata = 0.0D0
	Rdata = 0.0D0
	Wdata = 0.0D0


	!	set the seed number for random number generator

	call rnset(xseed)


	!	allocate space for cross section data generation

	allocate(Worker(Nindiv))
	allocate(ASkp(na,nx), HRkp(na,nx))


	!	intialize distribution of workers


	do indiv = 1, Nindiv
		Worker(indiv)%A = kss
		Worker(indiv)%H = hbar
		Worker(indiv)%X = mod(indiv, nx) + 1
	end do


	!	initial aggregate capital

	Kdata(1) = sum(Worker%A)/Nindiv


	!	start generating artificial time series data

	do time = 1, Nperiod


		!	interpolate decision rules at (K, p)

		do ix = 1, nx
		do ia = 1, na

			ASkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, AS(ia,ix,:,:))
			HRkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, HR(ia,ix,:,:))

		end do
		end do


		!	cross-section data for asset and hours given (Kapital, pt)

		do indiv = 1, Nindiv

			Worker(indiv)%H  = lininterp1(Worker(indiv)%A, agrid, HRkp(:,Worker(indiv)%X))
			Worker(indiv)%AP = lininterp1(Worker(indiv)%A, agrid, ASkp(:,Worker(indiv)%X))
			Worker(indiv)%AP = min(max(Worker(indiv)%AP, agrid(1)), agrid(na))
			
			Ldata(time) = Ldata(time) + exgrid(Worker(indiv)%X)*Worker(indiv)%H

		end do


		!	generate time series data

		Ldata(time) = Ldata(time)/Nindiv
		Rdata(time) = Pdata(time)*alpha*(Kdata(time)/Ldata(time))**(alpha-1.0D0)
		Wdata(time) = Pdata(time)*(1.0D0-alpha)*(Kdata(time)/Ldata(time))**(alpha)
		

		!	prepare for the simulation in the next period

		Kdata(time+1) = min(max(sum(Worker%AP)/Nindiv, kgrid(1)), kgrid(nk))
		Worker%A  = Worker%AP
		call NextIdiosyncraticShock()


		!	show how simulation is going on

		if (mod(time, 100) == 0) then
			write(*,'(A,I6)')   "time = ", time	
			write(*,'(A,F8.5)') "p = ", Pdata(time)
			write(*,'(A,F8.5)') "K = ", Kdata(time)
			write(*,'(A,F8.5)') "Kp = ", Kdata(time+1)
			write(*,'(A,F8.5)') "L = ", Ldata(time)
			write(*,'(A,F8.5)') "w = ", Wdata(time)
			write(*,'(A,F8.5)') "r = ", Rdata(time)-delta
			write(*,'(/)')
		end if

	end do


	!	deallocate memory

	deallocate(Worker)
	deallocate(ASkp, HRkp)


	!	save time sries data

	open(1, file='Output\Kdata.txt', status='replace')
	open(2, file='Output\Rdata.txt', status='replace')
	open(3, file='Output\Wdata.txt', status='replace')
	open(4, file='Output\Ldata.txt', status='replace')

	write(1, '(f12.6)') Kdata
	write(2, '(f12.6)') Rdata
	write(3, '(f12.6)') Wdata
	write(4, '(f12.6)') Ldata

	close(1)
	close(2)
	close(3)
	close(4)


end subroutine




!	----------------------------------------------------------------------
!	File name: FinalSimulation.f90
!
!	generate artificial data for all aggregate variables, which will be 
!	used in ModelStatistics().
!	----------------------------------------------------------------------


subroutine FinalSimulation()

	implicit none


	!	initialize time series data

	Kdata = 0.0D0
	Ldata = 0.0D0
	Rdata = 0.0D0
	Wdata = 0.0D0
	Edata = 0.0D0


	!	set the seed number for random number generator

	call rnset(xseed)


	!	allocate space for cross section data generation

	allocate(Worker(Nindiv))
	allocate(ASkp(na,nx), HRkp(na,nx))


	!	intialize distribution of workers


	do indiv = 1, Nindiv
		Worker(indiv)%A = kss
		Worker(indiv)%H = hbar
		Worker(indiv)%X = mod(indiv, nx) + 1
	end do


	!	initial aggregate capital

	Kdata(1) = sum(Worker%A)/Nindiv


	!	open a file for panel data set

	open(100, file='Output\Panel.txt', status='unknown')


	!	start generating artificial time series data

	do time = 1, Nperiod


		!	interpolate decision rules at (K, p)

		do ix = 1, nx
		do ia = 1, na

			ASkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, AS(ia,ix,:,:))
			HRkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, HR(ia,ix,:,:))

		end do
		end do


		!	cross-section data for asset and hours given (Kapital, pt)

		do indiv = 1, Nindiv

			Worker(indiv)%H  = lininterp1(Worker(indiv)%A, agrid, HRkp(:,Worker(indiv)%X))
			Worker(indiv)%AP = lininterp1(Worker(indiv)%A, agrid, ASkp(:,Worker(indiv)%X))
			Worker(indiv)%AP = min(max(Worker(indiv)%AP, agrid(1)), agrid(na))

			Ldata(time) = Ldata(time) + exgrid(Worker(indiv)%X)*Worker(indiv)%H
			Edata(time) = Edata(time) + Worker(indiv)%H

		end do


		!	generate time series data for aggregate macro variables

		Kdata(time+1) = min(max(sum(Worker%AP)/Nindiv, kgrid(1)), kgrid(nk))
		K2data(time) = dsqrt(sum(Worker%AP**2.0D0)/Nindiv - Kdata(time)**2.0D0)
!		K2data(time) = sum(Worker%AP**2.0D0)/Nindiv
		Ldata(time) = Ldata(time)/Nindiv
		Rdata(time) = Pdata(time)*alpha*(Kdata(time)/Ldata(time))**(alpha-1.0D0)
		Wdata(time) = Pdata(time)*(1.0D0-alpha)*(Kdata(time)/Ldata(time))**(alpha)
		Edata(time) = Edata(time)/Nindiv
		Ydata(time) = Pdata(time)*(Kdata(time)**alpha)*(Ldata(time)**(1.0D0-alpha))
		Idata(time) = Kdata(time+1) - (1.0D0-delta)*Kdata(time)
		Cdata(time) = Ydata(time) - Idata(time)
		Adata(time) = Ydata(time)/Edata(time)


		!	construct panel data if panel == .true.

		if (time > startpanel .and. panel == .true.) then
		
			do indiv = 1, Nindiv

				Worker(indiv)%W = Wdata(time)*exgrid(Worker(indiv)%X)*Worker(indiv)%H
				Worker(indiv)%C = (1.0D0+Rdata(time)-delta)*Worker(indiv)%A		&
									+ Worker(indiv)%W - Worker(indiv)%AP

				! write individual data on panel file

				write(100, '(3f12.6)') Worker(indiv)%A, Worker(indiv)%W, Worker(indiv)%C

			end do
	
		end if


		!	prepare for the simulation in the next period

		Worker%A  = Worker%AP
		call NextIdiosyncraticShock()


		!	show how simulation is going on

		if (mod(time, 100) == 0) then
			write(*,'(A,I6)')   "time = ", time	
			write(*,'(A,F8.5)') "p = ", Pdata(time)
			write(*,'(A,F8.5)') "K = ", Kdata(time)
			write(*,'(A,F8.5)') "Kp = ", Kdata(time+1)
			write(*,'(A,F8.5)') "L = ", Ldata(time)
			write(*,'(A,F8.5)') "w = ", Wdata(time)
			write(*,'(A,F8.5)') "r = ", Rdata(time)-delta
			write(*,'(/)')
		end if

	end do


	!	deallocate memory

	deallocate(Worker)
	deallocate(ASkp, HRkp)
	close(100)


end subroutine


!	----------------------------------------------------------------------
!	subroutine name: AggregateShockSeries()
!
!	generate time series of aggregate productivity shocks.
!	----------------------------------------------------------------------

subroutine AggregateShockSeries()

	implicit none

	integer, parameter:: pseed = 1


	! set the seed number for random number generator

	call rnset(pseed)


	!	normal random numbers for aggregate productivity shocks

	call drnnor(Nperiod, Pdata)
	call dscal(Nperiod, sigep, Pdata, 1)
		

	! generate Pdata

	do time = 2, Nperiod
		Pdata(time) = min(max(rhop*Pdata(time-1) + Pdata(time), pgrid(1)), pgrid(np))
	end do

	Pdata = dexp(Pdata)


	open(1, file='Output\Pdata.txt', status='replace')
	write(1, '(F12.6)') Pdata

end subroutine


!	----------------------------------------------------------------------
!	subroutine name: NextIdiosyncraticShock()
!
!	generate next period indiosyncratic productivity shocks.
!	----------------------------------------------------------------------


subroutine NextIdiosyncraticShock()

	implicit none

	integer j
	real(8) ushock


	do indiv = 1, Nindiv

		ushock = drnunf()

		do j = 1, nx

			if (ushock <= CtrX(Worker(indiv)%X,j)) then
				Worker(indiv)%X = j
				exit
			end if
		
		end do

	end do

end subroutine


end module
