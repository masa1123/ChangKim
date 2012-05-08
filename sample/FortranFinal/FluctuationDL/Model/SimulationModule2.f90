!	----------------------------------------------------------------------
!	Module name: SimulationModule.f90
!
!	collection of sub-programs to simulate artificial data
!	----------------------------------------------------------------------



Module SimulationModule2

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use Numerical_Libraries
	
	implicit none

	integer, parameter:: xseed = 2
	real(8), dimension(:,:), allocatable:: ASEkp, ASNkp, VFEkp, VFNkp
	real(8) vf_work, vf_home, PreviousHour
	type Individual
		real(8) :: A, AP, H
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
	allocate(ASEkp(na,nx), ASNkp(na,nx), VFEkp(na,nx), VFNkp(na,nx))


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


		! interpolate value function and asset decision rule at (Kapital, pt)

		do ix = 1, nx
		do ia = 1, na

			VFEkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	VFE(ia,ix,:,iPdata(time)))
			VFNkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	VFN(ia,ix,:,iPdata(time)))

			ASEkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	agrid(ASE(ia,ix,:,iPdata(time))))
			ASNkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	agrid(ASN(ia,ix,:,iPdata(time))))

		end do
		end do


		!	cross-section data for asset and hours given (Kapital, pt)

		do indiv = 1, Nindiv


			vf_work = lininterp1(Worker(indiv)%A, agrid, VFEkp(:,Worker(indiv)%X))
			vf_home = lininterp1(Worker(indiv)%A, agrid, VFNkp(:,Worker(indiv)%X))


			if (vf_work >= vf_home) then

				Worker(indiv)%H  = hbar
				Worker(indiv)%AP = lininterp1(Worker(indiv)%A, agrid, ASEkp(:, Worker(indiv)%X))
				Worker(indiv)%AP = min(max(Worker(indiv)%AP, agrid(1)), agrid(na))
				Ldata(time) = Ldata(time) + exgrid(Worker(indiv)%X)*Worker(indiv)%H

			else

				Worker(indiv)%H  = 0.0D0
				Worker(indiv)%AP = lininterp1(Worker(indiv)%A, agrid, ASNkp(:, Worker(indiv)%X))
				Worker(indiv)%AP = min(max(Worker(indiv)%AP, agrid(1)), agrid(na))

			end if


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
	deallocate(ASEkp, ASNkp, VFEkp, VFNkp)


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
	EINdata = 0.0D0
	EOUTdata = 0.0D0


	!	set the seed number for random number generator

	call rnset(xseed)


	!	allocate space for cross section data generation

	allocate(Worker(Nindiv))
	allocate(ASEkp(na,nx), ASNkp(na,nx), VFEkp(na,nx), VFNkp(na,nx))


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


		! interpolate value function and asset decision rule at (Kapital, pt)

		do ix = 1, nx
		do ia = 1, na

			VFEkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	VFE(ia,ix,:,iPdata(time)))
			VFNkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	VFN(ia,ix,:,iPdata(time)))

			ASEkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	agrid(ASE(ia,ix,:,iPdata(time))))
			ASNkp(ia,ix) = polyinterp1(Kdata(time), kgrid,	agrid(ASN(ia,ix,:,iPdata(time))))

		end do
		end do


		!	cross-section data for asset and hours given (Kapital, pt)

		do indiv = 1, Nindiv


			vf_work = lininterp1(Worker(indiv)%A, agrid, VFEkp(:,Worker(indiv)%X))
			vf_home = lininterp1(Worker(indiv)%A, agrid, VFNkp(:,Worker(indiv)%X))
			PreviousHour = Worker(indiv)%H


			if (vf_work >= vf_home) then

				Worker(indiv)%H  = hbar
				Worker(indiv)%AP = lininterp1(Worker(indiv)%A, agrid, ASEkp(:, Worker(indiv)%X))
				Worker(indiv)%AP = min(max(Worker(indiv)%AP, agrid(1)), agrid(na))
				Ldata(time) = Ldata(time) + exgrid(Worker(indiv)%X)*Worker(indiv)%H
				Edata(time) = Edata(time) + 1.0D0
				if (PreviousHOur == 0.0D0) EINdata(time) = EINdata(time) + 1

			else

				Worker(indiv)%H  = 0.0D0
				Worker(indiv)%AP = lininterp1(Worker(indiv)%A, agrid, ASNkp(:, Worker(indiv)%X))
				Worker(indiv)%AP = min(max(Worker(indiv)%AP, agrid(1)), agrid(na))
				if (PreviousHour == hbar) EOUTdata(time) = EOUTdata(time) + 1

			end if


		end do


		!	prepare for the simulation in the next period

		Kdata(time+1) = min(max(sum(Worker%AP)/Nindiv, kgrid(1)), kgrid(nk))
		Worker%A  = Worker%AP
		call NextIdiosyncraticShock()


		!	generate time series data for aggregate macro variables

		Ldata(time) = Ldata(time)/Nindiv
		Rdata(time) = Pdata(time)*alpha*(Kdata(time)/Ldata(time))**(alpha-1.0D0)
		Wdata(time) = Pdata(time)*(1.0D0-alpha)*(Kdata(time)/Ldata(time))**(alpha)
		Edata(time) = Edata(time)/Nindiv
		EINdata(time) = EINdata(time)/Nindiv
		EOUTdata(time) = EOUTdata(time)/Nindiv
		hEINdata(time) = EINdata(time)/Edata(time)
		hEOUTdata(time) = EOUTdata(time)/Edata(time)
		Ydata(time) = Pdata(time)*(Kdata(time)**alpha)*(Ldata(time)**(1.0D0-alpha))
		Idata(time) = Kdata(time+1) - (1.0D0-delta)*Kdata(time)
		Cdata(time) = Ydata(time) - Idata(time)
		Adata(time) = Ydata(time)/(hbar*Edata(time))


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
	deallocate(ASEkp, ASNkp, VFEkp, VFNkp)


end subroutine


!	----------------------------------------------------------------------
!	subroutine name: AggregateShockSeries()
!
!	generate time series of aggregate productivity shocks.
!	----------------------------------------------------------------------


subroutine AggregateShockSeries()

	implicit none

	integer, parameter:: pseed = 1
	integer jp, kp
	real(8) ushock


	! set the seed number for random number generator

	call rnset(pseed)


	! generate Zdata

	jp = np/2 + 1

	do time = 1, Nperiod

		ushock = drnunf()

		do kp = 1, np

			if (ushock <= CtrP(jp,kp)) then
				jp = kp
				exit
			end if

		end do

		iPdata(time) = jp
		Pdata(time) = epgrid(jp)

	end do

	open(1, file='Output\Pdata.txt', status='replace')
	write(1, '(f12.6)') Pdata


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
