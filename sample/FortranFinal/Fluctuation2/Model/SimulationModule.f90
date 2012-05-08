!	----------------------------------------------------------------------
!	Module name: SimulationModule.f90
!
!	collection of sub-programs to simulate artificial data
!	----------------------------------------------------------------------



Module SimulationModule

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use optimizationMod
	use Numerical_Libraries
	
	implicit none

	integer, parameter:: xseed = 2
	real(8), dimension(:,:,:), allocatable:: ASEax, ASNax
	real(8), dimension(:,:), allocatable:: ASEkp, ASNkp, VFEkp, VFNkp
	real(8):: vf_work, vf_home, PreviousHour
	real(8) breakE(nx), cscoefE(4,nx), breakN(nx), cscoefN(4,nx)
	integer:: Nindiv
	type Individual
		real(8) :: A, AP, H, W, C
		integer :: X
	end type Individual
	type(Individual), dimension(:), allocatable :: Worker

contains




!	----------------------------------------------------------------------
!	File name: SimulateData.f90
!
!	generate artificial data for all aggregate variables, which will be 
!	used in ModelStatistics().
!	----------------------------------------------------------------------


subroutine SimulateData()

	implicit none


	!	number of inidividuals for simulation

	if (final) then
		Nindiv = Nindiv_Fin_Sim
	else
		Nindiv = Nindiv_Sim
	end if


	!	initialize time series data

	Kdata = 0.0D0
	K2data = 0.0D0
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
	allocate(ASEax(nk,nk2,np), ASNax(nk,nk2,np), ASEkp(na,nx), ASNkp(na,nx), VFEkp(na,nx), VFNkp(na,nx))


	!	intialize distribution of workers

	do indiv = 1, Nindiv
		Worker(indiv)%A = kss
		Worker(indiv)%H = hbar
		Worker(indiv)%X = mod(indiv, nx) + 1
	end do


	!	initial aggregate capital

	Kdata(1) = sum(Worker%A)/Nindiv
	K2data(1) = sum(Worker%A**2.0D0)/Nindiv
	

	!	start generating artificial time series data

	do time = 1, Nperiod


		! interpolate value function and asset decision rule at (Kapital, pt)

		do ix = 1, nx
		do ia = 1, na

			!	interpolate value functions

			VFEkp(ia,ix) = polyinterp3(Kdata(time), K2data(time), dlog(Pdata(time)), kgrid, k2grid, pgrid, VFE(ia,ix,:,:,:))
			VFNkp(ia,ix) = polyinterp3(Kdata(time), K2data(time), dlog(Pdata(time)), kgrid, k2grid, pgrid, VFN(ia,ix,:,:,:))

			! temporary matrix for interpolation of decision rules for asset

			do ik = 1, nK
			do ik2 = 1, nK2
			do ip = 1, np
				ASEax(ik,ik2,ip) = agrid(ASE(ia,ix,ik,ik2,ip))
				ASNax(ik,ik2,ip) = agrid(ASN(ia,ix,ik,ik2,ip))
			end do
			end do
			end do

			!	interpolate decision rules

			ASEkp(ia,ix) = polyinterp3(Kdata(time), K2data(time), dlog(Pdata(time)), kgrid, k2grid, pgrid, ASEax)
			ASNkp(ia,ix) = polyinterp3(Kdata(time), K2data(time), dlog(Pdata(time)), kgrid, k2grid, pgrid, ASNax)

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


		!	generate time series data for aggregate macro variables

		Kdata(time+1) = min(max(sum(Worker%AP)/Nindiv, kgrid(1)), kgrid(nk))
		K2data(time+1) = min(max(sum(Worker%A**2.0D0)/Nindiv, k2grid(1)), k2grid(nk2))
!		K2data(time+1) = dsqrt(sum(Worker%A**2.0D0)/Nindiv - Kdata(time)**2.0D0)
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


		!	prepare for the simulation in the next period

		Worker%A = Worker%AP
		call NextIdiosyncraticShock()


		!	show how simulation is going on

		if (mod(time, 100) == 0) then
			write(*,'(A,I6)')   "time = ", time	
			write(*,'(A,F10.5)') "p =   ", Pdata(time)
			write(*,'(A,F10.5)') "K =   ", Kdata(time)
			write(*,'(A,F10.5)') "Kp =  ", Kdata(time+1)
			write(*,'(A,F10.5)') "K2 =  ", K2data(time)
			write(*,'(A,F10.5)') "K2p = ", K2data(time+1)
			write(*,'(A,F10.5)') "L =   ", Ldata(time)
			write(*,'(A,F10.5)') "w =   ", Wdata(time)
			write(*,'(A,F10.5)') "r =   ", Rdata(time)-delta
			write(*,'(/)')
		end if


	end do


	!	deallocate memory

	deallocate(Worker)
	deallocate(ASEax, ASNax, ASEkp, ASNkp, VFEkp, VFNkp)
	
	close(100)
	close(200)
	close(300)


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
