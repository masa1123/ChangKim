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
	real(8), dimension(:,:), allocatable:: ASEax, ASNax, ASEkp, ASNkp, VFEkp, VFNkp
	real(8):: vf_work, vf_home, PreviousHour
	real(8) ResX(na), ResW(na), mu_ax(na,nx)
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
!	generate artificial data for aggregate capital, interest rate, and 
!	wage rate through simulating individuals' decision rules. These data
!	will be used in RegressLOM().
!	----------------------------------------------------------------------


subroutine SimulateData()

	implicit none


	!	number of inidividuals for simulation

	Nindiv = Nindiv_Sim


	!	initialize time series data

	Kdata = 0.0D0
	Ldata = 0.0D0
	Rdata = 0.0D0
	Wdata = 0.0D0


	!	set the seed number for random number generator

	call rnset(xseed)


	!	allocate space for cross section data generation

	allocate(Worker(Nindiv))
	allocate(ASEax(nk,np), ASNax(nk,np), ASEkp(na,nx), ASNkp(na,nx), VFEkp(na,nx), VFNkp(na,nx))


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

			!	interpolate value functions

			VFEkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, VFE(ia,ix,:,:))
			VFNkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, VFN(ia,ix,:,:))

			! temporary matrix for interpolation of decision rules for asset

			do ik = 1, nK
			do ip = 1, np
				ASEax(ik,ip) = agrid(ASE(ia,ix,ik,ip))
				ASNax(ik,ip) = agrid(ASN(ia,ix,ik,ip))
			end do
			end do

			!	interpolate decision rules

			ASEkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, ASEax)
			ASNkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, ASNax)

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
	deallocate(ASEax, ASNax, ASEkp, ASNkp, VFEkp, VFNkp)


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


	!	number of inidividuals for simulation

	Nindiv = Nindiv_Fin_Sim


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
	allocate(ASEax(nk,np), ASNax(nk,np), ASEkp(na,nx), ASNkp(na,nx), VFEkp(na,nx), VFNkp(na,nx))


	!	intialize distribution of workers

	do indiv = 1, Nindiv
		Worker(indiv)%A = kss
		Worker(indiv)%H = hbar
		Worker(indiv)%X = mod(indiv, nx) + 1
	end do


	!	initial aggregate capital

	Kdata(1) = sum(Worker%A)/Nindiv
	

	!	open files

	open(100, file='Output\Panel.txt', status='unknown')
	open(200, file='Output\ResX.txt', status='unknown')
	open(300, file='Output\mu_ax.txt', status='unknown')


	!	start generating artificial time series data

	do time = 1, Nperiod


		! interpolate value function and asset decision rule at (Kapital, pt)

		do ix = 1, nx
		do ia = 1, na

			!	interpolate value functions

			VFEkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, VFE(ia,ix,:,:))
			VFNkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, VFN(ia,ix,:,:))

			! temporary matrix for interpolation of decision rules for asset

			do ik = 1, nK
			do ip = 1, np
				ASEax(ik,ip) = agrid(ASE(ia,ix,ik,ip))
				ASNax(ik,ip) = agrid(ASN(ia,ix,ik,ip))
			end do
			end do

			!	interpolate decision rules

			ASEkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, ASEax)
			ASNkp(ia,ix) = polyinterp2(Kdata(time), dlog(Pdata(time)), kgrid, pgrid, ASNax)

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
		K2data(time) = dsqrt(sum(Worker%A**2.0D0)/Nindiv - Kdata(time)**2.0D0)
!		K2data(time) = sum(Worker%A**2.0D0)/Nindiv
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


		!	construct panel data, measure and compute reservaton wages if panel == .true.

		if (((time >= pstart1 .and. time <= pend1) .or. (time >= pstart2 .and. time <= pend2)) .and. panel == .true.) then
		

			mu_ax = 0.0D0

			do indiv = 1, Nindiv

				!	panel

				Worker(indiv)%W = Wdata(time)*exgrid(Worker(indiv)%X)*Worker(indiv)%H
				Worker(indiv)%C = (1.0D0+Rdata(time)-delta)*Worker(indiv)%A		&
									+ Worker(indiv)%W - Worker(indiv)%AP
				write(100, '(F12.6,I6,2F12.6)') Worker(indiv)%A, Worker(indiv)%X, Worker(indiv)%H, Worker(indiv)%C


				!	measure

				if (Worker(indiv)%A < (agrid(1) + agrid(2))/2.0D0) then
					mu_ax(1,Worker(indiv)%X) = mu_ax(1,Worker(indiv)%X) + 1.0D0
				else if (Worker(indiv)%A >= (agrid(na-1) + agrid(na))/2.0D0) then
					mu_ax(na,Worker(indiv)%X) = mu_ax(na,Worker(indiv)%X) + 1.0D0
				else 

					do ia = 2, na-1
						if ((Worker(indiv)%A < (agrid(ia) + agrid(ia+1))/2.0D0) .and.	&
							(Worker(indiv)%A >= (agrid(ia-1) + agrid(ia))/2.0D0)) then
							mu_ax(ia,Worker(indiv)%X) = mu_ax(ia,Worker(indiv)%X) + 1.0D0
							exit
						end if
					end do

				end if

			end do
	

			!	reservation wages

			call FindReservationX()
			write(200, '(2F12.6)') (ResX(ia), ResW(ia), ia=1,na)


			!	write measure

			mu_ax = mu_ax/Nindiv*100.0D0
			write(300, '(<nx>F12.6)') ((mu_ax(ia,ix),ix=1,nx),ia=1,na)

		end if


		!	prepare for the simulation in the next period

		Worker%A = Worker%AP
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



!	----------------------------------------------------------------------
!	subroutine name: FindReservationX()
!
!	compute the reservation productivity.
!	----------------------------------------------------------------------


subroutine FindReservationX()

	implicit none


	!	compute the reservation productivity

	ResX = 0.0d0

	do ia = 1,na

		!	check if VFEkp and VFNkp intersect

		if (VFEkp(ia,nx) < VFNkp(ia,nx)) then
			ResX(ia) = xgrid(nx)

		elseif (VFEkp(ia,1) > VFNkp(ia,1)) then
			ResX(ia) = xgrid(1)

		else

			!	spline value functions VFEkp(ka,kx,:) and VFNkp(ka,kx,:)

			call dcsint(nx, xgrid, VFEkp(ia,1:nx), breakE, cscoefE)
			call dcsint(nx, xgrid, VFNkp(ia,1:nx), breakN, cscoefN)

			!	solve for x^* such that VFEkp(ia,x*) = VFNkp(ia,x*)

			call fmin(DiffValue, xgrid(1), xgrid(nx), ResX(ia))

		end if

	end do

	ResW = dexp(ResX)*Wdata(time)


end subroutine



real(8) function DiffValue(x)

	implicit none

	real(8), intent(in):: x


	DiffValue = dabs(dcsval(x, nx-1, breakE, cscoefE) - dcsval(x, nx-1, breakN, cscoefN))

end function



end module
