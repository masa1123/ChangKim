!	----------------------------------------------------------------------
!	Module name: MeasureModule.f90
!
!	Finds the invariant measure given the individual decision rule and 
!	the transition process of idiosyncratic shocks.
!	----------------------------------------------------------------------


module MeasureModule

	use Globals
	use GridsModule
	use LinInterpModule
	use Numerical_Libraries

	implicit none

	real(8), parameter:: tol_measure = 1.0D-8
	integer, parameter:: maxiter = 1E6
	real(8) Fmu(namu,nx), invAS(namu,nx)


contains



subroutine InitializeMeasure()

	implicit none

	do ix = 1, nx
	do ia = 1, namu

		cmu(ia,ix) = (agridmu(ia)-agridmu(1))/(agridmu(namu)-agridmu(1))*prX(ix)

	end do
	end do

	open(1, file='Output\cmu.txt', status='unknown')
	write(1, '(<nx>f12.6)') ((cmu(ia,ix),ix=1,nx),ia=1,namu)
	close(1)


end subroutine



subroutine InverseAssetDecision()

	implicit none


	do ix = 1, nx
	do ia = 1, namu

		if (agridmu(ia) <= AS(na,ix) .and. agridmu(ia) >= AS(1,ix)) then
			invAS(ia,ix) = lininterp1(agridmu(ia), AS(:,ix), agrid)
		else if (agridmu(ia) > AS(na,ix)) then
			invAS(ia,ix) = agridmu(namu)
		else
			invAS(ia,ix) = agridmu(1)
		end if

	end do
	end do

	open(1, file='Output\invAS.txt', status='unknown')
	write(1, '(<nx>f12.6)') ((invAS(ia,ix),ix=1,nx),ia=1,namu)


end subroutine



subroutine InvariantMeasure()

	implicit none

	real(8):: err
	integer:: iter, jx
	real(8):: break(namu,nx), cscoef(4,namu,nx)


	!	iconstruct asset grids for invariant measure

	call CapitalGridMu()


	!	initialize measure

	call InitializeMeasure()


	!	calcualte inverse function of asset decision rule

	call InverseAssetDecision()


	!	start iteration for distribution (cumulative measure)

	do iter = 1, maxiter


		!	cubic spline of cumulative measures for interpolation

!		do ix = 1, nx
!			call dcsint(namu, agridmu, cmu(:,ix), break(:,ix), cscoef(:,:,ix))
!		end do


		!	calculate measure before shocks change

		Fmu = 0.0D0

		do ix = 1, nx
		do ia = 1, namu

			do jx = 1, nx
				Fmu(ia,ix) = Fmu(ia,ix)	+ trX(jx,ix)	&
				*lininterp1(invAS(ia,jx), agridmu, cmu(:,jx))
				!*dcsval(invAS(ia,jx), namu-1, break(:,jx), cscoef(:,:,jx))
				
			end do

		end do
		end do


		!	calculate the differences between old and new measures

		err = maxval(dabs(Fmu-cmu))


		!	update measures 

		cmu = Fmu


		!	convergence test of the measures

		if (err < tol_measure) then
			print*, "Convergence achieved at InvariantMeasure iteration", iter
			print*, ""
			exit
		else if (mod(iter, 100) == 0) then
			print*, "InvariantMeasure iteration", iter
			print*, "err = ", err
			print*, ""
		end if

	end do


	!	shouldn't reach here

	if (iter > maxiter) then
		print*, "Convergence failed in InvariantMeasure...."
		stop
	end if


	!	measure (density)

	do ia = 1, namu-1

		mu(ia,:) = cmu(ia+1,:) - cmu(ia,:)

	end do
	

end subroutine




end module
