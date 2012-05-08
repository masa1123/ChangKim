!	----------------------------------------------------------------------
!	Module name: MeasureModule.f90
!
!	Finds the invariant measure given the individual decision rule and 
!	the transition process of idiosyncratic shocks.
!	----------------------------------------------------------------------


module MeasureModule

	use Globals

	implicit none


contains

	
subroutine InvariantMeasure()

	use Globals

	implicit none

	real(8), parameter:: tol_measure = 1.0D-8
	integer, parameter:: maxiter = 1E6
	real(8) Newmu(na,nx), Fmu(na,nx), err
	integer iter, i, k


	!	start iteration

	do iter = 1, maxiter


		!	calculate measure before shocks change

		Fmu = 0.0D0

		do k = 1, nx
		do i = 1, na

			Fmu(AS(i,k),k) = Fmu(AS(i,k),k) + mu(i,k)

		end do
		end do

		
		!	next period measures

		Newmu = 0.0D0
		Newmu = matmul(Fmu, trX)
		Newmu = Newmu/sum(Newmu)

		!	calculate the differences between old and new measures

		err = maxval(dabs(Newmu-mu))


		!	update measures 

		mu = Newmu


		!	convergence test of the measures

		if (err < tol_measure) then
			print*, "Convergence achieved at InvariantMeasure iteration", iter
			print*, "err = ", err
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


	!	cumulative measure

	cmu = 0.0D0
	
	do k = 1, nx
	
		cmu(1,k) = mu(1,k)
	
		do i = 2, na
			cmu(i,k) = cmu(i-1,k) + mu(i,k)
		end do

		if (k >= 2) then
			cmu(:,k) = cmu(:,k-1) + cmu(:,k)
		end if

	end do


end subroutine




subroutine AssetEarningsDistribution()

	implicit none


	!	joint density for asset and earnings 

	AsstEarnDensity = 0.0D0
	
	do ia = 1,na
	do ix = 1,nx
		
		if (HR(ia,ix) == 0.0D0) then
			AsstEarnDensity(ia,0) = AsstEarnDensity(ia,0) + mu(ia,ix)
		else
			AsstEarnDensity(ia,ix) = AsstEarnDensity(ia,ix) + mu(ia,ix)
		end if

	end do
	end do


	!	marginal density for asset

	AsstDensity = sum(AsstEarnDensity,dim=2)


	!	marginal density for earnings

	EarnDensity = sum(AsstEarnDensity,dim=1)


	!	joint distribution for asset and earnings

	AsstEarnDist = AsstEarnDensity

	do ix = 1,nx
		AsstEarnDist(:,ix) = AsstEarnDist(:,ix-1) + AsstEarnDensity(:,ix)
	end do

	do ia = 2,na
		AsstEarnDist(ia,:) = AsstEarnDist(ia-1,:) + AsstEarnDist(ia,:)
	end do


	!	marginal distribution for asset

	AsstDist = AsstDensity

	do ia = 2,na
		AsstDist(ia) = AsstDist(ia-1) + AsstDensity(ia)
	end do


	!	marginal distribution for earnings

	EarnDist = EarnDensity

	do ix = 1,nx
		EarnDist(ix) = EarnDist(ix-1) + EarnDensity(ix)
	end do


end subroutine




end module
