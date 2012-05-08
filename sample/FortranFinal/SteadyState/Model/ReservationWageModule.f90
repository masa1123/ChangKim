!	----------------------------------------------------------------------
!	File name: ReservationWageModule.f90
!
!	This routine solves for reservation wage for each asset level at 
!	which an individual with the asset is indifferent between working 
!	and nonworking.
!	----------------------------------------------------------------------


module ReservationWageModule


	use Numerical_Libraries
	use Globals
	use OptimizationMod

	implicit none

	real(8) breakE(nx), cscoefE(4,nx), breakN(nx), cscoefN(4,nx)


	contains


	subroutine FindReservationWages()

		ReservationX = 0.0d0

		do ia = 1,na

			!	check if VFE and VFN intersect

			if (HR(ia,nx) == 0.0d0) then
				ReservationX(ia) = xgrid(nx)

			elseif (HR(ia,1) > 0.0d0) then
				ReservationX(ia) = xgrid(1)
	
			else

				!	spline value functions VFE(ka,kx,:) and VFN(ka,kx,:)

				call dcsint(nx, xgrid, VFE(ia,1:nx), breakE, cscoefE)
				call dcsint(nx, xgrid, VFN(ia,1:nx), breakN, cscoefN)

				!	solve for x^* such that VFE(ia,x*) = VFN(ia,x*)

				call fmin(DiffValue, xgrid(1), xgrid(nx), ReservationX(ia))

			end if

		end do

		ReservationWage = dexp(ReservationX)*wage


	end subroutine



	real(8) function DiffValue(x)

		implicit none

		real(8), intent(in):: x


		DiffValue = dabs(dcsval(x, nx-1, breakE, cscoefE) - dcsval(x, nx-1, breakN, cscoefN))

	end function


	subroutine LaborSupply()

		real(8) break_cmu(na,nx), cscoef_cmu(na,4,nx)
		integer ja

		do ia = 1, na

			call dcsint(nx, xgrid, cmu(ia,:), break_cmu(ia,:), cscoef_cmu(ia,:,:))

		end do

		do ia = 1, na
		do ja = 1, na

			Reserv_cmu(ia,ja) = 100.0D0*dcsval(ReservationX(ia), nx-1, break_cmu(ja,:), cscoef_cmu(ja,:,:))

		end do
		end do

		do ia = 1, na

			LS(ia) = Reserv_cmu(ia,ia)

		end do


	end subroutine


end module
