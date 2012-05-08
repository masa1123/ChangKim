!	----------------------------------------------------------------------
!	Module name : ValueModule.f90
!	----------------------------------------------------------------------	*/


module ValueModule

	use Globals
	use LinInterpModule
	use OptimizationMod
	use Numerical_Libraries

	implicit none

	real(8):: break(na), cscoef(4,na)


contains


subroutine SolveValueFunction()

	implicit none

	integer, parameter:: MAXITER = 10000
	real(8), parameter:: tol_cnvrg = 1.0D-4
	real(8) NewVF(na,nx), NewAS(na,nx), NewHR(na,nx), errV, errA, errH,		&
			lbd, ubd, xmin
	integer iter


	! start iteration

	do iter = 1, MAXITER


		! evaluate expected value function
	
		PF = matmul(VF, transpose(trX))


		!	maximization procedure

		do ix = 1, nx

			call dcsint(na, agrid, PF(:,ix), break, cscoef)

			do ia = 1, na

				if (ia == 1) then
					if (ix == 1) then
						lbd = amin
					else
						lbd = AS(ia,ix-1)
					end if
				else
					lbd = AS(ia-1,ix)
				end if
					
				if (ia == na) then
					ubd = amax
				else
					ubd = AS(ia+1,ix)
				end if

				call fmin(NegValueF, lbd, ubd, xmin)

				NewHR(ia,ix) = hour
				NewAS(ia,ix) = xmin
				NewVF(ia,ix) = ValueF(NewAS(ia,ix))

			end do
	
		end do


		!	calculate errors

		errV = maxval(dabs(NewVF - VF))
		errA = maxval(dabs(NewAS - AS))
		errH = maxval(dabs(NewHR - HR))


		! update value function and decision rules

		VF  = NewVF
		AS  = NewAS
		HR  = NewHR


		! convergence check

		if (errV < tol_cnvrg) then
			write(*, '(A,I5,/)') "Individual Optimization is done at iteration " ,iter
			exit 
		else if (mod(iter,20) == 0) then
			write(*, '(A,I5)') "SolveValueFunction iteration ", iter
			write(*, '(A,F10.7,A,F10.7,A,F10.7,/)')		&
				"errV = ", errV, "   errA = ", errA, "   errH = ", errH
		end if

	end do


end subroutine



real(8) function SolveHour(xhour)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) cons, xhour

	cons = max((1.0D0+irate)*agrid(ia) + wage*xhour*exgrid(ix) - ap_hour, tiny)
	SolveHour = wage*exgrid(ix)/cons - B*xhour**(1.0D0/gama)

end function



real(8) function Utility(ap)

	implicit none

	real(8), parameter:: tiny = 1.0D-20, tol_hour = 1.0D-5
	real(8) ap, cm, ls, fval1, fval2, hour1, hour2
	integer i_bisc
	
	ap_hour = ap

	hour1 = HR(ia,ix)
	fval1 = SolveHour(hour1)

	if (dabs(fval1) <= tol_hour) then
	
		hour = HR(ia,ix)
	
	else

		do i_bisc = 1, 100

			if (fval1 > 0.0D0) then
				hour2 = 1.1D0*hour1
			else
				hour2 = 0.9D0*hour1
			end if

			fval2 = SolveHour(hour2)

			if (fval1*fval2 < 0.0D0) then
				hour = bisect(SolveHour, hour1, hour2, tol_hour)
				exit
			else
				hour1 = hour2
				fval1 = fval2
			end if

		end do

		if (i_bisc >= 100) then
			print *, "Hour wasn't found for (ia, ix) = ", ia, ix, HR(ia,ix) 
		end if

	end if

	cm = max((1.0D0+irate)*agrid(ia) + wage*hour*exgrid(ix) - ap, tiny)
	ls = B*hour**(1.0D0+1.0D0/gama)/(1.0D0+1.0D0/gama)

	Utility = log(cm) - ls

end function



real(8) function ValueF(ap)

	implicit none

	real(8) ap

	ValueF = Utility(ap) + bta*dcsval(ap, na-1, break, cscoef)

end function


real(8) function NegValueF(ap)

	implicit none

	real(8) ap

	NegValueF = -ValueF(ap)

end function


end module
