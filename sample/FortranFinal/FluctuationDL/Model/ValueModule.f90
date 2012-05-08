!	----------------------------------------------------------------------
!	module name: ValueModule.f90
!
!	collection of sub-programs to calculate value function and decision 
!	rule for the next period asset holdings.
!	----------------------------------------------------------------------


module ValueModule

	use Globals
	use LinInterpModule
	use PolyInterpModule
	use OptimizationMod
	use Numerical_Libraries

	implicit none

	real(8):: break(na), cscoef(4,na), Value


contains



!	----------------------------------------------------------------------
!	subroutine name: SolveValueFunction()
!	
!	calcualtes value functions and asset holdings for individuals.
!	----------------------------------------------------------------------


subroutine SolveValueFunction()

	implicit none

	real(8) NewVF(na,nx,nk,np), NewAS(na,nx,nk,np), NewHR(na,nx,nk,np),		&
			errV, errA, errH, lbd, ubd, xmin
	integer iter


	! start iteration

	do iter = 1, MAXITER


		!	evaluate expected value function

		call ExpectedValueFunction()


		!	maximization procedure

		do ip = 1, np
		do ik = 1, nk

			irate = dexp(Rcoef(1) + Rcoef(2)*dlog(kgrid(ik)) + Rcoef(3)*pgrid(ip)) - delta
			wage  = dexp(Wcoef(1) + Wcoef(2)*dlog(kgrid(ik)) + Wcoef(3)*pgrid(ip))

			do ix = 1, nx

				call dcsint(na, agrid, PF(:,ix,ik,ip), break, cscoef)

				do ia = 1, na

					if (ia == 1) then
						if (ix == 1) then
							lbd = amin
						else
							lbd = AS(ia,ix-1,ik,ip)
						end if
					else
						lbd = AS(ia-1,ix,ik,ip)
					end if
				
					if (ia == na) then
						ubd = amax
					else
						ubd = AS(ia+1,ix,ik,ip)
					end if

					call fmin(NegValueF, lbd, ubd, xmin)

					NewHR(ia,ix,ik,ip) = hour
					NewAS(ia,ix,ik,ip) = xmin
					NewVF(ia,ix,ik,ip) = Value

				end do

			end do

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

		if (errV < tol_value) then
			write(*, '(A,I5,/)') "Individual Optimization is done at iteration " ,iter
			exit 
		else if (mod(iter,showerr) == 0) then
			write(*, '(A,I5)') "SolveValueFunction iteration ", iter
			write(*, '(A,F10.7,A,F10.7,A,F10.7,/)')		&
				"errV = ", errV, "   errA = ", errA, "   errH = ", errH
		end if

	end do

	
end subroutine



!	----------------------------------------------------------------------
!	subroutine name: SolveHour(xhour)
!	
!	solve the optimal hours worked given a and a'
!	----------------------------------------------------------------------


real(8) function SolveHour(xhour)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) cons, xhour

	cons = max((1.0D0+irate)*agrid(ia) + wage*xhour*exgrid(ix) - ap_hour, tiny)
	SolveHour = wage*exgrid(ix)/cons - B*xhour**(1.0D0/gama)


end function




!	----------------------------------------------------------------------
!	subroutine name: Utility(ja)
!	
!	evaluates instantatneous utility at asset agrid(ja).
!	----------------------------------------------------------------------


real(8) function Utility(ap)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) ap, cm, ls, fval1, fval2, hour1, hour2
	integer i_bisc

	ap_hour = ap

	hour1 = HR(ia,ix,ik,ip)
	fval1 = SolveHour(hour1)

	if (dabs(fval1) <= tol_hour) then
	
		hour = HR(ia,ix,ik,ip)
	
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
			print *, "Hour wasn't found for (ia,ix,ik,ip) = ", ia, ix, ik, ip
		end if

	end if

	cm = max((1.0D0+irate)*agrid(ia) + wage*hour*exgrid(ix) - ap, tiny)
	ls = B*hour**(1.0D0+1.0D0/gama)/(1.0D0+1.0D0/gama)

	Utility = log(cm) - ls


end function



!	----------------------------------------------------------------------
!	subroutine name: ValueF(ap)
!	
!	evaluates value at asset agrid(ap).
!	----------------------------------------------------------------------


real(8) function ValueF(ap)

	implicit none

	real(8) ap

	Value = Utility(ap) + bta*dcsval(ap, na-1, break, cscoef)
	ValueF = Value


end function


real(8) function NegValueF(ap)

	implicit none

	real(8) ap

	NegValueF = -ValueF(ap)

end function



!	----------------------------------------------------------------------
!	subroutine name: ExpectedValueFunction()
!	
!	calcualtes expected value functions.
!	----------------------------------------------------------------------


subroutine ExpectedValueFunction()

	implicit none

	real(8) kp, VFP(na,nx,np)
	integer ipp, ixp


	PF = 0
	
	do ik = 1, nk
	do ip = 1, np


		!	next period aggregate capital

		kp = dexp(Kcoef(1) + Kcoef(2)*dlog(kgrid(ik)) + Kcoef(3)*pgrid(ip))


		do ipp = 1, np


			!	next period value function: VFP(ia,ixzp,ipp) with kp(ik,ip)

			do ixp = 1, nx
			do ia = 1, na

				VFP(ia,ixp,ipp) = polyinterp1(kp, kgrid, VF(ia,ixp,:,ipp))

			end do
			end do

		end do

		
		!	calculate expected value function: PF(ia,ixp,ik,ip)

		do ix = 1, nx
		do ia = 1, na

			do ipp = 1, np
			do ixp = 1, nx

				PF(ia,ix,ik,ip) = PF(ia,ix,ik,ip)		&
								 + VFP(ia,ixp,ipp)*trX(ix,ixp)*trP(ip,ipp)

			end do
			end do

		end do
		end do

	end do
	end do

end subroutine


end module
