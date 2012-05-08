!	----------------------------------------------------------------------
!	Module name : ValueModule.f90
!	----------------------------------------------------------------------	*/


module ValueModule

	use Globals
	use LinInterpModule

contains


subroutine SolveValueFunction()

	implicit none

	integer, parameter:: MAXITER = 10000
	real(8), parameter:: tol_cnvrg = 1.0D-3
	real(8) NewVF(na,nx), NewVFE(na,nx), NewVFN(na,nx), errV
	integer NewAS(na,nx), NewASE(na,nx), NewASN(na,nx), iter, errA


	! start iteration

	do iter = 1, MAXITER


		! evaluate expected value function
	
		PF = matmul(VF, transpose(trX))


		!	maximization procedure

		do ix = 1, nx
		do ia = 1, na

			hour = hbar
			income = (1.0+irate)*agrid(ia) + wage*hour*exgrid(ix)
			call DecisionRule(NewASE(ia,ix), NewVFE(ia,ix))

			hour = 0.0D0
			income = (1.0+irate)*agrid(ia)
			call DecisionRule(NewASN(ia,ix), NewVFN(ia,ix))

			if (NewVFE(ia,ix) >= NewVFN(ia,ix)) then
				HR(ia,ix) = hbar
				NewAS(ia,ix) = NewASE(ia,ix)
				NewVF(ia,ix) = NewVFE(ia,ix)
			else
				HR(ia,ix) = 0.0D0
				NewAS(ia,ix) = NewASN(ia,ix)
				NewVF(ia,ix) = NewVFN(ia,ix)
			end if

		end do
		end do


		!	calculate errors

		errV = maxval(dabs(NewVF - VF))
		errA = maxval(abs(NewAS - AS))


		! update value function and decision rules

		VFE = NewVFE
		VFN = NewVFN
		VF  = NewVF
		ASE = NewASE
		ASN = NewASN
		AS  = NewAS


		! convergence check

		if (errV < tol_cnvrg) then
			print*, "Individual Optimization is done at iteration " ,iter
			print*, " "
			exit 
		else if (mod(iter,50) == 0) then
			print*, "SolveValueFunction iteration ", iter
			print*, "errV = ", errV, "   errA = ", errA
			print*, " "
		end if

	end do


end subroutine



subroutine DecisionRule(ap2, vf2)


	implicit none

	real(8) vf1, vf2, vf3
	integer ap1, ap2, ap3, ubd
	logical ifstop


	ubd = locate(income, agrid)

	if (hour == hbar) then
		ap1 = max(ASE(ia,ix) - 1, 1)
		ap2 = ASE(ia,ix)
		ap3 = min(ASE(ia,ix) + 1, na)
	else
		ap1 = max(ASN(ia,ix) - 1, 1)
		ap2 = ASN(ia,ix)
		ap3 = min(ASN(ia,ix) + 1, na)
	end if

	if (ap3 > ubd) then
		ap3 = max(ubd, 1)
		ap2 = max(ap3 - 1, 1)
		ap1 = max(ap2 - 1, 1)
	end if

	ifstop = .true.


	do while (ifstop)

		! ap1, ap2 and ap2 are distinct

		if (ap1 < ap2 .and. ap2 < ap3) then

			vf1 = ValueF(ap1)
			vf2 = ValueF(ap2)
			vf3 = ValueF(ap3)

			if (vf1 < vf2 .and. vf2 > vf3) then
				ifstop = .false.
			else if (vf1 < vf2 .and. vf2 < vf3) then
				ap1 = ap2
				ap2 = min(ap1 + 1, ubd)
				ap3 = min(ap2 + 1, ubd)
			else if (vf1 > vf2 .and. vf2 > vf3) then
				ap3 = ap2
				ap2 = max(ap3 - 1, 1)
				ap1 = max(ap2 - 1, 1)
			else
				ap3 = ap2
				ap2 = max(ap3 - 1, 1)
				ap1 = max(ap2 - 1, 1)
			end if

		! ap1 = ap2 < ap2: in this case ap1 = ap2 = 1

		else if (ap1 == ap2 .and. ap2 < ap3) then

			vf2 = ValueF(ap2)
			vf3 = ValueF(ap3)
			
			if (vf2 > vf3) then
				ifstop = .false.
			else 
				ap2 = ap1 + 1
				ap3 = ap2 + 1
			end if

		! ap1 < ap2 = ap2: in this case ap2 = ap3 = ubd

		else if (ap1 < ap2 .and. ap2 == ap3) then

			vf1 = ValueF(ap1)
			vf2 = ValueF(ap2)
			
			if (vf1 < vf2) then
				ifstop = .false.
			else 
				ap2 = ap3 - 1
				ap1 = ap2 - 1
			end if

		! ap1 = ap2 = ap2: either 1 or ubd

		else if (ap1 == ap2 .and. ap2 == ap3) then

			vf2 = ValueF(ap2)
			ifstop = .false.

		end if

	end do


end subroutine



real(8) function Utility(ja)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) cm, uh
	integer ja

	cm = max(income - agrid(ja), tiny)
	uh = hour**(1.0D0+1.0D0/gama)/(1.0D0+1.0D0/gama)

	Utility = log(cm) - B*uh

end function


real(8) function ValueF(ja)

	implicit none

	integer ja

	ValueF = Utility(ja) + bta*PF(ja, ix)

end function


end module
