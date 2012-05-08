!	----------------------------------------------------------------------
!	module name: ValueModule.f90
!
!	collection of sub-programs to calculate value function and decision 
!	rule for the next period asset holdings.
!	----------------------------------------------------------------------


module ValueModule

	use Globals
	use PolyInterpModule
	use LinInterpModule

contains



!	----------------------------------------------------------------------
!	subroutine name: SolveValueFunction()
!	
!	calcualtes value functions and asset holdings for individuals.
!	----------------------------------------------------------------------

subroutine SolveValueFunction()

	implicit none

	real(8) errV
	integer errA
	real(8) NewVF(na,nx,nk,nk2,np), NewVFE(na,nx,nk,nk2,np), NewVFN(na,nx,nk,nk2,np)
	integer NewAS(na,nx,nk,nk2,np), NewASE(na,nx,nk,nk2,np), NewASN(na,nx,nk,nk2,np)



	! start iteration

	do iter = 1, MAXITER


		!	evaluate expected value function

		call ExpectedValueFunction()


		!	maximization procedure

		do ip = 1, np
		do ik2 = 1, nk2
		do ik = 1, nk

			irate = dexp(Rcoef(1) + Rcoef(2)*dlog(kgrid(ik)) + Rcoef(3)*dlog(k2grid(ik2)) + Rcoef(4)*pgrid(ip)) - delta
			wage  = dexp(Wcoef(1) + Wcoef(2)*dlog(kgrid(ik)) + Wcoef(3)*dlog(k2grid(ik2)) + Wcoef(4)*pgrid(ip))

			do ix = 1, nx
			do ia = 1, na

				hour = hbar
				income = (1.0+irate)*agrid(ia) + wage*hour*exgrid(ix)
				call DecisionRule(NewASE(ia,ix,ik,ik2,ip), NewVFE(ia,ix,ik,ik2,ip))

				hour = 0.0D0
				income = (1.0+irate)*agrid(ia)
				call DecisionRule(NewASN(ia,ix,ik,ik2,ip), NewVFN(ia,ix,ik,ik2,ip))

				if (NewVFE(ia,ix,ik,ik2,ip) >= NewVFN(ia,ix,ik,ik2,ip)) then
					HR(ia,ix,ik,ik2,ip) = hbar
					NewAS(ia,ix,ik,ik2,ip) = NewASE(ia,ix,ik,ik2,ip)
					NewVF(ia,ix,ik,ik2,ip) = NewVFE(ia,ix,ik,ik2,ip)
				else
					HR(ia,ix,ik,ik2,ip) = 0.0D0
					NewAS(ia,ix,ik,ik2,ip) = NewASN(ia,ix,ik,ik2,ip)
					NewVF(ia,ix,ik,ik2,ip) = NewVFN(ia,ix,ik,ik2,ip)
				end if

			end do
			end do

		end do
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

		if (errV < tol_value .and. errA == 0) then
			print*, "Individual Optimization is done at iteration " ,iter
			print*, " "
			exit 
		else if (mod(iter,showerr) == 0 .or. (iter <= 100 .and. iterLOM == 1)) then
			print*, "SolveValueFunction iteration ", iter
			print*, "errV = ", errV, "   errA = ", errA
			print*, " "
		end if

	end do

	
end subroutine


!	----------------------------------------------------------------------
!	subroutine name: DecisionRule.f90
!	
!	calcualtes individual's decision rule for asset accumulations.
!	----------------------------------------------------------------------

subroutine DecisionRule(ap2, vf2)

	implicit none

	real(8) vf1, vf2, vf3
	integer ap1, ap2, ap3, ubd
	logical ifstop


	ubd = locate(income, agrid)

	if (hour == hbar) then
		ap1 = max(ASE(ia,ix,ik,ik2,ip) - 1, 1)
		ap2 = ASE(ia,ix,ik,ik2,ip)
		ap3 = min(ASE(ia,ix,ik,ik2,ip) + 1, na)
	else
		ap1 = max(ASN(ia,ix,ik,ik2,ip) - 1, 1)
		ap2 = ASN(ia,ix,ik,ik2,ip)
		ap3 = min(ASN(ia,ix,ik,ik2,ip) + 1, na)
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


!	----------------------------------------------------------------------
!	subroutine name: Utility(ja)
!	
!	evaluates instantatneous utility at asset agrid(ja).
!	----------------------------------------------------------------------

real(8) function Utility(ja)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) cm, uh
	integer ja

	cm = max(income - agrid(ja), tiny)
	uh = hour**(1.0D0+1.0D0/gama)/(1.0D0+1.0D0/gama)

	Utility = log(cm) - B*uh

end function



!	----------------------------------------------------------------------
!	subroutine name: ValueF(ja)
!	
!	evaluates value at asset agrid(ja).
!	----------------------------------------------------------------------

real(8) function ValueF(ja)

	implicit none

	integer ja

	ValueF = Utility(ja) + bta*PF(ja,ix,ik,ik2,ip)

end function



!	----------------------------------------------------------------------
!	subroutine name: ExpectedValueFunction()
!	
!	calcualtes expected value functions.
!	----------------------------------------------------------------------


subroutine ExpectedValueFunction()

	implicit none

	real(8) kp, k2p, VFP(na,nx,np)
	integer ipp, ixp


	PF = 0
	
	do ik = 1, nk
	do ik2 = 1, nk2
	do ip = 1, np


		!	next period aggregate capital

		kp = dexp(Kcoef(1) + Kcoef(2)*dlog(kgrid(ik)) + Kcoef(3)*dlog(k2grid(ik2)) + Kcoef(4)*pgrid(ip))
		k2p = dexp(K2coef(1) + K2coef(2)*dlog(kgrid(ik)) + K2coef(3)*dlog(k2grid(ik2)) + K2coef(4)*pgrid(ip))


		do ipp = 1, np


			!	next period value function: VFP(ia,ixzp,ipp) with kp(ik,ip)

			do ixp = 1, nx
			do ia = 1, na

				VFP(ia,ixp,ipp) = polyinterp2(kp, k2p, kgrid, k2grid, VF(ia,ixp,:,:,ipp))

			end do
			end do

		end do

		
		!	calculate expected value function: PF(ia,ixp,ik,ip)

		do ix = 1, nx
		do ia = 1, na

			do ipp = 1, np
			do ixp = 1, nx

				PF(ia,ix,ik,ik2,ip) = PF(ia,ix,ik,ik2,ip) + VFP(ia,ixp,ipp)*trX(ix,ixp)*trP(ip,ipp)

			end do
			end do

		end do
		end do

	end do
	end do
	end do


end subroutine


end module
