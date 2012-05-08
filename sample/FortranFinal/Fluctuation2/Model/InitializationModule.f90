!	----------------------------------------------------------------------
!	Module name: InitializationModule.f90
!	----------------------------------------------------------------------


module InitializationModule

	use Globals
	use ValueModule

	implicit none


contains


subroutine InitializeParameters()

	implicit none

	open(1, file='Input\Parameters.txt', status = 'old')
	read(1, *) bta
	read(1, *) B
	read(1, *) gama
	read(1, *) avex
	read(1, *) rhox
	read(1, *) sigex
	read(1, *) avep
	read(1, *) rhop
	read(1, *) sigep

	sigx = sigex/dsqrt(1.0D0-rhox**2)
	sigp = sigep/dsqrt(1.0D0-rhop**2)


	!	write parameter values on log file

	open(0, file='Output\Logfile.txt', status='unknown')

	write(0, '(A/)') "Parameter Values:"
	write(0, '(A,F12.9)') "beta     = ", bta
	write(0, '(A,F12.4)') "B        = ", B
	write(0, '(A,F7.4)') "gamma    = ", gama
	write(0, '(A,F7.4)') "ave_x    = ", avex
	write(0, '(A,F7.4)') "rho_x    = ", rhox
	write(0, '(A,F7.4)') "sig_ex   = ", sigex
	write(0, '(A,F7.4)') "ave_p    = ", avep
	write(0, '(A,F7.4)') "rho_p    = ", rhop
	write(0, '(A,F7.4)') "sig_ep   = ", sigep


end subroutine



subroutine InitializeCoefficients()

	open(1, file='Input\CoefficientsIN.txt', status = 'old')
	read(1, *) Kcoef, K2coef, Rcoef, Wcoef


end subroutine



subroutine InitializeValue()

	implicit none


	! initialize value function and decision rules

	do ip = 1, np
	do ik2 = 1, nk2
	do ik = 1, nk

		irate = dexp(Rcoef(1) + Rcoef(2)*dlog(kgrid(ik)) + Rcoef(3)*dlog(k2grid(ik2)) + Rcoef(4)*pgrid(ip)) - delta
		wage  = dexp(Wcoef(1) + Wcoef(2)*dlog(kgrid(ik)) + Wcoef(3)*dlog(k2grid(ik2)) + Wcoef(4)*pgrid(ip))

		do ia = 1, na
		do ix = 1, nx

			hour = hbar
			income = (1.0+irate)*agrid(ia) + wage*hour*exgrid(ix)
			ASE(ia,ix,ik,ik2,ip) = ia
			VFE(ia,ix,ik,ik2,ip) = Utility(ia)

			hour = 0.0D0
			income = (1.0+irate)*agrid(ia)
			ASN(ia,ix,ik,ik2,ip) = ia
			VFN(ia,ix,ik,ik2,ip) = Utility(ia)

			AS(ia,ix,ik,ik2,ip) = ia
			VF(ia,ix,ik,ik2,ip) = max(VFE(ia,ix,ik,ik2,ip),VFN(ia,ix,ik,ik2,ip))

		end do
		end do

	end do
	end do
	end do

end subroutine


end module
