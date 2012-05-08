!	----------------------------------------------------------------------
!	File name: Initialize.f90
!	----------------------------------------------------------------------



subroutine InitializeParameters()

	use Globals

	implicit none

	open(1, file='Input\Parameters.txt', status = 'old')
	read(1, *) gama, kss, B


end subroutine


subroutine InitializeValue()

	use Globals
	use ValueModule

	implicit none


	! initialize value function and decision rules

	do ik = 1, nk
		do ip = 1, np
			KS(ik,ip) = kgrid(ik)
			VF(ik,ip) = Utility(KS(ik,ip))
		end do
	end do

end subroutine 



subroutine SteadyState()

	use Globals

	implicit none

	real(8) :: hss, wss, css


	hss = kss*(alpha/(irate+delta))**(1.0D0/(alpha-1.0D0))
	css = kss**alpha*hss**(1.0D0-alpha) - delta*kss
	wss = (1.0D0-alpha)*(kss/hss)**alpha


	B = wss/css*hss**(-1.0D0/gama)

	write(*, '(A,F12.6)') "B   = ", B
	write(*, '(A,F12.6)') "Hss = ", hss
	write(*, '(A,F12.6)') "Wss = ", Wss

end subroutine
