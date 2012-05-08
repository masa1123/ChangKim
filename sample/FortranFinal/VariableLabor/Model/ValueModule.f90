!	----------------------------------------------------------------------
!	File name : Value.cpp
!	----------------------------------------------------------------------	*/


module ValueModule

	use Numerical_Libraries
	use OptimizationMod
	use Globals

contains


real(8) function Prodf(k, h, p)

	implicit none

	real(8) k, h, p

	Prodf = p*k**alpha*h**(1.0D0-alpha)

end function



real(8) function SolveHour(xhour)

	implicit none

	real(8), parameter:: tiny = 1.0D-20
	real(8) wage, cons, xhour

	wage = epgrid(ip)*(1.0D0-alpha)*(kgrid(ik)/xhour)**alpha
	cons = max(Prodf(kgrid(ik), xhour, epgrid(ip))		&
				+ (1.0D0-delta)*kgrid(ik) - kp_hour, tiny)
	SolveHour = wage/cons - B*xhour**(1.0D0/gama)

end function



real(8) function Utility(kp)

	implicit none

	real(8), parameter:: tiny = 1.0D-20, tol_hour = 1.0D-6
	real(8) kp, cm, ls

	kp_hour = kp
	hour = bisect(SolveHour, 0.0D0, 0.9D0, tol_hour)

	cm = max(Prodf(kgrid(ik), hour, epgrid(ip))			&
				+ (1.0D0-delta)*kgrid(ik) - kp, tiny)
	ls = B*hour**(1.0D0+1.0D0/gama)/(1.0D0+1.0D0/gama)

	Utility = log(cm) - ls

end function


real(8) function ValueF(kp)

	implicit none

	real(8) kp

	Value = Utility(kp) + bta*dcsval(kp, nk-1, break, cscoef)
	ValueF = Value

end function


real(8) function NegValueF(kp)

	implicit none

	real(8) kp

	NegValueF = -ValueF(kp)

end function

end module
