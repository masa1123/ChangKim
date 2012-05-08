!	----------------------------------------------------------------------
!	Module name: IntegrandsModule.f90
!	
!	Collection of functions that are used in numerical integrations.
!	----------------------------------------------------------------------


module IntegrandsModule

	use Globals

	implicit none

contains


real(8) function normpdf(x)

	real(8), intent(in):: x
	
	normpdf = (1.0D0/(dsqrt(2.0D0*pi)*sigx))*dexp(-0.5D0*(x/sigx)**2.0D0)

end function



real(8) function expx_normpdf(x)

	real(8), intent(in):: x
	
	expx_normpdf = dexp(x)*normpdf(x)

end function



real(8) function flowfunc(x,eps)

	real(8), intent(in):: x, eps
	
	flowfunc = normpdf(x)*(1.0d0/(dsqrt(2.0d0*pi)*sigex))*dexp(-0.5d0*(eps/sigex)**2.0d0)

end function



real(8) function critical1S(x)

	real(8), intent(in):: x
	
	critical1S = xss - rhox*x

end function



real(8) function critical2(x)

	real(8), intent(in):: x
	
	critical2 = 10.0D0*sigex

end function



real(8) function critical3(x)

	real(8), intent(in):: x
	
	critical3 = -10.0D0*sigex

end function



real(8) function critical1F(x)

	real(8), intent(in):: x
	
	critical1F = Xdata(t) - rhox*x

end function



end module