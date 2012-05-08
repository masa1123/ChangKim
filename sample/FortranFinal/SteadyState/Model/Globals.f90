!	----------------------------------------------------------------------
!	File name: Globals.f90
!
!	Defines global variables.
!	----------------------------------------------------------------------


module Globals

	implicit none

	save


	!	variables for solving for an equilibrium of the model

	real(8), parameter:: amin = -2.0D0
!	real(8), parameter:: amax = 310.0D0	! rho=0.933, sigx=0.265
	real(8), parameter:: amax = 250.0D0	! rho=0.929, sigx=0.227


	!	grids for individual asset holdings

	real(8), parameter:: a0 = amin, a1 = 2.0D0, a2 = 5.0D0, a3 = 10.0D0, a4 = 20.0D0, a5 = 40.0D0, a6 = 85.0D0, a7 = amax
	real(8), parameter:: s1 = 0.02D0, s2 = 0.03D0, s3 = 0.05D0, s4 = 0.1D0, s5 = 0.2D0, s6 = 0.3D0, s7 = 0.4D0
	integer, parameter:: na1 = int((a1-a0)/s1) + 1, na2 = na1 + int((a2-a1)/s2),		&
						   na3 = na2 + int((a3-a2)/s3), na4 = na3 + int((a4-a3)/s4),	&
						   na5 = na4 + int((a5-a4)/s5), na6 = na5 + int((a6-a5)/s6),	&
						   na7 = na6 + int((a7-a6)/s7)
	integer, parameter:: na = na7, nx = 17


	!	model parameters

	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0, irate = 0.01d0, rd = irate + delta, hbar = 1.0D0/3.0D0
	

	!	model parameters to be read

	real(8):: bta, B, gama, avex, rhox, sigex


	!	auxiliary variables
	
	real(8):: kss, wage, income, sigx, hour
	

	!	storages

	real(8):: agrid(na), xgrid(nx), exgrid(nx), trX(nx,nx), CtrX(nx,nx),		&
			  prX(nx), VF(na,nx), VFE(na,nx), VFN(na,nx), HR(na,nx), PF(na,nx), &
			  mu(na,nx), cmu(na,nx)
	integer:: AS(na,nx), ASE(na,nx), ASN(na,nx), ia, ix


	!	timer

	character(8)::  starttime, finishtime
	character(10):: startdate, finishdate


	!	variables for distributions

	real(8):: AsstEarnDensity(na,0:nx), AsstEarnDist(na,0:nx), AsstDensity(na),	&
			  AsstDist(na), EarnDensity(0:nx), EarnDist(0:nx)


	!	asset quintile
	
!	integer, parameter:: asst20 = 60, asst40 = 247, asst60 = 401, asst80 = 514	!	sigx = 0.225, rho = 0.95
	integer, parameter:: asst20 = 85, asst40 = 274, asst60 = 414, asst80 = 526	!	sigx = 0.287, rho = 0.939


	!	variables for reservation wages and labor supply

	real(8):: ReservationWage(na), ReservationX(na), Reserv_cmu(na,na), LS(na)


	contains

	subroutine GlobalVariables()

		kss = (alpha/rd)**(1.0D0/(1.0D0-alpha))
		wage = (1.0D0-alpha)*(kss**alpha)
		sigx = sigex/dsqrt(1.0D0-rhox**2)

	end subroutine

end module
