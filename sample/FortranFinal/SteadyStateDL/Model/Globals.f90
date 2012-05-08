!	----------------------------------------------------------------------
!	File name: Globals.f90
!
!	Defines global variables.
!	----------------------------------------------------------------------


module Globals

	implicit none

	save


	!	variables for solving for an equilibrium of the model

!	real(8), parameter:: amin = -2.0D0,	amax = 280.0D0	! rho=0.933, sigx=0.265
	real(8), parameter:: amin = -2.0D0,	amax = 200.0D0	! rho=0.929, sigx=0.227


	!	grids for individual asset holdings

	real(8), parameter:: a0 = amin, a1 = -1.5D0, a2 = 0.0D0, a3 = 5.0D0,	&
						 a4 = 15.0D0, a5 = 40.0D0, a6 = 80.0D0, a7 = amax

	real(8), parameter:: s1 = 0.02D0, s2 = 0.05D0, s3 = 0.1D0, s4 = 0.2D0,	&
						 s5 = 0.5D0, s6 = 1.0D0, s7 = 2.0D0

	integer, parameter:: na1 = int((a1-a0)/s1) + 1, na2 = na1 + int((a2-a1)/s2),	&
						 na3 = na2 + int((a3-a2)/s3), na4 = na3 + int((a4-a3)/s4),	&
						 na5 = na4 + int((a5-a4)/s5), na6 = na5 + int((a6-a5)/s6),	&
						 na7 = na6 + int((a7-a6)/s7)

	integer, parameter:: na = na7, nx = 17


	!	asset grids that will be used for the invariant measure

	real(8), parameter:: smu1 = 0.01D0, smu2 = 0.025D0, smu3 = 0.04D0,		&
						 smu4 = 0.08D0,	smu5 = 0.1D0, smu6 = 0.2D0, smu7 = 0.5D0

	integer, parameter:: namu1 = int((a1-a0)/smu1) + 1,		&
						 namu2 = namu1 + int((a2-a1)/smu2),	&
						 namu3 = namu2 + int((a3-a2)/smu3), &
						 namu4 = namu3 + int((a4-a3)/smu4),	&
						 namu5 = namu4 + int((a5-a4)/smu5), &
						 namu6 = namu5 + int((a6-a5)/smu6),	&
						 namu7 = namu6 + int((a7-a6)/smu7)

	integer, parameter:: namu = namu7


	!	model parameters

	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0, irate = 0.01d0,		&
						 rd = irate + delta, hbar = 0.2D0
	
	!	model parameters to be read

	real(8):: bta, B, gama, avex, rhox, sigex


	!	auxiliary variables
	
	real(8):: kss, wage, income, sigx, hour, ap_hour
	

	!	storages

	real(8):: agrid(na), xgrid(nx), exgrid(nx), trX(nx,nx), CtrX(nx,nx),		&
			  prX(nx), CprX(nx), AS(na,nx), VF(na,nx), HR(na,nx), PF(na,nx),	&
			  agridmu(namu), mu(namu,nx), cmu(namu,nx), HRmu(namu,nx),	EFmu(namu,nx)
	integer:: ia, ix


	!	timer

	character(8)::  starttime, finishtime
	character(10):: startdate, finishdate


	contains

	subroutine GlobalVariables()

		kss = (alpha/rd)**(1.0D0/(1.0D0-alpha))
		wage = (1.0D0-alpha)*(kss**alpha)
		sigx = sigex/dsqrt(1.0D0-rhox**2)

	end subroutine

end module
