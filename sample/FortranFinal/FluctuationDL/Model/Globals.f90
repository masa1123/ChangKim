!	----------------------------------------------------------------------
!	File name: Globals.f90
!
!	Defines global variables.
!	----------------------------------------------------------------------


module Globals

	implicit none

	save
	

	!	bounds for individual and aggregate capital 

!	real(8), parameter:: amin = -2.0D0,	amax = 220.0D0, kss = 11.3745522464D0	! (rho=0.933, sigx=0.265)
	real(8), parameter:: amin = -2.0D0,	amax = 200.0D0, kss = 10.2647026623D0	! (rho=0.929, sigx=0.227)
	real(8), parameter:: kbnd = 0.1D0, pbnd = 3.0D0
	real(8), parameter:: kmax = kss*(1.0D0 + kbnd), kmin = kss*(1.0D0 - kbnd)
						 

	!	grids for individual asset holdings

	real(8), parameter:: a0 = amin, a1 = -1.5D0, a2 = 0.0D0, a3 = 5.0D0, a4 = 15.0D0, a5 = 40.0D0, a6 = 80.0D0, a7 = amax
	real(8), parameter:: s1 = 0.02D0, s2 = 0.05D0, s3 = 0.1D0, s4 = 0.2D0, s5 = 0.5D0, s6 = 1.0D0, s7 = 2.0D0

	integer, parameter:: na1 = int((a1-a0)/s1) + 1, na2 = na1 + int((a2-a1)/s2),		&
						   na3 = na2 + int((a3-a2)/s3), na4 = na3 + int((a4-a3)/s4),	&
						   na5 = na4 + int((a5-a4)/s5), na6 = na5 + int((a6-a5)/s6),	&
						   na7 = na6 + int((a7-a6)/s7)
	integer, parameter:: na = na7, nx = 17, nk = 11, np = 9


	!	maximum number of iterations

	integer, parameter:: maxiterLOM = 40, MAXITER = 10000, showerr = 10


	!	tolerence for converce check

	real(8), parameter:: tol_value = 1.0D-4, tol_LOM = 1.0D-4, tol_hour = 1.0D-4


	!	number of periods, individuals for simulation

	integer, parameter:: Nperiod = 3500, Nindiv = 50000, Nskip = 500, startpanel = 3281


	!	adjustment factor

	real(8), parameter:: sfac = 0.4D0

	
	!	model parameters

	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0, hbar = 1.0D0/3.0D0
	real(8) bta, B, gama, avex, rhox, sigex, avep, rhop, sigep, sigx, sigp


	!	auxiliary variables

	real(8) irate, wage, income, hour, ap_hour


	!	grid vectors

	real(8) agrid(na), xgrid(nx), exgrid(nx), kgrid(nk), pgrid(np), epgrid(np)


	!	transition probability

	real(8) trX(nx,nx), CtrX(nx,nx), prX(nx), trP(np,np), CtrP(np,np)


	!	simulated time series data

	real(8) Kdata(Nperiod+1), Ldata(Nperiod), Rdata(Nperiod), Wdata(Nperiod), Pdata(Nperiod), Ydata(Nperiod),	&
			 K2data(Nperiod), Idata(Nperiod), Cdata(Nperiod), Adata(Nperiod), Edata(Nperiod)
	integer iPdata(Nperiod)


	!	coefficients of equations

	real(8) Kcoef(3), Rcoef(3), Wcoef(3), NewKcoef(3), NewRcoef(3), NewWcoef(3)


	!	indexing variables
		
	integer ia, ix, ik, ip, indiv, time


	!	spaces for value functions and decision rules

	real(8) VF(na,nx,nk,np), AS(na,nx,nk,np), HR(na,nx,nk,np), PF(na,nx,nk,np)


	!	construct panel data

	logical panel

end module
