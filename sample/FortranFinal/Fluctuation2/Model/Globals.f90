!	----------------------------------------------------------------------
!	File name: Globals.f90
!
!	Defines global variables.
!	----------------------------------------------------------------------


module Globals

	implicit none

	save
	
	!	bounds for individual and aggregate capital

	!	sigx=0.227 rhox=0.929

	real(8), parameter:: amin = -2.0D0,	amax = 250.0D0
	
	real(8), parameter:: kbnd = 0.10D0, kss = 11.6001418987D0
	real(8), parameter:: kmax = kss*(1.0D0 + kbnd), kmin = kss*(1.0D0 - kbnd)

	real(8), parameter:: k2bnd = 0.15D0, k2ss = 348.9237959719D0
	real(8), parameter:: k2max = k2ss*(1.0D0 + kbnd), k2min = k2ss*(1.0D0 - kbnd)
	
	real(8), parameter:: pbnd = 2.9D0
	

	!	number of grids

	real(8), parameter:: a0 = amin, a1 = 2.0D0, a2 = 5.0D0, a3 = 10.0D0, a4 = 20.0D0, a5 = 40.0D0, a6 = 85.0D0, a7 = amax
	real(8), parameter:: s1 = 0.02D0, s2 = 0.03D0, s3 = 0.05D0, s4 = 0.1D0, s5 = 0.2D0, s6 = 0.3D0, s7 = 0.4D0
	integer, parameter:: na1 = int((a1-a0)/s1) + 1, na2 = na1 + int((a2-a1)/s2), na3 = na2 + int((a3-a2)/s3),		&
						   na4 = na3 + int((a4-a3)/s4),	na5 = na4 + int((a5-a4)/s5), na6 = na5 + int((a6-a5)/s6),	&
						   na7 = na6 + int((a7-a6)/s7)
	integer, parameter:: na = na7, nx = 17, nk = 7, nk2 = 7, np = 9


	!	maximum number of iterations

	integer, parameter:: maxiterLOM = 40, MAXITER = 10000, showerr = 10


	!	tolerence for converce check

	real(8), parameter:: tol_value = 1.0D-4, tol_LOM = 1.0D-4


	!	number of periods, individuals for simulation

	integer, parameter:: Nperiod = 3500, Nskip = 500, Nindiv_Sim = 50000, Nindiv_Fin_Sim = 200000
						 

	!	adjustment factor

	real(8), parameter:: sfac = 0.5D0

	
	!	model parameters

	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0, hbar = 1.0D0/3.0D0
	real(8) bta, B, gama, avex, rhox, sigex, avep, rhop, sigep, sigx, sigp


	!	auxiliary variables

	real(8) irate, wage, income, hour


	!	grid vectors

	real(8) agrid(na), xgrid(nx), exgrid(nx), kgrid(nk), k2grid(nk2), pgrid(np), epgrid(np)


	!	transition probability

	real(8) trX(nx,nx), CtrX(nx,nx), prX(nx), trP(np,np), CtrP(np,np)


	!	simulated time series data

	real(8) Kdata(Nperiod+1), K2data(Nperiod+1), Ldata(Nperiod), Rdata(Nperiod), Wdata(Nperiod), Pdata(Nperiod),	&
			 Ydata(Nperiod), Idata(Nperiod), Cdata(Nperiod), Adata(Nperiod), Edata(Nperiod), EINdata(Nperiod),		&
			 EOUTdata(Nperiod), hEINdata(Nperiod), hEOUTdata(Nperiod)


	!	coefficients of equations

	integer, parameter:: ncoef = 4
	real(8) Kcoef(4), K2coef(4), Rcoef(4), Wcoef(4), NewKcoef(4), NewK2coef(4), NewRcoef(4), NewWcoef(4)


	!	indexing variables
		
	integer ia, ix, ik, ik2, ip, indiv, time, iter, iterLOM


	!	final simulation
	
	logical final
	

	!	spaces for value functions and decision rules

	integer AS(na,nx,nk,nk2,np), ASE(na,nx,nk,nk2,np), ASN(na,nx,nk,nk2,np)
	real(8) VF(na,nx,nk,nk2,np), HR(na,nx,nk,nk2,np), PF(na,nx,nk,nk2,np),	VFE(na,nx,nk,nk2,np), VFN(na,nx,nk,nk2,np)


end module
