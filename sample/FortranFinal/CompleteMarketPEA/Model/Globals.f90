!	----------------------------------------------------------------------
!	File name: Globals.f90
!
!	Defines global variables.
!	----------------------------------------------------------------------


module Globals

	implicit none

	save
	 
	!	fixed numbers

	integer, parameter:: nT = 3500, neq = 1, np = 3, Max_Iter = 10000,		&
						 Max_Iter_LS = 100, t_LS = 500
	real(8), parameter:: tol_Iter = 1.0D-6, tol_LS = 1.0D-6, sfac = 0.25D0


	!	model parameters

	real(8), parameter:: pi = 3.14159265358979D0
	real(8), parameter:: alpha = 0.36D0, delta = 0.025D0, irate = 0.01D0,	&
						 bta = 1.0D0/(1.0D0+irate),	Ess = 0.6D0, hbar = 1.0D0/3.0D0
	real(8) B, gamac, gamal, rhop, sigep, rhox, sigex, sigx, sigp
	
	
	!	steady state variables
	
	real(8) Kss, Hss, Lss, Xss, Yss, EINss, EOUTss


	!	storage of similation results

	real(8) Kdata(nT+1), logKdata(nT+1), Pdata(nT), logPdata(nT), Cdata(nT),	&
			Xdata(nT), Ldata(nT), Ydata(nT), Idata(nT), Wdata(nT), Rdata(nT),	&
			Adata(nT), Edata(nT), EINdata(nT), EOUTdata(nT), hEINdata(nT),		&
			hEOUTdata(nT), phiQ(nT), psiQ(nT)

	
	!	parameters of the conditional expectation function

	real(8) coef(np), coefNew(np)


	!	control for final simulation

	logical final


	!	variables for numerical integration

	real(8):: errest
	real(8):: errabs, errrel
	integer:: irule
	

	!	integer for indexing

	integer t


end module
