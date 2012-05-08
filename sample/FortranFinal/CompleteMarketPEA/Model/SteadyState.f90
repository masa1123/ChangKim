!	----------------------------------------------------------------------
!	File name: SteadyState.f90
!
!	Solve for capital, efficiency units, B in the steady state.
!	----------------------------------------------------------------------


subroutine SteadyState()

	use Globals
	use IntegrandsModule
	use qdag_int
	use twodq_int
	use rnunf_int
	use rnnor_int
	use anorin_int
	use imsl_libraries

	implicit none


	!   find the reservation productivity that satisfies E = 60% in the steady state

	xss = d_anorin(1.0D0-Ess)
	xss = sigx*xss


	!   calculate total efficiency units (L) with xss in the steady state

	call d_qdag(normpdf, xss, 10.0D0*sigx, Hss, errabs=errabs, errrel=errrel, irule=irule, errest=errest)
	call d_qdag(expx_normpdf, xss, 10.0D0*sigx, Lss, errabs=errabs, errrel=errrel, irule=irule, errest=errest)

	Hss = hbar*Hss
	Lss = hbar*Lss

	
	!   calculate capital in the steady state by solving the euler equation 
	!	in the steady state
	!
	!   1 = bta*{F_K(K,L)+1-delta}

	Kss = (alpha/(irate+delta))**(1.0D0/(1.0D0-alpha))*Lss

	!   calculate aggreagete variables in the steady state: K/L, Output, Worker Flows

	Yss = Kss**alpha*Lss**(1.0D0-alpha)
	call d_twodq(flowfunc, -10.0D0*sigx, xss, critical1S, critical2, EINss, errabs=errabs, errrel=errrel, irule=irule, errest=errest)
	call d_twodq(flowfunc, xss, 10.0D0*sigx, critical3, critical1S, EOUTss, errabs=errabs, errrel=errrel, irule=irule, errest=errest)


	!   find B that is consistent with xss derived from the f.o.c.
	!
	!    F(K,L)-delta*K     1-(1-hbar)^(1-gama)
	!   ----------------*B*--------------------- = exp(xss)
	!     hbar*F_L(K,L)              1-gama

	B = dexp(xss)*hbar/(Yss-delta*Kss)*(1.0D0-alpha)*(Kss/Lss)**alpha		&
		/(hbar**(1.0D0+1.0D0/gamal)/(1.0D0+1.0D0/gamal))


	!	save steady state values

	open(1, file='Output\SteadyState.txt', status='unknown')

	write(1, '(A/)') "Parameter Values:"
	write(1, '(A,F13.10)') "beta     = ", bta
	write(1, '(A,F10.4)') "B        = ", B
	write(1, '(A,F10.4)') "gamma    = ", gamal
	write(1, '(A,F10.4)') "rho_x    = ", rhox
	write(1, '(A,F10.4)') "sig_x    = ", sigex
	write(1, '(/)') 
	write(1, '(A/)') "Aggregate Variables:"
	write(1, '(A,F15.10)') "Aggregate Capital                     = ", Kss
	write(1, '(A,F15.10)') "Aggregate Hours Worked                = ", Hss
	write(1, '(A,F15.10)') "Aggregate Effeciency Units of Hours   = ", Lss
	write(1, '(A,F15.10)') "Aggregate Capital per Effeciency Unit = ", Kss/Lss
	write(1, '(A,F15.10)') "Aggregate Output                      = ", Yss
	write(1, '(/)') 
	write(1, '(A,F15.10)') "Employment/Population                 = ", Ess
	write(1, '(A,F15.10," (",F8.5,")")') "Inflow to Employment                  = ", EINss, EINss/Ess
	write(1, '(A,F15.10," (",F8.5,")")') "Outflow from Employment               = ", EOUTss, EOUTss/Ess


end subroutine



