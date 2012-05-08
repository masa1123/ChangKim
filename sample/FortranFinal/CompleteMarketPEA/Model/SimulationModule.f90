!	----------------------------------------------------------------------
!	Module name: SimulationModule.f90
!
!	collection of subroutines for similating parameterized expectation 
!	and artificial data.
!	----------------------------------------------------------------------


module SimulationModule

	use Globals
	use IntegrandsModule
	use qdag_int
	use twodq_int
	use rnunf_int
	use rnnor_int
	use imsl_libraries

	implicit none

contains


subroutine Simulation()


	!	initial capital

	Kdata(1) = Kss
	logKdata(1) = dlog(Kss)


	!	generate time series of consumption and capital

	do t = 1, nT

		!	consumption using parametric expectation psi

		psiQ(t) = coef(1)*exp(coef(2)*logKdata(t) + coef(3)*logPdata(t))
		Cdata(t) = (bta*psiQ(t))**(-1.0D0/gamac)

		!	efficiency units of labor

		Xdata(t) = xss
		call FindReserv(Xdata(t))
		call d_qdag(expx_normpdf, Xdata(t), 10.0D0*sigx, Ldata(t), errabs=errabs, errrel=errrel, irule=irule, errest=errest)
		Ldata(t) = hbar*Ldata(t)

		!	next period capital

		Kdata(t+1) = Pdata(t)*Kdata(t)**alpha*Ldata(t)**(1.0D0-alpha)	&
					+ (1.0D0-delta)*Kdata(t) - Cdata(t)

		if (Kdata(t+1) < 0.0D0) then
			print*, 'consumption too large'
			stop
		endif

		logKdata(t+1) = log(Kdata(t+1))

	end do

		
	!	calculate phi: the expression inside the conditional expectation
	!	which will be used to update psi

	do t = 1, nT-1
	
		phiQ(t) = Cdata(t+1)**(-gamac)*(alpha*Pdata(t+1)*(Kdata(t+1)/Ldata(t+1))**(alpha-1.0D0) + 1.0D0-delta)
	
	end do


	!	additional variables after convergence of the coefficients

	if (final) then

		do t = 1, nT
			
			call d_qdag(normpdf, Xdata(t), 10.0D0*sigx, Edata(t), errabs=errabs, errrel=errrel, irule=irule, errest=errest)
			Ydata(t) = Pdata(t)*Kdata(t)**alpha*Ldata(t)**(1.0D0-alpha)
			Idata(t) = Ydata(t) - Cdata(t)
			Wdata(t) = (1.0D0-alpha)*Pdata(t)*(Kdata(t)/Ldata(t))**alpha
			Rdata(t) = alpha*Pdata(t)*(Kdata(t)/Ldata(t))**(alpha-1.0D0)
			Adata(t) = Ydata(t)/(hbar*Edata(t))
			if (t >= 2) then
				call d_twodq(flowfunc, -10.0D0*sigx, Xdata(t-1), critical1F, critical2, EINdata(t), errabs=errabs, errrel=errrel, irule=irule, errest=errest)
				call d_twodq(flowfunc, Xdata(t-1), 10.0D0*sigx, critical3, critical1F, EOUTdata(t), errabs=errabs, errrel=errrel, irule=irule, errest=errest)
				hEINdata(t) = EINdata(t)/Edata(t)
				hEOUTdata(t) = EOUTdata(t)/Edata(t)
			end if

		end do

	end if

end subroutine





subroutine FindReserv(xR)

	real(8), parameter:: tol_xR = 1.0D-6
	real(8) fval, fder, xR
	integer iter_xR


	!	iteration to find hour

	do iter_xR = 1, 200
	
		call EvaluateF(xR, fval, fder)
		xR = xR - fval/fder

		if(dabs(fval) < tol_xR) then
			exit
		end if

	end do

	if (iter_xR >= 200) print *, "Can't find xR"

end subroutine




subroutine EvaluateF(xR, fval, fder)

	real(8) fval, fder, xR, xlabor, mpl

	call d_qdag(expx_normpdf, xR, 10.0D0*sigx, xlabor, errabs=errabs, errrel=errrel, irule=irule, errest=errest)
	xlabor = hbar*xlabor

	mpl = (1.0D0-alpha)*Pdata(t)*(kdata(t)/xlabor)**alpha
	fval = (Cdata(t)**(-gamac))*hbar*mpl*dexp(xR)				&
			- B*hbar**(1.0D0+1.0D0/gamal)/(1.0D0+1.0D0/gamal)
	
	fder = (Cdata(t)**(-gamac))*hbar*mpl*dexp(xR)*				&
			(1.0D0 + alpha*hbar*normpdf(xR)*dexp(xR)/xlabor)

end subroutine




subroutine GenerateShocks()

	use Globals

	integer, parameter:: pseed = 1
	real(8) nshock


	!	set the seed number for random number generator

	call rnset(pseed)


	!	normal random numbers for aggregate productivity shocks

	call d_rnnor(logPdata)
	call dscal(nT, sigep, logPdata, 1)


	! generate Pdata

	do t = 2, nT
		logPdata(t) = min(max(rhop*logPdata(t-1) + logPdata(t), -3.0*sigp), 3.0*sigp)
	end do

	Pdata = dexp(logPdata)

end subroutine


end module
