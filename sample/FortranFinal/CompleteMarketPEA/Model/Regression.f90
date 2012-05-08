Subroutine Regression()

	use Globals
	use rlse_int

	implicit none
	
	integer iter_LS
	real(8) err, coef0(3), psiQ0(nT), dpsiQ0(nT,3), dep0(nT), sst, sse


	!	nonlinear least square regression: linearized model

	iter_LS = 0
	err = 1000.D0
	coefNew = coef
	
	do while (err > tol_LS .and. iter_LS < Max_Iter_LS)

		iter_LS = iter_LS + 1

		coef0 = coefNew

		do t = t_LS, nT-1

			!	psi(beta0)

			psiQ0(t) = coef0(1)*exp(coef0(2)*logKdata(t) + coef0(3)*logPdata(t))

			!	regressors

			dpsiQ0(t,1) = psiQ0(t)/coef0(1)
			dpsiQ0(t,2) = psiQ0(t)*logKdata(t)
			dpsiQ0(t,3) = psiQ0(t)*logPdata(t)


			!	dependent variable

			dep0(t) = phiQ(t) - psiQ0(t) + dot_product(dpsiQ0(t,:), coef0)

		end do

		
		!	regression: OLS

		call d_rlse(dep0(t_LS:nT-1), dpsiQ0(t_LS:nT-1,:), coefNew, intcep = 0, sst=sst, sse=sse)


		!	convergence check

		err = maxval(dabs(coefNew - coef0))

	end do


end subroutine
