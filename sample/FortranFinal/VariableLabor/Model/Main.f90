!	----------------------------------------------------------------------
!	File name : VariableLabor.f90
!
!	coordinates all routines that solve the representative agent model 
!	with variable labor input.
!	----------------------------------------------------------------------


program VariableLabor

    include 'link_f90_dll.h'

	use Numerical_Libraries
	use Globals
	use ValueModule
	use SimulationModule
	use TauchenMod
	use OptimizationMod

	implicit none

	integer, parameter:: maxiter = 10000
	integer:: iter, i, j
	real(8):: errV, errK, errH, err
	real(8):: lbd, ubd, xmin, PF(nk,np), NewVF(nk,np), NewKS(nk,np), NewHR(nk,np)


	!	initialize parameters

	call InitializeParameters()


	! steady state capital and corresponding discount factor

	bta = 1.0D0/(1.0D0+irate)
	sigp = sig/dsqrt(1.0D0-rhop**2)


	! initialize capital grids

	call CapitalGrid()
	call ShockGrid()


	!	calculate Tauchen's transition probability matrix 

	call TauchenTr(avep, sigp, rhop, pgrid, trP)

	open(3, file='Output\trP.txt', status='unknown')
	write(3, '(<np>f12.6)') ((trP(i,j),j=1,np),i=1,np)


	! initialize value function and decision rules

	call InitializeValue()
	

	! steady state

	call SteadyState()


	! iterate value function

	do iter = 1, maxiter

		! conditional expectation of value function

		PF = matmul(VF, transpose(trP))


		! maximization procedure


		do ip = 1, np

			call dcsint(nk, kgrid, PF(:,ip), break, cscoef)
	
			do ik = 1, nk

				lbd = 0.9*kgrid(ik)
				ubd = 1.1*kgrid(ik)
				call fmin(NegValueF, lbd, ubd, xmin)
	
				NewKS(ik,ip) = xmin
				NewVF(ik,ip) = Value
				NewHR(ik,ip) = hour

			end do
		end do


		! calculate the errors

		errV = maxval(dabs(NewVF - VF))
		errK = maxval(dabs(NewKS - KS))
		errH = maxval(dabs(NewHR - HR))
		err  = max(errV, errK, errH)


		! update value function and decision rules

		VF = NewVF
		KS = NewKS
		HR = NewHR


		! convergence check

		if (errV < tol_cnvrg) then
			write(*, '(A,I6)') "Convergence achieved at iteration  ", iter
			exit 
		else if (mod(iter,100) == 0) then
			write(*, '(A,I6)') "iteration  ", iter
			write(*, '(A,F10.7,A,F10.7,A,F10.7,/)')		&
				"errV = ", errV, "   errK = ", errK, "   errH = ", errH
		end if

	end do


	! save value function and decision rules

	open(1, file='Output\VF.txt', status='unknown')
	open(2, file='Output\KS.txt', status='unknown')
	open(3, file='Output\HR.txt', status='unknown')

	write(1, '(<np>F12.6)') ((VF(ik,ip),ip=1,np),ik=1,nk)
	write(2, '(<np>F12.6)') ((KS(ik,ip),ip=1,np),ik=1,nk)
	write(3, '(<np>F12.6)') ((HR(ik,ip),ip=1,np),ik=1,nk)


	! generate artificial data for statistics and regression

	call SimulateData()


	! generate artificial data for impulse response

!	call ImpulseResponse()


end program VariableLabor
