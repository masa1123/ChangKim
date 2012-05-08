!	----------------------------------------------------------------------
!	Module name: TransitionModule.f90
!
!	calculates transition probabilities for idiosyncratic shocks and 
!	aggregate productivity shock.
!	----------------------------------------------------------------------


module TransitionModule

	use Globals
	use TauchenMod

	implicit none

	integer i, j

contains


subroutine IdiosyncraticTransition()

	implicit none
	
	real(8) prXp(nx), err


	! transition probability matrix for idiosyncratic shocks

	call TauchenTr(avex, sigx, rhox, xgrid, trX)
	
	open(1, file='Output\trX.txt', status='unknown')
	write(1, '(<nx>f12.6)') ((trX(i,j),j=1,nx),i=1,nx)


	!	cumulative transition matrix of idiosyncratic shocks for random drawing

	CtrX(:,1) = trX(:,1)
	do ix = 2, nx
		CtrX(:,ix) = CtrX(:,ix-1) + trX(:,ix)
	end do

	open(2, file='Output\CtrX.txt', status='unknown')
	write(2, '(<nx>f12.6)') ((CtrX(i,j),j=1,nx),i=1,nx)


	! invariant distribution of shocks for idiosyncratic shocks

	call TauchenNew(avex, sigx, xgrid, prX)

	do while (err < 1.0e-10)

		prXp = matmul(transpose(trX), prX)
		err = maxval(dabs(prXp-prX))
		prX = prXp

	end do


	if (dabs(1.0D0 - sum(prX)) <= 1.0e-6) then
		prX = prX/sum(prX)
	else
		print*, "Fatal Error in InvariantProbability"
		stop
	end if


	! save transition probability matrix and invariant probability

	open(3, file='Output\prX.txt', status='unknown')
	write(3, '(f20.16)') prX


end subroutine



subroutine AggregateTransition()

	implicit none
	

	! transition probability matrix for aggregate shock: more than 2 states case

	call TauchenTr(avep, sigp, rhop, pgrid, trP)

	open(1, file='Output\trP.txt', status='unknown')
	write(1, '(<np>f12.6)') ((trP(i,j),j=1,np),i=1,np)


	!	cumulative transition matrix of aggregate shock for random drawing

	CtrP(:,1) = trP(:,1)
	do ip = 2, np
		CtrP(:,ip) = CtrP(:,ip-1) + trP(:,ip)
	end do

	open(2, file='Output\CtrP.txt', status='unknown')
	write(2, '(<np>f12.6)') ((CtrP(i,j),j=1,np),i=1,np)


end subroutine


subroutine AggregateTransition2()

	implicit none
	

	! transition probability matrix for aggregate shock: 2 states case

	trP(1,1) = 0.95D0
	trP(1,2) = 0.05D0
	trP(2,1) = 0.05D0
	trP(2,2) = 0.95D0

	open(1, file='Output\trP.txt', status='unknown')
	write(1, '(<np>f12.6)') ((trP(i,j),j=1,np),i=1,np)


	!	cumulative transition matrix of aggregate shock for random drawing

	CtrP(:,1) = trP(:,1)
	CtrP(:,2) = CtrP(:,1) + trP(:,2)

	open(2, file='Output\CtrP.txt', status='unknown')
	write(2, '(<np>f12.6)') ((CtrP(i,j),j=1,np),i=1,np)


end subroutine


end module

