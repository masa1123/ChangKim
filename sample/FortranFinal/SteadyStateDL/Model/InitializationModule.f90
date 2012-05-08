!	----------------------------------------------------------------------
!	Module name: InitializationModule.f90
!	----------------------------------------------------------------------


module InitializationModule

	use Globals
	use TauchenMod
	use ValueModule

	implicit none


contains


subroutine InitializeParameters()

	open(1, file='Input\Parameters.txt', status = 'old')
	read(1, *) bta
	read(1, *) B
	read(1, *) gama
	read(1, *) avex
	read(1, *) rhox
	read(1, *) sigex

end subroutine



subroutine InitializeTransition()

	implicit none
	
	integer i, j
	real(8) prXp(nx), err


	! transition probability matrix

	call TauchenTr(avex, sigx, rhox, xgrid, trX)
	
	open(1, file='Output\trX.txt', status='unknown')
	write(1, '(<nx>f12.6)') ((trX(i,j),j=1,nx),i=1,nx)


	!	cumulative transition probability matrix for random drawing

	CtrX(:,1) = trX(:,1)
	do ix = 2, nx
		CtrX(:,ix) = CtrX(:,ix-1) + trX(:,ix)
	end do

	open(2, file='Output\CtrX.txt', status='unknown')
	write(2, '(<nx>f12.6)') ((CtrX(i,j),j=1,nx),i=1,nx)


	! invariant distribution of shocks

	call TauchenNew(avex, sigx, xgrid, prX)

	do while (err < 1.0e-10)

		prXp = matmul(transpose(trX), prX)
		err = maxval(dabs(prXp-prX))
		prX = prXp
	
	end do


	if (dabs(1.0D0 - sum(prX)) <= 1.0e-6) then
		prX = prX/sum(prX)
		CprX(1) = prX(1) 
		do ix = 2, nx
			CprX(ix) = CprX(ix-1) + prX(ix)
		end do
	else
		print*, "Fatal Error in InvariantProbability"
		stop
	end if


	! save transition probability matrix and invariant probability

	open(3, file='Output\prX.txt', status='unknown')
	open(4, file='Output\CprX.txt', status='unknown')

	write(3, '(f20.16)') prX
	write(4, '(f20.16)') CprX

end subroutine



subroutine InitializeValue()

	implicit none


	! initialize value function and decision rules given (x,z)

	do ia = 1, na
	do ix = 1, nx

		AS(ia,ix) = agrid(ia)
		HR(ia,ix) = hbar
		VF(ia,ix) = Utility(AS(ia,ix))

	end do
	end do


end subroutine



subroutine ReadValue()

	implicit none


	! initialize value function and decision rules using previous results

	open(1, file='Output\VF.txt', status='unknown')
	open(2, file='Output\AS.txt', status='unknown')
	read(1, *) ((VF(ia,ix),ix=1,nx),ia=1,na)
	read(2, *) ((AS(ia,ix),ix=1,nx),ia=1,na)
	close(1)
	close(2)


end subroutine



end module
