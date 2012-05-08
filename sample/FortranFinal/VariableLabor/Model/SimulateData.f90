module SimulationModule

	use Globals
	use ValueModule
	use PolyInterpModule
	use Numerical_Libraries

contains


subroutine SimulateData()

	implicit none

	integer t
	real(8) Ydata(NT), Kdata(NT+1), Edata(NT), Cdata(NT), Idata(NT), Pdata(NT), &
			Rdata(NT), Wdata(NT), Adata(NT), KSp(nk), HRp(nk)


	! generate aggregate shock series

	call GenerateShock(Pdata)


	! start generating data

	Kdata(1) = kss

	do t = 1, NT

		do ik = 1, nk
			KSp(ik) = polyinterp1(dlog(Pdata(t)), pgrid, KS(ik,:))
			HRp(ik) = polyinterp1(dlog(Pdata(t)), pgrid, HR(ik,:))
		end do

		call dcsint(nk, kgrid, KSp, kbreak, kcscoef)
		call dcsint(nk, kgrid, HRp, hbreak, hcscoef)

		Kdata(t+1) = dcsval(Kdata(t), nk-1, kbreak, kcscoef)
		Edata(t) = dcsval(Kdata(t), nk-1, hbreak, hcscoef)
		Ydata(t) = Prodf(Kdata(t), Edata(t), Pdata(t))
		Idata(t) = Kdata(t+1) - (1.0D0-delta)*Kdata(t)
		Cdata(t) = Ydata(t) - Idata(t)
		Rdata(t) = Pdata(t)*alpha*(Kdata(t)/Edata(t))**(alpha-1.0D0)
		Wdata(t) = Pdata(t)*(1.0D0-alpha)*(Kdata(t)/Edata(t))**alpha
		Adata(t) = Ydata(t)/Edata(t)
	
	end do


	! save artificial data

	open(1, file='Output\Ydata.txt', status='unknown')
	open(2, file='Output\Kdata.txt', status='unknown')
	open(3, file='Output\Cdata.txt', status='unknown')
	open(4, file='Output\Idata.txt', status='unknown')
	open(5, file='Output\Pdata.txt', status='unknown')
	open(6, file='Output\Edata.txt', status='unknown')
	open(7, file='Output\Rdata.txt', status='unknown')
	open(8, file='Output\Wdata.txt', status='unknown')
	open(9, file='Output\Adata.txt', status='unknown')
	open(0, file='Output\Ldata.txt', status='unknown')

	write(1, '(f12.6)') Ydata
	write(2, '(f12.6)') Kdata
	write(3, '(f12.6)') Cdata
	write(4, '(f12.6)') Idata
	write(5, '(f12.6)') Pdata
	write(6, '(f12.6)') Edata
	write(7, '(f12.6)') Rdata
	write(8, '(f12.6)') Wdata
	write(9, '(f12.6)') Adata
	write(0, '(f12.6)') Edata


end subroutine




subroutine GenerateShock(Pdata)

	implicit none

	integer, parameter:: aseed = 1
	integer i, j, k
	real(8) ushock, trprob, Pdata(NT)


	! set the seed number for random number generator

	call rnset(aseed)


	!	normal random numbers for aggregate productivity shocks

	call drnnor(NT, Pdata)
	call dscal(NT, sig, Pdata, 1)
		

	! generate Pdata

	do i = 2, NT
		Pdata(i) = min(max(rhop*Pdata(i-1) + Pdata(i), pgrid(1)), pgrid(np))
	end do

	Pdata = dexp(Pdata)


end subroutine


end module
