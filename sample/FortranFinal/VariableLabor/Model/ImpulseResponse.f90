!	----------------------------------------------------------------------
!	File name: SimulateData.cpp
!	----------------------------------------------------------------------


subroutine ImpulseResponse()

	use Globals
	use ValueModule

	implicit none

	integer t
	real(8) Ydata(NT), Kdata(NT+1), Hdata(NT), Cdata(NT), Idata(NT), Zdata(NT),	&
			Rdata(NT), Wdata(NT), &
			Zkbreak(nk,2), Zkcscoef(4,nk,2), Zhbreak(nk,2), Zhcscoef(4,nk,2)


	! spline interpolation of decision rule for a normal state
	
	call dcsint(nk, kgrid, KS(:,5), kbreak, kcscoef)
	Zkcscoef(:,:,1) = kcscoef
	Zkbreak(:,1) = kbreak

	call dcsint(nk, kgrid, HR(:,5), hbreak, hcscoef)
	Zhcscoef(:,:,1) = hcscoef
	Zhbreak(:,1) = hbreak


	! spline interpolation of decision rule for a 1.68% increase in productivity
	
	call dcsint(nk, kgrid, KS(:,6), kbreak, kcscoef)
	Zkcscoef(:,:,2) = kcscoef
	Zkbreak(:,2) = kbreak

	call dcsint(nk, kgrid, HR(:,6), hbreak, hcscoef)
	Zhcscoef(:,:,2) = hcscoef
	Zhbreak(:,2) = hbreak


	! aggregate shock series

	Zdata = ezgrid(5)
	Zdata(101:NT) = ezgrid(6)


	! start generating data
	
	Kdata(1) = kss

	kcscoef = Zkcscoef(:,:,1)
	kbreak  = Zkbreak(:,1)
	hcscoef = Zhcscoef(:,:,1)
	hbreak  = Zhbreak(:,1)


	do t = 1, NT

		if (t == 101) then
			kcscoef = Zkcscoef(:,:,2)
			kbreak  = Zkbreak(:,2)
			hcscoef = Zhcscoef(:,:,2)
			hbreak  = Zhbreak(:,2)
!		else if (t >= 102) then 
!			kcscoef = Zkcscoef(:,:,1)
!			kbreak  = Zkbreak(:,1)
!			hcscoef = Zhcscoef(:,:,1)
!			hbreak  = Zhbreak(:,1)
		end if

		Kdata(t+1) = dcsval(Kdata(t), nk-1, kbreak, kcscoef)
		Hdata(t) = dcsval(Kdata(t), nk-1, hbreak, hcscoef)
		Ydata(t) = Prodf(Kdata(t), Hdata(t), Zdata(t))
		Idata(t) = Kdata(t+1) - (1.0D0-delta)*Kdata(t)
		Cdata(t) = Ydata(t) - Idata(t)
		Rdata(t) = Zdata(t)*alpha*(Kdata(t)/Hdata(t))**(alpha-1.0D0)
		Wdata(t) = Zdata(t)*(1.0D0-alpha)*(Kdata(t)/Hdata(t))**alpha
	
	end do


	! save artificial data

	open(1, file='Output\Ydata2.txt', status='unknown')
	open(2, file='Output\Kdata2.txt', status='unknown')
	open(3, file='Output\Cdata2.txt', status='unknown')
	open(4, file='Output\Idata2.txt', status='unknown')
	open(5, file='Output\Zdata2.txt', status='unknown')
	open(6, file='Output\Hdata2.txt', status='unknown')
	open(7, file='Output\Rdata2.txt', status='unknown')
	open(8, file='Output\Wdata2.txt', status='unknown')

	write(1, '(f12.6)') Ydata
	write(2, '(f12.6)') Kdata
	write(3, '(f12.6)') Cdata
	write(4, '(f12.6)') Idata
	write(5, '(f12.6)') Zdata
	write(6, '(f12.6)') Hdata
	write(7, '(f12.6)') Rdata
	write(8, '(f12.6)') Wdata

end subroutine
