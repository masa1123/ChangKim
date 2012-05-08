subroutine SaveTimeSeries()

	use Globals

	implicit none


	!	save simulated time series data


	open(0, file='Input\ExpectationParameters.txt', status='unknown')
	open(1, file='Output\phi.txt', status='unknown')
	open(2, file='Output\psi.txt', status='unknown')
	open(3, file='Output\Kdata.txt', status='unknown')
	open(4, file='Output\Pdata.txt', status='unknown')
	open(5, file='Output\Cdata.txt', status='unknown')
	open(6, file='Output\Ldata.txt', status='unknown')
	open(7, file='Output\Ydata.txt', status='unknown')
	open(8, file='Output\Idata.txt', status='unknown')
	open(9, file='Output\Wdata.txt', status='unknown')
	open(10, file='Output\Rdata.txt', status='unknown')
	open(11, file='Output\Xdata.txt', status='unknown')
	open(12, file='Output\Adata.txt', status='unknown')
	open(13, file='Output\Edata.txt', status='unknown')
	open(14, file='Output\EINdata.txt', status='unknown')
	open(15, file='Output\EOUTdata.txt', status='unknown')
	open(16, file='Output\hEINdata.txt', status='unknown')
	open(17, file='Output\hEOUTdata.txt', status='unknown')

	write(0, '(f20.15)') coef
	write(1, '(f12.6)') (phiQ(t), t=t_LS,nT-1)
	write(2, '(f12.6)') (psiQ(t), t=t_LS,nT-1)
	write(3, '(f12.6)') Kdata
	write(4, '(f12.6)') Pdata
	write(5, '(f12.6)') Cdata
	write(6, '(f12.6)') Ldata
	write(7, '(f12.6)') Ydata
	write(8, '(f12.6)') Idata
	write(9, '(f12.6)') Wdata
	write(10, '(f12.6)') Rdata
	write(11, '(f12.6)') Xdata
	write(12, '(f12.6)') Adata
	write(13, '(f12.6)') Edata
	write(14, '(f12.6)') EINdata
	write(15, '(f12.6)') EOUTdata
	write(16, '(f12.6)') hEINdata
	write(17, '(f12.6)') hEOUTdata


end subroutine
