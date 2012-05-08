!	----------------------------------------------------------------------
!	File name: SaveValueFunction.f90
!	----------------------------------------------------------------------
	
	

subroutine SaveCoefficients()

	use Globals

	open(1, file='Output\CoefficientsOUT.txt', status='unknown')
	write(1, '(3F10.6)') Kcoef(1), Kcoef(2), Kcoef(3)
	write(1, '(3F10.6)') Rcoef(1), Rcoef(2), Rcoef(3)
	write(1, '(3F10.6)') Wcoef(1), Wcoef(2), Wcoef(3)

end subroutine



subroutine SaveSimulatedTimeSeriesData()

	use Globals

	open(1, file='Output\Kdata.txt', status='unknown')
	open(2, file='Output\Rdata.txt', status='unknown')
	open(3, file='Output\Wdata.txt', status='unknown')
	open(4, file='Output\Ldata.txt', status='unknown')
	open(5, file='Output\Edata.txt', status='unknown')
	open(8, file='Output\Ydata.txt', status='unknown')
	open(9, file='Output\Idata.txt', status='unknown')
	open(10, file='Output\Cdata.txt', status='unknown')
	open(11, file='Output\Adata.txt', status='unknown')
	open(12, file='Output\K2data.txt', status='unknown')

	write(1, '(F12.6)') Kdata
	write(2, '(F12.6)') Rdata
	write(3, '(F12.6)') Wdata
	write(4, '(F12.6)') Ldata
	write(5, '(F12.6)') Edata
	write(8, '(F12.6)') Ydata
	write(9, '(F12.6)') Idata
	write(10, '(F12.6)') Cdata
	write(11, '(F12.6)') Adata
	write(12, '(F12.6)') K2data


end subroutine



subroutine SaveValueFunction()

	use Globals


	open(1, file='Output\VF.txt', status='unknown')
	open(2, file='Output\AS.txt', status='unknown')
	open(3, file='Output\HR.txt', status='unknown')
	

	do ip = 1, np
	do ik = 1, nk
			
		write(1, '(A,I2,A,I2,/)') "ik = ", ik, "     ip = ", ip
		write(1, '(<nx>F15.6)') ((VF(ia,ix,ik,ip),ix=1,nx),ia=1,na)
		write(1, '(A,/)')

		write(2, '(A,I2,A,I2,/)') "ik = ", ik, "     ip = ", ip
		write(2, '(<nx>I6)') ((AS(ia,ix,ik,ip),ix=1,nx),ia=1,na)
		write(2, '(A,/)')

		write(3, '(A,I2,A,I2,/)') "ik = ", ik, "     ip = ", ip
		write(3, '(<nx>F8.3)') ((HR(ia,ix,ik,ip),ix=1,nx),ia=1,na)
		write(3, '(A,/)')

	end do
	end do


end subroutine

