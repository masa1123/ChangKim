!	----------------------------------------------------------------------
!	File name: SaveValueFunction.f90
!	----------------------------------------------------------------------
	
	

subroutine SaveCoefficients()

	use Globals

	open(1, file='Output\CoefficientsOUT.txt', status='unknown')
	write(1, '(<ncoef>F10.6)') Kcoef
	write(1, '(<ncoef>F10.6)') K2coef
	write(1, '(<ncoef>F10.6)') Rcoef
	write(1, '(<ncoef>F10.6)') Wcoef
	close(1)
	
	
end subroutine



subroutine SaveSimulatedTimeSeriesData()

	use Globals

	open(1, file='Output\Kdata.txt', status='unknown')
	open(2, file='Output\Rdata.txt', status='unknown')
	open(3, file='Output\Wdata.txt', status='unknown')
	open(4, file='Output\Ldata.txt', status='unknown')
	open(5, file='Output\Edata.txt', status='unknown')
	open(6, file='Output\EINdata.txt', status='unknown')
	open(7, file='Output\EOUTdata.txt', status='unknown')
	open(8, file='Output\Ydata.txt', status='unknown')
	open(9, file='Output\Idata.txt', status='unknown')
	open(10, file='Output\Cdata.txt', status='unknown')
	open(11, file='Output\Adata.txt', status='unknown')
	open(12, file='Output\hEINdata.txt', status='unknown')
	open(13, file='Output\hEOUTdata.txt', status='unknown')
	open(14, file='Output\Pdata.txt', status='replace')
	open(15, file='Output\K2data.txt', status='unknown')

	write(1, '(F12.6)') Kdata
	write(2, '(F12.6)') Rdata
	write(3, '(F12.6)') Wdata
	write(4, '(F12.6)') Ldata
	write(5, '(F12.6)') Edata
	write(6, '(F12.6)') EINdata
	write(7, '(F12.6)') EOUTdata
	write(8, '(F12.6)') Ydata
	write(9, '(F12.6)') Idata
	write(10, '(F12.6)') Cdata
	write(11, '(F12.6)') Adata
	write(12, '(F12.6)') hEINdata
	write(13, '(F12.6)') hEOUTdata
	write(14, '(f12.6)') Pdata
	write(15, '(F12.6)') K2data

	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	close(6)
	close(7)
	close(8)
	close(9)
	close(10)
	close(11)
	close(12)
	close(13)
	close(14)
	close(15)


end subroutine



subroutine SaveValueFunction2()

	use Globals


	open(1, file='Output\VF.txt', status='unknown')
	open(2, file='Output\AS.txt', status='unknown')
	open(3, file='Output\HR.txt', status='unknown')
	open(4, file='Output\VFE.txt', status='unknown')
	open(5, file='Output\VFN.txt', status='unknown')
	open(6, file='Output\ASE.txt', status='unknown')
	open(7, file='Output\ASN.txt', status='unknown')
	

	write(1, '(<nx>F15.6)') (((((VF(ia,ix,ik,ik2,ip),ix=1,nx),ia=1,na),ik=1,nk),ik2=1,nk2),ip=1,np)
	write(2, '(<nx>I6)')    (((((AS(ia,ix,ik,ik2,ip),ix=1,nx),ia=1,na),ik=1,nk),ik2=1,nk2),ip=1,np)
	write(3, '(<nx>F8.3)')  (((((HR(ia,ix,ik,ik2,ip),ix=1,nx),ia=1,na),ik=1,nk),ik2=1,nk2),ip=1,np)
	write(4, '(<nx>F15.6)') (((((VFE(ia,ix,ik,ik2,ip),ix=1,nx),ia=1,na),ik=1,nk),ik2=1,nk2),ip=1,np)
	write(5, '(<nx>F15.6)') (((((VFN(ia,ix,ik,ik2,ip),ix=1,nx),ia=1,na),ik=1,nk),ik2=1,nk2),ip=1,np)
	write(6, '(<nx>I6)')    (((((ASE(ia,ix,ik,ik2,ip),ix=1,nx),ia=1,na),ik=1,nk),ik2=1,nk2),ip=1,np)
	write(7, '(<nx>I6)')    (((((ASN(ia,ix,ik,ik2,ip),ix=1,nx),ia=1,na),ik=1,nk),ik2=1,nk2),ip=1,np)


end subroutine
