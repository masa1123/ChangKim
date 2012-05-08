!	----------------------------------------------------------------------
!	Module name: SaveResultsModule.f90
!
!	Save all results
!	----------------------------------------------------------------------


module SaveResultsModule

	use Globals

	implicit none


contains


subroutine SaveEquilibrium()

	implicit none


	!	save value function and decision rules

	open(1, file='Output\VF.txt', status='unknown')
	open(2, file='Output\AS.txt', status='unknown')
	open(3, file='Output\HR.txt', status='unknown')
	open(4, file='Output\VFE.txt', status='unknown')
	open(5, file='Output\ASE.txt', status='unknown')
	open(6, file='Output\VFN.txt', status='unknown')
	open(7, file='Output\ASN.txt', status='unknown')

	write(1, '(<nx>f15.6)') ((VF(ia,ix),ix=1,nx),ia=1,na)
	write(2, '(<nx>f8.3)') ((agrid(AS(ia,ix)),ix=1,nx),ia=1,na)
	write(3, '(<nx>f8.3)') ((HR(ia,ix),ix=1,nx),ia=1,na)
	write(4, '(<nx>f15.6)') ((VFE(ia,ix),ix=1,nx),ia=1,na)
	write(5, '(<nx>f8.3)') ((agrid(ASE(ia,ix)),ix=1,nx),ia=1,na)
	write(6, '(<nx>f15.6)') ((VFN(ia,ix),ix=1,nx),ia=1,na)
	write(7, '(<nx>f8.3)') ((agrid(ASN(ia,ix)),ix=1,nx),ia=1,na)


	!	save the measure 

	mu = 100.0D0*mu
	cmu = 100.0D0*cmu

	open(8, file='Output\mu.txt', status='unknown')
	open(9, file='Output\cmu.txt', status='unknown')
	write(8, '(<nx>f12.6)') ((mu(ia,ix),ix=1,nx),ia=1,na)
	write(9, '(<nx>f12.6)') ((cmu(ia,ix),ix=1,nx),ia=1,na)


end subroutine


subroutine SaveDistribution()

	open(1, file='Output\AsstEarnDensity.txt', status='unknown')
	open(2, file='Output\AsstDensity.txt', status='unknown')
	open(3, file='Output\EarnDensity.txt', status='unknown')
	open(4, file='Output\AsstEarnDist.txt', status='unknown')
	open(5, file='Output\AsstDist.txt', status='unknown')
	open(6, file='Output\EarnDist.txt', status='unknown')

	write(1, '(<nx+1>F12.6)') ((AsstEarnDensity(ia,ix),ix=0,nx),ia=1,na)
	write(2, '(F12.6)') AsstDensity
	write(3, '(F12.6)') EarnDensity
	write(4, '(<nx+1>F12.6)') ((AsstEarnDist(ia,ix),ix=0,nx),ia=1,na)
	write(5, '(F12.6)') AsstDist
	write(6, '(F12.6)') EarnDist


end subroutine


subroutine SaveReservationWages()

	open(1, file='Output\ReservationX.txt', status='unknown')
	open(2, file='Output\ReservationWages.txt', status='unknown')

	write(1, '(F12.6)') ReservationX
	write(2, '(F12.6)') ReservationWage

end subroutine



subroutine SaveLaborSupply()

	open(1, file='Output\LS.txt', status='unknown')
	write(1, '(F12.6)') LS

end subroutine


end module
