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
	open(4, file='Output\HRmu.txt', status='unknown')
	open(5, file='Output\EFmu.txt', status='unknown')

	write(1, '(<nx>f15.6)') ((VF(ia,ix),ix=1,nx),ia=1,na)
	write(2, '(<nx>f12.6)') ((AS(ia,ix),ix=1,nx),ia=1,na)
	write(3, '(<nx>f10.6)') ((HR(ia,ix),ix=1,nx),ia=1,na)
	write(4, '(<nx>f10.6)') ((HRmu(ia,ix),ix=1,nx),ia=1,namu)
	write(5, '(<nx>f10.6)') ((EFmu(ia,ix),ix=1,nx),ia=1,namu)


	!	save the measure 

	mu = 100.0D0*mu
	cmu = 100.0D0*cmu

	open(8, file='Output\mu.txt', status='unknown')
	open(9, file='Output\cmu.txt', status='unknown')
	write(8, '(<nx>f12.6)') ((mu(ia,ix),ix=1,nx),ia=1,namu)
	write(9, '(<nx>f12.6)') ((cmu(ia,ix),ix=1,nx),ia=1,namu)


end subroutine



end module
