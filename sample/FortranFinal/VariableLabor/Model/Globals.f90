!	----------------------------------------------------------------------
!	File name: Globals.f90
!
!	Defines global variables.
!	----------------------------------------------------------------------


module Globals

	implicit none

	integer, parameter:: nk = 31, np = 9, NT = 3500
	real(8), parameter:: alpha = 0.36d0, delta = 0.025d0, irate = 0.01d0,	&
						 avep = 0.0D0, rhop = 0.95D0, sig = 0.007D0
	real(8), parameter:: rd = irate + delta
	real(8), parameter:: tol_cnvrg = 1.0D-6
	real(8) B, gama, kss, bta, sigp
	real(8) kgrid(nk), VF(nk,np), KS(nk,np), HR(nk,np), pgrid(np),			&
			epgrid(np), trP(np,np), hour, kp_hour, break(nk), cscoef(4,nk), &
			kbreak(nk), kcscoef(4,nk), hbreak(nk), hcscoef(4,nk), Value
	integer ik, ip

end module



!	when changing np, make appropriate changes in 
!	1. ShockGrid()
!	2. part that generates transition matrix in Main()