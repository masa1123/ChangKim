!	----------------------------------------------------------------------
!	File name: GenerateShock.cpp
!	----------------------------------------------------------------------


subroutine GenerateShock(Zdata)

	use Numerical_Libraries
	use Globals

	integer, parameter:: aseed = 1
	integer i, j, k
	real(8) ushock, trprob, Zdata(NT)


	! set the seed number for random number generator

	call rnset(aseed)


	! generate Zdata

	j = 1

	do i = 1, NT

		ushock = drnunf()
		trprob = 0.0D0

		do k = 1, nz

			trprob = trprob + trZ(j,k)
			if (ushock <= trprob) then
				j = k
				exit
			end if

		end do

		Zdata(i) = ezgrid(j)
	
	end do

end subroutine

