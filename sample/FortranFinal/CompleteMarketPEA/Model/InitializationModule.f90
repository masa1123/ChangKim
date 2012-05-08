module InitializationModule

	use Globals

contains


subroutine InitializeModelParameters()

	implicit none

	!	read model parameters

	open(1, file='Input\ModelParameters.txt', status = 'old')
	read(1, *) gamac
	read(1, *) gamal
	read(1, *) rhox
	read(1, *) sigex
	read(1, *) rhop
	read(1, *) sigep
	close(1)

	sigx = sigex/dsqrt(1.0D0-rhox**2.0)
	sigp = sigep/dsqrt(1.0D0-rhop**2.0)

end subroutine



subroutine InitializeExpectationParameters()

	implicit none


	!	read initial coefficients of parameterized expectation

	open(1, file='Input\ExpectationParameters.txt', status = 'old')
	read(1, *) coef
	close(1)

end subroutine


end module
