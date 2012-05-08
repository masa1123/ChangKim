!	----------------------------------------------------------------------
!	Module name: RegressionModule.f90
!
!	regresses simulated data to get new law of motion for aggregate
!	capital and price equations.
!	----------------------------------------------------------------------


module RegressionModule

	use Globals
	use rlse_int
	use rstat_int

	implicit none

	integer, parameter:: nobs = Nperiod-Nskip
	real(8) RHS(nobs,2)
	
contains



subroutine RegressLOM()

	implicit none

	real(8) sst, sse

	
	!	construct common regressors

	RHS(:,1) = dlog(Kdata(Nskip+1:Nperiod))
	RHS(:,2) = dlog(Pdata(Nskip+1:Nperiod))


	!	law of motion for aggregate capital

	call d_rlse(dlog(Kdata(Nskip+2:Nperiod+1)), RHS, NewKcoef, sst=sst, sse=sse)


	!	interest rate eqeuation

	call d_rlse(dlog(Rdata(Nskip+1:Nperiod)), RHS, NewRcoef, sst=sst, sse=sse)


	!	law of motion for aggregate capital

	call d_rlse(dlog(Wdata(Nskip+1:Nperiod)), RHS, NewWcoef, sst=sst, sse=sse)


end subroutine




!	----------------------------------------------------------------------
!	File name: FinalRegressLOM.f90
!
!	regresses simulated data to get new law of motion for aggregate
!	capital and price equations and summary statistics.
!	----------------------------------------------------------------------


subroutine FinalRegressLOM()

	implicit none

	real(8) sst, sse, dfe, R(3,3), WK(16), AOV(15), COEF(3,5), COVB(3,3), SQSS(2,4)
	integer nrmiss, irbef(1)
	
	
	!	construct common regressors

	RHS(:,1) = dlog(Kdata(Nskip+1:Nperiod))
	RHS(:,2) = dlog(Pdata(Nskip+1:Nperiod))


	!	open the log file

	open(0, file='Output\Logfile.txt', status='unknown')


	!	final regressions with summary statistics

	!	capital

	call dr2se(nobs, dlog(Kdata(Nskip+2:Nperiod+1)), 2, RHS, nobs, 1, NewKcoef, sst, sse, R, 3, dfe, nrmiss, WK)
	call d_rstat(irbef, NewKcoef, R, dfe, sse, AOV, SQSS, COEF, COVB, ief = -2, print='A')

	write(0, '(//A/)')    "Summary statistics for Converged Capital Equation"
	write(0, '(A,F20.4)') "Degrees of freedom for regression:               ", AOV(1)
	write(0, '(A,F20.4)') "Degrees of freedom for error:                    ", AOV(2)
	write(0, '(A,F20.4)') "Total degrees of freedom:                        ", AOV(3)
	write(0, '(A,F20.4)') "Sum of squares for regression:                   ", AOV(4)
	write(0, '(A,F20.4)') "Sum of squares for error:                        ", AOV(5)
	write(0, '(A,F20.4)') "Total sum of squares:                            ", AOV(6)
	write(0, '(A,F20.4)') "Regression mean square:                          ", AOV(7)
	write(0, '(A,F20.4)') "Error mean square:                               ", AOV(8)
	write(0, '(A,F20.4)') "F-statistic:                                     ", AOV(9)
	write(0, '(A,F20.4)') "p-value:                                         ", AOV(10)
	write(0, '(A,F20.4)') "R^2 (in percent):                                ", AOV(11)
	write(0, '(A,F20.4)') "Adjusted R^2 (in percent):                       ", AOV(12)
	write(0, '(A,F20.4)') "Estimated standard deviation of the model error: ", AOV(13)
	write(0, '(A,F20.4)') "Mean of the response (dependent) variable:       ", AOV(14)
	write(0, '(A,F20.4)') "Coefficient of variation (in percent):           ", AOV(15)

	write(0, '(/A/)')    "Inference on Coefficients"
	write(0, '(A)') "Varialbe       Estimate        St. err         t-stat          p-val"
    write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " Const  ", COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4)
	write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " K      ", COEF(2,1), COEF(2,2), COEF(2,3), COEF(2,4)
	write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " p      ", COEF(3,1), COEF(3,2), COEF(3,3), COEF(3,4)


	!	interest rate

	call dr2se(nobs, dlog(Rdata(Nskip+1:Nperiod)), 2, RHS, nobs, 1, NewRcoef, sst, sse, R, 3, dfe, nrmiss, WK)
	call d_rstat(irbef, NewRcoef, R, dfe, sse, AOV, SQSS, COEF, COVB, ief = -2, print='A')

	write(0, '(//A/)')    "Summary statistics for Converged Interest Rate Equation"
	write(0, '(A,F20.4)') "Degrees of freedom for regression:               ", AOV(1)
	write(0, '(A,F20.4)') "Degrees of freedom for error:                    ", AOV(2)
	write(0, '(A,F20.4)') "Total degrees of freedom:                        ", AOV(3)
	write(0, '(A,F20.4)') "Sum of squares for regression:                   ", AOV(4)
	write(0, '(A,F20.4)') "Sum of squares for error:                        ", AOV(5)
	write(0, '(A,F20.4)') "Total sum of squares:                            ", AOV(6)
	write(0, '(A,F20.4)') "Regression mean square:                          ", AOV(7)
	write(0, '(A,F20.4)') "Error mean square:                               ", AOV(8)
	write(0, '(A,F20.4)') "F-statistic:                                     ", AOV(9)
	write(0, '(A,F20.4)') "p-value:                                         ", AOV(10)
	write(0, '(A,F20.4)') "R^2 (in percent):                                ", AOV(11)
	write(0, '(A,F20.4)') "Adjusted R^2 (in percent):                       ", AOV(12)
	write(0, '(A,F20.4)') "Estimated standard deviation of the model error: ", AOV(13)
	write(0, '(A,F20.4)') "Mean of the response (dependent) variable:       ", AOV(14)
	write(0, '(A,F20.4)') "Coefficient of variation (in percent):           ", AOV(15)

	write(0, '(/A/)')    "Inference on Coefficients"
	write(0, '(A)') "Varialbe       Estimate        St. err         t-stat          p-val"
    write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " Const  ", COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4)
	write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " K      ", COEF(2,1), COEF(2,2), COEF(2,3), COEF(2,4)
	write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " p      ", COEF(3,1), COEF(3,2), COEF(3,3), COEF(3,4)


	!	wage

	call dr2se(nobs, dlog(Wdata(Nskip+1:Nperiod)), 2, RHS, nobs, 1, NewWcoef, sst, sse, R, 3, dfe, nrmiss, WK)
	call d_rstat(irbef, NewWcoef, R, dfe, sse, AOV, SQSS, COEF, COVB, ief = -2, print='A')

	write(0, '(//A/)')    "Summary statistics for Converged Wage Equation"
	write(0, '(A,F20.4)') "Degrees of freedom for regression:               ", AOV(1)
	write(0, '(A,F20.4)') "Degrees of freedom for error:                    ", AOV(2)
	write(0, '(A,F20.4)') "Total degrees of freedom:                        ", AOV(3)
	write(0, '(A,F20.4)') "Sum of squares for regression:                   ", AOV(4)
	write(0, '(A,F20.4)') "Sum of squares for error:                        ", AOV(5)
	write(0, '(A,F20.4)') "Total sum of squares:                            ", AOV(6)
	write(0, '(A,F20.4)') "Regression mean square:                          ", AOV(7)
	write(0, '(A,F20.4)') "Error mean square:                               ", AOV(8)
	write(0, '(A,F20.4)') "F-statistic:                                     ", AOV(9)
	write(0, '(A,F20.4)') "p-value:                                         ", AOV(10)
	write(0, '(A,F20.4)') "R^2 (in percent):                                ", AOV(11)
	write(0, '(A,F20.4)') "Adjusted R^2 (in percent):                       ", AOV(12)
	write(0, '(A,F20.4)') "Estimated standard deviation of the model error: ", AOV(13)
	write(0, '(A,F20.4)') "Mean of the response (dependent) variable:       ", AOV(14)
	write(0, '(A,F20.4)') "Coefficient of variation (in percent):           ", AOV(15)

	write(0, '(/A/)')    "Inference on Coefficients"
	write(0, '(A)') "Varialbe       Estimate        St. err         t-stat          p-val"
    write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " Const  ", COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4)
	write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " K      ", COEF(2,1), COEF(2,2), COEF(2,3), COEF(2,4)
	write(0, '(A, F15.6, F15.6, F15.6, F15.6)') " p      ", COEF(3,1), COEF(3,2), COEF(3,3), COEF(3,4)


end subroutine


!	----------------------------------------------------------------------
!	File name: FinalRegressLOM2.f90
!	----------------------------------------------------------------------


subroutine FinalRegressLOM2()

	implicit none

	real(8) sst, sse, dfe, R(4,4), WK(21), AOV(15), COEF(4,5), COVB(4,4), SQSS(3,4)
	integer nrmiss, irbef(1)
	real(8) RHS(nobs,3), K2coef(4), R2coef(4), W2coef(4)
	

	!	construct common regressors

	RHS(:,1) = dlog(Kdata(Nskip+1:Nperiod))
	RHS(:,2) = dlog(K2data(Nskip+1:Nperiod))
	RHS(:,3) = dlog(Pdata(Nskip+1:Nperiod))


	!	open the log file

	open(0, file='Output\Logfile.txt', status='unknown')


	!	final regressions with summary statistics

	!	capital

	call dr2se(nobs, dlog(Kdata(Nskip+2:Nperiod+1)), 3, RHS, nobs, 1, K2coef, sst, sse, R, 4, dfe, nrmiss, WK)
	call d_rstat(irbef, K2coef, R, dfe, sse, AOV, SQSS, COEF, COVB, ief = -3, print='A')

	write(0, '(//A/)')    "Summary statistics for Converged Capital Equation"
	write(0, '(A,F20.4)') "Degrees of freedom for regression:               ", AOV(1)
	write(0, '(A,F20.4)') "Degrees of freedom for error:                    ", AOV(2)
	write(0, '(A,F20.4)') "Total degrees of freedom:                        ", AOV(3)
	write(0, '(A,F20.4)') "Sum of squares for regression:                   ", AOV(4)
	write(0, '(A,F20.4)') "Sum of squares for error:                        ", AOV(5)
	write(0, '(A,F20.4)') "Total sum of squares:                            ", AOV(6)
	write(0, '(A,F20.4)') "Regression mean square:                          ", AOV(7)
	write(0, '(A,F20.4)') "Error mean square:                               ", AOV(8)
	write(0, '(A,F20.4)') "F-statistic:                                     ", AOV(9)
	write(0, '(A,F20.4)') "p-value:                                         ", AOV(10)
	write(0, '(A,F20.4)') "R^2 (in percent):                                ", AOV(11)
	write(0, '(A,F20.4)') "Adjusted R^2 (in percent):                       ", AOV(12)
	write(0, '(A,F20.4)') "Estimated standard deviation of the model error: ", AOV(13)
	write(0, '(A,F20.4)') "Mean of the response (dependent) variable:       ", AOV(14)
	write(0, '(A,F20.4)') "Coefficient of variation (in percent):           ", AOV(15)

	write(0, '(/A/)')    "Inference on Coefficients"
	write(0, '(A)') "Varialbe       Estimate        St. err         t-stat          p-val"
    write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " Const  ", COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " K      ", COEF(2,1), COEF(2,2), COEF(2,3), COEF(2,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " K2     ", COEF(3,1), COEF(3,2), COEF(3,3), COEF(3,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " p      ", COEF(4,1), COEF(4,2), COEF(4,3), COEF(4,4)


	!	interest rate

	call dr2se(nobs, dlog(Rdata(Nskip+1:Nperiod)), 3, RHS, nobs, 1, R2coef, sst, sse, R, 4, dfe, nrmiss, WK)
	call d_rstat(irbef, R2coef, R, dfe, sse, AOV, SQSS, COEF, COVB, ief = -3, print='A')

	write(0, '(//A/)')    "Summary statistics for Converged Interest Rate Equation"
	write(0, '(A,F20.4)') "Degrees of freedom for regression:               ", AOV(1)
	write(0, '(A,F20.4)') "Degrees of freedom for error:                    ", AOV(2)
	write(0, '(A,F20.4)') "Total degrees of freedom:                        ", AOV(3)
	write(0, '(A,F20.4)') "Sum of squares for regression:                   ", AOV(4)
	write(0, '(A,F20.4)') "Sum of squares for error:                        ", AOV(5)
	write(0, '(A,F20.4)') "Total sum of squares:                            ", AOV(6)
	write(0, '(A,F20.4)') "Regression mean square:                          ", AOV(7)
	write(0, '(A,F20.4)') "Error mean square:                               ", AOV(8)
	write(0, '(A,F20.4)') "F-statistic:                                     ", AOV(9)
	write(0, '(A,F20.4)') "p-value:                                         ", AOV(10)
	write(0, '(A,F20.4)') "R^2 (in percent):                                ", AOV(11)
	write(0, '(A,F20.4)') "Adjusted R^2 (in percent):                       ", AOV(12)
	write(0, '(A,F20.4)') "Estimated standard deviation of the model error: ", AOV(13)
	write(0, '(A,F20.4)') "Mean of the response (dependent) variable:       ", AOV(14)
	write(0, '(A,F20.4)') "Coefficient of variation (in percent):           ", AOV(15)

	write(0, '(/A/)')    "Inference on Coefficients"
	write(0, '(A)') "Varialbe       Estimate        St. err         t-stat          p-val"
    write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " Const  ", COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " K      ", COEF(2,1), COEF(2,2), COEF(2,3), COEF(2,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " K2     ", COEF(3,1), COEF(3,2), COEF(3,3), COEF(3,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " p      ", COEF(4,1), COEF(4,2), COEF(4,3), COEF(4,4)


	!	wage

	call dr2se(nobs, dlog(Wdata(Nskip+1:Nperiod)), 3, RHS, nobs, 1, W2coef, sst, sse, R, 4, dfe, nrmiss, WK)
	call d_rstat(irbef, W2coef, R, dfe, sse, AOV, SQSS, COEF, COVB, ief = -3, print='A')

	write(0, '(//A/)')    "Summary statistics for Converged Wage Equation"
	write(0, '(A,F20.4)') "Degrees of freedom for regression:               ", AOV(1)
	write(0, '(A,F20.4)') "Degrees of freedom for error:                    ", AOV(2)
	write(0, '(A,F20.4)') "Total degrees of freedom:                        ", AOV(3)
	write(0, '(A,F20.4)') "Sum of squares for regression:                   ", AOV(4)
	write(0, '(A,F20.4)') "Sum of squares for error:                        ", AOV(5)
	write(0, '(A,F20.4)') "Total sum of squares:                            ", AOV(6)
	write(0, '(A,F20.4)') "Regression mean square:                          ", AOV(7)
	write(0, '(A,F20.4)') "Error mean square:                               ", AOV(8)
	write(0, '(A,F20.4)') "F-statistic:                                     ", AOV(9)
	write(0, '(A,F20.4)') "p-value:                                         ", AOV(10)
	write(0, '(A,F20.4)') "R^2 (in percent):                                ", AOV(11)
	write(0, '(A,F20.4)') "Adjusted R^2 (in percent):                       ", AOV(12)
	write(0, '(A,F20.4)') "Estimated standard deviation of the model error: ", AOV(13)
	write(0, '(A,F20.4)') "Mean of the response (dependent) variable:       ", AOV(14)
	write(0, '(A,F20.4)') "Coefficient of variation (in percent):           ", AOV(15)

	write(0, '(/A/)')    "Inference on Coefficients"
	write(0, '(A)') "Varialbe       Estimate        St. err         t-stat          p-val"
    write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " Const  ", COEF(1,1), COEF(1,2), COEF(1,3), COEF(1,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " K      ", COEF(2,1), COEF(2,2), COEF(2,3), COEF(2,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " K2     ", COEF(3,1), COEF(3,2), COEF(3,3), COEF(3,4)
	write(0, '(A, F15.10, F15.6, F15.6, F15.6)') " p      ", COEF(4,1), COEF(4,2), COEF(4,3), COEF(4,4)


end subroutine



end module