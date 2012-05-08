!	----------------------------------------------------------------------
!	File name: AggregateCapital.cpp
!
!	Calculates the aggregate capital by numerically integrating 
!	the measure, and then interest rate.
!	----------------------------------------------------------------------


subroutine AggregateVariables()

	use Globals
	use dfport

	implicit none

	real(8) AK, AH, AE, AKE, irateb, wageb, Prod, Empl, EIN, EOUT,		&
			EmplQuintile(5), EINQuintile(5), EOUTQuintile(5), WealthMobility(5,5)
	integer i, k, kk


	!	aggregate capital, hours

	AK = dot_product(agrid, sum(mu, dim=2))
	AH = sum(HR*mu)
	AE = dot_product(exgrid, sum(HR*mu, dim=1))
	AKE = AK/AE

	
	!	real interest rate, wage

	irateb = alpha*AKE**(alpha-1.0D0) - delta
	wageb  = (1.0D0-alpha)*AKE**alpha


	!	output

	Prod = AK**alpha*AE**(1.0D0-alpha)


	!	employment rates

	Empl = sum(HR*mu/hbar)

	EmplQuintile(1) = sum(HR(1:asst20,:)*mu(1:asst20,:)/hbar)/sum(mu(1:asst20,:))
	EmplQuintile(2) = sum(HR(asst20+1:asst40,:)*mu(asst20+1:asst40,:)/hbar)/sum(mu(asst20+1:asst40,:))
	EmplQuintile(3) = sum(HR(asst40+1:asst60,:)*mu(asst40+1:asst60,:)/hbar)/sum(mu(asst40+1:asst60,:))
	EmplQuintile(4) = sum(HR(asst60+1:asst80,:)*mu(asst60+1:asst80,:)/hbar)/sum(mu(asst60+1:asst80,:))
	EmplQuintile(5) = sum(HR(asst80+1:na,:)*mu(asst80+1:na,:)/hbar)/sum(mu(asst80+1:na,:))


	!	worker flows in and out of employment

	EIN = 0.0D0
	EOUT = 0.0D0
	EINQuintile = 0.0D0
	EOUTQuintile = 0.0D0

	do k = 1, nx
	do i = 1, na

		do kk = 1, nx
			if (HR(i,k) == hbar .and. HR(AS(i,k),kk) == 0.0D0) then

				EOUT = EOUT + mu(i,k)*trX(k,kk)

				if (i <= asst20) then
					EOUTQuintile(1) = EOUTQuintile(1) + mu(i,k)*trX(k,kk)
				else if (i > asst20 .and. i <= asst40) then
					EOUTQuintile(2) = EOUTQuintile(2) + mu(i,k)*trX(k,kk)
				else if (i > asst40 .and. i <= asst60) then
					EOUTQuintile(3) = EOUTQuintile(3) + mu(i,k)*trX(k,kk)
				else if (i > asst60 .and. i <= asst80) then
					EOUTQuintile(4) = EOUTQuintile(4) + mu(i,k)*trX(k,kk)
				else 
					EOUTQuintile(5) = EOUTQuintile(5) + mu(i,k)*trX(k,kk)
				end if

			else if (HR(i,k) == 0.0D0 .and. HR(AS(i,k),kk) == hbar) then

				EIN = EIN + mu(i,k)*trX(k,kk)

				if (i <= asst20) then
					EINQuintile(1) = EINQuintile(1) + mu(i,k)*trX(k,kk)
				else if (i > asst20 .and. i <= asst40) then
					EINQuintile(2) = EINQuintile(2) + mu(i,k)*trX(k,kk)
				else if (i > asst40 .and. i <= asst60) then
					EINQuintile(3) = EINQuintile(3) + mu(i,k)*trX(k,kk)
				else if (i > asst60 .and. i <= asst80) then
					EINQuintile(4) = EINQuintile(4) + mu(i,k)*trX(k,kk)
				else 
					EINQuintile(5) = EINQuintile(5) + mu(i,k)*trX(k,kk)
				end if

			end if
		end do

	end do
	end do


	!	wealth mobility

	do k = 1, nx
	do i = 1, na

		if (i <= asst20 .and. AS(i,k) <= asst20) then
			WealthMobility(1,1) = WealthMobility(1,1) + mu(i,k)
		else if (i <= asst20 .and. AS(i,k) <= asst40) then
			WealthMobility(1,2) = WealthMobility(1,2) + mu(i,k)
		else if (i <= asst20 .and. AS(i,k) <= asst60) then
			WealthMobility(1,3) = WealthMobility(1,3) + mu(i,k)
		else if (i <= asst20 .and. AS(i,k) <= asst80) then
			WealthMobility(1,4) = WealthMobility(1,4) + mu(i,k)
		else if (i <= asst20 .and. AS(i,k) <= na) then
			WealthMobility(1,5) = WealthMobility(1,5) + mu(i,k)

		else if (i <= asst40 .and. AS(i,k) <= asst20) then
			WealthMobility(2,1) = WealthMobility(2,1) + mu(i,k)
		else if (i <= asst40 .and. AS(i,k) <= asst40) then
			WealthMobility(2,2) = WealthMobility(2,2) + mu(i,k)
		else if (i <= asst40 .and. AS(i,k) <= asst60) then
			WealthMobility(2,3) = WealthMobility(2,3) + mu(i,k)
		else if (i <= asst40 .and. AS(i,k) <= asst80) then
			WealthMobility(2,4) = WealthMobility(2,4) + mu(i,k)
		else if (i <= asst40 .and. AS(i,k) <= na) then
			WealthMobility(2,5) = WealthMobility(2,5) + mu(i,k)

		else if (i <= asst60 .and. AS(i,k) <= asst20) then
			WealthMobility(3,1) = WealthMobility(3,1) + mu(i,k)
		else if (i <= asst60 .and. AS(i,k) <= asst40) then
			WealthMobility(3,2) = WealthMobility(3,2) + mu(i,k)
		else if (i <= asst60 .and. AS(i,k) <= asst60) then
			WealthMobility(3,3) = WealthMobility(3,3) + mu(i,k)
		else if (i <= asst60 .and. AS(i,k) <= asst80) then
			WealthMobility(3,4) = WealthMobility(3,4) + mu(i,k)
		else if (i <= asst60 .and. AS(i,k) <= na) then
			WealthMobility(3,5) = WealthMobility(3,5) + mu(i,k)

		else if (i <= asst80 .and. AS(i,k) <= asst20) then
			WealthMobility(4,1) = WealthMobility(4,1) + mu(i,k)
		else if (i <= asst80 .and. AS(i,k) <= asst40) then
			WealthMobility(4,2) = WealthMobility(4,2) + mu(i,k)
		else if (i <= asst80 .and. AS(i,k) <= asst60) then
			WealthMobility(4,3) = WealthMobility(4,3) + mu(i,k)
		else if (i <= asst80 .and. AS(i,k) <= asst80) then
			WealthMobility(4,4) = WealthMobility(4,4) + mu(i,k)
		else if (i <= asst80 .and. AS(i,k) <= na) then
			WealthMobility(4,5) = WealthMobility(4,5) + mu(i,k)

		else if (i <= na .and. AS(i,k) <= asst20) then
			WealthMobility(5,1) = WealthMobility(5,1) + mu(i,k)
		else if (i <= na .and. AS(i,k) <= asst40) then
			WealthMobility(5,2) = WealthMobility(5,2) + mu(i,k)
		else if (i <= na .and. AS(i,k) <= asst60) then
			WealthMobility(5,3) = WealthMobility(5,3) + mu(i,k)
		else if (i <= na .and. AS(i,k) <= asst80) then
			WealthMobility(5,4) = WealthMobility(5,4) + mu(i,k)
		else if (i <= na .and. AS(i,k) <= na) then
			WealthMobility(5,5) = WealthMobility(5,5) + mu(i,k)
		end if

	end do
	end do

	WealthMobility(1,:) = WealthMobility(1,:)/sum(mu(1:asst20,:))
	WealthMobility(2,:) = WealthMobility(2,:)/sum(mu(asst20+1:asst40,:))
	WealthMobility(3,:) = WealthMobility(3,:)/sum(mu(asst40+1:asst60,:))
	WealthMobility(4,:) = WealthMobility(4,:)/sum(mu(asst60+1:asst80,:))
	WealthMobility(5,:) = WealthMobility(5,:)/sum(mu(asst80+1:na,:))


	!	print results on screen

	write(*, '(A,F13.10)') "beta     = ", bta
	write(*, '(A,F10.4)') "B        = ", B
	write(*, '(A,F10.4)') "ave_x    = ", avex
	write(*, '(A,F10.4)') "rho_x    = ", rhox
	write(*, '(A,F10.4)') "sig_x    = ", sigex
	write(*, *) " "
	print*, "Aggregate Capital                     = ", AK
	print*, "Aggregate Hours Worked                = ", AH
	print*, "Aggregate Effeciency Units of Hours   = ", AE
	print*, "Aggregate Capital per Effeciency Unit = ", AKE
	print*, "Aggregate Output                      = ", Prod
	print*, ""
	print*, "Discount Factor                       = ", bta
	print*, "Real Interest Rate                    = ", irateb
	print*, "Wage Rate per Effeciency Unit         = ", wageb
	print*, ""
	print*, "Employment/Population                 = ", Empl
	print*, "Inflow to Employment                  = ", EIN
	print*, "Outflow from Employment               = ", EOUT
	
	
	!	save results on file

	open(1, file='Output\Aggregates.txt', status='unknown')

	write(1, '(A/)') "Parameter Values:"
	write(1, '(A,F13.10)') "beta     = ", bta
	write(1, '(A,F10.4)') "B        = ", B
	write(1, '(A,F10.4)') "gamma    = ", gama
	write(1, '(A,F10.4)') "ave_x    = ", avex
	write(1, '(A,F10.4)') "rho_x    = ", rhox
	write(1, '(A,F10.4)') "sig_x    = ", sigex
	write(1, '(/)') 
	write(1, '(A/)') "Simulated Aggregate Variables:"
	write(1, '(A,F15.10)') "Aggregate Capital                     = ", AK
	write(1, '(A,F15.10)') "Aggregate Hours Worked                = ", AH
	write(1, '(A,F15.10)') "Aggregate Effeciency Units of Hours   = ", AE
	write(1, '(A,F15.10)') "Aggregate Capital per Effeciency Unit = ", AKE
	write(1, '(A,F15.10)') "Aggregate Output                      = ", Prod
	write(1, '(/)') 
	write(1, '(A,F15.10)') "Real Interest Rate                    = ", irateb
	write(1, '(A,F15.10)') "Wage Rate per Effeciency Unit         = ", wageb
	write(1, '(/)') 
	write(1, '(A,F15.10)') "Employment/Population                 = ", Empl
	write(1, '(A,F15.10," (",F8.5,")")') "Inflow to Employment                  = ", EIN, EIN/Empl
	write(1, '(A,F15.10," (",F8.5,")")') "Outflow from Employment               = ", EOUT, EOUT/Empl
	write(1, '(/)') 
	write(1, '(A/)') "Employment Rates by Asset Quintile"
	write(1, '(A,F15.10)') "Employment Rate of 1st Quintile          = ", EmplQuintile(1)
	write(1, '(A,F15.10)') "Employment Rate of 2nd Quintile          = ", EmplQuintile(2)
	write(1, '(A,F15.10)') "Employment Rate of 3rd Quintile          = ", EmplQuintile(3)
	write(1, '(A,F15.10)') "Employment Rate of 4th Quintile          = ", EmplQuintile(4)
	write(1, '(A,F15.10)') "Employment Rate of 5th Quintile          = ", EmplQuintile(5)
	write(1, '(/)') 
	write(1, '(A/)') "Worker Flows by Asset Quintile"
	write(1, '(A,F15.10," (",F8.5,")")') "Inflows to Employment of 1st Quintile    = ", EINQuintile(1), EINQuintile(1)/sum(HR(1:asst20,:)*mu(1:asst20,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Inflows to Employment of 2nd Quintile    = ", EINQuintile(2), EINQuintile(2)/sum(HR(asst20+1:asst40,:)*mu(asst20+1:asst40,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Inflows to Employment of 3rd Quintile    = ", EINQuintile(3), EINQuintile(3)/sum(HR(asst40+1:asst60,:)*mu(asst40+1:asst60,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Inflows to Employment of 4th Quintile    = ", EINQuintile(4), EINQuintile(4)/sum(HR(asst60+1:asst80,:)*mu(asst60+1:asst80,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Inflows to Employment of 5th Quintile    = ", EINQuintile(5), EINQuintile(5)/sum(HR(asst80+1:na,:)*mu(asst80+1:na,:)/hbar)
	write(1, *) " "
	write(1, '(A,F15.10," (",F8.5,")")') "Outflows from Employment of 1st Quintile = ", EOUTQuintile(1), EOUTQuintile(1)/sum(HR(1:asst20,:)*mu(1:asst20,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Outflows from Employment of 2nd Quintile = ", EOUTQuintile(2), EOUTQuintile(2)/sum(HR(asst20+1:asst40,:)*mu(asst20+1:asst40,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Outflows from Employment of 3rd Quintile = ", EOUTQuintile(3), EOUTQuintile(3)/sum(HR(asst40+1:asst60,:)*mu(asst40+1:asst60,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Outflows from Employment of 4th Quintile = ", EOUTQuintile(4), EOUTQuintile(4)/sum(HR(asst60+1:asst80,:)*mu(asst60+1:asst80,:)/hbar)
	write(1, '(A,F15.10," (",F8.5,")")') "Outflows from Employment of 5th Quintile = ", EOUTQuintile(5), EOUTQuintile(5)/sum(HR(asst80+1:na,:)*mu(asst80+1:na,:)/hbar)
	write(1, '(/)') 
	write(1, '(A/)') "Wealth Mobility"
	write(1, '(5F10.5)') ((WealthMobility(i,k),k=1,5),i=1,5)


	!	set the finishing time

	call date(finishdate)
	finishtime = clock()


	!	print the starting and finishing times

	write(1, '(/)') 
	write(1, '(A,A,A)') "The program started at  ", startdate, starttime
	write(1, '(A,A,A)') "The program finished at ", finishdate, finishtime

end subroutine
