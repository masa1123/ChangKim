!	----------------------------------------------------------------------
!	File name: SimulateData.f90
!
!	generate artificial data for aggregate capital, interest rate, and 
!	wage rate through simulating individuals' decision rules. These data
!	will be used in RegressLOM().
!	----------------------------------------------------------------------


subroutine SimulatePanel()

	use Globals
	use LinInterpModule
	use Numerical_Libraries
	
	implicit none

	integer, parameter:: xseed = 2
	integer, parameter:: Nindiv = 10000, Ntime = 120, Nskip = 500
	real(8) PanelA(Nindiv,Ntime), PanelH(Nindiv,Ntime), PanelC(Nindiv,Ntime),	&
			PanelW(Nindiv,Ntime)
	real(8) YrPanelA(Nindiv,Ntime/4), YrPanelH(Nindiv,Ntime/4),		&
			YrPanelC(Nindiv,Ntime/4),	YrPanelW(Nindiv,Ntime/4)
	real(8) TsK(Ntime), TsE(Ntime), TsL(Ntime), TsW(Ntime), TsC(Ntime)
	real(8) CrA(Nindiv,2), CrH(Nindiv,2)
	integer CrX(Nindiv,Ntime)
	real(8) vf_work, vf_home, Kaggr, Eaggr
	integer indiv, time

	
	!	set the seed number for random number generator

	call rnset(xseed)


	!	intialize distribution of workers

	do indiv = 1, Nindiv

		CrA(indiv,1) = kss
		CrX(indiv,1) = mod(indiv, nx) + 1
		CrH(indiv,1) = hbar

	end do


	!	start generating artificial panel data

	do time = -Nskip, Ntime


		!	cross-section data for asset and hours

		do indiv = 1, Nindiv

			vf_work = lininterp1(CrA(indiv,1), agrid, VFE(:,CrX(indiv,1)))
			vf_home = lininterp1(CrA(indiv,1), agrid, VFN(:,CrX(indiv,1)))

			if (vf_work >= vf_home) then

				CrH(indiv,2) = hbar
				CrA(indiv,2) = lininterp1(CrA(indiv,1), agrid, agrid(ASE(:, CrX(indiv,1))))
				CrA(indiv,2) = min(max(CrA(indiv,2), agrid(1)), agrid(na))

			else

				CrH(indiv,2) = 0.0D0
				CrA(indiv,2) = lininterp1(CrA(indiv,1), agrid, agrid(ASN(:, CrX(indiv,1))))
				CrA(indiv,2) = min(max(CrA(indiv,2), agrid(1)), agrid(na))

			end if

		end do


		!	record panel data (quarterly)

		if (time >= 1) then

			PanelA(:,time) = CrA(:,1)
			PanelH(:,time) = CrH(:,2)
			PanelW(:,time) = wage*CrH(:,2)*exgrid(CrX(:,1))/hbar
			PanelC(:,time) = wage*CrH(:,2)*exgrid(CrX(:,1)) + (1.0D0+irate)*CrA(:,1) - CrA(:,2)


			!	construct yearly panel data

			if (mod(time, 4) == 0.0) then

				YrPanelA(:,time/4) = sum(PanelA(:,time-3:time),dim=2)/4
				YrPanelH(:,time/4) = sum(PanelH(:,time-3:time),dim=2)
				do indiv = 1, Nindiv
					if (YrPanelH(indiv,time/4) > 0) then
						YrPanelW(indiv,time/4) &
						= sum(PanelW(indiv,time-3:time))*hbar/YrPanelH(indiv,time/4)
					else
						YrPanelW(indiv,time/4) = 0.0D0
					end if
				end do
				YrPanelC(:,time/4) = sum(PanelC(:,time-3:time),dim=2)/4

			end if


			!	record time-series data (quarterly)

			TsK(time) = sum(PanelA(:,time))/Nindiv
			TsE(time) = sum(PanelH(:,time)*3.0D0)/Nindiv
			TsL(time) = sum(CrH(:,2)*exgrid(CrX(:,1)))/Nindiv
			TsW(time) = (1.0D0-alpha)*(TsK(time)/TsL(time))**(alpha)
			TsC(time) = sum(PanelC(:,time))/Nindiv


		end if


		!	show aggregate variables 

		if (mod(time,10) == 0) then

			Kaggr = sum(CrA(:,2))/Nindiv
			Eaggr = sum(CrH(:,2)*3.0D0)/Nindiv
			write(*, '(A,I4,A,F8.4,A,F8.4)')	&
			"time = ", time, "   Kaggr = ", Kaggr, "   Eaggr = ", Eaggr

		end if


		!	generate next period idiosyncratic shocks

		call NextIdiosyncraticShock(CrX, Nindiv)


		!	update cross section data one period

		CrA(:,1) = CrA(:,2)
		CrH(:,1) = CrH(:,2)
		CrX(:,1) = CrX(:,2)


	end do


	!	save time sries data

!	open(1, file='Output\PanelA.txt', status='unknown')
!	open(2, file='Output\PanelH.txt', status='unknown')
!	open(3, file='Output\PanelW.txt', status='unknown')
!	open(4, file='Output\PanelC.txt', status='unknown')
!	open(5, file='Output\YrPanelA.txt', status='unknown')
!	open(6, file='Output\YrPanelH.txt', status='unknown')
!	open(7, file='Output\YrPanelW.txt', status='unknown')
!	open(8, file='Output\YrPanelC.txt', status='unknown')
!	open(9, file='Output\TsK.txt', status='unknown')
!	open(10, file='Output\TsE.txt', status='unknown')
!	open(11, file='Output\TsW.txt', status='unknown')
!	open(12, file='Output\TsC.txt', status='unknown')
!	open(13, file='Output\TsL.txt', status='unknown')

!	write(1, '(<Ntime>f12.6)') ((PanelA(indiv,time), time=1,Ntime),indiv=1,Nindiv)
!	write(2, '(<Ntime>f12.6)') ((PanelH(indiv,time), time=1,Ntime),indiv=1,Nindiv)
!	write(3, '(<Ntime>f12.6)') ((PanelW(indiv,time), time=1,Ntime),indiv=1,Nindiv)
!	write(4, '(<Ntime>f12.6)') ((PanelC(indiv,time), time=1,Ntime),indiv=1,Nindiv)
!	write(5, '(<Ntime/4>f12.6)') ((YrPanelA(indiv,time), time=1,Ntime/4),indiv=1,Nindiv)
!	write(6, '(<Ntime/4>f12.6)') ((YrPanelH(indiv,time), time=1,Ntime/4),indiv=1,Nindiv)
!	write(7, '(<Ntime/4>f12.6)') ((YrPanelW(indiv,time), time=1,Ntime/4),indiv=1,Nindiv)
!	write(8, '(<Ntime/4>f12.6)') ((YrPanelC(indiv,time), time=1,Ntime/4),indiv=1,Nindiv)
!	write(9, '(f12.6)') TsK
!	write(10, '(f12.6)') TsE
!	write(11, '(f12.6)') TsW
!	write(12, '(f12.6)') TsC
!	write(13, '(f12.6)') TsL


	!	data set for stata

	open(15, file='Output\sYrPanel.txt', status='unknown')
	write(15, '(4f12.6)') ((YrPanelA(indiv,time), YrPanelH(indiv,time), YrPanelW(indiv,time), YrPanelC(indiv,time), time=1,Ntime/4),indiv=1,Nindiv)
	open(16, file='Output\sQrTS.txt', status='unknown')
	write(16, '(4f12.6)') (TsK(time), TsE(time), TsW(time), TsC(time), time=1,Ntime)


end subroutine




subroutine NextIdiosyncraticShock(CrX, Nindiv)

	use Numerical_Libraries
	use Globals

	integer j, Nindiv, CrX(Nindiv,2)
	real(8) ushock


	do indiv = 1, Nindiv

		ushock = drnunf()

		do j = 1, nx

			if (ushock <= CtrX(CrX(indiv,1),j)) then
				CrX(indiv,2) = j
				exit
			end if
		
		end do

	end do


end subroutine
