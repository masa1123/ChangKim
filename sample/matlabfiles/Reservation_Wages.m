%	This creates Figure 3 Reservation Wage Schedule


clear all;
clc;


%   move to the directory th at contains results

cd ..\FortranRevision\SteadyState\Output;



%	global variables

alpha = 0.36;
rd = 0.01+0.025;
hbar = 1/3;

kss = (alpha/rd)^(1/(1-alpha));
wage = (1-alpha)*(kss^alpha);

load AsstEarnDensity.txt;
load AsstDensity.txt;
load EarnDensity.txt;
load AsstEarnDist.txt;
load AsstDist.txt;
load EarnDist.txt;
load ReservationX.txt;
load ReservationWages.txt;

load agrid.txt;
load xgrid.txt;
load exgrid.txt;

A1=agrid;
AE1=AsstEarnDensity;
Earn1=EarnDensity;
Asst1=AsstDensity;

earnings = [0; hbar*wage*exgrid];


%	asset levels of interest (earlier version) % sigx=0.287 rho=0.939
%asst20 = 85;
%asst40 = 274;
%asst60 = 414;
%asst80 = 526;


%	asset levels of interest (new version) % sigx=0.227 rho=0.929
asst20 = 93;
asst40 = 273;
asst60 = 404;
asst80 = 506;

%	ratio of mean asset of each subgroup to sample mean asset

MeanAsst = agrid'*AsstDensity;
%unit= 60524/MeanAsst  % all group
unit= 102744/MeanAsst  % high school graduate 35<=age<=55 


%	share of asset holdings of each subgroup out of total asset

ShareAsst0020 = 100*agrid(1:asst20)'*AsstDensity(1:asst20)/MeanAsst;
ShareAsst2040 = 100*agrid(asst20+1:asst40)'*AsstDensity(asst20+1:asst40)/MeanAsst;
ShareAsst4060 = 100*agrid(asst40+1:asst60)'*AsstDensity(asst40+1:asst60)/MeanAsst;
ShareAsst6080 = 100*agrid(asst60+1:asst80)'*AsstDensity(asst60+1:asst80)/MeanAsst;
ShareAsst80100 = 100*agrid(asst80+1:end)'*AsstDensity(asst80+1:end)/MeanAsst;


%	share of earnings of each subgroup out of total earnings

TotalEarnings = dot(earnings,sum(AsstEarnDensity));

ShareEarn0020 = 100*dot(earnings,sum(AsstEarnDensity(1:asst20,:)))/TotalEarnings;
ShareEarn2040 = 100*dot(earnings,sum(AsstEarnDensity(asst20+1:asst40,:)))/TotalEarnings;
ShareEarn4060 = 100*dot(earnings,sum(AsstEarnDensity(asst40+1:asst60,:)))/TotalEarnings;
ShareEarn6080 = 100*dot(earnings,sum(AsstEarnDensity(asst60+1:asst80,:)))/TotalEarnings;
ShareEarn80100 = 100*dot(earnings,sum(AsstEarnDensity(asst80+1:end,:)))/TotalEarnings;


%	Lorenz curve

AsstRange = [0 20 40 60 80 100];

LorenzAsst(1) = 0;
LorenzAsst(2) = ShareAsst0020;
LorenzAsst(3) = LorenzAsst(2) + ShareAsst2040;
LorenzAsst(4) = LorenzAsst(3) + ShareAsst4060;
LorenzAsst(5) = LorenzAsst(4) + ShareAsst6080;
LorenzAsst(6) = LorenzAsst(5) + ShareAsst80100;

EarnRange = [0; 100*EarnDist];
MeanEarn = earnings'*EarnDensity;

LorenzEarn(1) = 0;
LorenzEarn(2) = 100*earnings(1:1)'*EarnDensity(1:1)/MeanEarn;
LorenzEarn(3) = 100*earnings(1:2)'*EarnDensity(1:2)/MeanEarn;
LorenzEarn(4) = 100*earnings(1:3)'*EarnDensity(1:3)/MeanEarn;
LorenzEarn(5) = 100*earnings(1:4)'*EarnDensity(1:4)/MeanEarn;
LorenzEarn(6) = 100*earnings(1:5)'*EarnDensity(1:5)/MeanEarn;
LorenzEarn(7) = 100*earnings(1:6)'*EarnDensity(1:6)/MeanEarn;
LorenzEarn(8) = 100*earnings(1:7)'*EarnDensity(1:7)/MeanEarn;
LorenzEarn(9) = 100*earnings(1:8)'*EarnDensity(1:8)/MeanEarn;
LorenzEarn(10) = 100*earnings(1:9)'*EarnDensity(1:9)/MeanEarn;
LorenzEarn(11) = 100*earnings(1:10)'*EarnDensity(1:10)/MeanEarn;
LorenzEarn(12) = 100*earnings(1:11)'*EarnDensity(1:11)/MeanEarn;
LorenzEarn(13) = 100*earnings(1:12)'*EarnDensity(1:12)/MeanEarn;
LorenzEarn(14) = 100*earnings(1:13)'*EarnDensity(1:13)/MeanEarn;
LorenzEarn(15) = 100*earnings(1:14)'*EarnDensity(1:14)/MeanEarn;
LorenzEarn(16) = 100*earnings(1:15)'*EarnDensity(1:15)/MeanEarn;
LorenzEarn(17) = 100*earnings(1:16)'*EarnDensity(1:16)/MeanEarn;
LorenzEarn(18) = 100*earnings(1:17)'*EarnDensity(1:17)/MeanEarn;
LorenzEarn(19) = 100*earnings(1:18)'*EarnDensity(1:18)/MeanEarn;


%   Gini coefficient

disp('Gini-Coefficient for Wealth and Earnings: sigx=0.227 rhox=0.929')
GiniAsst = 0;
for i = 2:size(AsstRange,2)
   GiniAsst = GiniAsst + 0.5*(AsstRange(i) - AsstRange(i-1))*(LorenzAsst(i) + LorenzAsst(i-1));
end
GiniAsst = (5000 - GiniAsst)/5000

GiniEarn = 0;
for i = 2:size(EarnRange,1)
   GiniEarn = GiniEarn + 0.5*(EarnRange(i) - EarnRange(i-1))*(LorenzEarn(i) + LorenzEarn(i-1));
end
GiniEarn = (5000 - GiniEarn)/5000


% For US, from PSID 1984 high school and 35<=age<=55 

Range = [0 1 5 10 20 40 60 80 90 95 99 100];
%        1 2 3 4  5  6  7  8  9  10 11  12
LorenzAsstD(1) = 0;
LorenzAsstD(2) = 0;
LorenzAsstD(3) = 0;
LorenzAsstD(4) = 0;
LorenzAsstD(5) = 1.03;
LorenzAsstD(6) = LorenzAsstD(5) + 7.07;
LorenzAsstD(7) = LorenzAsstD(6) + 13.01;
LorenzAsstD(8) = LorenzAsstD(7) + 21.10;
LorenzAsstD(9) = LorenzAsstD(8) + 17.98;
LorenzAsstD(10) = LorenzAsstD(9) + 13.06;
LorenzAsstD(11) = LorenzAsstD(10) + 15.54;
LorenzAsstD(12) = LorenzAsstD(11) + 11.17;

LorenzEarnD(1) = 0;
LorenzEarnD(2) = 0;
LorenzEarnD(3) = 0;
LorenzEarnD(4) = 0;
LorenzEarnD(5) = 4.77;
LorenzEarnD(6) = LorenzEarnD(5) + 12.83;
LorenzEarnD(7) = LorenzEarnD(6) + 18.46;
LorenzEarnD(8) = LorenzEarnD(7) + 24.85;
LorenzEarnD(9) = LorenzEarnD(8) + 16.72;
LorenzEarnD(10) = LorenzEarnD(9) + 10.07;
LorenzEarnD(11) = LorenzEarnD(10) + 10.22;
LorenzEarnD(12) = LorenzEarnD(11) + 2.31;


% keep variables for an integrated graph

LorenzAsst_S287R939 = LorenzAsst;
LorenzEarn_S287R939 = LorenzEarn;


% plot figures

cd ..\..\..\Matlab;

%   Lorenz curves: wealth

figure(1)
plot(0:10:100, 0:10:100, ':', AsstRange, LorenzAsst_S287R939, '--', Range, LorenzAsstD, '-', 'LineWidth', 2);
axis([0 100 0 100]);
legend('45^o','Model', 'U.S.', 0)';
print Lorenzw -depsc
%title('Lorenz Curve for Wealth', 'FontSize',14);


%   Lorenz curves: earnings

figure(2)
plot(0:10:100, 0:10:100, ':', EarnRange, LorenzEarn_S287R939, '--', Range, LorenzEarnD, '-', 'LineWidth', 2);
axis([0 100 0 100]);
legend('45^o','Model', 'U.S.', 0)';
print Lorenze -depsc
%title('Lorenz Curve for Earnings', 'FontSize',14);


%   plot reservation wages for the benchmark: Figure 5 in the paper

MeanAsst = agrid'*AsstDensity;
unit= 102744/MeanAsst;
hbar=1/3;

MedianAsst = agrid(340);
nasub = 514; %agrid where unit*nasub=$200,000
vertical_line = 0:50:100000;
atmeanasset = unit*MeanAsst*ones(size(vertical_line,1),1);
atmedianasset = unit*MedianAsst*ones(size(vertical_line,1),1);


amax = 940; %maximal asset level plotted.
figure(3)
subplot(211)
area(unit*agrid(1:nasub), hbar*unit*ReservationWages(amax+20)*ones(nasub,1), 'faceColor','y')
hold on
plot(unit*agrid(1:amax), hbar*unit*ReservationWages(1:amax),'-', 'LineWidth', 2)
axis(unit*[agrid(1) agrid(amax) 0 hbar*ReservationWages(amax+20)])
ylabel('Reservation Wages')
title('A: All Asset Levels', 'FontSize', 12)
hold off
subplot(212)
plot(unit*agrid(1:nasub), hbar*unit*ReservationWages(1:nasub),'-', 'LineWidth', 2)
axis([unit*agrid(1) unit*agrid(nasub) 0 1.1e4])
xlabel('Assets in 1983 dollars')
ylabel('Reservation Wages')
title('B: Assets less than $200,000', 'FontSize', 12)
hold on
plot(atmeanasset, vertical_line, 'r-', 'LineWidth', 1)
plot(atmedianasset, vertical_line, 'r-', 'LineWidth', 1)
text(MeanAsst*unit-2000, 4000, 'Mean Wealth')
text(MedianAsst*unit-2000, 3000, 'Median Wealth')
hold off
print Figure3 -depsc
print Reservation -dpsc
subplot



