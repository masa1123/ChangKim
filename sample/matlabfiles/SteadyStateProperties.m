%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: SteadyStateProperties
% calculate properties of distributions and plot them
%   this code generates Table 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;


%	global variables

alpha = 0.36;
rd = 0.01+0.025;
hbar = 1/3;

kss = (alpha/rd)^(1/(1-alpha));
wage = (1-alpha)*(kss^alpha);


% Model 1: sig
cd ..\FortranRevision\SteadyState\Output;

load AsstEarnDensity.txt;
load AsstDensity.txt;
load EarnDensity.txt;
load AsstEarnDist.txt;
load AsstDist.txt;
load EarnDist.txt;
load agrid.txt;
load exgrid.txt;

Earnings = [0; hbar*wage*exgrid];

%	asset levels of interesting

asst20 = 93; 
asst40 = 273;
asst60 = 404;
asst80 = 506;
asst90 = 557;
asst95 = 607;
asst98 = 647;
asst99 = 678;


%	ratio of mean asset of each subgroup to sample mean asset

MeanAsst = agrid'*AsstDensity;
MeanAsst0020 = agrid(1:asst20)'*AsstDensity(1:asst20)/sum(AsstDensity(1:asst20));
MeanAsst2040 = agrid(asst20+1:asst40)'*AsstDensity(asst20+1:asst40)/sum(AsstDensity(asst20+1:asst40));
MeanAsst4060 = agrid(asst40+1:asst60)'*AsstDensity(asst40+1:asst60)/sum(AsstDensity(asst40+1:asst60));
MeanAsst6080 = agrid(asst60+1:asst80)'*AsstDensity(asst60+1:asst80)/sum(AsstDensity(asst60+1:asst80));
MeanAsst80100 = agrid(asst80+1:end)'*AsstDensity(asst80+1:end)/sum(AsstDensity(asst80+1:end));
MeanAsst9095 = agrid(asst90+1:asst95)'*AsstDensity(asst90+1:asst95)/sum(AsstDensity(asst90+1:asst95));
MeanAsst9599 = agrid(asst95+1:asst99)'*AsstDensity(asst95+1:asst99)/sum(AsstDensity(asst95+1:asst99));
MeanAsst99100 = agrid(asst99+1:end)'*AsstDensity(asst99+1:end)/sum(AsstDensity(asst99+1:end));

RatioAsst0020 = MeanAsst0020/MeanAsst;
RatioAsst2040 = MeanAsst2040/MeanAsst;
RatioAsst4060 = MeanAsst4060/MeanAsst;
RatioAsst6080 = MeanAsst6080/MeanAsst;
RatioAsst80100 = MeanAsst80100/MeanAsst;
RatioAsst9095 = MeanAsst9095/MeanAsst;
RatioAsst9599 = MeanAsst9599/MeanAsst;
RatioAsst99100 = MeanAsst99100/MeanAsst;


%	share of asset holdings of each subgroup out of total asset

ShareAsst0020 = 100*agrid(1:asst20)'*AsstDensity(1:asst20)/MeanAsst;
ShareAsst2040 = 100*agrid(asst20+1:asst40)'*AsstDensity(asst20+1:asst40)/MeanAsst;
ShareAsst4060 = 100*agrid(asst40+1:asst60)'*AsstDensity(asst40+1:asst60)/MeanAsst;
ShareAsst6080 = 100*agrid(asst60+1:asst80)'*AsstDensity(asst60+1:asst80)/MeanAsst;
ShareAsst80100 = 100*agrid(asst80+1:end)'*AsstDensity(asst80+1:end)/MeanAsst;
ShareAsst8090 = 100*agrid(asst80+1:asst90)'*AsstDensity(asst80+1:asst90)/MeanAsst;
ShareAsst9095 = 100*agrid(asst90+1:asst95)'*AsstDensity(asst90+1:asst95)/MeanAsst;
ShareAsst9599 = 100*agrid(asst95+1:asst99)'*AsstDensity(asst95+1:asst99)/MeanAsst;
ShareAsst99100 = 100*agrid(asst99+1:end)'*AsstDensity(asst99+1:end)/MeanAsst;


%	share of earnings of each subgroup out of total earnings

TotalEarnings = dot(Earnings,sum(AsstEarnDensity));

ShareEarn0020 = 100*dot(Earnings,sum(AsstEarnDensity(1:asst20,:)))/TotalEarnings;
ShareEarn2040 = 100*dot(Earnings,sum(AsstEarnDensity(asst20+1:asst40,:)))/TotalEarnings;
ShareEarn4060 = 100*dot(Earnings,sum(AsstEarnDensity(asst40+1:asst60,:)))/TotalEarnings;
ShareEarn6080 = 100*dot(Earnings,sum(AsstEarnDensity(asst60+1:asst80,:)))/TotalEarnings;
ShareEarn80100 = 100*dot(Earnings,sum(AsstEarnDensity(asst80+1:end,:)))/TotalEarnings;
ShareEarn9095 = 100*dot(Earnings,sum(AsstEarnDensity(asst90+1:asst95,:)))/TotalEarnings;
ShareEarn9599 = 100*dot(Earnings,sum(AsstEarnDensity(asst95+1:asst99,:)))/TotalEarnings;
ShareEarn99100 = 100*dot(Earnings,sum(AsstEarnDensity(asst99+1:end,:)))/TotalEarnings;


%	Lorenz curve

AsstRange = [0 20 40 60 80 90 95 99 100];

LorenzAsst(1) = 0;
LorenzAsst(2) = ShareAsst0020;
LorenzAsst(3) = LorenzAsst(2) + ShareAsst2040;
LorenzAsst(4) = LorenzAsst(3) + ShareAsst4060;
LorenzAsst(5) = LorenzAsst(4) + ShareAsst6080;
LorenzAsst(6) = LorenzAsst(5) + ShareAsst8090;
LorenzAsst(7) = LorenzAsst(6) + ShareAsst9095;
LorenzAsst(8) = LorenzAsst(7) + ShareAsst9599;
LorenzAsst(9) = LorenzAsst(8) + ShareAsst99100;

EarnRange = [0; 100*EarnDist];
MeanEarn = Earnings'*EarnDensity;

NonZeroEarnings = Earnings(2:end);
NonZeroEarnDensity = EarnDensity(2:end)/(1-EarnDensity(1));
NonZeroEarnDist = (EarnDist(2:end)-EarnDist(1))/(1-EarnDensity(1));

NonZeroEarnRange = [0; 100*NonZeroEarnDist];
NonZeroMeanEarn = NonZeroEarnings'*NonZeroEarnDensity;

LorenzEarn(1) = 0;
for i=1:size(Earnings,1)
    LorenzEarn(i+1) = 100*Earnings(1:i)'*EarnDensity(1:i)/MeanEarn;
end
LorenzNonZeroEarn(1) = 0;
for i=1:size(NonZeroEarnings,1)
    LorenzNonZeroEarn(i+1) = 100*NonZeroEarnings(1:i)'*NonZeroEarnDensity(1:i)/NonZeroMeanEarn;
end


%   Gini coefficient

GiniAsst = 0;
for i = 2:size(AsstRange,2)
   GiniAsst = GiniAsst + 0.5*(AsstRange(i) - AsstRange(i-1))*(LorenzAsst(i) + LorenzAsst(i-1));
end
GiniAsst = (5000 - GiniAsst)/5000;

GiniEarn = 0;
for i = 2:size(EarnRange,1)
   GiniEarn = GiniEarn + 0.5*(EarnRange(i) - EarnRange(i-1))*(LorenzEarn(i) + LorenzEarn(i-1));
end
GiniNonZeroEarn = 0;
for i = 2:size(NonZeroEarnRange,1)
   GiniNonZeroEarn = GiniNonZeroEarn + 0.5*(NonZeroEarnRange(i) - NonZeroEarnRange(i-1))*(LorenzNonZeroEarn(i) + LorenzNonZeroEarn(i-1));
end
GiniEarn = (5000 - GiniEarn)/5000;
GiniNonZeroEarn = (5000 - GiniNonZeroEarn)/5000;


%	output to a log file

out = fopen('Properties of Distribution.txt','W');
fprintf(out, 'Properties of Asset-Earning Distribution');
fprintf(out, '\n\n');
fprintf(out, 'sigx = 0.227, rhox = 0.929');
fprintf(out, '\n\n\n');

fprintf(out, '                    Households in Wealth Quintile          The Wealth Rich\n');
fprintf(out, '                 -----------------------------------    --------------------\n');
fprintf(out, '            ');
fprintf(out, '%8s', '1st', '2nd', '3rd', '4th', '5th', '10-5', '5-1', '1');
fprintf(out, '\n\n');

fprintf(out, 'ShareAsst   ');
fprintf(out, '%8.2f', ShareAsst0020, ShareAsst2040, ShareAsst4060, ShareAsst6080, ShareAsst80100, ShareAsst9095, ShareAsst9599, ShareAsst99100);
fprintf(out, '\n');

fprintf(out, 'RatioAsst   ');
fprintf(out, '%8.2f', RatioAsst0020, RatioAsst2040, RatioAsst4060, RatioAsst6080, RatioAsst80100, RatioAsst9095, RatioAsst9599, RatioAsst99100);
fprintf(out, '\n');

fprintf(out, 'ShareEarn   ');
fprintf(out, '%8.2f', ShareEarn0020, ShareEarn2040, ShareEarn4060, ShareEarn6080, ShareEarn80100, ShareEarn9095, ShareEarn9599, ShareEarn99100);
fprintf(out, '\n\n\n');

Top01Bottom40 = ShareAsst99100/(ShareAsst0020+ShareAsst2040);
Top01Bottom20 = ShareAsst99100/ShareAsst0020;
Top05Bottom20 = (ShareAsst9599+ShareAsst99100)/ShareAsst0020;
Top10Bottom20 = (ShareAsst9095+ShareAsst9599+ShareAsst99100)/ShareAsst0020;

fprintf(out, '\n');
fprintf(out, 'Gini Coefficient for Asset = %7.4f', GiniAsst);
fprintf(out, '\n');
fprintf(out, 'Gini Coefficient for Earnings = %7.4f', GiniEarn);
fprintf(out, '\n');
fprintf(out, 'Gini Coefficient for NonZeroEarnings = %7.4f', GiniNonZeroEarn);
fprintf(out, '\n\n');
fprintf(out, 'Ratio of Top 1 percent to Bottom 40 percent = %7.4f', Top01Bottom40);
fprintf(out, '\n');
fprintf(out, 'Ratio of Top 1 percent to Bottom 20 percent = %7.4f', Top01Bottom20);
fprintf(out, '\n');
fprintf(out, 'Ratio of Top 5 percent to Bottom 20 percent = %7.4f', Top05Bottom20);
fprintf(out, '\n');
fprintf(out, 'Ratio of Top 10 percent to Bottom 20 percent = %7.4f', Top10Bottom20);
fprintf(out, '\n\n\n\n');
fclose(out);


% For US, from Diaz-Ginmenez et al.

%Range = [0 1 5 10 20 40 60 80 90 95 99 100];

%LorenzAsstD(1) = 0;
%LorenzAsstD(2) = 0;
%LorenzAsstD(3) = 0;
%LorenzAsstD(4) = 0;
%LorenzAsstD(5) = 0;
%LorenzAsstD(6) = LorenzAsstD(5) + 1.74;
%LorenzAsstD(7) = LorenzAsstD(6) + 5.72;
%LorenzAsstD(8) = LorenzAsstD(7) + 13.43;
%LorenzAsstD(9) = LorenzAsstD(8) + 79.49-12.62-23.95-29.55;
%LorenzAsstD(10) = LorenzAsstD(9) + 12.62;
%LorenzAsstD(11) = LorenzAsstD(10) + 23.95;
%LorenzAsstD(12) = LorenzAsstD(11) + 29.55;

%LorenzEarnD(1) = 0;
%LorenzEarnD(2) = 0;
%LorenzEarnD(3) = 0;
%LorenzEarnD(4) = 0;
%LorenzEarnD(5) = 0;
%LorenzEarnD(6) = LorenzEarnD(5) + 3.19;
%LorenzEarnD(7) = LorenzEarnD(6) + 12.49;
%LorenzEarnD(8) = LorenzEarnD(7) + 23.33;
%LorenzEarnD(9) = LorenzEarnD(8) + 61.39-12.38-16.37-14.76;
%LorenzEarnD(10) = LorenzEarnD(9) + 12.38;
%LorenzEarnD(11) = LorenzEarnD(10) + 16.37;
%LorenzEarnD(12) = LorenzEarnD(11) + 14.76;


% For US, from PSID 1984

Range = [0 1 5 10 20 40 60 80 90 95 99 100];

LorenzAsstD(1) = 0;
LorenzAsstD(2) = 0;
LorenzAsstD(3) = 0;
LorenzAsstD(4) = 0;
LorenzAsstD(5) = 0;
LorenzAsstD(6) = LorenzAsstD(5) + 0.5;
LorenzAsstD(7) = LorenzAsstD(6) + 5.06;
LorenzAsstD(8) = LorenzAsstD(7) + 18.74;
LorenzAsstD(9) = LorenzAsstD(8) + 18.04;
LorenzAsstD(10) = LorenzAsstD(9) + 15.06;
LorenzAsstD(11) = LorenzAsstD(10) + 20.49;
LorenzAsstD(12) = LorenzAsstD(11) + 22.62;

LorenzEarnD(1) = 0;
LorenzEarnD(2) = 0;
LorenzEarnD(3) = 0;
LorenzEarnD(4) = 0;
LorenzEarnD(5) = 0;
LorenzEarnD(6) = LorenzEarnD(5) + 5.25;
LorenzEarnD(7) = LorenzEarnD(6) + 14.94;
LorenzEarnD(8) = LorenzEarnD(7) + 27.42;
LorenzEarnD(9) = LorenzEarnD(8) + 20.84;
LorenzEarnD(10) = LorenzEarnD(9) + 12.34;
LorenzEarnD(11) = LorenzEarnD(10) + 13.36;
LorenzEarnD(12) = LorenzEarnD(11) + 5.83;


% plot Lorenz curves

figure(1)
plot(0:10:100, 0:10:100, ':', AsstRange, LorenzAsst, '-', EarnRange, LorenzEarn, '--', 'LineWidth', 2);
axis([0 100 0 100]);
legend('45^o line', 'Asset Holdings', 'Earnings', 0)';
title('Lorenz Curves of the Model', 'FontSize',14);


figure(2)
plot(0:10:100, 0:10:100, ':', Range, LorenzAsstD, '-', Range, LorenzEarnD, '--', 'LineWidth', 2);
axis([0 100 0 100]);
legend('45^o line', 'Wealth', 'Earnings', 0)';
title('Lorenz Curves of U.S. Economy', 'FontSize',14);

figure(3)
plot(0:10:100, 0:10:100, ':', AsstRange, LorenzAsst, '-', Range, LorenzAsstD, '--', 'LineWidth', 2);
axis([0 100 0 100]);
legend('45^o line', 'Model I', 'U.S.', 0)';
print -dpsc LorenzAsst_M1
title('Lorenz Curve for Wealth', 'FontSize',14);

figure(4)
plot(0:10:100, 0:10:100, ':', EarnRange, LorenzEarn, '-', Range, LorenzEarnD, '--', 'LineWidth', 2);
axis([0 100 0 100]);
legend('45^o line', 'Model I', 'U.S.', 0)';
print -dpsc LorenzEarn_M1
title('Lorenz Curve for Earnings', 'FontSize',14);


EarnRange_M1 = EarnRange;
LorenzAsst_M1 = LorenzAsst;
LorenzEarn_M1 = LorenzEarn;

cd ..\..\..\Matlab;

