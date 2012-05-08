
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: DataStatistic.m
% Compute Labor market wedge and statistics for U.S. economy
% Last Revised: 12-25 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

%   load data

% US time series
% Data_Selection = 1 old version
% Data_Selection = 2 new version

data_selection = 2;

switch data_selection
    case 1

load ..\usdata\preference_wage.txt;  
% This data set contains wage, hours, services consumption, non-durable consumption, population over 16 for 1964:1-2003:2
load ..\usdata\preference.dat;  
% This data set contains output, consumption, hours, labor productivity, TFP, investment for 1948:1-2003:2

data_wage=preference_wage;

[m n]= size(data_wage);

W = data_wage(:,1);
H = data_wage(:,2);
Cn = data_wage(:,3);
Cs = data_wage(:,4);
P = data_wage(:,5);
C = Cn+Cs;
W = log(W);

data=preference;

data=data(end-m+1:end,:);

Y = data(:,1); % output
C = data(:,2);
H = data(:,3);
A = data(:,4); % average labor productivity
S = data(:,5);
I = data(:,6);

nobs = size(H,1);
time = 1964+zeros(nobs,1);
for i = 2:nobs;
    time(i) = time(i-1)+.25;
end;

clear data data_wage

case 2
[data, name] = xlsread ('d:\Research\Preferences\usdata\us_labor.xls'); 
% dataset has lbip lbmn lbout lbcp7 gcq gifq gdpq p16 for 1948:1-2002:2
%   y    h   y/h   w     c    i    y   pop
%   1    2    3    4     5    6    7    8    
data = log(data);


skip = 10 % skip the first 10 years

data = data(4*skip:end,:);

[m n]= size(data);
Y = data(:,1); % output
H = data(:,2);
A = data(:,3); % average labor productivity
W = data(:,4);
C = data(:,5);
I = data(:,6);

nobs = size(H,1);
time = 1948+skip+zeros(nobs,1);
for i = 2:nobs;
    time(i) = time(i-1)+.25;
end;

end


Y = Y-mean(Y);
I = I-mean(I);
C = C-mean(C);
A = A-mean(A);
H = H-mean(H);
W = W-mean(W);

% labor supply elasticity
gamma=1.5;

B = -(A-C-(1/gamma)*H);
B2 = -(W-C-(1/gamma)*H);
MRS = (1/gamma)*H + C;




% Put things together
Nvar=9;
tsdata = zeros(nobs,Nvar);
tsdata = [Y C I A W MRS B B2 H]; 

for i = 1:Nvar
	tsdata(:,i) = hpfilter(tsdata(:,i), 1600);
end

%	statistics

stdev = zeros(Nvar,1);
rstdev = zeros(Nvar,1);
corr = zeros(Nvar,Nvar);
lag = 4;
autocorr = zeros(Nvar,lag);
crosscorr = zeros(Nvar*(Nvar-1)/2,2*lag+1);

for i = 1:Nvar
   stdev(i) = std(tsdata(:,i));
   rstdev(i) = stdev(i)/stdev(1);
end

corr = corrcoef(tsdata);

for i=1:Nvar
   for j=1:lag
      corrtmp = corrcoef(tsdata(1:end-j,i), tsdata(j+1:end,i));
      autocorr(i,j) = corrtmp(2,1);
   end
end


%	cross correlation

m = 0;
for i=1:Nvar-1
   for j=i+1:Nvar
      m = m + 1;
      for k=-lag:lag
         if k >= 0 
            corrtmp = corrcoef(tsdata(1:end-k,i), tsdata(k+1:end,j));
         else
            corrtmp = corrcoef(tsdata(1-k:end,i), tsdata(1:end+k,j));
         end
         crosscorr(m,lag+k+1) = corrtmp(2,1);
      end
   end
end

outS = fopen('DataStats.txt','W');

fprintf(outS, '%7s', 'Y', 'C', 'I', 'Y/H', 'W', 'MRS', 'B', 'B2', 'H');
fprintf(outS, '\n\n\n');

fprintf(outS, 'Standard deviations (percentage) \n\n');
fprintf(outS, '%7.2f', 100*stdev');
fprintf(outS, '\n\n\n');

fprintf(outS, 'Standard deviations relative to Output \n\n');
fprintf(outS, '%7.2f', rstdev');
fprintf(outS, '\n\n\n');

fprintf(outS, 'Correlation Coefficients \n\n');
for i=1:Nvar
	fprintf(outS, '%7.2f', corr(i,:));
   fprintf(outS, '\n');
end
fprintf(outS, '\n\n');

fprintf(outS, 'Auto-Correlation Coefficients \n\n');
for i=1:Nvar
   fprintf(outS, '%7.2f', autocorr(i,:));
   fprintf(outS, '\n');
end
fprintf(outS, '\n\n');

fprintf(outS, 'Cross-Correlation Coefficients \n\n');
for m=1:size(crosscorr,1)
   fprintf(outS, '%7.2f', crosscorr(m,:));
   fprintf(outS, '\n');
   end
fprintf(outS, '\n\n');

fclose(outS)


hY = tsdata(:,1);
hC = tsdata(:,2);
hI = tsdata(:,3);
hA = tsdata(:,4);
hW = tsdata(:,5);
hMRS = tsdata(:,6);
hB = tsdata(:,7);
hB2 = tsdata(:,8);
hH = tsdata(:,9);


figure(1)
plot(time, hH, '-', time, hA, 'r--','linewidth',2)
%axis(axis_default);
legend(' ln H ', ' ln Y/H')
print -dpsc DataH_A

figure(2)
plot(time, hH, '-',time, hW,'r--','linewidth',2)
legend(' ln H ','ln W')
print -dpsc DataH_W

figure(3)
plot(time, hH, '-', time, hB,'r--','linewidth',2);
legend(' ln H ', ' ln Wedge')
print -dpsc DataH_B
print -depsc Figure2

%title('Hours and Wages');

figure(4)
plot(time, hH, '-', time, hB2,'--','linewidth',2);
legend(' ln H ', 'ln (Wedge2)')
print -dpsc DataH_B2

figure(5)
plot(time, hMRS, '-', time, hA,'r--','linewidth',2);
legend( 'ln MRS','ln Y/H')
print -dpsc DataMRS_A
print -depsc Figure1

figure(6)
plot(time, hMRS, '-', time, hW,'r--','linewidth',2);
legend( 'ln MRS','ln W')
print -dpsc DataMRS_W



%figure(7)
%plot(time, Solow_trend, '-', time, S, 'm--')
%axis(axis_default);
%legend('Solow', 'Trend')

%figure(8)
%plot(time, Solow_detrend, '-', time, -hB, 'r--', time, recession, 'c:')
%legend('Solow', '-B')

%figure(9)
%plot(time, hS, '-', time, -hB, 'r--')
%axis(axis_default);
%legend('ln TFP', ' -ln B');
%title('TFP and Preference Shifts: Data (HP filtered)');

%output1 = fopen('SolowHP.dat','w');
%fprintf(output1,'%8.4f\n', hS);
%fclose(output1);

%output2 = fopen('Solow_detrend.dat','w');
%fprintf(output2,'%8.4f\n' , Solow_detrend);
%fclose(output2);

%figure(11)
%plot(time, H-mean(H), '-', time, WtoC-mean(WtoC),'--');
%legend('ln H', 'ln W/C')
%print DataHWC_raw
%title('Data');

%figure(12)
%plot(time, H-mean(H),'-', time, -(B-mean(B)),'--');
%legend('ln H', ' - ln B')
%print DataHB_raw
%title('Data');

%figure(14)
%plot(time, WtoC-mean(WtoC), '-', time, B-mean(B),'--');
%legend('ln W/C', 'ln B')
%print DataWCB_raw
%title('Data');

%figure(31)
%plot(time, hH, '-', time, rB,'--');
%legend('ln H', 'estimated B')

%corWtoC_H = corrcoef(WtoC, H)
%corB_H = corrcoef(B, H)
%corWtoC_B = corrcoef(WtoC, B)
