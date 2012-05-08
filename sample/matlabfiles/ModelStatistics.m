clear;
clc;

%%%%%% This File generates statistics for Table 2, 3, 4

%   choose the model
%   1: Fluctuation 
%   2: Complete Market
%   3: FluctuationDL 
%   4: Fluctuation2 Second moments 


model = 4;

switch model
    case 1
        cd ..\FortranRevision\Fluctuation\Output;
        Nskip = 500;
    case 2
        cd ..\FortranRevision\CompleteMarketPEA\Output;
        Nskip = 500;
    case 3
        cd ..\FortranRevision\FluctuationDL\Output;
        Nskip = 500;
    case 4
        cd ..\FortranRevision\Fluctuation2\Output;
        Nskip = 500;
end


%   load data

outS = fopen('ModelStats.txt','W');
outE = fopen('LS Elasticity.txt','W');


load Kdata.txt
load Rdata.txt
load Wdata.txt
load Ldata.txt
load Edata.txt
load Ydata.txt
load Idata.txt
load Cdata.txt
load Adata.txt
load Pdata.txt

cd ..\..\..\Matlab

gama=1.5;
MRS = (Edata.^(1/gama)).*Cdata;
MPL = Ydata./Edata;
B  = MRS./MPL;


%   number of variables and observations

Nperiod = size(Ydata, 1);
Nvar = 11;
nobs = Nperiod - Nskip;


%   put things together

tsdata = zeros(nobs,Nvar);
tsdata(:,1)  = Ydata(Nskip+1:Nperiod);
tsdata(:,2)  = Cdata(Nskip+1:Nperiod);
tsdata(:,3)  = Idata(Nskip+1:Nperiod);
tsdata(:,4)  = Adata(Nskip+1:Nperiod);
tsdata(:,5)  = Wdata(Nskip+1:Nperiod);
tsdata(:,6)  = Rdata(Nskip+1:Nperiod) - 0.025;
tsdata(:,7)  = Kdata(Nskip+1:Nperiod);
tsdata(:,8)  = Ldata(Nskip+1:Nperiod);
tsdata(:,9)  = Edata(Nskip+1:Nperiod);
tsdata(:,10) = MRS(Nskip+1:Nperiod);
tsdata(:,11) = B(Nskip+1:Nperiod);

ave = mean(tsdata);

%	take log

tsdata = log(tsdata);


%	HP filter

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


fprintf(outS, '%7s', 'Y', 'C', 'I','Y/H', 'W', 'R','K','L','E','MRS','Wedge');
fprintf(outS, '\n\n\n');

fprintf(outS, 'Averages \n\n');
fprintf(outS, '%7.2f', ave');
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




%   =================================================================   %
%               Regression labor supply elasticity                      %
%   =================================================================   %


%   drop first 500 observations and take log

Kp = log(Kdata(Nskip+2:end));
K  = log(Kdata(Nskip+1:end-1));
W  = log(Wdata(Nskip+1:end));
R  = log(Rdata(Nskip+1:end));
P  = log(Pdata(Nskip+1:end));
L  = log(Ldata(Nskip+1:end));
E  = log(Edata(Nskip+1:end));
A  = log(Adata(Nskip+1:end));
C  = log(Cdata(Nskip+1:end));


%   regressions

fprintf(outE, 'Regressions of (Unfiltered) MPL, APL, Employment');
fprintf(outE, '\n\n\n');


[coef_, coefint_, resid_, residint_, stats_] = regress(E, [ones(size(E,1),1) A], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;

fprintf(outE, 'Regression of Employment on Average Productivity\n');
fprintf(outE, 'coefficients: \n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, '\n');
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');

[coef_, coefint_, resid_, residint_, stats_] = regress(E, [ones(size(E,1),1) W], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;
fprintf(outE, 'Regression of Employment on Wage\n');
fprintf(outE, 'coefficients:\n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, '\n');
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');

[coef_, coefint_, resid_, residint_, stats_] = regress(E, [ones(size(E,1),1) A C], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;
fprintf(outE, 'Regression of Employment on Average Productivity and Consumption\n');
fprintf(outE, 'coefficients:\n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');

[coef_, coefint_, resid_, residint_, stats_] = regress(E, [ones(size(E,1),1) W C], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;
fprintf(outE, 'Regression of Employment on Wage and Consumption\n');
fprintf(outE, 'coefficients:\n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, '\n');
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');

[coef_, coefint_, resid_, residint_, stats_] = regress(E, [ones(size(E,1),1) A-C], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;
fprintf(outE, 'Regression of Employment on A/C\n');
fprintf(outE, 'coefficients:\n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, '\n');
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');

[coef_, coefint_, resid_, residint_, stats_] = regress(E, [ones(size(E,1),1) W-C], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;
fprintf(outE, 'Regression of Employment on W/C\n');
fprintf(outE, 'coefficients:\n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, '\n');
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');

[coef_, coefint_, resid_, residint_, stats_] = regress(L, [ones(size(L,1),1) A-C], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;
fprintf(outE, 'Regression of Efficiency on A/C\n');
fprintf(outE, 'coefficients:\n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, '\n');
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');

[coef_, coefint_, resid_, residint_, stats_] = regress(L, [ones(size(L,1),1) W-C], 0.05);
sd_ = (coef_-coefint_(:,1))/2;
tratio_ = coef_./sd_;
fprintf(outE, 'Regression of Efficiency on W/C\n');
fprintf(outE, 'coefficients:\n');
fprintf(outE, '%12.5f', coef_);
fprintf(outE, '\n');
fprintf(outE, 't ratio: \n');
fprintf(outE, '%12.5f', tratio_);
fprintf(outE, '\n');
fprintf(outE, 'R-squared: \n');
fprintf(outE, '%12.5f', stats_(1));
fprintf(outE, '\n');
fprintf(outE, '===================\n');
fprintf(outE, '\n\n\n');

fclose(outE);


