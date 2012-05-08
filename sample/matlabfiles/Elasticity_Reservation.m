clear;
clc;

% This file computes the labor supply elasticity based on the reservation
% wage distribution

cd ..\FortranRevision\SteadyState\Output;


%	open an output file

out = fopen('Elasticity.txt','W');

fprintf(out, 'Labor Supply Elasticity (Steady State)');
fprintf(out, '\n\n\n');


load ReservationWages.txt;
load LS.txt;

%   labor supply elasticity from participation-reservation wage plot

L_95_225 = log(LS);
W_95_225 = log(ReservationWages);
Elas_95_225 = gradient(L_95_225)./gradient(W_95_225);
L_95_225 = LS;
LSE_95_225 = [L_95_225 Elas_95_225];

[xx,x58] = min(abs(LS-58));
[xx,x60] = min(abs(LS-60));
[xx,x62] = min(abs(LS-62));


%   labor supply elasticity using regressed reservation wages

[xx,n1] = min(abs(LS-50));
[xx,n2] = min(abs(LS-70));

ResW2 = ReservationWages(n1:n2);
LS2   = LS(n1:n2);

n = size(LS2,1);

X = [ones(n,1) LS2 LS2.^2 LS2.^3];
[B,BINT,R,RINT,STATS] = regress(ResW2, X);
STATS(1)
ResW2 = ReservationWages(n1:n2) - R;

Elas = gradient(log(LS2))./gradient(log(ResW2));
LSE_95_225_R = [LS2 Elas];

[xx,x58r] = min(abs(LS2-58));
[xx,x60r] = min(abs(LS2-60));
[xx,x62r] = min(abs(LS2-62));


%   write the labor supply elasticities

fprintf(out, 'Labor Supply Elasticity');
fprintf(out, '\n\n');
fprintf(out, '    Empl         Elasticities\n');
fprintf(out, '%10s', ' ', '(I)', '(II)');
fprintf(out, '\n');
fprintf(out, '    ---------------------------\n');
fprintf(out, '%10.3f', LSE_95_225(x58,1), LSE_95_225(x58,2), LSE_95_225_R(x58r,2));
fprintf(out, '\n');
fprintf(out, '%10.3f', LSE_95_225(x60,1), LSE_95_225(x60,2), LSE_95_225_R(x60r,2));
fprintf(out, '\n');
fprintf(out, '%10.3f', LSE_95_225(x62,1), LSE_95_225(x62,2), LSE_95_225_R(x62r,2));
fprintf(out, '\n\n\n');

fclose(out);

cd ..\..\..\Matlab;


