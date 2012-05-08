clear all;
clc;

% This file creates Figure 4 and 5 from simulated aggregate time series
%   labor supply elasticity

Nskip = 500;
Nend = 3500;


% ====== Heterogeneity + Incomplete Market + Indivisible Labor ========== %

cd ..\FortranRevision\Fluctuation\Output;

load Kdata.txt;
load Rdata.txt;
load Wdata.txt;
load Ldata.txt;
load Edata.txt;
load Ydata.txt;
load Cdata.txt;
load Adata.txt;
load Pdata.txt;
load Idata.txt;

load lwedge.dat;
load lwedge2.dat;
load w_error.dat;
load r_error.dat;

cd ..\..\..\Matlab

A  = (Adata(Nskip+1:Nend));
C  = (Cdata(Nskip+1:Nend));
Y  = (Ydata(Nskip+1:Nend));
E  = (Edata(Nskip+1:Nend));
L  = (Ldata(Nskip+1:Nend));
W  = (Wdata(Nskip+1:Nend));
P  = (Pdata(Nskip+1:Nend));
K  = (Kdata(Nskip+1:Nend));
I  = (Idata(Nskip+1:Nend));
AtoC = A./C;
WtoC = W./C;

gama = 1.5;
%B = AtoC.*E.^(-1/gama); % Preference residuals: B  = A./C.*(E.^(-1));
MRS = C.*E.^(1/gama);
Prod = Y./E;
B = MRS./Prod;

A_HII = A./mean(A)-1;
B_HII = B./mean(B)-1;
C_HII = C./mean(C)-1;
W_HII = W./mean(W)-1;
P_HII = P./mean(P)-1;
Y_HII = Y./mean(Y)-1;
E_HII = E./mean(E)-1;
L_HII = L./mean(L)-1;
K_HII = K./mean(K)-1;
I_HII = I./mean(I)-1;
AtoC_HII = AtoC./mean(AtoC)-1;
WtoC_HII = WtoC./mean(WtoC)-1;
MRS = MRS./mean(MRS)-1;

% ====== Heterogeneity + Complete Market + Indivisible Labor ========== %

cd ..\FortranRevision\CompleteMarketPEA\Output;

load Kdata.txt;
load Rdata.txt;
load Wdata.txt;
load Ldata.txt;
load Edata.txt;
load Ydata.txt;
load Cdata.txt;
load Adata.txt;
load Pdata.txt;
load Idata.txt;

cd ..\..\..\Matlab

A  = (Adata(Nskip+1:Nend));
C  = (Cdata(Nskip+1:Nend));
Y  = (Ydata(Nskip+1:Nend));
E  = (Edata(Nskip+1:Nend));
L  = (Ldata(Nskip+1:Nend));
W  = (Wdata(Nskip+1:Nend));
P  = (Pdata(Nskip+1:Nend));
K  = (Kdata(Nskip+1:Nend));
I  = (Idata(Nskip+1:Nend));
AtoC = A./C;
WtoC = W./C;

gama = 1.5;
MRS = C.*E.^(1/gama);
Prod = Y./E;
B = MRS./Prod;

A_HCI = A./mean(A)-1;
B_HCI = B./mean(B)-1;
C_HCI = C./mean(C)-1;
W_HCI = W./mean(W)-1;
P_HCI = P./mean(P)-1;
Y_HCI = Y./mean(Y)-1;
E_HCI = E./mean(E)-1;
L_HCI = L./mean(L)-1;
K_HCI = K./mean(K)-1;
I_HCI = I./mean(I)-1;
AtoC_HCI = AtoC./mean(AtoC)-1;
WtoC_HCI = WtoC./mean(WtoC)-1;


% ====== Heterogeneity + Incomplete Market + Divisible Labor ========== %

cd ..\FortranRevision\FluctuationDL\Output;

load Kdata.txt;
load Rdata.txt;
load Wdata.txt;
load Ldata.txt;
load Edata.txt;
load Ydata.txt;
load Cdata.txt;
load Adata.txt;
load Pdata.txt;
load Idata.txt;

cd ..\..\..\Matlab

A  = (Adata(Nskip+1:Nend));
C  = (Cdata(Nskip+1:Nend));
Y  = (Ydata(Nskip+1:Nend));
E  = (Edata(Nskip+1:Nend));
L  = (Ldata(Nskip+1:Nend));
W  = (Wdata(Nskip+1:Nend));
P  = (Pdata(Nskip+1:Nend));
K  = (Kdata(Nskip+1:Nend));
I  = (Idata(Nskip+1:Nend));
AtoC = A./C;
WtoC = W./C;

gama = 0.4;
MRS = C.*E.^(1/gama);
Prod = Y./E;
B = MRS./Prod;

A_HID = A./mean(A)-1;
B_HID = B./mean(B)-1;
C_HID = C./mean(C)-1;
W_HID = W./mean(W)-1;
P_HID = P./mean(P)-1;
Y_HID = Y./mean(Y)-1;
E_HID = E./mean(E)-1;
L_HID = L./mean(L)-1;
K_HID = K./mean(K)-1;
I_HID = I./mean(I)-1;
AtoC_HID = AtoC./mean(AtoC)-1;
WtoC_HID = WtoC./mean(WtoC)-1;


% ====== Representative agent model  ========== %

cd ..\FortranRevision\VariableLabor\Output;

load Kdata.txt;
load Rdata.txt;
load Wdata.txt;
load Ldata.txt;
load Edata.txt;
load Ydata.txt;
load Cdata.txt;
load Adata.txt;
load Pdata.txt;
load Idata.txt;

cd ..\..\..\Matlab

A  = (Adata(Nskip+1:Nend));
C  = (Cdata(Nskip+1:Nend));
Y  = (Ydata(Nskip+1:Nend));
E  = (Edata(Nskip+1:Nend));
L  = (Ldata(Nskip+1:Nend));
W  = (Wdata(Nskip+1:Nend));
P  = (Pdata(Nskip+1:Nend));
K  = (Kdata(Nskip+1:Nend));
I  = (Idata(Nskip+1:Nend));
AtoC = A./C;
WtoC = W./C;

gama = 0.4;
MRS = C.*E.^(1/gama);
Prod = Y./E;
B = MRS./Prod;

A_RA = A./mean(A)-1;
B_RA = B./mean(B)-1;
C_RA = C./mean(C)-1;
W_RA = W./mean(W)-1;
P_RA = P./mean(P)-1;
Y_RA = Y./mean(Y)-1;
E_RA = E./mean(E)-1;
L_RA = L./mean(L)-1;
K_RA = K./mean(K)-1;
I_RA = I./mean(I)-1;
AtoC_RA = AtoC./mean(AtoC)-1;
WtoC_RA = WtoC./mean(WtoC)-1;


% ====== 2nd Moments Heterogeneity + Incomplete Market + Indivisible Labor ========== %

cd ..\FortranRevision\Fluctuation2\Output;

load Kdata.txt;
load Rdata.txt;
load Wdata.txt;
load Ldata.txt;
load Edata.txt;
load Ydata.txt;
load Cdata.txt;
load Adata.txt;
load Pdata.txt;
load Idata.txt;

load lwedge.dat;
load lwedge2.dat;
load w_error.dat;
load r_error.dat;

cd ..\..\..\Matlab

A  = (Adata(Nskip+1:Nend));
C  = (Cdata(Nskip+1:Nend));
Y  = (Ydata(Nskip+1:Nend));
E  = (Edata(Nskip+1:Nend));
L  = (Ldata(Nskip+1:Nend));
W  = (Wdata(Nskip+1:Nend));
P  = (Pdata(Nskip+1:Nend));
K  = (Kdata(Nskip+1:Nend));
I  = (Idata(Nskip+1:Nend));
AtoC = A./C;
WtoC = W./C;

gama = 1.5;
%B = AtoC.*E.^(-1/gama); % Preference residuals: B  = A./C.*(E.^(-1));
MRS = C.*E.^(1/gama);
Prod = Y./E;
B = MRS./Prod;

A_HII2 = A./mean(A)-1;
B_HII2 = B./mean(B)-1;
C_HII2 = C./mean(C)-1;
W_HII2 = W./mean(W)-1;
P_HII2 = P./mean(P)-1;
Y_HII2 = Y./mean(Y)-1;
E_HII2 = E./mean(E)-1;
L_HII2 = L./mean(L)-1;
K_HII2 = K./mean(K)-1;
I_HII2 = I./mean(I)-1;
AtoC_HII2 = AtoC./mean(AtoC)-1;
WtoC_HII2 = WtoC./mean(WtoC)-1;
MRS2 = MRS./mean(MRS)-1;



% ========= Plot series from the models together ============= %

%cd temp;

tstart = 3400;
tend = 3450;


range = tstart-Nskip+1:tend-Nskip;
time = 1:1:50';


figure(1)
plot(time, A_HII(range), '-', time, A_HCI(range), '--', time, A_HID(range), '-.', time, A_RA(range), ':', 'linewidth', 2);
%title('Labor Productivity (A)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc A;

figure(2)
plot(time, B_HII(range), '-', time, B_HCI(range), '--', time, B_HID(range), '-.', time, B_RA(range), ':', 'linewidth', 2);
%title('Preference Residuals (B)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc B;
print -depsc Figure5;

figure(3)
plot(time, C_HII(range), '-', time, C_HCI(range), '--', time, C_HID(range), '-.', time, C_RA(range), ':', 'linewidth', 2);
%title('Consumption (C)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc C;

figure(4)
plot(time, W_HII(range), '-', time, W_HCI(range), '--', time, W_HID(range), '-.', time, W_RA(range), ':', 'linewidth', 2);
%title('Wage (W)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc W;

figure(6)
plot(time, Y_HII(range), '-', time, Y_HCI(range), '--', time, Y_HID(range), '-.', time, Y_RA(range), ':', 'linewidth', 2);
%title('Output (Y)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc Y;

figure(7)
plot(time, E_HII(range), '-', time, E_HCI(range), '--', time, E_HID(range), '-.', time, E_RA(range), ':', 'linewidth', 2);
%title('Employment (E)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc E;

figure(8)
plot(time, L_HII(range), '-', time, L_HCI(range), '--', time, L_HID(range), '-.', time, L_RA(range), ':', 'linewidth', 2);
%title('Labor Supply in Efficiency Units (L)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc L;

figure(10)
plot(time, AtoC_HII(range), '-', time, AtoC_HCI(range), '--', time, AtoC_HID(range), '-.', time, AtoC_RA(range), ':', 'linewidth', 2);
%title('Productivity over Consumption (A/C)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc AtoC;

figure(11)
plot(time, WtoC_HII(range), '-', time, WtoC_HCI(range), '--', time, WtoC_HID(range), '-.', time, WtoC_RA(range), ':', 'linewidth', 2);
%title('Wage over Consumption (W/C)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc WtoC;

figure(12)
plot(time, K_HII(range), '-', time, K_HCI(range), '--', time, K_HID(range), '-.', time, K_RA(range), ':', 'linewidth', 2);
%title('Capital (K)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc K;

figure(13)
plot(time, I_HII(range), '-', time, I_HCI(range), '--', time, I_HID(range), '-.', time, I_RA(range), ':', 'linewidth', 2);
%title('Investment (I)', 'FontSize',14)
legend('Heterogeneity + Incomplete Market + Indivisible Labor', ...
       'Heterogeneity + Complete Market + Indivisible Labor', ...
       'Heterogeneity + Incomplete Market + Divisible Labor', ...
       'Representative Agent Model' ,0);
print -dpsc I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(15)
plot(time, E_HII(range), '-', time, A_HII(range), '--', time, K_HII(range), '-.', time, AtoC_HII(range), '-.', 'linewidth', 2);
legend('Employment', 'Aggregate Wage', 'Aggregate Capital', 'Aggregate Wage over Consumption', 0);
grid on
print -dpsc EAKP_HII;

figure(16)
plot(time, E_HID(range), '-', time, A_HID(range), '--', time, K_HID(range), '-.', time, AtoC_HID(range), '-.', 'linewidth', 2);
legend('Employment', 'Aggregate Wage', 'Aggregate Capital', 'Aggregate Wage over Consumption', 0);
grid on
print -dpsc EAKP_HID;

figure(17)
plot(time, E_HCI(range), '-', time, A_HCI(range), '--', time, K_HCI(range), '-.', time, AtoC_HCI(range), '-.', 'linewidth', 2);
legend('Employment', 'Aggregate Wage', 'Aggregate Capital', 'Aggregate Wage over Consumption', 0);
grid on
print -dpsc EAKP_HCI;

figure(18)
plot(time, E_RA(range), '-', time, A_RA(range), '--', time, K_RA(range), '-.', time, AtoC_RA(range), '-.', 'linewidth', 2);
legend('Employment', 'Aggregate Wage', 'Aggregate Capital', 'Aggregate Wage over Consumption', 0);
grid on
print -dpsc EAKP_RA;

figure(61)
plot(time, E_HII(range), 'b-', time, B_HII(range), 'r--', 'linewidth', 2);
legend('ln H', 'ln Wedge', 0);
print -dpsc ModelH_B;
print -depsc Figure4;

figure(62)
plot(time, B_HII(range), '-', time, w_error(range), ':', time, r_error(range), '--','linewidth', 2);
legend('ln Wedge','Forecast Error in Wage' , 'Forecast Error in interest rate',0);
print -dpsc forecast_error;

figure(63)
plot(time,lwedge(range),'-',time,lwedge2(range),':','linewidth', 2);
legend('ln Wedge','ln Wedge_o',0);
print -dpsc wedge_lowerbound;



figure(91)
plot(time, E_HII(range), '-', time, E_HII2(range), '--', time, E_RA(range), ':','linewidth', 2);
title('Employment', 'FontSize',14)
legend('First Moment', ' Second Moments', 'Representative',0);

figure(92)
plot(time, E_HII2(range), 'b-', time, B_HII2(range), 'r--', 'linewidth', 2);
legend('ln H', 'ln Wedge', 0);
title('Second Moment')

figure(93)
plot(time, B_HII(range), '-', time, B_HII2(range), '--','linewidth', 2);
title('Wedge', 'FontSize',14)
legend('First Moment', ' Second Moments', 0);


%cd ..;
