clear;
clc;
close all;

N = 6; % total beam spot number 60
K = 2; % beam spot number per slot 25
M = 3; % time slot number 20
N_UT = 100; % total user number 200
radius_cov = 300e3; % total covered radius 300e3(NJ-SH),5947e3(China-Ukraine)
sigma = 1e-2; % noise variance sigma -126.47(dB）

%% Sum rate vesus Ptotal(observation:N) 
figure(1)
load('FP_N6_fixdb_BF11.mat')
RSUM_OP_r4=RSUM_OP1;
load('FP_N6_r3_fixdb_BF1.mat')
RSUM_OP_r3=RSUM_OP1;
load('ite4_N6_fixdb_P.mat')

% P_tot = [40:5:60];
P_tot = [45:5:70];
plot(P_tot,1e2*RSUM_NON(2:1:7),'-ob','LineWidth',1.6,'MarkerSize',8);
hold on;
plot(P_tot,1e2*RSUM_OP_r4(2:1:7),'-or','LineWidth',1.6,'MarkerSize',8);
hold on;
plot(P_tot,1e2*RSUM_OP_r3(1:1:6),'-og','LineWidth',1.6,'MarkerSize',8);
hold on;
% plot(P_tot,1e2*RSUM_OP(2:1:7),'-or','LineWidth',1.6,'MarkerSize',8);
% hold on;

xlabel('$P_{\rm{total}}$ (W)','Interpreter','latex','FontName','Times New Roman','FontSize',14)%$$Electrical Power (w)$$\
ylabel('$R_{\rm{sum}}$ (Mbits/sec/Hz)' ,'Interpreter','latex','FontName','Times New Roman','FontSize',14)
% legend('Non-BF','FP','CCCP')
legend('Non-BF','FP(r=0.4)','FP(r=0.3)')
set(legend,'fontSize',14,'FontName','Times New Roman','interpreter','latex');
grid on%'Uniform ,N=3','TG ,N=3','ABG ,N=3',
