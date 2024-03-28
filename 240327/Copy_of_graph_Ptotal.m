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
load('FP_N6_r3_BigP_BF1.mat')
load('ite5_N6_BigP_fixdb.mat')
load('FP_N6_80120_r3_BF1.mat')

% P_tot = [40:5:60];
% P_tot = [60:20:120];
plot(P_tot0,1e2*RSUM_NON1,'-ob','LineWidth',1.6,'MarkerSize',8);
hold on;
plot(P_tot0,1e2*RSUM_OP1,'-or','LineWidth',1.6,'MarkerSize',8);
hold on;


xlabel('$P_{\rm{total}}$ (W)','Interpreter','latex','FontName','Times New Roman','FontSize',14)%$$Electrical Power (w)$$\
ylabel('$R_{\rm{sum}}$ (Mbits/sec/Hz)' ,'Interpreter','latex','FontName','Times New Roman','FontSize',14)
% legend('Non-BF','FP','CCCP')
legend('Non-BF','FP(r=0.3)')
set(legend,'fontSize',14,'FontName','Times New Roman','interpreter','latex');
grid on%'Uniform ,N=3','TG ,N=3','ABG ,N=3',
