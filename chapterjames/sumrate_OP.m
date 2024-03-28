 clear;
 clc;
 close all;
% beam positin(BP) 
%% system parameters
N = 6; % total beam spot number 60
K = 2; % beam spot number per slot 25
M = 3; % time slot number 20
N_UT = 50; % total user number 200
radius_cov = 300e3; % total covered radius 300e3(NJ-SH),5947e3(China-Ukraine)

%% beam position division
dist = rand(N_UT,1)*radius_cov; % UT's module of original point
azimuth = rand(N_UT,1)*2*pi; % UT's azimuth
UT_x = dist.*cos(azimuth); % UT's X-axis
UT_y = dist.*sin(azimuth); % UT's Y-axis
UT = [UT_x,UT_y]; % UT's coordinate
UT_clust_cell = cell(N,1); % UT coordinate cell of  beam center n
radius_max = zeros(N,1); % beamwidth of beam center n
dist_bb = zeros(N,N); % distance bewteen beam i and beam j
% 调用Kmeans进行波束划分
% X N*P的数据矩阵
% Idx N*1的向量,存储的是每个点的聚类标号
% Ctrs K*P的矩阵,存储的是K个聚类质心位置
% SumD 1*K的和向量,存储的是类间所有点与该类质心点距离之和
% D N*K的矩阵，存储的是每个点与所有质心的距离
% opts = statset('Display','final');
% [Idx,Ctrs,SumD,D] = kmeans(UT,N,'Replicates',50,'Options',opts);
% for n = 1:N
%     UT_clust_cell{n} = UT(Idx==n,:);
%     UT_clust_n1 = cell2mat(UT_clust_cell(n));
%     dist_UT_C = zeros(size(UT_clust_n1,1),1); % distance between UT and beam center
%     for n1 = 1:size(UT_clust_n1,1)
%         dist_UT_C(n1,1) = pdist([Ctrs(n,:);UT_clust_n1(n1,:)]);
%     end
%     radius_max(n,1) = max(dist_UT_C); % beamwidth of beam center k
%     for n2 = 1:N
%         dist_bb(n,n2) = pdist([Ctrs(n,:);Ctrs(n2,:)]);
%     end
%     plot(UT(Idx==n,1),UT(Idx==n,2),'r.','MarkerSize',14)
%     hold on
%     plot(Ctrs(n,1),Ctrs(n,2),'kx','MarkerSize',14,'LineWidth',2)
%     hold on
%     pos = [Ctrs(n,1)-radius_max(n,1),Ctrs(n,2)-radius_max(n,1),2*radius_max(n,1),2*radius_max(n,1)]; 
%     rectangle('Position',pos,'Curvature',[1 1])
%     hold on
%     axis equal;
% end
% save('BP_div_N6.mat')
load('BP_div_N6.mat')


%% channel model parameters
Phi = zeros(N,N); % uniform distribute phaze
lamda = 299792.458e3/20e9; % wavelength: light velocity 299792.458e3/carrier frequenze 20(GHz) 
d0 = 1000e3; % distance between satllite to BP
Gr = 10; % receive antenna gain 40(dB)
Bolzamann = 1.38e-23; % Bolzamann constant 1.38*10^(-23)(J/m)
B = 100e6; % bandwidth 240(MHz)
% B = 1;
Temprt = 293; % noise temperature 293(K)
sigma = 1e-2; % noise variance sigma -126.47(dB）or 1.63(dB）
% for n1 = 1:N
%     for n2 = 1:N
%         Phi(n1,n2) = unifrnd(1,2*pi);
%     end
% end
% save('Phi_N6.mat')
load('Phi_N6.mat')

%% channel model
theta = zeros(N,N); % off-axis angle bewteen beam i and beam j
theta_3dB = zeros(N,1); % 3dB gain angle
u = zeros(N,N); % coefficient of theta
Gp = zeros(N,N); % peak beam gain
Gt = zeros(N,N); % transmit beam gain
for n1 = 1:N
    for n2 = 1:N
%         theta_3dB(n1,n2) = unifrnd(2.3,4)*pi/180;       
        theta(n1,n2) = 0;
%         theta(n1,n2) = 0;
%         theta_3dB(n1,1) = atan(radius_max(n1,1)/d0);
%         theta_3dB(n1,1) = pi/3;
        u(n1,n2) = 2.07123*sin(theta(n1,n2))/sin(theta_3dB(n1,1));
        Gp(n1,n2) = 0.65*(2^2)*(pi^2)/((theta_3dB(n1,1))^2); % peak beam gain
        if theta(n1,n2) == 0
           Gt(n1,n2) = 10;
        else
           Gt(n1,n2) = Gp(n1,n2)*((besselj(1,u(n1,n2)))/(2*u(n1,n2))+(36*besselj(3,u(n1,n2)))/((u(n1,n2))^3))^2;
        end
    end
end
H = zeros(N,N); % channel matrix
d1 = zeros(N,N); % distance between satellite and BP
for n1 = 1:N
    for n2 = 1:N
        d1(n1,n2) = sqrt(d0^2+radius_max(n1,1)^2);
        H(n1,n2) = sqrt(((lamda/(4*pi*d1(n1,n2)))^2)*Gr*Gt(n1,n2)/(Bolzamann*B*Temprt))*exp(1i*Phi(n1,n2));
    end
end

%% sum rate
Numx0 = 100;
rate_thred_non = [0];
rate_thred = [0];
P_tot0 = [70];
% P_tot0 = [40];
intfr_thred = [4];
observe_p = length(P_tot0);
Omega = zeros(N,N);
influ = zeros(N,N);
X_cell = cell(Numx0,1);
P0_cell = cell(M,1);
% parfor i_Dc = 1:observe_L
    %% illumination pattern 
% for x0 = 1:Numx0 
%     X0 = zeros(N,M);
%     [index_X] = illumination(N,M,K);
%     for t =1:M
%         for k = 1:K
%             X0(index_X(k,t),t) = 1;
%         end
%     end
%     X0_cell{x0} = X0;
% end
% save('X0_cell6.mat')
load('X0_cell6.mat')
X_cell = X0_cell;

for iii = 1:observe_p
    P_tot=P_tot0(iii);
    Rsum_ini = 1; 
    Rsum_op = 2;
    i_iter = 0;
% while abs(trace(Rsum_op-Rsum_ini)) >= 1e-1
%     Rsum_ini = Rsum_op; 
    i_iter = i_iter + 1;
    if i_iter == 1
       [Rsum_NON,R_non,loc_non_x,Numx2,MAX_X_non]=NONOP(P_tot,i_iter,H,X_cell,Numx0,rate_thred_non,intfr_thred,sigma);
       rate_thredR = min(R_non(R_non~=0));
       rate_thred = 0.3;
       rate_thred_ite(i_iter) = rate_thredR;
       MAX_X = MAX_X_non;
       loc_x = loc_non_x;
       Rsum_NON0(i_iter) = Rsum_NON;
       R_NON_cell{i_iter} = R_non;
    end
    [P0_cell,Rsum_op,R_op,MAX_Xop]=FP(N,K,M,H,MAX_X,P_tot,rate_thred,sigma);
    Rsum_OP0(i_iter) = Rsum_op;
    loc_x_itr(i_iter) = loc_x;
    MAX_X_iter{i_iter}=MAX_Xop;
    R_OP_cell{i_iter} = R_op;
    fprintf('Main|Ptotal=%d|iter=%d|rsum=%d\n\n',P_tot,i_iter,Rsum_op);
    for t = 1:M
        pt1=[];
        for n = 1:N
            P1 = cell2mat(P0_cell(n,t));
            [V1,D1] = eigs(P1);
            [value1,num1] = max(diag(D1));
            pt2 = abs(sqrt(value1)*V1(:,num1));
            pt1=[pt1,pt2];
        end
        P0_cell{t}=pt1;
    end
%     [SUMR,R0,flag]=verify(P0_cell,H,MAX_Xop,rate_thred);
%    for t = 1:M
%        for n = 1:N
%            P = cell2mat(P0_cell(n,t));
%            [V1,D1] = eigs(P);
%            [value1,num1] = max(diag(D1));
%            pt1 = abs(sqrt(value1)*V1(:,num1));
%            P_cell{n,t} = pt1*pt1';
%        end
%    end
% end
    rate_thred_p(iii) = rate_thred;
    RSUM_NON1(iii) = Rsum_NON;
    R_NON_cell{iii} = R_non;
    RSUM_OP1(iii) = Rsum_op;
    R_OP_cell{iii} = R_op;
end

save('CCCP3_N6_70_r3_BF1.mat')
