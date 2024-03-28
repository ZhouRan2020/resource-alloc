 clear;
 clc;
 close all;
% beam positin(BP) 
%% system parameters
N = 8; % total beam spot number 60
K = 2; % beam spot number per slot 25
M = 4; % time slot number 20
N_UT = 100; % total user number 200
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
% save('BP_div_N8.mat')
load('BP_div_N8.mat')
% load('Phi_N8.mat')

%% channel model parameters
Phi = zeros(N,N); % uniform distribute phaze
lamda = 299792.458e3/20e9; % wavelength: light velocity 299792.458e3/carrier frequenze 20(GHz) 
d0 = 1000e3; % distance between satllite to BP
Gr = 10; % receive antenna gain 40(dB)
Bolzamann = 1.38e-23; % Bolzamann constant 1.38*10^(-23)(J/m)
B = 100e6; % bandwidth 240(MHz)
% B = 1;
Temprt = 293; % noise temperature 293(K)
sigma = 1e-6; % noise variance sigma -126.47(dB）or 1.63(dB）
% for n1 = 1:N
%     for n2 = 1:N
%         Phi(n1,n2) = unifrnd(1,2*pi);
%     end
% end
load('Phi_N8.mat')

%% channel model
theta = zeros(N,N); % off-axis angle bewteen beam i and beam j
theta_3dB = zeros(N,1); % 3dB gain angle
u = zeros(N,N); % coefficient of theta
Gp = zeros(N,N); % peak beam gain
Gt = zeros(N,N); % transmit beam gain
for n1 = 1:N
    for n2 = 1:N
%         theta_3dB(n1,n2) = unifrnd(2.3,4)*pi/180;       
%         theta(n1,n2) = 0;
    theta(n1,n2) = atan(dist_bb(n1,n2)/d0);
    theta_3dB(n1,1) = atan(radius_max(n1,1)/d0);
    u(n1,n2) = 2.07123*sin(theta(n1,n2))/sin(theta_3dB(n1,1));
    Gp(n1,n2) = 0.65*(4^2)*(pi^2)/((theta_3dB(n1,1))^2); % peak beam gain
    if theta(n1,n2) == 0
       Gt(n1,n2) = Gp(n1,n2);
    else
       Gt(n1,n2) = Gp(n1,n2)*((besselj(1,u(n1,n2)))/(2*u(n1,n2))+(36*besselj(3,u(n1,n2)))/((u(n1,n2))^3))^2;
    end
    end
end
H = zeros(N,N); % channel matrix
d1 = zeros(N,N); % distance between satellite and BP
for n1 = 1:N
    for n2 = 1:N
        d1(n1,n2) = sqrt(d0^2+dist_bb(n1,n2)^2);
        H(n1,n2) = sqrt(((lamda/(4*pi*d1(n1,n2)))^2)*Gr*Gt(n1,n2)/(Bolzamann*B*Temprt))*exp(1i*Phi(n1,n2));
    end
end

%% sum rate
unit0 = 1e8;
Nnumx0 = 50;
rate_thred = [0.2:0.4:1.8];
P_tot = [80:20:160];
penalty_thresh = [8.6:0.4:10.2];
observe_p = length(P_tot);
i_Pen = 2;
X = zeros(N,M);
Omega = zeros(N,N);
influ = zeros(N,N);
X0_cell = cell(Nnumx0,1);

% parfor i_Dc = 1:observe_L
ii = 1;
for x0 = 1:Nnumx0    
    %% illumination pattern 
%     X0 = zeros(N,M);
%     [index_X] = illumination(N,M,K);
%     for t =1:M
%         for k = 1:K
%             X0(index_X(k,t),t) = 1;
%         end
%     end
%     X0_cell{x0} = X0;
    load('X0_cell.mat')
    X0 = cell2mat(X0_cell(x0));
    x_limit = X0(:);
    for n1 = 1:N
        for n2 = 1:N
           Omega(n1,n2) = abs((H(:,n1)')*H(:,n2))/(norm(H(:,n1))*norm(H(:,n2)));
        end
    end
    A = kron(eye(M),Omega);
    penalty = x_limit'*A*x_limit;
    if penalty <= penalty_thresh(i_Pen)
       X_cell{ii} = X0;
       ii = ii+1;
    end
end
Nnumx = ii-1;
Rsum_non0 = zeros(Nnumx,observe_p);
R_non0_cell = cell(Nnumx,observe_p);
Rsum_non = zeros(Nnumx,observe_p);
R_non_cell = cell(Nnumx,observe_p);
Rsum_op1 = zeros(Nnumx,observe_p);
R_op_cell = cell(Nnumx,observe_p);
SUMR_NON0 = zeros(observe_p,1);
SUMR_NON = zeros(observe_p,1);
SUMR_OP = zeros(observe_p,1);
for i_Dc = 1:observe_p
    for x = 1:Nnumx
        X = cell2mat(X_cell(x));
        %% sum rate of non-BF-non-Prd case
        Temp_intfr0 = zeros(N,M);
        intfr0 = zeros(N,M);    
            for t = 1:M
                W_ini0 = ones(N,N);
                P_ini0 = sqrt(P_tot(i_Dc)/(M*K))*diag(X(:,t));
                for n1 = 1:N
                    for n2 = 1:N 
                        if n2 ~= n1
                           Temp_intfr0(n2,t) = abs(H(:,n1)'*(P_ini0(:,n2).*W_ini0(:,n2)))^2;
                        else
                           Temp_intfr0(n2,t) = 0;
                        end
                    end
                    intfr0(n1,t) = sum(Temp_intfr0(:,t));
                    R_non0(n1,t) = B*log2(1 + X(n1,t)*(abs(H(:,n1)'*(P_ini0(:,n1).*W_ini0(:,n1)))^2)/(X(n2,t)*intfr0(n1,t)+sigma));  
                end
            end
            Rsum_non0(x,i_Dc) = sum(R_non0(:))/unit0;
            R_non0_cell{x,i_Dc} = R_non0/unit0;

       %% sum rate of non-BF-Prd case
        Temp_intfr = zeros(N,M);
        intfr = zeros(N,M);    
            for t = 1:M
                W_ini1 = ones(N,N);
%                 W_ini1 = (H')*inv(H*H'+sigma/sqrt(P_tot(i_Dc))/N*eye(N));
                P_ini = sqrt(P_tot(i_Dc)/(M*K))*diag(X(:,t));
                for n1 = 1:N
                    for n2 = 1:N 
                        if n2 ~= n1
                           Temp_intfr(n2,t) = abs(H(:,n1)'*(P_ini(:,n2).*W_ini1(:,n2)))^2;
                        else
                           Temp_intfr(n2,t) = 0;
                        end
                    end
                    intfr(n1,t) = sum(Temp_intfr(:,t));
                    R_non(n1,t) = B*log2(1 + X(n1,t)*(abs(H(:,n1)'*(P_ini(:,n1).*W_ini1(:,n1)))^2)/(X(n2,t)*intfr(n1,t)+sigma));  
                end
            end
            Rsum_non(x,i_Dc) = sum(R_non(:))/unit0;
            R_non_cell{x,i_Dc} = R_non/unit0;

    %% sum rate of beamforming case  
            W_ini2 = ones(N*N,N*N);
%             W_ini0 = (H')*inv(H*H'+sigma/(P_tot(i_Dc)/N/M)*eye(N));
%             w = W_ini1(:);
%             W_ini2 = w*w';
            P_ini1 = diag(ones(N*N,1)); 
            P_ini2 = diag(ones(N*N,1));
            P_ini3 = diag(ones(N*N,1));
            P_ini4 = diag(ones(N*N,1));
%             P_ini1 = 1e0*ones(N*N,N*N); 
%             P_ini2 = 1e0*ones(N*N,N*N);
%             P_ini3 = 1e0*ones(N*N,N*N);
%             P_ini4 = 1e0*ones(N*N,N*N);
%             P_ini5 = diag(ones(N*N,1));
            Rsum_ini = 1;
    %         P_ini5 = ones(N*N,N*N);
            status = 1;
            i_iteration = 1;
            for itr = 1                
                [P1,P2,P3,P4,Rsum_op,R_op,status] = beamforming(N,M,H,W_ini2,X,P_tot(i_Dc),P_ini1,P_ini2,P_ini3,P_ini4,rate_thred(1));
                fprintf('case=%d|Ptotal=%d|iteration=%d|rsum=%d\n',x,P_tot(i_Dc),i_iteration,Rsum_op);
                if status == 0                  
                   break
                end
            end
            while abs(trace(Rsum_op-Rsum_ini)) >= 1e-1
                   P_ini1 = P1;
                   P_ini2 = P2;
                   P_ini3 = P3;
                   P_ini4 = P4;
%                    P_ini5 = P5;
                   Rsum_ini = Rsum_op;
                   [P1,P2,P3,P4,Rsum_op,R_op,status] = beamforming(N,M,H,W_ini2,X,P_tot(i_Dc),P_ini1,P_ini2,P_ini3,P_ini4,rate_thred(1));
                   i_iteration = i_iteration +1;
                   fprintf('case=%d|Ptotal=%d|iteration=%d|rsum=%d\n',x,P_tot(i_Dc),i_iteration,Rsum_op);
%                    Rsum_op_itr(i_iteration) = Rsum_op;
                   if status == 0
                      break
                   end  
            end
            Rsum_op1(x,i_Dc) = sum(R_op(:));
            R_op_cell{x,i_Dc} = R_op;
    end
    [SUMR_NON0(i_Dc),loc_non0(i_Dc)] = max(Rsum_non0(:,i_Dc));
    [SUMR_NON(i_Dc),loc_non(i_Dc)] = max(Rsum_non(:,i_Dc));
    [SUMR_OP(i_Dc),loc_op(i_Dc)] = max(Rsum_op1(:,i_Dc));
end

save('NonPrd_all_N8_penlty2.mat')
