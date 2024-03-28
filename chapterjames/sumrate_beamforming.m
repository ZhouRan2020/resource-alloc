clear all;
clc;
close all;
% beam positin(BP)
%% system parameters
N = 6; % total beam spot number 60 6个波位
K = 2; % beam spot number per slot 25 每个时隙最多同时点亮2个波位
M = 3; % time slot number 20 一共有3个时隙
N_UT = 100; % total user number 200 一共200个用户终端
radius_cov = 300e3; % total covered radius 300e3(NJ-SH),5947e3(China-Ukraine)
%每个波位的覆盖半径300km
%% beam position division
dist = rand(N_UT,1)*radius_cov; % UT's module of original point
%每个用户到波位中心点距离的模
azimuth = rand(N_UT,1)*2*pi; % UT's azimuth
%每个用户到波位中心的角度
UT_x = dist.*cos(azimuth); % UT's X-axis
%用户相对波位中心的x坐标
UT_y = dist.*sin(azimuth); % UT's Y-axis
%用户相对波位中心的y坐标
UT = [UT_x,UT_y]; % UT's coordinate
%两个tall矩阵的并联
UT_clust_cell = cell(N,1); % UT coordinate cell of  beam center n
%为每个波位存储一个矩阵，表示。。。
radius_max = zeros(N,1); % beamwidth of beam center n
%为每个波位存储一个数，表示最大半径
dist_bb = zeros(N); % distance bewteen beam i and beam j
%波位与波位之间的距离
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
load('BP_div_N6.mat')
% load('Phi_N8.mat')

%% channel model parameters
Phi = zeros(N); % uniform distribute phaze 均匀分布的相位
lamda = 299792.458e3/20e9; % wavelength: light velocity 299792.458e3/carrier frequenze 20(GHz)
%波长=c/频率
d0 = 1000e3; % distance between satllite to BP 卫星到波位的高度
Gr = 10; % receive antenna gain 40(dB) 接收天线增益
Bolzamann = 1.38e-23; % Bolzamann constant 1.38*10^(-23)(J/m)
B = 100e6; % bandwidth 240(MHz) 带宽
% B = 1;
Temprt = 293; % noise temperature 293(K) 噪声温度
sigma = 1e-6; % noise variance sigma -126.47(dB）or 1.63(dB）噪声功率
for n1 = 1:N
    for n2 = 1:N
        Phi(n1,n2) = unifrnd(1,2*pi); %生成均匀分布的相位矩阵
    end
end

%% channel model
theta = zeros(N); % off-axis angle bewteen beam i and beam j
%波束之间的离轴角
theta_3dB = zeros(N,1); % 3dB gain angle
%3dB增益角
u = zeros(N); % coefficient of theta
%角度相关系数
Gp = zeros(N); % peak beam gain峰值波束增益
Gt = zeros(N); % transmit beam gain传输波束增益
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
%至此，生成了H信道矩阵
%% sum rate
unit0 = 1e8;
Nnumx = 60;
rate_thred = 0.2;%速率门限，速率的下限
P_tot = 80;%总的功率
observe_L = length(P_tot);%测试点的个数
X = zeros(N,M);%某一个窗口内的点亮情况
Omega = zeros(N);%波束间的干扰系数
influ = zeros(N);%
Rsum_non0 = zeros(Nnumx,observe_L);%为每次实验存储一个数据行向量
R_non0_cell = cell(Nnumx,observe_L);
%为每个实验点存储一个矩阵，记录了一次实验中N*M空间内波位的流量分布

Rsum_non = zeros(Nnumx,observe_L);
R_non_cell = cell(Nnumx,observe_L);

Rsum_op1 = zeros(Nnumx,observe_L);
R_op_cell = cell(Nnumx,observe_L);

SUMR_NON0 = zeros(observe_L,1);%每个测试点一个数字，表示Nnumx次实验加权得来的总值
SUMR_NON = zeros(observe_L,1);
SUMR_OP = zeros(observe_L,1);

% for i_Dc = 1:observe_L
%  parfor x = 1:Nnumx
for x = 1:Nnumx %开始一次独立试验
    %% illumination pattern
    X = zeros(N,M);
    index_X = illumination(N,M,K);%获得一个N*M随机的点亮模式
    for t =1:M
        for k = 1:K
            X(index_X(k,t),t) = 1;%将点亮模式序号转化为0-1的矩阵形式
        end
    end
    x_limit = X(:);%一列，规模是NM
    for n1 = 1:N
        for n2 = 1:N
            %获得干扰矩阵，size是N*N
            Omega(n1,n2) = abs((H(:,n1)')*H(:,n2))/(norm(H(:,n1))*norm(H(:,n2)));
        end
    end
    %克罗内克积：MN*MN
    A = kron(eye(M),Omega);
    % penalty = x_limit'*A*x_limit;
    c_thresh = 0.1;
    %         load('X.mat')
    %         %% sum rate of non-BF-non-Prd case
    %         Temp_intfr0 = zeros(N,M);
    %         intfr0 = zeros(N,M);
    %             for t = 1:M
    %                 W_ini0 = ones(N,N);
    %                 P_ini0 = sqrt(P_tot(i_Dc)/(M*K))*diag(X(:,t));
    %                 for n1 = 1:N
    %                     for n2 = 1:N
    %                         if n2 ~= n1
    %                            Temp_intfr0(n2,t) = abs(H(:,n1)'*(P_ini0(:,n2).*W_ini0(:,n2)))^2;
    %                         else
    %                            Temp_intfr0(n2,t) = 0;
    %                         end
    %                     end
    %                     intfr0(n1,t) = sum(Temp_intfr0(:,t));
    %                     R_non0(n1,t) = B*log2(1 + X(n1,t)*(abs(H(:,n1)'*(P_ini0(:,n1).*W_ini0(:,n1)))^2)/(X(n2,t)*intfr0(n1,t)+sigma));
    %                 end
    %             end
    %             Rsum_non0(x,i_Dc) = sum(R_non0(:))/unit0;
    %             R_non0_cell{x,i_Dc} = R_non0/unit0;
    %
    %        %% sum rate of non-BF-Prd case
    %         Temp_intfr = zeros(N,M);
    %         intfr = zeros(N,M);
    %             for t = 1:M
    % %                 W_ini = ones(N,N);
    %                 W_ini1 = (H')*inv(H*H'+sigma/sqrt(P_tot(i_Dc))/N*eye(N));
    %                 P_ini = sqrt(P_tot(i_Dc)/(M*K))*diag(X(:,t));
    %                 for n1 = 1:N
    %                     for n2 = 1:N
    %                         if n2 ~= n1
    %                            Temp_intfr(n2,t) = abs(H(:,n1)'*(P_ini(:,n2).*W_ini1(:,n2)))^2;
    %                         else
    %                            Temp_intfr(n2,t) = 0;
    %                         end
    %                     end
    %                     intfr(n1,t) = sum(Temp_intfr(:,t));
    %                     R_non(n1,t) = B*log2(1 + X(n1,t)*(abs(H(:,n1)'*(P_ini(:,n1).*W_ini1(:,n1)))^2)/(X(n2,t)*intfr(n1,t)+sigma));
    %                 end
    %             end
    %             Rsum_non(x,i_Dc) = sum(R_non(:))/unit0;
    %             R_non_cell{x,i_Dc} = R_non/unit0;

    %% sum rate of beamforming case
    %             W_ini = ones(N*N,N*N);
    %             W_ini0 = (H')*inv(H*H'+sigma/(P_tot(i_Dc)/N/M)*eye(N));
    %规模：N*N
    W_ini1 = (H')/(H*H'+sigma/sqrt(P_tot)/N*eye(N));
    w = W_ini1(:);%列向量，长度NN
    W_ini2 = w*w';%维度是NN*NN
    P_ini1 = diag(ones(N*N,1));%维度是NN*NN
    P_ini2 = diag(ones(N*N,1));
    P_ini3 = diag(ones(N*N,1));
    %             P_ini4 = diag(ones(N*N,1));
    %             P_ini1 = 1e0*ones(N*N,N*N);
    %             P_ini2 = 1e0*ones(N*N,N*N);
    %             P_ini3 = 1e0*ones(N*N,N*N);
    %             P_ini4 = 1e0*ones(N*N,N*N);
    %             P_ini5 = diag(ones(N*N,1));
    Rsum_ini = 1;
    %             P_ini5 = ones(N*N,N*N);
    %迭代次数
    status = 1;
    i_iteration = 1;
    for itr = 1
        %进入子函数，优化了三个时隙内的波束赋形矩阵，分别为P1，P2，P3
        [P1,P2,P3,Rsum_op,R_op,status] = beamformingM3(N,M,H,X,P_tot,P_ini1,P_ini2,P_ini3,rate_thred);
        fprintf('case=%d|Ptotal=%d|iteration=%d|rsum=%d\n',x,P_tot,i_iteration,Rsum_op);
        if status == 0
            break
        end
    end
    %矩阵的迹是对角线上的元素之和
    %希望Rsum_op和Rsum_ini的差距不能过大
    while abs(trace(Rsum_op-Rsum_ini)) >= 1e-1
        P_ini1 = P1;
        P_ini2 = P2;
        P_ini3 = P3;
        %                    P_ini4 = P4;
        %                    P_ini5 = P5;
        Rsum_ini = Rsum_op;
        [P1,P2,P3,Rsum_op,R_op,status] = beamformingM3(N,M,H,X,P_tot,P_ini1,P_ini2,P_ini3,rate_thred);
        i_iteration = i_iteration +1;
        fprintf('case=%d|Ptotal=%d|iteration=%d|rsum=%d\n',x,P_tot,i_iteration,Rsum_op);
        %                    Rsum_op_itr(i_iteration) = Rsum_op;
        if status == 0
            break
        end
    end
    Rsum_op1(x,1) = sum(R_op(:));
    R_op_cell{x,1} = R_op;
end
%     [SUMR_NON0(i_Dc),loc_non0(i_Dc)] = max(Rsum_non0(:,i_Dc));
%     [SUMR_NON(i_Dc),loc_non(i_Dc)] = max(Rsum_non(:,i_Dc));
%     [SUMR_OP(i_Dc),loc_op(i_Dc)] = max(Rsum_op1(:,i_Dc));
% end

% save('Prd_nonop_N8_enum60_P1.mat')
