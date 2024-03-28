 clear;
 clc;
 close all;
 
%% system parameters
N_UT = 100; % total user number 200
N = 10; % total beam spot number 60
radius_cov = 300e3; % total covered radius 300e3(NJ-SH),5947e3(China-Ukraine)

%% beam position division
dist = rand(N_UT,1)*radius_cov; % UT's module to original point
azimuth = rand(N_UT,1)*2*pi; % UT's azimuth
UT_x = dist.*cos(azimuth); % UT's X-axis
UT_y = dist.*sin(azimuth); % UT's Y-axis
UT = [UT_x,UT_y]; % UT's coordinate

%% 调用Kmeans函数
%X N*P的数据矩阵
%Idx N*1的向量,存储的是每个点的聚类标号
%Ctrs K*P的矩阵,存储的是K个聚类质心位置
%SumD 1*K的和向量,存储的是类间所有点与该类质心点距离之和
%D N*K的矩阵，存储的是每个点与所有质心的距离
opts = statset('Display','final');
[Idx,Ctrs,SumD,D] = kmeans(UT,N,'Replicates',50,'Options',opts);

UT_clust_cell = cell(N,1);
radius_max = zeros(N,1);

figure(1); % divided beam position by kmeans
for n = 1:N
    UT_clust_cell{n} = UT(Idx==n,:);
    UT_clust_k = cell2mat(UT_clust_cell(n));
    dist_UT_C = zeros(size(UT_clust_k,1),1);
    for kk = 1:size(UT_clust_k,1)
        dist_UT_C(kk,1) = pdist([Ctrs(n,:);UT_clust_k(kk,:)]);
    end
    radius_max(n,1) = max(dist_UT_C); % beamwidth of beam center k
    plot(UT(Idx==n,1),UT(Idx==n,2),'r.','MarkerSize',14)
    hold on
    plot(Ctrs(n,1),Ctrs(n,2),'kx','MarkerSize',14,'LineWidth',2)
    hold on
    pos = [Ctrs(n,1)-radius_max(n,1),Ctrs(n,2)-radius_max(n,1),2*radius_max(n,1),2*radius_max(n,1)]; 
    rectangle('Position',pos,'Curvature',[1 1])
    hold on
    axis equal;
end

figure(2); % UT's distribution 
plot(UT_x,UT_y,'sr','LineWidth',1.6,'MarkerSize',8)
set(gca,'XLim',[-radius_cov,radius_cov]);
set(gca,'YLim',[-radius_cov,radius_cov]);
legend('UT')
set(legend,'fontSize',14,'FontName','Times New Roman','interpreter','latex');
hold on;
