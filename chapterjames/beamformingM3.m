function [P1,P2,P3,obj_sr,rate,status] = beamformingM3(N,M,H,X,P_tot,P_ini1,P_ini2,P_ini3,rate_thred)
sigma = 1e-2; % noise variance sigma -126.47(dB）or 1.63(dB）
B = 1;
alpha = 1e0; 
beta = 1e-2;% scale factor
%% matrix definition
%填充了一些自定义的矩阵
for t = 1:M
    for n = 1:N
        for n2 = 1:N
            if n2 ~= n
                Temp_X1(n2,t) = 0;
                Temp_X2(n2,t) = X(n2,t);
            else
                Temp_X1(n2,t) = X(n2,t);
                Temp_X2(n2,t) = 0;
            end
        end
        H_hat = H(:,n)*H(:,n)';
        G_cell{n,t} = kron(diag(X(:,t)),H_hat);
        G_hat_cell{n,t} = kron(diag(Temp_X2(:,t)),H_hat);
    end
end

%% sum-rate maximization problem (CCCP method)
cvx_quiet(0)
cvx_begin SDP quiet
% cvx_begin SDP
% cvx_solver Mosek
cvx_solver SeDuMi
% cvx_solver sdpt3
cvx_precision best
variable P1(N*N,N*N) complex symmetric
variable P2(N*N,N*N) complex symmetric
variable P3(N*N,N*N) complex symmetric
variable rate(N,M)
variable Dsum(N,1)
variable Rsum
maximize Rsum
subject to
% t = 1 time slot
for n = 1:N
    G = cell2mat(G_cell(n,1));
    G_hat = cell2mat(G_hat_cell(n,1));
    Temp_sr(n,1) = B*1/log(2) * log(sigma + real(trace(P1*G) ))- B*1/log(2) * log(sigma + real(trace(P_ini1*G_hat)))...
        - B*real(trace(G_hat*((P1 - P_ini1)))/(log(2)*(sigma + real(trace(P_ini1*G_hat)))));
    alpha*Temp_sr(n,1) >= alpha*rate(n,1);
end
% t = 2 time slot
for n = 1:N
    G = cell2mat(G_cell(n,2));
    G_hat = cell2mat(G_hat_cell(n,2));
    Temp_sr(n,2) = B*1/log(2) * log(sigma + real(trace(P2*G) ))- B*1/log(2) * log(sigma + real(trace((P_ini2)*G_hat)))...
        - B*real(trace(G_hat*((P2 - P_ini2)))/(log(2)*(sigma + real(trace(P_ini2*G_hat)))));
    alpha*Temp_sr(n,2) >= alpha*rate(n,2);
end
% t = 3 time slot
for n = 1:N
    G = cell2mat(G_cell(n,3));
    G_hat = cell2mat(G_hat_cell(n,3));
    Temp_sr(n,3) = B*1/log(2) * log(sigma + real(trace(P3*G) ))- B*1/log(2) * log(sigma + real(trace((P_ini3)*G_hat)))...
        - B*real(trace(G_hat*((P3 - P_ini3)))/(log(2)*(sigma + real(trace(P_ini3*G_hat)))));
    alpha*Temp_sr(n,3) >= alpha*rate(n,3);
end

alpha*rate(:) >= 0;

for n = 1:N
    Dsum(n) == sum(rate(n,:));
    Dsum(n) >= rate_thred;
end
Rsum == sum(rate(:));

beta*real(trace(P1)) <= beta*P_tot/M;
beta*real(trace(P2)) <= beta*P_tot/M;
beta*real(trace(P3)) <= beta*P_tot/M;

beta*P1 ==  hermitian_semidefinite(N*N);
beta*P2 ==  hermitian_semidefinite(N*N);
beta*P3 ==  hermitian_semidefinite(N*N);

cvx_end

if strcmp(cvx_status, 'Solved')||strcmp(cvx_status, 'Inaccurate/Solved')
    %求解成功
    status = 1;
    P1 = P1;
    P2 = P2;
    P3 = P3;
    obj_sr = cvx_optval;
    rate = rate;
else
    %求解失败
    status = 0;
    P1 = NaN;
    P2 = NaN;
    P3 = NaN;
    obj_sr = NaN;
    rate = NaN;
end

end
