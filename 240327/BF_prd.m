function [P1,P2,obj_sr,rate,status] = BF_prd(N,M,H,B,W_cell,X_t_precode,P_ini1,P_ini2,rate_thred)
sigma = 2^(-126.47/10); % noise variance sigma -126.47(dB）or 1.63(dB）
B = 1; % unit bandwidth
alpha = 10^(0); beta = 10^(0);% scale factor
%% matrix definition
for t = 1:M 
    W = cell2mat(W_cell(t));
    Temp_W = [];
    for n1 = 1:N            
        Temp_W = [Temp_W;W(:,n1)];
    end
    W_hat_cell{t} = Temp_W*Temp_W';
end
for n = 1:N    
    for t = 1:M 
        for n2 = 1:N 
            if n2 ~= n
               Temp_BD1{n2,t} = 0;
               Temp_BD2{n2,t} = X_t_precode(n2,t);
            else
               Temp_BD1{n2,t} = X_t_precode(n2,t);
               Temp_BD2{n2,t} = 0;
            end
        end               
        BD1_cell{n,t} = blkdiag(Temp_BD1{:,t});
        BD2_cell{n,t} = blkdiag(Temp_BD2{:,t});
        H_hat_cell{n} = H(:,n)*H(:,n)';
        BD1 = cell2mat(BD1_cell(n,t));
        BD2 = cell2mat(BD2_cell(n,t));
        H_hat = cell2mat(H_hat_cell(n));
        G_cell{n,t} = kron(BD1,H_hat);
        G_hat_cell{n,t} = kron(BD2,H_hat);
    end
end

%% sum-rate maximization problem (CCCP method)
cvx_quiet(0)
% cvx_solver SeDuMi
cvx_solver sdpt3 
cvx_begin SDP
cvx_precision best
variable P1(N*N,N*N)
variable P2(N*N,N*N)
variable rate(N,M)
variable Rsum
maximize Rsum
subject to
% t = 1 time slot
for n = 1:N
    G = cell2mat(G_cell(n,1));
    G_hat = cell2mat(G_hat_cell(n,1));
    W_hat = cell2mat(W_hat_cell(1));
    Temp_sr(n,1) = B*1/log(2) * log(sigma + real(trace((P1.*W_hat)*G) + trace((P1.*W_hat)*G_hat)))- B*1/log(2) * log(sigma + real(trace((P_ini1.*W_hat)*G_hat)))...
    - B*real(trace(G_hat*((P1 - P_ini1).*W_hat))/(log(2)*(sigma + real(trace((P_ini1.*W_hat)*G_hat)))));
    alpha*Temp_sr(n,1) >= alpha*rate(n,1);
    if X_t_precode(n,1)~=0
       rate(n,1) >= rate_thred;
    end
end
% t = 2 time slot
for n = 1:N
    G = cell2mat(G_cell(n,2));
    G_hat = cell2mat(G_hat_cell(n,2));
    W_hat = cell2mat(W_hat_cell(2));
    Temp_sr(n,2) = B*1/log(2) * log(sigma + real(trace((P2.*W_hat)*G) + trace((P2.*W_hat)*G_hat)))- B*1/log(2) * log(sigma + real(trace((P_ini2.*W_hat)*G_hat)))...
    - B*real(trace(G_hat*((P2 - P_ini2).*W_hat))/(log(2)*(sigma + real(trace((P_ini2.*W_hat)*G_hat)))));
    alpha*Temp_sr(n,2) >= alpha*rate(n,2);
    if X_t_precode(n,2)~=0
       rate(n,2) >= rate_thred;
    end
end
alpha*Rsum == alpha*sum(rate(:));
beta*P1(:) >= 0;
beta*P2(:) >= 0;
beta*rate(:) >= 0;
alpha*(trace(P1)+trace(P2)) <= alpha*100;
alpha*P1 ==  semidefinite(N*N);
alpha*P2 ==  semidefinite(N*N);

cvx_end

if strcmp(cvx_status, 'Solved')||strcmp(cvx_status, 'Inaccurate/Solved')
   status = 1;
   P1 = P1;
   P2 = P2;
   obj_sr = cvx_optval;
   rate = rate;
else
   status = 0;
   P1 = NaN;
   P2 = NaN;
   obj_sr = NaN;
   rate = NaN;
end

end