function [SUMR_OP,R_op0,loc_x,Numx2,MAX_X]=selectX_OP(H,X_cell,Numx0,P_tot,P0_cell_op,rate_thred,intfr_thred,i_iter)
%% system parameters
N = 6; % total beam spot number 60
K = 2; % beam spot number per slot 25
M = 3; % time slot number 20
B = 1; % bandwidth 240(MHz)
sigma = 1e-2; % noise variance sigma -126.47(dB
%% 
unit0 = 1e8;
i_penlt = length(intfr_thred);
flag = zeros(N,1);
% Rsum_op0 = zeros(Numx0,1);
% R_op0_cell = cell(Numx0,1);
ii1 = 1;
ii2 = 1;
%%
for nx = 1:Numx0
    X = cell2mat(X_cell(nx));
%% constrants of non-adjacency beam
    xv_t = X(:);
    Ht = H';
    for n1 = 1:N
        for n2 = 1:N
           Omega(n1,n2) = abs((H(n1,:))*Ht(n2,:)')/(norm(H(n1,:))*norm(Ht(n2,:)));
        end
    end
    A = kron(eye(M),Omega);
    intfr = xv_t'*A*xv_t;
    if intfr <= intfr_thred(i_penlt)
       selectX0_cell{ii1} = X;
       ii1 = ii1+1;
    end
    Numx1 = ii1-1;
end
%% constraint of rate threshold
% P_ini = (P_tot/M)*ones(N*N);
for nx = 1:Numx1
    X1 = cell2mat(selectX0_cell(nx)); 
    for t0 = 1:M
     for n1 = 1:N
       for n2 = 1:N 
            if n2 ~= n1
                Temp_X1(n2,t0) = 0;
                Temp_X2(n2,t0) = X1(n2,t0);
            else
                Temp_X1(n2,t0) = X1(n2,t0);
                Temp_X2(n2,t0) = 0;
            end
        end 
        H_hat = H(:,n1)*H(:,n1)';
        G_t = kron(diag(X1(:,t0)),H_hat);
        G_hat_t = kron(diag(Temp_X2(:,t0)),H_hat);
        if i_iter == 1
            for t1 = 1:M
%                 p_ini = ones(N,1).*X1(:,t1);
                p_ini = ones(N,1);
                p_ini_hat = [p_ini',p_ini',p_ini',p_ini',p_ini',p_ini']';
                P_ini = p_ini_hat*p_ini_hat';
                P0_cell{t1} = (P_tot/M)*P_ini;
                rate_thred =0;
            end
        else
            P0_cell = P0_cell_op;
        end
%         P0_cell = P0_cell_op;
        P0 = cell2mat(P0_cell(t0));
        sr(n1,t0) = (B*1/log(2) * log(sigma + real(trace(P0*G_t)))...
                     - B*1/log(2) * log(sigma + real(trace(P0*G_hat_t))));
%         sr1(n1,t0) = B*log2(1 + real(trace(G_t*P0))/(real(trace(G_hat_t*P0))+sigma));
        
     end
    end
    for n1 = 1:N
        Dsum(n1) = sum(sr(n1,:));
        if Dsum(n1) >= B*rate_thred
           flag(n1) = 1;
        end
    end
    if all(flag)
        loc_selcetx(ii2) = nx;
        sr_slct_cell{ii2} = sr;
        selectX_cell{ii2} = X1;
        ii2 = ii2+1;
    end
    Numx2 = ii2-1;
end
%% maximum sumrate
for nx_s = 1:Numx2
    sr_slct = cell2mat(sr_slct_cell(nx_s));
    Rsum_op0(nx_s) = sum(sr_slct(:));
    R_op0_cell{nx_s} = sr_slct;
end
[SUMR_OP,loc_x] = max(Rsum_op0(:));
R_op0 = cell2mat(R_op0_cell(loc_x));
MAX_X = cell2mat(selectX_cell(loc_x));
end  