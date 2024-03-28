function [SUMR_NON,R_non0,loc_x,Numx2,MAX_X ]=NONOP1(P0_cell,H,X_cell,Numx0,rate_thred,intfr_thred)
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
Rsum_op0 = zeros(Numx0,1);
R_non0_cell = cell(Numx0,1);
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
for nx = 1:Numx1
    X1 = cell2mat(selectX0_cell(nx));
    %% case1
Temp_intfr0 = zeros(N,M);
    intfr0 = zeros(N,M); 
    for t = 1:M
        for n1 = 1:N
            P = cell2mat(P0_cell(n1,t));
            H_hat = H(:,n1)*H(:,n1)';
            for n2 = 1:N
                P0 = cell2mat(P0_cell(n2,t));
                if n2 ~= n1
                   Temp_intfr0(n2,t) = X1(n2,t)*real(trace(P0*H_hat));
                else
                   Temp_intfr0(n2,t) = 0;
                end
            end
            intfr0(n1,t) = sum(Temp_intfr0(:,t));
            R_non0(n1,t) = B*log2(1 + X1(n1,t)*real(trace(P*H_hat))/(intfr0(n1,t)+sigma));  
        end
    end

    for n1 = 1:N
        Dsum(n1) = sum(R_non0(n1,:));
        if Dsum(n1) >= B*rate_thred
           flag(n1) = 1;
        end
    end
    if all(flag)
        loc_selcetx(ii2) = nx;
        r_slct_cell{ii2} = R_non0;
        selectX_cell{ii2} = X1;
        ii2 = ii2+1;
    end
    Numx2 = ii2-1;
end
%% maximum sumrate
for nx_s = 1:Numx2
    sr_slct = cell2mat(r_slct_cell(nx_s));
    Rsum_non0(nx_s) = sum(sr_slct(:));
    R_non0_cell{nx_s} = sr_slct;
end
[SUMR_NON,loc_x] = max(Rsum_non0(:));
R_non0 = cell2mat(R_non0_cell(loc_x));
MAX_X = cell2mat(selectX_cell(loc_x));

end

%% case2
%     Temp_intfr0 = zeros(N,M);
%     intfr0 = zeros(N,M); 
%     P0{1} = P11;
%     P0{2} = P22;
%     P0{3} = P33;
% %     P0{4} = P44;
%     for t = 1:M
%         P = cell2mat(P0(t));
% %         P = sqrt(P_tot/(M*K))*diag(X(:,t));
%         for n1 = 1:N
%             for n2 = 1:N 
%                 if n2 ~= n1
%                    Temp_intfr0(n2,t) = abs(H(:,n1)'*(P(:,n2)))^2;
%                 else
%                    Temp_intfr0(n2,t) = 0;
%                 end
%             end
%             intfr0(n1,t) = sum(Temp_intfr0(:,t));
%             R_non0(n1,t) = B*log2(1 + X(n1,t)*(abs(H(:,n1)'*(P(:,n1)))^2)/(X(n2,t)*intfr0(n1,t)+sigma));  
%         end
%     end
%     Rsum_non0(x,1) = sum(R_non0(:))/unit0;
%     R_non0_cell{x,1} = R_non0/unit0;
%     [SUMR_OP,loc_non0] = max(Rsum_non0(:,1));  
    