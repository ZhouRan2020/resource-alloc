function [SUMR,R0,flag]=verify(P0_cell,H,MAX_Xop,rate_thred)
%% system parameters
N = 6; % total beam spot number 60
K = 2; % beam spot number per slot 25
M = 3; % time slot number 20
B = 1; % bandwidth 240(MHz)
sigma = 1e-2; % noise variance sigma -126.47(dB
%% 
flag = zeros(N,1);
%% 
X1 = MAX_Xop;
    %% case1
Temp_intfr0 = zeros(N,M);
intfr0 = zeros(N,M); 
for t = 1:M
    P = cell2mat(P0_cell(t));
    for n1 = 1:N
        for n2 = 1:N 
            if n2 ~= n1
               Temp_intfr0(n2,t) = X1(n2,t)*(abs((H(:,n1)'*P(:,n2)))^2);
            else
               Temp_intfr0(n2,t) = 0;
            end
        end
        intfr0(n1,t) = sum(Temp_intfr0(:,t));
        R0(n1,t) = B*log2(1 + X1(n1,t)*((abs(H(:,n1)'*P(:,n1)))^2)/(intfr0(n1,t)+sigma));  
    end
end
for n1 = 1:N
    Dsum(n1) = sum(R0(n1,:));
    if Dsum(n1) >= B*rate_thred
       flag(n1) = 1;
    end
end
SUMR = sum(R0(:));
end

    