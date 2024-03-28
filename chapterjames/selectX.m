function [SUMR_NON0,loc_non0]=selectX(B,sigma,N,M,H,X_cell,Nnumx,P11,P22,P33)
unit0 = 1e8;
Rsum_non0 = zeros(Nnumx,1);
R_non0_cell = cell(Nnumx,1);

for x = 1:Nnumx
    X = cell2mat(X_cell(x));
    Temp_intfr0 = zeros(N,M);
    intfr0 = zeros(N,M); 
    P0{1} = P11;
    P0{2} = P22;
    P0{3} = P33;
%     P0{4} = P44;
    for t = 1:M
        P = cell2mat(P0(t));
%         P = sqrt(P_tot/(M*K))*diag(X(:,t));
        for n1 = 1:N
            for n2 = 1:N 
                if n2 ~= n1
                   Temp_intfr0(n2,t) = abs(H(:,n1)'*(P(:,n2)))^2;
                else
                   Temp_intfr0(n2,t) = 0;
                end
            end
            intfr0(n1,t) = sum(Temp_intfr0(:,t));
            R_non0(n1,t) = B*log2(1 + X(n1,t)*(abs(H(:,n1)'*(P(:,n1)))^2)/(X(n2,t)*intfr0(n1,t)+sigma));  
        end
    end
    Rsum_non0(x,1) = sum(R_non0(:))/unit0;
    R_non0_cell{x,1} = R_non0/unit0;
    [SUMR_NON0,loc_non0] = max(Rsum_non0(:,1));
end
  