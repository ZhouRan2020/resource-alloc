function [X,x,Rsum_op,R_op,status]=FP_Xop(N,K,M,H,P0_cell,v,r_thred,sigma)
alpha = 1e0;
beta = 1e0;
%% matrix definition
for t = 1:M
    for n = 1:N
        e_t = zeros(M,1);
        e_n = zeros(N,1);
        H_hat = H(:,n)*H(:,n)';
        P = cell2mat(P0_cell(n,t));
        e_t(t) = 1; e_n(n) = 1;
        vv_cell{n,t} = 2*v(n,t)*sqrt(real(trace(H_hat*P)))*kron(e_t,e_n);
        d_hat = real(trace(H_hat*P))*ones(N,1);
        d_hat(n) = 0;
        d_cell{n,t} = kron(e_t,d_hat);
    end
end

%% sum-rate maximization (FP_Xop method)
    cvx_quiet(0)
   cvx_begin SDP quiet
%     cvx_begin SDP
    cvx_solver Mosek
%     cvx_solver SeDuMi
%    cvx_solver sdpt3 
%     cvx_precision best
    variable x(N*M,1) binary
%     variable x(N*M,1)
%     variable X(N*M,N*M)
    variable rate(N,M)
    maximize sum(rate(:))
%     find sum(rate(:))
    subject to
    for t = 1:M
        for n = 1:N 
            vv = cell2mat(vv_cell(n,t));
            d = cell2mat(d_cell(n,t));
            1/log(2) * log(1+x'*vv-v(n,t)^2*(x'*d+sigma)) >= rate(n,t);
        end
    end    
%     for n = 1:N
%         sum(rate(n,:)) >= r_thred;
%     end
%     for t = 1:M
%          for n = 1:N
%              rate(n,t) >= 0;
%          end
%     end
    %   
    kron(eye(M),ones(1,N))*x <= K*ones(M,1);
    X = reshape(x,N,M);
    for n = 1:N
        X(n,:)*ones(M,1) >= 1;
    end
%     alpha*[X,x;x',1] >= 0;
%     for i = N*M
%        alpha*X(i,i) == alpha*x(i);
%        beta*x(i) <= beta*1;
%        beta*x(i) >= 0;
%     end
%     alpha*X == semidefinite(N*M);
    cvx_end
    
    if strcmp(cvx_status, 'Solved')||strcmp(cvx_status, 'Inaccurate/Solved')
       status = 1;
       X = reshape(abs(x),N,M);
       x = x;
       Rsum_op = cvx_optval;
       R_op = rate;
    else
       status = 0;
       X = NaN;
       x = NaN;
       Rsum_op = NaN;
       R_op = NaN;
    end