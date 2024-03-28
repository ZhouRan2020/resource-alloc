function [P0_cell,Rsum_op,R_op,MAX_Xop]=FP(N,K,M,H,MAX_X,P_tot,r_thred,sigma)
% X = MAX_X;
% rand1=[0;0;1;1;0;0];
% X = [rand1,rand1,rand1];
X = ones(N,M);
Rsum_ini = 1;
Rsum_op = 2;
i_iter = 1;
for t = 1:M
    for n = 1:N
        P0_cell{n,t} = 1e2*ones(N,N); % Initialize P0
        P1_cell{n,t} = 1e0*ones(N,N); % Initialize P1
    end
end
% load('P_cell00.mat')
while abs(trace(Rsum_op-Rsum_ini)) >= 1e-1
      Rsum_ini = Rsum_op;
      for t = 1:M
          for n1 = 1:N
              H_hat = H(:,n1)*H(:,n1)';
              for n = 1:N
                  P0 = cell2mat(P0_cell(n,t));
                  Temp_D(n,t) = X(n,t)*real(trace(P0*H_hat));
              end
              Deno(n1,t) = sum(Temp_D(:,t));
              P0 = cell2mat(P0_cell(n1,t));
              for n2 = 1:N
                  p02 = cell2mat(P0_cell(n2,t));
                  if n2 ~= n1
                     Temp_intfr(n2,t) = X(n2,t)*real(trace(p02*H_hat));
                  else
                     Temp_intfr(n2,t) = 0;
                  end
               end
               intfr(n1,t) = sum(Temp_intfr(:,t));
               v(n1,t) = X(n1,t)*real(trace(P0*H_hat))/(intfr(n1,t)+sigma);
%                eta(n1,t) = sqrt((1+v(n1,t))*X(n1,t)*real(trace(P1*H_hat)))/(Deno(n1,t)+sigma);
               eta(n1,t) = sqrt(X(n1,t)*real(trace(P0*H_hat)))/(intfr(n1,t)+sigma);
         end
      end
       [MAX_Xop,MAX_xop,Rsum_op0,R_op0,status0]=FP_Xop(N,K,M,H,P0_cell,eta,0,sigma);
       MAX_XOP{i_iter} = MAX_Xop;
       MAX_xOP{i_iter} = MAX_xop;
       Rsum_NON_cell(i_iter) = Rsum_op0;
       R_NON_cell{i_iter} = R_op0;
       fprintf('Ptotal=%d|iteration=%d|rsum0=%d\n',P_tot,i_iter,Rsum_op0);
       if status0 == 0
          break
       end
%        MAX_Xop = MAX_X;
 %%    
        [P1,P2,P3,Rsum_op,R_op]=CCCP(N,M,H,MAX_Xop,P_tot,r_thred);
        [p01] = GRrecover(P1,N);
        [p02] = GRrecover(P2,N);
        [p03] = GRrecover(P3,N);
        for n = 1:N
            P01 = reshape(p01',N,N);
            P02 = reshape(p02',N,N);
            P03 = reshape(p03',N,N);
            P0_cell{n,1} = P01(:,n)*P01(:,n)';
            P0_cell{n,2} = P02(:,n)*P02(:,n)';
            P0_cell{n,3} = P03(:,n)*P03(:,n)';
        end
       i_iter = i_iter +1;
       fprintf('Ptotal=%d|iteration=%d|rsum=%d\n',P_tot,i_iter,Rsum_op); 
end