function [P1,P2,P3,Rsum_op,R_op]=CCCP(N,M,H,MAX_X,P_tot,rate_thred)
% X = cell2mat(X_cell(loc_non0));
X = MAX_X;
P1 = 1e0*ones(N*N)+0.1i*ones(N*N); 
P2 = 1e0*ones(N*N)+0.1i*ones(N*N);
P3 = 1e0*ones(N*N)+0.1i*ones(N*N);
Rsum_ini = 1;
Rsum_op = 2;
status = 1;
i_iteration = 1;
% for itr = 1                
%     [P1,P2,P3,Rsum_op,R_op,status] = beamformingM3(N,M,H,X,P_tot,P_ini1,P_ini2,P_ini3,rate_thred);
%     fprintf('Ptotal=%d|iteration=%d|rsum=%d\n',P_tot,i_iteration,Rsum_op);
%     if status == 0                  
%        break
%     end
% end
while abs(trace(Rsum_op-Rsum_ini)) >= 1e-1
       P_ini1 = P1;
       P_ini2 = P2;
       P_ini3 = P3;
       Rsum_ini = Rsum_op;
       [P1,P2,P3,Rsum_op,R_op,status] = beamformingM3(N,M,H,X,P_tot,P_ini1,P_ini2,P_ini3,rate_thred);
       i_iteration = i_iteration +1;
       fprintf('Ptotal=%d|iteration=%d|rsum=%d\n',P_tot,i_iteration,Rsum_op);
%      Rsum_op_itr(i_iteration) = Rsum_op;
       if status == 0
          break
       end  
end
