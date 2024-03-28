function [W_hat,P1,P2,Rsum_non,Rsum_op1] = prd(N,M,H_cell,X_exhau,beamwidth)
B = 1; % unit bandwidth
%% generate influence factor
for k =1:length(beamwidth)
    H = cell2mat(H_cell(k));
   for n1 = 1:N
       for n2 = 1:N
%            influ(n1,n2) = abs((H(:,n1)')*H(:,n2))/(norm(H(:,n1))*norm(H(:,n2)));
           influ(n1,n2) = abs((H(:,n1)')*H(:,n2))/(norm(H(:,n1))*norm(H(:,n2)));
       end
   end    
%  influ =  [1,0.221617628315299,0.353144912337155,0.268325293108536;0.221617628315299,1,0.345303017432870,0.638870175953735;0.353144912337155,0.345303017432870,1,0.217999740678818;0.268325293108536,0.638870175953735,0.217999740678818,1];
   c_thresh = 0.1;
   
%% generate cluster and precoding matrix
%  for j1 = 182:length(Num_act) % 182(M=2),251(M=2)
 for j1 = 1 % 182(M=2),251(M=2)
     X_t = cell2mat(X_exhau(j1));
    % -----------------------------------------------------------
    % -----------------------------------------------------------     
     for t = 1:M
         num_act(t) = 0;
         for n = 1:N           
             if (X_t(n,t) == 1)
                num_act(t) = num_act(t)+1;
                index_act{num_act(t),1} = n;
                index_clnum(num_act(t),1) = n;
             else
                num_act(t) = num_act(t); 
             end
             W_hat0 = zeros(N-num_act(t)); % the precoding matrix is 0 for non-illuminated beam 
         end
       % clustering(non-repeat step)-----------------------------------------------------------------------
         num_size = [];
         n1 = 1;  
         while n1 <= length(index_clnum)
               del=[];
               size_clus=1;
               if n1 ~= length(index_clnum)
                  for n2 = n1+1 : length(index_clnum)
                      index_n1 = cell2mat(index_act(n1));
                      index_n2 = cell2mat(index_act(n2));
                      if ~(influ(index_n2,index_clnum(n1)) < c_thresh)
                          del = [del;n2];
                          index_act{n1} = [index_n1;index_n2];
                          size_clus = size_clus+1;
                      end
                  end
               end
               index_clnum(del) = [];
               index_act(del) = [];
               num_size(n1) = size_clus;
%                index_loc = sum(num_size(:)==size_clus);
%                clus_list{index_loc,size_clus} = index_act{n1};
               n1=n1+1;
         end
         clus_list{t} = index_act; % cluster list in t time slot
       % precoding ------------------------------------------------------------
         sigma = 2^(-126.47/10); % noise variance sigma -126.47(dB）
         P_tot = 100; % power total 100(W)
         L(t) = length(index_act); % number of cluster in t time slot
         for i1=1:L(t)
             index_act_i1 = cell2mat(index_act(i1));
             for ii1 = 1:length(index_act_i1)
                 for ii2 = 1:length(index_act_i1)
                     Hii(ii1,ii2) = H(index_act_i1(ii1),index_act_i1(ii2));                    
                 end
             end
             Hi{i1,t} = Hii;
             if (length(index_act_i1) ~= 1)
                 W_hat{i1,t} = (Hii')*inv(Hii*Hii'+sigma/(P_tot/N)*eye(length(index_act_i1)));
                 Hii = [];
             else
                 W_hat{i1,t} = 1;
                 Hii = [];
             end
         end        
         W_cell{t} = blkdiag(W_hat{:,t},W_hat0); % BDiag(W_hat)       
     end
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    %% the sum rate of non-beamforming case
    for t = 1:M
        P_ini0_cell{t} = (P_tot/(M*num_act(t)))*diag([ones(num_act(t),1);zeros((N-num_act(t)),1)]); % initialize Power matrix
        W_ini = cell2mat(W_cell(t));
        P_ini0 = cell2mat(P_ini0_cell(t));
        for n1 = 1:N
            for n2 = 1:N 
                if n2 ~= n1
                   Temp_intfr(n2,t) = norm(H(:,n2)*H(:,n2)'*(P_ini0(:,n2).*W_ini(:,n2)))^2;
                else
                   Temp_intfr(n2,t) = 0;
                end
            end
            intfr(n1,t) = sum(Temp_intfr(:,t));
            R_non(n1,t) = B*log2(1 + (norm(H(:,n1)*H(:,n1)'*(P_ini0(:,n1).*W_ini(:,n1)))^2)/(intfr(n1,t)+sigma));  
        end
    end
    Rsum_non(1,k) = sum(R_non(:));
    R_non_cell{1,k} = R_non;
    demand1 = min(R_non(find(R_non~=0)));
    demand2 = mean(R_non(find(R_non~=0)));
    demand3 = max(R_non(find(R_non~=0)));
    rate_thred = [demand1:0.1:demand2];
    Rsum_non1 = Rsum_non(1,1)*ones(1,length(rate_thred));
    %% the sum rate of beamforming case
 for ii = 1:length(rate_thred)   
    for t = 1:M
        for n = 1:N
            X_t_precode(n,t) = (R_non(n,t)~=0);
        end
    end   
    P_ini1 = (P_tot/(M*num_act(1)^2))*diag([ones(N*N-(N-num_act(1))^2,1);zeros((N-num_act(1))^2,1)]); % initialize Power matrix
    P_ini2 = (P_tot/(M*num_act(2)^2))*diag([ones(N*N-(N-num_act(2))^2,1);zeros((N-num_act(2))^2,1)]); % initialize Power matrix
    status = 1;
    i_iteration = 1;
    for itr = 1
        fprintf('case=%d|iteration=%d\n',j1,i_iteration);
        [P1,P2,Rsum_op,R_op,status] = beamforming(N,M,H,B,W_cell,X_t_precode,P_ini1,P_ini2,rate_thred(ii));
        Rsum_op_itr(i_iteration) = Rsum_op;
        if status == 0                  
           break
        end
    end
    while  abs(trace(P1-P_ini1)) >= 0.001 && abs(trace(P2-P_ini2)) >= 0.001
           P_ini1 = P1;
           P_ini2 = P2;
           [P1,P2,Rsum_op,R_op,status] = beamforming(N,M,H,B,W_cell,X_t_precode,P_ini1,P_ini2,rate_thred(ii));          
           i_iteration = i_iteration + 1;
           Rsum_op_itr(i_iteration) = Rsum_op;
           fprintf('case=%d|iteration=%d\n',j1,i_iteration);
           if status == 0
              break
           end  
    end
    Rsum_op1(1,ii) = sum(R_op(:));
    R_op_cell{1,ii} = R_op;
 end 
 end
end

%  save('test_rate_threshd_beamwidth1.mat')
%% figure
figure(1)
plot(rate_thred,Rsum_non1,'-ob','LineWidth',1.6,'MarkerSize',8);
hold on;
plot(rate_thred,Rsum_op1,'-or','LineWidth',1.6,'MarkerSize',8);
hold on;
title('Sum Rate  vesus Rate Requirements','Interpreter','latex','FontName','Times New Roman','FontSize',20);
xlabel('${\bar R}$(bits/sec/Hz)','Interpreter','latex','FontName','Times New Roman','FontSize',14)%$$Electrical Power (w)$$\
ylabel('$R_{\rm{sum}}$(bits/sec/Hz)' ,'Interpreter','latex','FontName','Times New Roman','FontSize',14)
legend('Benchmark','Beamforming')
set(legend,'fontSize',14,'FontName','Times New Roman','interpreter','latex');
grid on%'Uniform ,N=3','TG ,N=3','ABG ,N=3',

figure(2)
plot([1:1:i_iteration],Rsum_op_itr,'-xb','LineWidth',1.6,'MarkerSize',8);
hold on;
% title('CCCP convergence','Interpreter','latex','FontName','Times New Roman','FontSize',20);
xlabel('Number of iteration','Interpreter','latex','FontName','Times New Roman','FontSize',14)%$$Electrical Power (w)$$\
ylabel('$R_{\rm{sum}}$(bits/sec/Hz)' ,'Interpreter','latex','FontName','Times New Roman','FontSize',14)
legend('CCCP convergence')
set(legend,'fontSize',14,'FontName','Times New Roman','interpreter','latex');
grid on%'Uniform ,N=3','TG ,N=3','ABG ,N=3',


end