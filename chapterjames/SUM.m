 clear;
 clc;
 close all;

load('NonPrd_op_N12_penlty160.mat')

P_tot = [80];

observe_L = length(P_tot);
SUMR_OP = zeros(observe_L,1);
for i_Dc = 1:observe_L
%     [SUMR_NON(i_Dc),loc_non] = max(Rsum_non(:,i_Dc));
    [SUMR_OP(i_Dc),loc_op] = max(Rsum_op1(:,i_Dc));
end
save('NonPrd_op_N12_penlty160_SUM.mat')