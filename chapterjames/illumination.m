function [index_X] = illumination(N,M,K)
round1 = zeros(K,M);
index_BP0 = find(zeros(N,1)==0);
for t = 1:M       
    R = randperm(length(index_BP0));
    choosekN_K = [];
for k = 1: K
    choosekN_K = [choosekN_K;index_BP0(R(k))];
end
    round1(:,t) = choosekN_K;
    round2 = index_BP0;
    tempt_N_K = setdiff(round2,choosekN_K);
    if length(round2) >= K
       index_BP0 = tempt_N_K;
    else
       break
    end
end
index_X = round1;
end