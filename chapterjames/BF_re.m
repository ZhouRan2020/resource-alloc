function [P_cell,Rsum_op,R_op,status]=BF_re(N,M,H,X,P_tot,v,eta,r_thred,sigma)
%% matrix definition
for t = 1:M
    for n = 1:N
        G_cell{n,t} = H(:,n)*H(:,n)';
    end
end   

%% sum-rate maximization problem (FP method)
    cvx_quiet(0)
    cvx_begin SDP quiet
%     cvx_begin SDP
%     cvx_solver Mosek
%      cvx_solver SeDuMi
    cvx_solver sdpt3 
    cvx_precision best
    variable P11(N,N) complex
    variable P21(N,N) complex
    variable P31(N,N) complex
    variable P41(N,N) complex
    variable P51(N,N) complex
    variable P61(N,N) complex
    variable P12(N,N) complex
    variable P22(N,N) complex
    variable P32(N,N) complex
    variable P42(N,N) complex
    variable P52(N,N) complex
    variable P62(N,N) complex
    variable P13(N,N) complex
    variable P23(N,N) complex
    variable P33(N,N) complex
    variable P43(N,N) complex
    variable P53(N,N) complex
    variable P63(N,N) complex
    variables z(N,M)
    variable rate(N,M)
    variable Dsum(N,1)
    variable Rsum
    maximize Rsum
    subject to
%     P_cell{1,1}=P11;P_cell{2,1}=P21;P_cell{3,1}=P31;P_cell{4,1}=P41;P_cell{5,1}=P51;P_cell{6,1}=P61;
    % t = 1 time slot
        G11 = cell2mat(G_cell(1,1));
        G21 = cell2mat(G_cell(2,1));
        G31 = cell2mat(G_cell(3,1));
        G41 = cell2mat(G_cell(4,1));
        G51 = cell2mat(G_cell(5,1));
        G61 = cell2mat(G_cell(6,1));
        intf11 = real(X(2,1)*trace(P21*G11)+X(3,1)*trace(P31*G11)...
            +X(4,1)*trace(P41*G11)+X(5,1)*trace(P51*G11)+X(6,1)*trace(P61*G11));
        f(1,1) = 1/log(2) * log(1+2*eta(1,1)*z(1,1)-(intf11+sigma)*eta(1,1)^2);
        f(1,1) >= rate(1,1);
        {z(1,1),X(1,1),trace(G11*P11)} == rotated_lorentz(1);
        %
        intf21 = real(X(1,1)*trace(P11*G21)+X(3,1)*trace(P31*G21)...
            +X(4,1)*trace(P41*G21)+X(5,1)*trace(P51*G21)+X(6,1)*trace(P61*G21));        
        f(2,1) = 1/log(2) * log(1+2*eta(2,1)*z(2,1)-(intf21+sigma)*eta(2,1)^2);
        f(2,1) >= rate(2,1);
        {z(2,1),X(2,1),trace(G21*P21)} == rotated_lorentz(1);
        %
        intf31 = real(X(1,1)*trace(P11*G31)+X(2,1)*trace(P21*G31)...
            +X(4,1)*trace(P41*G31)+X(5,1)*trace(P51*G31)+X(6,1)*trace(P61*G31));        
        f(3,1) = 1/log(2) * log(1+2*eta(3,1)*z(3,1)-(intf31+sigma)*eta(3,1)^2);
        f(3,1) >= rate(3,1);
        {z(3,1),X(3,1),trace(G31*P31)} == rotated_lorentz(1);
        %
        intf41 = real(X(1,1)*trace(P11*G41)+X(2,1)*trace(P21*G41)+X(3,1)*trace(P31*G41)...
            +X(5,1)*trace(P51*G41)+X(6,1)*trace(P61*G41));        
        f(4,1) = 1/log(2) * log(1+2*eta(4,1)*z(4,1)-(intf41+sigma)*eta(4,1)^2);
        f(4,1) >= rate(4,1);
        {z(4,1),X(4,1),trace(G41*P41)} == rotated_lorentz(1);
        %
        intf51 = real(X(1,1)*trace(P11*G51)+X(2,1)*trace(P21*G51)+X(3,1)*trace(P31*G51)...
            +X(4,1)*trace(P41*G51)+X(6,1)*trace(P61*G51));        
        f(5,1) = 1/log(2) * log(1+2*eta(5,1)*z(5,1)-(intf51+sigma)*eta(5,1)^2);
        f(5,1) >= rate(5,1);
        {z(5,1),X(5,1),trace(G51*P51)} == rotated_lorentz(1);
        %
        intf61 = real(X(1,1)*trace(P11*G61)+X(2,1)*trace(P21*G61)+X(3,1)*trace(P31*G61)...
            +X(4,1)*trace(P41*G61)+X(5,1)*trace(P51*G61));        
        f(6,1) = 1/log(2) * log(1+2*eta(6,1)*z(6,1)-(intf61+sigma)*eta(6,1)^2);
        f(6,1) >= rate(6,1);
        {z(6,1),X(6,1),trace(G61*P61)} == rotated_lorentz(1);
    % t = 2 time slot
        G12 = cell2mat(G_cell(1,2));
        G22 = cell2mat(G_cell(2,2));
        G32 = cell2mat(G_cell(3,2));
        G42 = cell2mat(G_cell(4,2));
        G52 = cell2mat(G_cell(5,2));
        G62 = cell2mat(G_cell(6,2));
        intf12 = real(X(2,2)*trace(P22*G12)+X(3,2)*trace(P32*G12)...
            +X(4,2)*trace(P42*G12)+X(5,2)*trace(P52*G12)+X(6,2)*trace(P62*G12));       
        f(1,2) = 1/log(2) * log(1+2*eta(1,2)*z(1,2)-(intf12+sigma)*eta(1,2)^2);
        f(1,2) >= rate(1,2);
        {z(1,2),X(1,2),trace(G12*P12)} == rotated_lorentz(1);
        %
        intf22 = real(X(1,2)*trace(P12*G22)+X(3,2)*trace(P32*G22)...
            +X(4,2)*trace(P42*G22)+X(5,2)*trace(P52*G22)+X(6,2)*trace(P62*G22));
        f(2,2) = 1/log(2) * log(1+2*eta(2,2)*z(2,2)-(intf22+sigma)*eta(2,2)^2);
        f(2,2) >= rate(2,2);
        {z(2,2),X(2,2),trace(G22*P22)} == rotated_lorentz(1);
        %
        intf32 = real(X(1,2)*trace(P12*G32)+X(2,2)*trace(P22*G32)...
            +X(4,2)*trace(P42*G32)+X(5,2)*trace(P52*G32)+X(6,2)*trace(P62*G32));
        f(3,2) = 1/log(2) * log(1+2*eta(3,2)*z(3,2)-(intf32+sigma)*eta(3,2)^2);
        f(3,2) >= rate(3,2);
        {z(3,2),X(3,2),trace(G32*P32)} == rotated_lorentz(1);
        %
        intf42 = real(X(1,2)*trace(P12*G42)+X(2,2)*trace(P22*G42)+X(3,2)*trace(P32*G42)...
            +X(5,2)*trace(P52*G42)+X(6,2)*trace(P62*G42));
        f(4,2) = 1/log(2) * log(1+2*eta(4,2)*z(4,2)-(intf42+sigma)*eta(4,2)^2);
        f(4,2) >= rate(4,2);
        {z(4,2),X(4,2),trace(G42*P42)} == rotated_lorentz(1);
        %
        intf52 = real(X(1,2)*trace(P12*G52)+X(2,2)*trace(P22*G52)+X(3,2)*trace(P32*G52)...
            +X(4,2)*trace(P42*G52)+X(6,2)*trace(P62*G52));
        f(5,2) = 1/log(2) * log(1+2*eta(5,2)*z(5,2)-(intf52+sigma)*eta(5,2)^2);
        f(5,2) >= rate(5,2);
        {z(5,2),X(5,2),trace(G52*P52)} == rotated_lorentz(1);
        %
        intf62 = real(X(1,2)*trace(P12*G62)+X(2,2)*trace(P22*G62)+X(3,2)*trace(P32*G62)...
            +X(4,2)*trace(P42*G62)+X(5,2)*trace(P52*G62));
        f(6,2) = 1/log(2) * log(1+2*eta(6,2)*z(6,2)-(intf62+sigma)*eta(6,2)^2);
        f(6,2) >= rate(6,2);
        {z(6,2),X(6,2),trace(G62*P62)} == rotated_lorentz(1);
    % t = 3 time slot
        G13 = cell2mat(G_cell(1,3));
        G23 = cell2mat(G_cell(2,3));
        G33 = cell2mat(G_cell(3,3));
        G43 = cell2mat(G_cell(4,3));
        G53 = cell2mat(G_cell(5,3));
        G63 = cell2mat(G_cell(6,3));
        intf13 = real(X(2,3)*trace(P23*G13)+X(3,3)*trace(P33*G13)...
            +X(4,3)*trace(P43*G13)+X(5,3)*trace(P53*G13)+X(6,3)*trace(P63*G13));        
        f(1,3) = 1/log(2) * log(1+2*eta(1,3)*z(1,3)-(intf13+sigma)*eta(1,3)^2);
        f(1,3) >= rate(1,3);
        {z(1,3),X(1,3),trace(G13*P13)} == rotated_lorentz(1);
        %
        intf23 = real(X(1,3)*trace(P13*G23)+X(3,3)*trace(P33*G23)...
            +X(4,3)*trace(P43*G23)+X(5,3)*trace(P53*G23)+X(6,3)*trace(P63*G23));
        f(2,3) = 1/log(2) * log(1+2*eta(2,3)*z(2,3)-(intf23+sigma)*eta(2,3)^2);
        f(2,3) >= rate(2,3);
        {z(2,3),X(2,3),trace(G23*P23)} == rotated_lorentz(1);
        %
        intf33 = real(X(1,3)*trace(P13*G33)+X(2,3)*trace(P23*G33)...
            +X(4,3)*trace(P43*G33)+X(5,3)*trace(P53*G33)+X(6,3)*trace(P63*G33));
        f(3,3) = 1/log(2) * log(1+2*eta(3,3)*z(3,3)-(intf33+sigma)*eta(3,3)^2);
        f(3,3) >= rate(3,3);
        {z(3,3),X(3,3),trace(G33*P33)} == rotated_lorentz(1);
        %
        intf43 = real(X(1,3)*trace(P13*G43)+X(2,3)*trace(P23*G43)+X(3,3)*trace(P33*G43)...
            +X(5,3)*trace(P53*G43)+X(6,3)*trace(P63*G43));
        f(4,3) = 1/log(2) * log(1+2*eta(4,3)*z(4,3)-(intf43+sigma)*eta(4,3)^2);
        f(4,3) >= rate(4,3);
        {z(4,3),X(4,3),trace(G43*P43)} == rotated_lorentz(1);
        %
        intf53 = real(X(1,3)*trace(P13*G53)+X(2,3)*trace(P23*G53)+X(3,3)*trace(P33*G53)...
            +X(4,3)*trace(P43*G53)+X(6,3)*trace(P63*G53));
        f(5,3) = 1/log(2) * log(1+2*eta(5,3)*z(5,3)-(intf53+sigma)*eta(5,3)^2);
        f(5,3) >= rate(5,3);
        {z(5,3),X(5,3),trace(G53*P53)} == rotated_lorentz(1);
        %
        intf63 = real(X(1,3)*trace(P13*G63)+X(2,3)*trace(P23*G63)+X(3,3)*trace(P33*G63)...
            +X(4,3)*trace(P43*G63)+X(5,3)*trace(P53*G63));
        f(6,3) = 1/log(2) * log(1+2*eta(6,3)*z(6,3)-(intf63+sigma)*eta(6,3)^2);
        f(6,3) >= rate(6,3);
        {z(6,3),X(6,3),trace(G63*P63)} == rotated_lorentz(1);
    for t = 1:M
         for n = 1:N
             rate(n,t) >= 0;
         end
    end     
    for n = 1:N
        Dsum(n) == sum(rate(n,:));
        Dsum(n) >= r_thred;
    end
    Rsum == sum(rate(:));
    real(trace(P11)+trace(P21)+trace(P31)+trace(P41)+trace(P51)+trace(P61)) <= P_tot/M;
    real(trace(P12)+trace(P22)+trace(P32)+trace(P42)+trace(P52)+trace(P62)) <= P_tot/M;
    real(trace(P13)+trace(P23)+trace(P33)+trace(P43)+trace(P53)+trace(P63)) <= P_tot/M;
    P11 ==  hermitian_semidefinite(N);
    P21 ==  hermitian_semidefinite(N);
    P31 ==  hermitian_semidefinite(N);
    P41 ==  hermitian_semidefinite(N);
    P51 ==  hermitian_semidefinite(N);
    P61 ==  hermitian_semidefinite(N);
    P12 ==  hermitian_semidefinite(N);
    P22 ==  hermitian_semidefinite(N);
    P32 ==  hermitian_semidefinite(N);
    P42 ==  hermitian_semidefinite(N);
    P52 ==  hermitian_semidefinite(N);
    P62 ==  hermitian_semidefinite(N);
    P13 ==  hermitian_semidefinite(N);
    P23 ==  hermitian_semidefinite(N);
    P33 ==  hermitian_semidefinite(N);
    P43 ==  hermitian_semidefinite(N);
    P53 ==  hermitian_semidefinite(N);
    P63 ==  hermitian_semidefinite(N);    
    cvx_end
    if strcmp(cvx_status, 'Solved')||strcmp(cvx_status, 'Inaccurate/Solved')
       status = 1;
       P_cell{1,1} = P11;P_cell{2,1} = P21;P_cell{3,1} = P31;P_cell{4,1} = P41;P_cell{5,1} = P51;P_cell{6,1} = P61;
       P_cell{1,2} = P12;P_cell{2,2} = P22;P_cell{3,2} = P32;P_cell{4,2} = P42;P_cell{5,2} = P52;P_cell{6,2} = P62;
       P_cell{1,3} = P13;P_cell{2,3} = P23;P_cell{3,3} = P33;P_cell{4,3} = P43;P_cell{5,3} = P53;P_cell{6,3} = P63;
       Rsum_op = cvx_optval;
       R_op = rate;
    else
       status = 0;
       for t = 1:M
           for n = 1:N
               P_cell{n,t} = NaN;
           end
       end       
       Rsum_op = NaN;
       R_op = NaN;
    end