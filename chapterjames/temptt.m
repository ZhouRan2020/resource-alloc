    % t = 1 time slot
        G11 = cell2mat(G_cell(1,1));
        f(1,1) = 1/log(2) * log(1+v(1,1))-v(1,1)+2*eta(1,1)*z1(1)-(X(1,1)*trace(P11*G11)+sigma)*eta(1,1)^2;
        f(n,1) >= rate(n,1);
        {z1(1),(1+v(1,1))*X(1,t),trace(G11*P11)} == rotated_lorentz(1);
        %
        G21 = cell2mat(G_cell(2,1));
        f(1,1) = 1/log(2) * log(1+v(1,1))-v(1,1)+2*eta(1,1)*z1(1)-(X(1,1)*trace(P21*G21)+sigma)*eta(1,1)^2;
        f(n,1) >= rate(n,1);
        {z1(1),(1+v(1,1))*X(2,t),trace(G21*P21)} == rotated_lorentz(1);
        %
        G11 = cell2mat(G_cell(1,1));
        f(1,1) = 1/log(2) * log(1+v(1,1))-v(1,1)+2*eta(1,1)*z1(1)-(X(1,1)*trace(P11*G11)+sigma)*eta(1,1)^2;
        f(n,1) >= rate(n,1);
        {z1(1),(1+v(1,1))*X(1,t),trace(G11*P11)} == rotated_lorentz(1);
        %
        G11 = cell2mat(G_cell(1,1));
        f(1,1) = 1/log(2) * log(1+v(1,1))-v(1,1)+2*eta(1,1)*z1(1)-(X(1,1)*trace(P11*G11)+sigma)*eta(1,1)^2;
        f(n,1) >= rate(n,1);
        {z1(1),(1+v(1,1))*X(1,t),trace(G11*P11)} == rotated_lorentz(1);
        %
        G11 = cell2mat(G_cell(1,1));
        f(1,1) = 1/log(2) * log(1+v(1,1))-v(1,1)+2*eta(1,1)*z1(1)-(X(1,1)*trace(P11*G11)+sigma)*eta(1,1)^2;
        f(n,1) >= rate(n,1);
        {z1(1),(1+v(1,1))*X(1,t),trace(G11*P11)} == rotated_lorentz(1);
        %
        G11 = cell2mat(G_cell(1,1));
        f(1,1) = 1/log(2) * log(1+v(1,1))-v(1,1)+2*eta(1,1)*z1(1)-(X(1,1)*trace(P11*G11)+sigma)*eta(1,1)^2;
        f(n,1) >= rate(n,1);
        {z1(1),(1+v(1,1))*X(1,t),trace(G11*P11)} == rotated_lorentz(1);