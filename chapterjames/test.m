    cvx_quiet(0)
%    cvx_begin SDP quiet
    cvx_begin SDP
    cvx_solver Mosek
%     cvx_solver SeDuMi
%    cvx_solver sdpt3 
    cvx_precision best
    variable x(4,1) binary
%     variable x(N*M,1)
    variable W
    maximize W
    subject to
W == x'*[1,-2,3,-4]';

    cvx_end
