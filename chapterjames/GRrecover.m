function [p0] = GRrecover(P,N)
L = 100; % number of Gaussian randomizations
G = sqrt(2) / 2 * (randn(N^2, N^2) + 1j * randn(N^2, N^2));
hr = sqrt(2) / 2 * (randn(N^2, 1) + 1j * randn(N^2, 1));
phi = diag(hr') * G;

max_F = 0;
max_v = 0;
[U, Sigma] = eig(P);
for l = 1 : L
    r = sqrt(2) / 2 * (randn(N^2,1) + 1j * randn(N^2,1));
    v = U * Sigma^(0.5) * r;
    v = exp(1j * angle(v / v(end)));
    v = v(1 : N^2);
    if v' * phi * phi' * v > max_F
        max_v = v;
        max_F = v' * phi * phi' * v;
    end
end
% max_v' * phi * phi' * max_v;
p0 = max_v;