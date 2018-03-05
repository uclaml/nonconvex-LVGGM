function [S0,Z0] = initialization(sigma,s,r)
sigma = adjust_sigma(sigma);
omega = eye(size(sigma,1))/sigma;
S0 = HardThresh(omega,s);
[U,D] = mexeig(omega-S0);
D_sq = sqrtM(full(D));
Z0 = U*D_sq(:,1:r);
end

