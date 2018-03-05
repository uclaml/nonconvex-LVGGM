function out_ATGD = ATGD(S_star,L_star,sigma,S0,Z0,maxIt,eta1,eta2,s,stoptol)
err_S = zeros(maxIt,1);  err_L = err_S;  err_omega = err_S;  obj = err_S;

err_S(1) = norm(S0-S_star,'fro');
obj(1) = lossF(sigma,S0,Z0*Z0');

incoh_L = zeros(maxIt,1);

S = S0;
err_L(1) = norm(Z0*Z0'-L_star,'fro');
Z= Z0;
L=Z*Z';

err_L_min =err_L(1);
err_S_min =err_S(1);
for t = 2:maxIt
    S = S - eta1 * nablaF(sigma,S,L);%/(t^(1/5));
    S = HardThresh(S,s);
    
    Z = Z - eta2*( 2*nablaF(sigma,S,L)*Z );
    L = Z*Z';
    
    err_S(t) = norm(S-S_star,'fro');  err_L(t) = norm(L-L_star,'fro');
    err_omega(t) =norm(L+S-S_star-L_star,'fro');
    obj(t) = lossF(sigma,S,L);
    if err_L(t)-err_L_min >stoptol     
        break;      
    end
    err_S_min = err_S(t);       err_L_min = err_L(t);
    if abs(err_L(t)-err_L(t-1))<stoptol %successive iterations
        break
    end
end

if t == maxIt
    fprintf('%s\n', 'Maximum number of iteration reached, ATGD may not converge.');
end

out_ATGD.S =S;                 out_ATGD.Z = Z;                   out_ATGD.obj = obj(1:t);
out_ATGD.err_S = err_S(1:t);   out_ATGD.err_L=err_L(1:t);        out_ATGD.incoh = incoh_L(1:t);
out_ATGD.L = Z*Z';             out_ATGD.err_omega=err_omega(1:t);

end

function [G] = nablaF(sigma,S,L)
d = size(sigma,1);
omega = S+L;

G = sigma - omega\eye(d);
end

function [V] = lossF(sigma,S,L)
d = size(sigma,1);

omega = S + L;
V = sum(sum(sigma.*omega)) - log(det(omega));
end

function mu = incoh(X)
mu = max(sqrt(sum(X.^2)));
end