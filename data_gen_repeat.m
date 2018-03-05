%generate data for LVGGM
%model: (X_O,X_L)~N(0,Sigma), omega^* = (Sigma_OO)^{-1}
clc
clear

%%%dimension--setting------------------------------------------------------
d = 100;n = 1000;r = 2;density = 0.05;signal = 'small';model='model';
%--------------------------------------------------------------------------

%reciprocal conditional number
rc = 0.5;

% number of repetions
exp_num = 10;
for k =1:exp_num   
    d_all = d+r;
    omega_tilde = sprandsym(d_all,density,rc,1);
    S_star = full(omega_tilde(1:d,1:d));
    L_star = - omega_tilde(1:d,(d+1):d_all)/(omega_tilde((d+1):d_all,(d+1):d_all))*omega_tilde((d+1):d_all,1:d);
    L_star = full(L_star);
    rank(L_star)
    L_star = L_star/d;
    
    omega_star = S_star + L_star;
    sigma_star = eye(d)/omega_star;
    sigma_tilde = eye(d_all)/omega_tilde;
    
    %check model
    norm(sigma_tilde(1:d,1:d)-sigma_star,'fro')
       
    density_real = sum(sum(S_star~=0))/(d^2);
    X = mvnrnd(zeros(n,d),sigma_star);
    hsigma = X'*X/n;
    filename = ['./data/',model,'_',signal,'_','n',num2str(n),'_d',num2str(d),'_r',num2str(r),'_rep',num2str(k),'.mat'];
    save(filename,'hsigma','sigma_star','omega_star','S_star','L_star','density_real','r');
end


