clc
clear

% number of repetitions
expr = 10;

% data setting                                                                                      
d = 100;n=1000;r=2;model='model';signal='small';maxIt = 500;

%step sizes
eta1_para= 0; eta2_para= 0.3; 

% thresholding parameter
thresh_idx = -1.4;

method = 'ATGD';

time = zeros(expr,1); err_S = time; err_L = time; err_omega = time;
for k = 1:expr
    filename = ['./data/',model,'_',signal,'_n',num2str(n),'_d',num2str(d),'_r',num2str(r),'_rep',num2str(k),'.mat'];
    load(filename)
    
    nu = max(eigs(omega_star));
    tol = 4*nu^2/(16*nu^4+1);
    
    linsp_thresh = linspace(0,1,10);
    eta1 = tol*2^(eta1_para); 
    eta2 = tol*2^(eta2_para);
    thresh_ratio = exp(thresh_idx);
    s = round(thresh_ratio*density_real *d^2);    
    stoptol = 1e-7;
    
    [S0,Z0] = initialization(hsigma,s,r);
    timein = cputime;
    out_ATGD = ATGD(S_star,L_star,hsigma,S0,Z0,maxIt,eta1,eta2,s,stoptol);
    timeout = cputime - timein;
    time(k) = timeout;
    
    % error plot
    err_S_all = out_ATGD.err_S;            err_L_all = out_ATGD.err_L;
    err_omega_all = out_ATGD.err_omega;    objective = out_ATGD.obj;
    S_hat = out_ATGD.S;                    L_hat = out_ATGD.L;
    
    iter = length(err_S_all);
    
    err_S(k) = err_S_all(iter);            err_L(k) = err_L_all(iter);
    err_omega(k) = err_omega_all(iter);    
end

err_S_mean = mean(err_S);          err_S_sd = std(err_S);
err_L_mean = mean(err_L);          err_L_sd = std(err_L);
err_omega_mean = mean(err_omega);  err_omega_sd = std(err_omega);
time_mean = mean(time);

fprintf(['Average results on %2.0f experiments\n' ...
         'Error of S: %10.4f' char(177) '%0.4f\n' ...
         'Error of L: %10.4f' char(177) '%0.4f\n' ...
         'Error of Precision Matrix: %10.4f' char(177) '%0.4f\n' ...
         'Average Time:  %10.4f\n'],...
    [expr,err_S_mean,err_S_sd,err_L_mean,err_L_sd,err_omega_mean,err_omega_sd,time_mean]);
