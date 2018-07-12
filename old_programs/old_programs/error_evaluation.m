% simulated data for transformation of ellipsoid -> sphere

rng('shuffle');

%% objects
r0 = 30e-9; % final radius
r1 = 10e-9; % initial semiminor axis
r2 = 20e-9; % initial semimajor axis
start1 = 1; end1 = 3; num1 = 3; % values over which to vary tau1
start2 = 1; end2 = 3; num2 = 3; % values over which to vary tau2 
tau1vals = linspace(start1,end1,num1); % relaxation time for r1
tau2vals = linspace(start2,end2,num2); % relaxation time for r2
n = 300; % # of spatial leasurements
m = 50; % # of temporal measurements
err1_var_tau = zeros(length(tau1vals),length(tau2vals));
err2_var_tau = zeros(length(tau1vals),length(tau2vals));
errtot_var_tau = zeros(length(tau1vals),length(tau2vals));
errel1_var_tau = zeros(length(tau1vals),length(tau2vals));
errel2_var_tau = zeros(length(tau1vals),length(tau2vals));
erreltot_var_tau = zeros(length(tau1vals),length(tau2vals));

% create the grid
qi = logspace(8,10,n);
t = logspace(-1,3,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';

%% varying tau
count1 = 1;

for tau1_iter = tau1vals 
    count2 = 1;
        
    for tau2_iter = tau2vals

        % clean data
        Xclean = zeros(n,m);
        for ii = 1:n
            for jj = 1:m
                Xclean(ii,jj) = elipse_p_q(qgrid(ii,jj),r_t(r1,r0,T(ii,jj),tau1_iter),r_t(r2,r0,T(ii,jj),tau2_iter));
            end
        end
        
        % add noise
        sigma = 5e-5;
        G = sigma .* randn(n,m);
        X = Xclean+G;
        
        % number of singular values to keep
        [U, S, V] = svd(X, 'econ');
        
        % with noise
        coeff = optimal_SVHT_coef(m/n,1); % Gavish & Donoho
        x = diag(S);
        threshold = sqrt(n) * sigma * coeff;
        r = sum(x >= threshold);
        
        % DMD
        imode = 1;
        [w,e,b] = optdmd(X,t,r,imode); % w: dmd modes, e:eigenvalues, b:coefficients
        
        % only care about real eigenvalues
        e(abs(imag(e))>0.05) = 0;
        e = real(e);       
       
        eig_tau = -1./ e; % tau values corresponding to the eigenvalues
        
        [err1,ind1] = min(abs(eig_tau - tau1_iter));
        [err2,ind2] = min(abs(eig_tau - tau2_iter));
        
        % Can't use the same eigenvalue for both unless tau1 = tau2
        if tau1_iter ~= tau2_iter            
            if ind1 == ind2
                eig_tau(ind1) = inf; % make sure this value isnèt closest to tau so the next closest is used 
                if err1 > err2 % 
                    [err1,ind1] = min(abs(eig_tau - tau1_iter));
                else
                    [err2,ind2] = min(abs(eig_tau - tau2_iter));
                end
            end
        end
                
        errtot = sqrt(err1^2 + err2^2); 
        
        err1_var_tau(count1,count2) = err1;
        err2_var_tau(count1,count2) = err2;
        errtot_var_tau(count1,count2) = errtot;
        errel1_var_tau(count1,count2) = err1/tau1_iter;
        errel2_var_tau(count1,count2) = err2/tau2_iter;
        erreltot_var_tau(count1,count2) = sqrt((err1/tau1_iter)^2 + (err2/tau2_iter)^2);
        
        fprintf('iterations completed varying tau: %d / d\n',(20*(count1-1)),count2);
        count2 = count2 + 1;   
        
    end
    
    count1  = count1 + 1;

end

% meshgrid(X,Y), tau2 changes as you go along each row -> x, 
% tau1 changes as you go down each column -> y
% element (i,j) of error matrix corresponds to the error for tau1 =
% TAU1(i,j) and tau2 = TAU2(i,j)
[TAU2,TAU1] = meshgrid(tau2vals,tau1vals); 

% for a 2d plot of tau1/tau2 vs error
TAU1_TAU2 = TAU1./TAU2;
tau1_tau2 = reshape(TAU1_TAU2,1,[]);
errtot_var_tau_vec = reshape(errtot_var_tau,1,[]);
erreltot_var_tau_vec = reshape(erreltot_var_tau,1,[]);