% simulated data for transformation of ellipsoid -> sphere

rng('shuffle');

%% constants
r0 = 30e-9; % final radius
r1 = 10e-9; % initial semiminor axis
r2 = 20e-9; % initial semimajor axis
tau1 = 4; % relaxation time for r1
tau2 = 10; % relaxation time for r2
n = 300; % # of spatial leasurements
m = 30; % # of temporal measurements
rvals = 1:20;
numtrials = 100;
tau1_mat = zeros(numtrials,length(rvals));
tau2_mat = zeros(numtrials,length(rvals));


%% creating the data matrix

% grid
qi = logspace(7,9.1761,n);
t = logspace(-1,2,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';

% clean data
Xclean = zeros(n,m);
for ii = 1:n
    for jj = 1:m
        Xclean(ii,jj) = elipse_p_q(qgrid(ii,jj),r_t(r1,r0,T(ii,jj),tau1),r_t(r2,r0,T(ii,jj),tau2));
    end
end

%% loop 
count1 = 1;
for ii = 1:numtrials
    
    count2 = 1;
       
    % add noise proportional to intensity
    sigma = 0.03;
    G = sigma .* randn(n,m);
    X = Xclean .* (1+G);
    
    for r = rvals
        
        % DMD
        imode = 1;
       
        % sometimes gives errors "Input to SVD must not contain NaN or Inf." 
        % catch these errors so the rest of the program still executes
        try
            [w,e,b] = optdmd(X,t,r,imode); % w: dmd modes, e:eigenvalues, b:coefficients
        catch
            warning('Problem using optdmd.  Assigning value of 0.');
            e = zeros(r,1); % this will give an error in tau of Inf later on
            b = zeros(r,1);
            w = zeros(n,r);
        end
        
        % error in reconstruction
        Xdmd = w * diag(b) * exp(e*t);              
        
        % error in eigenvalues
        % only care about real eigenvalues
        e(abs(imag(e))>0.5) = 0;
        e = real(e);
        
        eig_tau = -1./ e;
        eig_tau_original = eig_tau;
        
        [erre1,ind1] = min(abs(eig_tau - tau1));
        [erre2,ind2] = min(abs(eig_tau - tau2));
        
        % Can't use the same eigenvalue for both unless tau1 = tau2
        if tau1 ~= tau2            
            if ind1 == ind2
                eig_tau(ind1) = inf; % make sure this value isn't closest to tau so the next closest is used 
                if erre1 > erre2 % 
                    [erre1,ind1] = min(abs(eig_tau - tau1));
                else
                    [erre2,ind2] = min(abs(eig_tau - tau2));
                end
            end
        end
        
        tau1_calc = eig_tau_original(ind1);
        tau2_calc = eig_tau_original(ind2);
        
        tau1_mat(count1,count2) = tau1_calc;
        tau2_mat(count1,count2) = tau2_calc;
        
        fprintf('iterations completed: %d / %d\n',(length(rvals)*(count1-1))+count2,numtrials*length(rvals));
        count2 = count2+1;
        
    end
    
    
    count1 = count1+1;
end

% vector containing the mean error across all the trials for each rank
meanerrerel = zeros(1,length(rvals));

tau1_mat(tau1_mat>10*tau1) = 0;
tau2_mat(tau2_mat>10*tau2) = 0;
tau1_mat(tau1_mat<0) = 0;
tau2_mat(tau2_mat<0) = 0;

for jj = 1:length(rvals) 
    col1 = tau1_mat(:,jj);
    col2 = tau2_mat(:,jj);
    tau1mean = sum(col1) ./ sum(col1~=0); % calculate mean not taking into account the zero elements
    tau2mean = sum(col2) ./ sum(col2~=0);
    meanerrerel(jj) = sqrt(((tau1mean-tau1)/tau1)^2 + ((tau2mean-tau2)/tau2)^2);
end

%% plot
plot(rvals,meanerrerel);
set(gca,'YScale','log');
