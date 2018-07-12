%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script which calculates total scattered intensity for 2 populations of 
% small spheres and big spheres, with n small spheres forming 1 big sphere,
% their concentrations of the form exp(-t/tau) and 1 - exp(-t/tau), 
% and adds multiplicative noise. 
% Then uses a modified optdmd algorithm to find the dynamic modes and a calculated
% relaxation time (also a denoised system reconstruction if wanted).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');

%% constants, variables, objects
numtrials = 100; %number of times we add noise and calculate the dmd 
tempvec = zeros(1,numtrials); % vector to store result of each trial
evalsvec = zeros(1,numtrials);

nell = 10; % number of small ellipsoids to make 1 big ellipse
% 1: big ellipsoid, 0: small ellipsoid
ratio0 = 1.5;
ratio1 = 1.5;
rb0 = 10e-9;  
ra0 = rb0 * ratio0;
rb1 = rb0 * (nell*(ratio0/ratio1))^(1/3); 
ra1 = rb1 * ratio1;

M0 = 180000;
M1 = M0 * nell;

c0 = 1e-9; % initial concentration (mol/cm^3) 

tau = 10;

r = 3; % rank

n = 300; % # of spatial measurements
m = 80; % # of temporal measurements

sigma = 0.03; % multiplicative error

% the grid of co-ordinates for the data
qi = logspace(7.7,9.1761,n);
t = logspace(-1,3,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% main program

% clean data
Xclean = zeros(n,m);
for ii = 1:n
    for jj = 1:m
        P0 = elipse_p_q(qgrid(ii,jj),ra0,rb0);
        P1 = elipse_p_q(qgrid(ii,jj),ra1,rb1);
        [C0,C1] = c_t(c0,T(ii,jj),tau,nell);
        Xclean(ii,jj) = M0^2 * P0* C0 + M1^2 * P1 * C1;
    end
end

% lpoop    
for kk = 1:numtrials
    
    % add noise
    G = sigma .* randn(n,m);
    X = Xclean.*(1+G);
    
%     % decide rank
%     [U, S, V] = svd(X, 'econ');
%     semilogy(diag(S)); % put breakpoint
    
    % DMD
    imode = 1;
    [w,e,b,final_error,evals] = multioptdmd3(X,t,r,0.035);
    e(abs(imag(e))>0.1) = 0;
    e = real(e);
    eig_tau = -1./ e;
    [~,ind] = min(abs(eig_tau - tau));
    tautemp = eig_tau(ind);
    tempvec(kk) = tautemp;
    evalsvec(kk) = evals;
    
    fprintf('iterations completed: %d / %d\n',kk,numtrials);
end

hold on
scatter(1:numtrials,tempvec)
plot(1:numtrials,tau.*ones(1,numtrials),'k','LineWidth',2)

tau_calc = mean(tempvec,'omitnan')
STD = std(tempvec,'omitnan');
error = STD/sqrt(numtrials);
