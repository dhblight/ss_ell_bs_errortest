%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script which calculates total scattered intensity for a system of small
% spheres -> big ellipsoids -> big spheres. Each transformation is considered 
% instantaneous and happens with a certain probability -> the concentrations 
% of each species are the solutions to the bateman equation and are given 
% by c_t_3species.m
%
% Uses a modified optdmd algorithm to find the dynamic modes and calculated
% relaxation times for the 2 processes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% seed random number generation based on time
rng('shuffle');

%% constants, variables, objects
numtrials = 50; %number of times we add noise and calculate the dmd 
tempmat = zeros(2,numtrials); % vector to store result of each trial
evalsvec = zeros(1,numtrials);
errvec = zeros(1,numtrials);

nsph = 10; % number of small spheres to make big ellipsoid
r0 = 10e-9; % radius of small spheres
r1 = r0 * nsph^(1/3); % radius of final big spheres
ratio = 3; % ra/rb for the ellipse, ra > rb for prolate (needle-like) ellipsoid 
rb = r1 * ratio^(-1/3);
ra = rb * ratio;

Mss = 180000; % molar mass of small spheres (g/mol)
M = Mss * nsph;

c0 = 1e-9; % initial concentration (mol/cm^3) 

tau1 = 1; % relaxation time for small spheres -> ellipsoid
tau2 = 200; % relaxation time for ellipsoid -> big sphere

r = 6; % rank

n = 300; % # of spatial measurements
m = 50; % # of temporal measurements

sigma = 0.03; % multiplicative error

% the grid of co-ordinates for the data
qi = logspace(7.7,9.1761,n);
t = logspace(-1,3,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% main program

% generate clean data
        
% form factors
Pss = sphere_p_q(qgrid,r0);
Pell = elipse_p_q(qgrid,ra,rb);
Pbs = sphere_p_q(qgrid,r1);

% concentrations
[Css,Cell,Cbs] = c_t_3species(c0,T,tau1,tau2,nsph);

% intentisty
Xclean = Mss^2 .* Pss .* Css + Mbs^2 .* Pbs .* Cbs + ...
    Mell^2 .* Pell .* Cell;


% loop    
for kk = 1:numtrials
    
    % add noise
    G = sigma .* randn(n,m);
    X = Xclean.*(1+G);
    
%     % decide rank, only need the first time
%     [U, S, V] = svd(X, 'econ');
%     semilogy(diag(S)); % put breakpoint
    
    % DMD
    imode = 1; % don't project to POD modes
    [w,e,b,final_error,evals] = multioptdmd3(X,t,r,0.05);
    e(abs(imag(e))>0.1) = 0; % only care about real eigenvalues
    e = real(e);
    eig_tau = -1./ e; % tau values corresponding to the eigenvalues
    eig_tau_original = eig_tau; % eig_tau will be modified but we still want the values
    
    % find indices of eigenvalues and their errors
    [err1,ind1] = min(abs(eig_tau - tau1));
    [err2,ind2] = min(abs(eig_tau - tau2));
    
    % Can't use the same eigenvalue for both unless tau1 = tau2
    if tau1 ~= tau2
        if ind1 == ind2
            eig_tau(ind1) = inf; % make sure this value isn't closest to tau so the next closest is used
            if err1 > err2 % calculate the next closest value
                [err1,ind1] = min(abs(eig_tau - tau1)); % for tau1 if it was closer to tau2 originally
            else
                [err2,ind2] = min(abs(eig_tau - tau2)); % or vice versa
            end
        end
    end
    
    % calculated tau value
    tau1_temp = eig_tau_original(ind1);
    tau2_temp = eig_tau_original(ind2);
    
    % store info
    tempmat(1,kk) = tau1_temp;
    tempmat(2,kk) = tau2_temp;
    evalsvec(kk) = evals;
    errvec(kk) = final_error;
    
    fprintf('iterations completed: %d / %d\n',kk,numtrials);
end

% hold on
% scatter(1:numtrials,tempvec)
% plot(1:numtrials,tau.*ones(1,numtrials),'k','LineWidth',2)

tau1_calc = mean(tempmat(1,:),'omitnan')
tau2_calc = mean(tempmat(2,:),'omitnan')
% STD = std(tempvec,'omitnan');
% error = STD/sqrt(numtrials)

% % saves data as "prefix_dd-mmm-yyyy_HH-MM.mat"
formatout = 'dd-mmm-yyyy_HH-MM';
save(strcat('ss_ell_bs','_',datestr(datetime('now'),formatout),'.mat'));
