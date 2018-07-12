% simulated data for transformation of ellipsoid -> sphere

rng('shuffle');

%% constants
r0 = 30e-9; % final radius
r1 = 10e-9; % initial semiminor axis
r2 = 20e-9; % initial semimajor axis
tau1 = 2; % relaxation time for r1
tau2 = 2; % relaxation time for r2
n = 300; % # of spatial leasurements
m = 30; % # of temporal measurements

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
        Xclean(ii,jj) = elipse_p_q(qgrid(ii,jj),r_t(r1,r0,T(ii,jj),tau1),r_t(r2,r0,T(ii,jj),tau2))* qgrid(ii,jj)^2;
    end
end

% % multipy by q^2? !! remember to adjust sigma accordingly !! 
% Xclean = Xclean .* qgrid .^2;

% add noise
sigma = 0.03;
G = sigma .* randn(n,m);
X = Xclean .* (1+G);

% % no noise
% X = Xclean;

%% number of singular values to keep
[U, S, V] = svd(X, 'econ');

% % with noise
% coeff = optimal_SVHT_coef(m/n,1); % Gavish & Donoho 
% x = diag(S);
% threshold = sqrt(n) * sigma * coeff;
% r = sum(x >= threshold);

% without noise
coeff = optimal_SVHT_coef(m/n,0); % Gavish & Donoho 
x = diag(S);
threshold = median(x) * coeff;
r = sum(x >= threshold);

%% DMD
% normal
imode = 1;
[w,e,b] = optdmd(X,t,5,imode); % w: dmd modes, e:eigenvalues, b:coefficients
Xdmd = w * diag(b) * exp(e*t); % reconstruction of data using dmd modes

% 
% % projected
% maxiter = 30;
% tol = sigma/10;
% opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',1e-9);
% imode = 2; % projected version
% [wproj,eproj,bproj,atilde] = optdmd(X,t,r,imode,opts,[],U);
% Xdmdproj = wproj * diag(bproj) * exp(eproj*t);
% pod1 = wproj(:,1) .* bproj(1) .* exp(eproj(1)*t); 

%% plots of data
%subplot(2,2,1)
surf(real(X))
set(gca,'YScale','log');
set(gca,'ZScale','log');

% subplot(2,2,2)
surf(real(Xclean))
set(gca,'YScale','log');
set(gca,'ZScale','log');
% 
% subplot(2,2,3)
% surf(real(Xdmd))
 set(gca,'YScale','log');
 set(gca,'ZScale','log');
% 
% subplot(2,2,4)
% surf(real(Xdmdproj))
% set(gca,'YScale','log');
% set(gca,'ZScale','log');
