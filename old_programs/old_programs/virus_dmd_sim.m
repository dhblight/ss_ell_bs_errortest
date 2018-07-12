% simulated data ellipsoid -> sphere

% constants
r0 = 30e-9;
r1 = 10e-9;
r2 = 20e-9;
tau = 3.7;
n = 300;
m = 200;

% creating the data matrix
qi = logspace(-3,0,n);
t = logspace(-1,2,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';


X = zeros(n,m);
for ii = 1:n
    for jj = 1:m
        X(ii,jj) = elipse_p_q(qgrid(ii,jj),r_t(r1,r0,T(ii,jj),tau),r_t(r2,r0,T(ii,jj),tau));
    end
end

% % number of singular values to keep
% [U, S, V] = svd(X, 'econ');
% coeff = 1.6089; % Gavish & Donoho for beta = 0.1
% threshold = median(diag(S)) * coeff;
% r = sum(diag(S) >= threshold);

% DMD
imode = 1;
[w,e,b] = optdmd(X,t,5,imode);
Xdmd = w * diag(b) * exp(e*t);
X1 = w(:,1) .* b(1) .* exp(e(1)*t);
X2 = w(:,2) .* b(2) .* exp(e(2)*t);
X3 = w(:,3) .* b(3) .* exp(e(3)*t);
X4 = w(:,4) .* b(4) .* exp(e(4)*t);
X5 = w(:,5) .* b(5) .* exp(e(5)*t);

subplot(2,2,1)
surf(real(X1))

subplot(2,2,2)
surf(real(X2))

subplot(2,2,3)
surf(real(X))

subplot(2,2,4)
surf(real(Xdmd))

% pi = zeros(1,300);
% compt = 1;
% for ii = qi
%     pi(compt) = elipse_p_q(ii,r1,r2);
%     compt = compt + 1;
% end
% 
% % pi = elipse_p_q(qi,r1,r2);
% 
% loglog(qi,pi)
