% constants
r0 = 4;
r1 = 1;
tau = 4;
n = 300;
m = 30;

% creating the data matrix
qi = logspace(0,3,n);
t = logspace(-1,2,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';


X = zeros(n,m);
for ii = 1:n
    for jj = 1:m
        X(ii,jj) = sphere_p_q(qgrid(ii,jj),r_t(r1,r0,T(ii,jj),tau));
    end
end
[U,S,V] = svd(X,'econ');

% DMD
imode = 1;
[w,e,b] = optdmd(X,t,8,imode);
Xdmd = w * diag(b) * exp(e*t);
X1 = w(:,1) .* b(1) .* exp(e(1)*t);
X2 = w(:,2) .* b(2) .* exp(e(2)*t);

subplot(1,2,1)
surf(real(X))
%set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'ZScale','log');

subplot(1,2,2)
surf(real(Xdmd))
%set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'ZScale','log');

