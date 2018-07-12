% creating the data matrix
xi = linspace(-10,10,400);
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
[Xgrid,T] = meshgrid(xi,t);

f1 = sech(Xgrid+3) .* (1*exp(1j*2.3*T));
f2 = (sech(Xgrid) .* tanh(Xgrid)) .* (2*exp(1j*2.8*T));
fclean = f1 + f2;
Xclean = fclean.';
sigma = 0.01;
g = randn(200,400);
f = f1 + f2 + sigma.* g;
X = f.'; 

lambda = (optimal_SVHT_coef(0.5,1) * sqrt(200) * sigma);

[U,S,V] = svd(X, 'econ');
x = diag(S); 
x( x < lambda ) = 0; 
num = sum(x >= lambda);
Xhat = U * diag(x) * V';

r = 2;
imode = 1;
[w,e,b] = optdmd(X,t,r,imode);
Xdmd = w*diag(b)*exp(e*t);

subplot(2,1,1)
surfl(real(X'));
shading interp; view(-20,60);

subplot(2,1,2)
surfl(real(Xhat'));
shading interp; view(-20,60);

errdata = norm(X-Xclean,'fro')/norm(Xclean,'fro')
errdenoised = norm(Xhat-Xclean,'fro')/norm(Xclean,'fro')
errdmd = norm(Xdmd-Xclean,'fro')/norm(Xclean,'fro')