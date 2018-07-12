% Practise using optdmd.m
% Using the example from section 1.4 of "Dynamic Mode Decomposition" 

% creating the data matrix
xi = linspace(-10,10,400);
t = linspace(0,4*pi,200);
dt = t(2) - t(1);
[Xgrid,T] = meshgrid(xi,t);

f1 = sech(Xgrid+3) .* (1*exp(1j*2.3*T));
f2 = (sech(Xgrid) .* tanh(Xgrid)) .* (2*exp(1j*2.8*T));
f = f1 + f2;
X = f.'; % data matrix

% [U,S,V] = svd(X, 'econ');
% can plot a semilog graph of diag(S) to check the required rank


% rank & mode
r = 2; imode = 1;

% evaluation of the optimised dmd
[w,e,b] = optdmd(X,t,r,imode);

% reconstructed values
X1 = w(:,1) .* b(1) .* exp(e(1)*t);
X2 = w(:,2) .* b(2) .* exp(e(2)*t);
Xdmd = w * diag(b) * exp(e*t);
%[Xdmd2,err] = best_reconstruction(X,e,t); % using function of Askham &
%Kutz,

% plots
colormap winter

subplot(2,4,[1,2])
surfl(real(X1'));
shading interp; view(-20,60);
set(gca, 'YTick', numel(t)/4 * (0:4));
set(gca, 'Yticklabel', {'0', '\pi', '2\pi', '3\pi', '4\pi'});
set(gca, 'XTick', linspace(1, numel(xi),3));
set(gca, 'Xticklabel', {'-10', '0', '10'});
title('$f_{1}$','Interpreter','latex');

subplot(2,4,[3,4])
surfl(real(X2'));
shading interp; view(-20,60);
set(gca, 'YTick', numel(t)/4 * (0:4));
set(gca, 'Yticklabel', {'0', '\pi', '2\pi', '3\pi', '4\pi'});
set(gca, 'XTick', linspace(1, numel(xi),3));
set(gca, 'Xticklabel', {'-10', '0', '10'});
title('$f_{2}$','Interpreter','latex');

subplot(2,4,[6,7])
surfl(real(Xdmd'));
shading interp; view(-20,60);
set(gca, 'YTick', numel(t)/4 * (0:4));
set(gca, 'Yticklabel', {'0', '\pi', '2\pi', '3\pi', '4\pi'});
set(gca, 'XTick', linspace(1, numel(xi),3));
set(gca, 'Xticklabel', {'-10', '0', '10'});
title('$f_{1}+f_{2}$','Interpreter','latex');
