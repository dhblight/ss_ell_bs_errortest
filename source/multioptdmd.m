function [w,e,b,err] = multioptdmd2(X,t,r,imode,num_vals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification to optdmd.m allowing for multiple starting values of alpha,
% decreasing the likelihood of an erroneous result caused by the
% initialisation algoithm finding a local minimizer of eq(47) in the 
% optimized DMD article instead of the global minimizer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');

[rows,cols] = size(X); 
errvals = zeros(30,num_vals); % 30 is max iterations in varpro opts
evals = zeros(r,num_vals);
wvals = zeros(rows,r,num_vals);
bvals = zeros(r,num_vals);
lasterr = zeros(1,num_vals);

% initial guess at alpha_init which we will perturb for e
[u,~,~] = svd(X,'econ');
ux1 = u'*X;
ux2 = ux1(:,2:end);
ux1 = ux1(:,1:end-1);
t1 = t(1:end-1);
t2 = t(2:end);
dx = (ux2-ux1)*diag(1./(t2-t1));
xin = (ux1+ux2)/2;
[u1,s1,v1] = svd(xin,'econ');
u1 = u1(:,1:r);
v1 = v1(:,1:r);
s1 = s1(1:r,1:r);
atilde = u1'*dx*v1/s1;
alpha_init0 = eig(atilde);
clear ux1 ux2 atilde t1 t2 dx xin

% fit to all of data
m = length(t);
[is,~] = size(X);
ia = r;
n = r;

factors = logspace(-2,0,num_vals-1);

for ii = 1:num_vals
    
    if ii == 1
        alpha_init = alpha_init0;
    else
        alpha_init = alpha_init0 .* factors(ii-1) .* randn(size(alpha_init0));
    end
    [wtemp,etemp,niter,err,imode,alphas] = varpro2(transpose(X),t, ...
        @varpro2expfun,@varpro2dexpfun,m,n,is,ia,alpha_init,varpro_opts_edit());
    wtemp = transpose(wtemp);        
    btemp = sqrt(sum(abs(wtemp).^2,1))'; % normalize
    wvals(:,:,ii) = wtemp*diag(1./btemp);
    evals(:,ii) = etemp;
    bvals(:,ii) = btemp;
    errvals(:,ii) = err;
    if errvals(:,ii) ~= zeros(30,1)
        lasterr(ii) = errvals(30,ii);
    else
        lasterr(ii) = inf;
    end
end

[~,ind] = min(lasterr);
w = wvals(:,:,ind);
e = evals(:,ind);
b = bvals(:,ind);
err = errvals(:,ind);

    

