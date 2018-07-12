function [w,e,b,final_error,evals] = multioptdmd3(X,t,r,maxerr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification to optdmd.m allowing for multiple starting values of alpha,
% decreasing the likelihood of an erroneous result caused by the
% initialisation algoithm finding a local minimizer of eq(47) in the 
% optimized DMD article instead of the global minimizer
%
% Input:
%
% X: Data Matrix
% t: time series
% r: rank to use 
% maxerr: maximum error given b varpro 2 which we accept (if unsure run
% multiple trials with a high maxerr and store the final_error values, then
% choose a sensible number based off these values (~ %error?? seems to work)
% tol: max difference between 2 successive values which are less than maxerr
% before we exit the loop
%
% Output:
%
% w: each w(:,ii) is the DMD mode
% e: e(ii) is the eigenvalue corresponding to w(:,ii)
% b: b(ii) is the coefficient corresponding to the weight of each mode
% final_error:
% evals: the number of times varpro2.m was called
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');

% initial guess at alpha_init which we will perturb 
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

% inputs for varpro2
m = length(t);
[is,~] = size(X);
ia = r;
n = r;

count = 1;
evals = 0;

while true % repeat until a satisfying error is obtained
    
    if evals == 0
        alpha_init = alpha_init0; % 1st try is our initial guess
    else
        alpha_init = alpha_init0 .* 10^(3*rand-2.5); % after that we multiply it by a random amount
    end
    
    % execution of varpro2.m to find w,e,b
    [wtemp,e,niter,err,imode,alphas] = varpro2edit(transpose(X),t, ...
        @varpro2expfun,@varpro2dexpfun,m,n,is,ia,alpha_init,varpro_opts_edit());
    wtemp = transpose(wtemp);
    b = sqrt(sum(abs(wtemp).^2,1))'; 
    w = wtemp*diag(1./b); % normalise w  
    
    % it was observed that when varpro2 lowers the error to zero this leads to less reliable results
    if all(err) > 0 % so avoid this
       
        lasterr(count) = err(30); 
        
        % break if the algorithm repeatedly finds the same error and it is below our maximum allowed value
        if count > 1 
            movmean = (lasterr(count) + lasterr(count-1))/2;
            if movmean < maxerr
                if abs(lasterr(count) - lasterr(count-1))/mean([lasterr(count),lasterr(count-1)])...
                        <= 0.01
                    final_error = lasterr(count);
                    evals = evals + 1;
                    break
                end
            end
        end
        count = count+1;
    end
    evals = evals + 1;
end