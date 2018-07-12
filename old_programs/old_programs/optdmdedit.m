function [w,e,b,err] = optdmdedit(X,t,r,imode,alpha_init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wrapper for computing the optimized DMD of data 
%
% Makes use of varpro2 routines
%
% Input:
%
% X - data matrix, length(t) columns with each a snapshot
%   of some time-varying system
% t - times corresponding to each snapshot
% r - rank of fit, i.e. number of exponentials to fit
% imode - flag, determines type of computation
%   imode = 1, fit full data, slower
%   imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})
% varargin{1} - options structure for varpro2. see varpro_opts.m
%   for details
% varargin{2} - initial guess. if not provided, one is computed
%   using trapezoidal rule approximation
% varargin{3} - orthogonal basis for projection (if POD modes precomputed
%   or a different basis desired)
%
% Output:
%
% w - each column is a DMD mode
% e - each entry e(i) is an eigenvalue corresponding to w(:,i)
% b - the best fit coefficient of each DMD mode
% varargout{1} - return projected system matrix 
%   A =  (u'*w)*diag(e)*(pinv(w)*u) where u is the first
%   r POD modes or varargin{3}.
% varargout{2} - return basis for projection
% varargout{3} - return full system matrix A = w*diag(e)*pinv(w)

% 
% X should be approximated by
%
% X ~ w*diag(b)*exp(e*t')
%
% if t is a column vector.
%
% Travis Askham 2017
%
% MIT License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = varpro_opts_edit();


% fit to all of data
m = length(t);
[is,~] = size(X);
ia = r;
n = r;
[w,e,niter,err,imode,alphas] = varpro2(transpose(X),t, ...
    @varpro2expfun,@varpro2dexpfun,m,n,is,ia,alpha_init,opts);

w = transpose(w);

% normalize
b = sqrt(sum(abs(w).^2,1))';
w = w*diag(1./b);
    
end


        

