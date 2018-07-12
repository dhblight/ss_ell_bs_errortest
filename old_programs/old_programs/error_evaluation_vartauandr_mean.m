%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for calculating the error in the relaxation times for an elipsoid
% changing into a sphere calculated by optdmd. Uses various configurations of 
% initial radii and relaxation times.
% 
% The function which calculates I(q) (elipse_p_q.m) takes a long time to
% evaluate leading to a long running time depending on the number of values
% of r and tau used (~40s per grid co-ordinate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');

%% constants, variables, objects
numtrials = 200; %number of times we add noise qnd calculate the dmd 

r0 = 40e-9; % radius of final sphere
r1 = 5e-9; % initial semiminor axis (value used when varying tau)
r2 = 20e-9; % initial semimajor axis (value used when varying tau)
r1vals = linspace(1e-9,30e-9,2); % values over which to vary r1
r2vals = linspace(1e-9,30e-9,2); % values over which to vary r2

tau1 = 2; % relaxation time of r1 -> r0 (value used when varying r)
tau2 = 2; % relaxation time of r2 -> r0 (value used when varying r)
tau1vals = linspace(1,20,20); % values over which to vary tau1 
tau2vals = linspace(1,20,20); % values iver which to vary tau2

n = 300; % # of spatial measurements
m = 30; % # of temporal measurements

errel1_var_tau = zeros(length(tau1vals),length(tau2vals)); % storing errors for different values of tau
errel2_var_tau = zeros(length(tau1vals),length(tau2vals));
erreltot_var_tau = zeros(length(tau1vals),length(tau2vals));

errel1_var_r = zeros(length(r1vals),length(r2vals)); % storing errors for different values of r
errel2_var_r = zeros(length(r1vals),length(r2vals));
erreltot_var_r = zeros(length(r1vals),length(r2vals));

rankusedtau = zeros(length(tau1vals),length(tau2vals)); % store the rank used
rankusedr = zeros(length(r1vals),length(r2vals));

tempmat = zeros(2,numtrials); % for storing the erros at each trial
tempvec = zeros(2,1);

%sigma = 5e-5; % for adding noise (additive)
sigma = 0.03;

% the grid of co-ordinates for the data
qi = logspace(7,9.1761,n);
t = logspace(-1,2,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% varying tau

count1 = 1; % tracks the index of the tau1 value being used as opposed to its actual value

for tau1_iter = tau1vals 
    count2 = 1; % same for tau2
        
    for tau2_iter = tau2vals

        % clean data
        Xclean = zeros(n,m);
        for ii = 1:n
            for jj = 1:m
                Xclean(ii,jj) = elipse_p_q(qgrid(ii,jj),r_t(r1,r0,T(ii,jj),tau1_iter),r_t(r2,r0,T(ii,jj),tau2_iter));
            end
        end
                
        % loop for number of trials
        for ii = 1:numtrials
            
            % add noise
            G = sigma .* randn(n,m);
            %X = Xclean+G;
            X = Xclean.*(1+G);
                                   
            % number of singular values to keep
            [U, S, V] = svd(X, 'econ');
                        
            %coeff = optimal_SVHT_coef(m/n,1); % Gavish & Donoho
            coeff = optimal_SVHT_coef(m/n,0);
            x = diag(S);           
            %threshold = sqrt(n) * sigma * coeff;
            threshold = coeff * median(x);
            r = sum(x >= threshold);
                                            
            % DMD
            imode = 1;
            %maxiter = 30;
            %tol = sigma;            
            %opts = varpro_opts('maxiter',maxiter,'tol',tol,'eps_stall',1e-9);            
            %imode = 2; % projected version
            
            % sometimes gives errors "Input to SVD must not contain NaN or Inf." for a certain rank 
            % catch these errors so the rest of the program still executes
            try 
                [w,e,b] = optdmd(X,t,r,imode); % w: dmd modes, e:eigenvalues, b:coefficients
                %[w,e,b,atilde] = optdmd(X,t,r,imode,opts,[],U);
            catch 
                % if error            
                r = round(m*(3/5)); % set a high rank
                count3 = 1;
                while count3 > 0
                    try
                        [w,e,b] = optdmd(X,t,r,imode);
                        %[w,e,b,atilde] = optdmd(X,t,r,imode,opts,[],U);
                    catch
                        count3 = count3+1;  
                        r = r - 1;
                    end
                    count3 = count3-1;
                    % if optdmd ran without error count3 = 0 -> exit while
                    % loop, if not count3 = 1 -> try again with a rank of
                    % one less until it succeeds.
                end
            end
                                   
            rankusedtau(count1,count2) = r; % store the rank
                        
            % only care about real eigenvalues
            e(abs(imag(e))>0.1) = 0;
            e = real(e);
                        
            
            eig_tau = -1./ e; % tau values corresponding to the eigenvalues
            eig_tau_original = eig_tau; % eig_tau will be modified but we still want the values
            
            % find indices of eigenvalues and their errors
            [err1,ind1] = min(abs(eig_tau - tau1_iter));
            [err2,ind2] = min(abs(eig_tau - tau2_iter));
            
%             % Can't use the same eigenvalue for both unless tau1 ~ tau2
%             if max(tau1_iter/tau2_iter,tau2_iter/tau1_iter) > 1.1
%                 if ind1 == ind2
%                     eig_tau(ind1) = inf; % make sure this value isn't closest to tau so the next closest is used
%                     if err1 > err2 %
%                         [err1,ind1] = min(abs(eig_tau - tau1_iter));
%                     else
%                         [err2,ind2] = min(abs(eig_tau - tau2_iter));
%                     end
%                 end
%             end
            
            % calculated tau value
            tau1_calc = eig_tau_original(ind1);
            tau2_calc = eig_tau_original(ind2);
            
            % store the calculated  tau values for each trial in a temporary matrix
            tempmat(1,ii) = tau1_calc;
            tempmat(2,ii) = tau2_calc;
            
        end
        
        % calculate the mean  of the calculated tau values over the trials
        tempmat(tempmat > max([tau1_iter,tau2_iter])*10) = 0; % remove obviously erroneous results
        for jj = 1:2
            row = tempmat(jj,:);
            tempvec(jj) = sum(row) ./ sum(row~=0); % calculate mean not taking into account the zero elements
        end
        
        % store error in matrices
        errel1_var_tau(count1,count2) = abs(tempvec(1)-tau1_iter)/tau1_iter;
        errel2_var_tau(count1,count2) = abs(tempvec(2)-tau2_iter)/tau2_iter;
        erreltot_var_tau(count1,count2) = sqrt((abs(tempvec(1)-tau1_iter)/tau1_iter)^2 + (abs(tempvec(2)-tau2_iter)/tau2_iter)^2);
        
        % display progress
        fprintf('iterations completed varying tau: %d / %d\n',(length(tau2vals)*(count1-1))+count2,length(tau1vals)*length(tau2vals));
        count2 = count2 + 1;              
        
    end
    
    count1  = count1 + 1;

end

% meshgrid(X,Y), tau2 changes as you go along each row -> x, 
% tau1 changes as you go down each column -> y
% element (i,j) of error matrix corresponds to the error for tau1 =
% TAU1(i,j) and tau2 = TAU2(i,j)
[TAU2,TAU1] = meshgrid(tau2vals,tau1vals); 

% for a 2d plot of tau1/tau2 vs error
TAU1_TAU2 = TAU1./TAU2;
tau1_tau2 = reshape(TAU1_TAU2,1,[]);
erreltot_var_tau_vec = reshape(erreltot_var_tau,1,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% varying radii
% count1 = 1;
% 
% for r1_iter = r1vals 
%     count2 = 1;
%         
%     for r2_iter = r2vals
% 
%         Xclean = zeros(n,m);
%         for ii = 1:n
%             for jj = 1:m
%                 Xclean(ii,jj) = elipse_p_q(qgrid(ii,jj),r_t(r1_iter,r0,T(ii,jj),tau1),r_t(r2_iter,r0,T(ii,jj),tau2));
%             end
%         end
%         
%         for ii = 1:numtrials
%             G = sigma .* randn(n,m);
%             %X = Xclean+G;
%             X = Xclean.*(1+G);
%             
%             [U, S, V] = svd(X, 'econ');            
%             %coeff = optimal_SVHT_coef(m/n,1);
%             coeff = optimal_SVHT_coef(m/n,0);
%             x = diag(S);
%             %threshold = sqrt(n) * sigma * coeff;
%             threshold = median(x) * coeff;            
%             r = sum(x >= threshold);
%             
%             imode = 1;
%             try 
%                 [w,e,b] = optdmd(X,t,r,imode); 
%             catch 
%                 r = round(m*(3/5)); 
%                 count3 = 1;
%                 while count3 > 0
%                     try
%                         [w,e,b] = optdmd(X,t,r,imode);
%                     catch
%                         count3 = count3+1;  
%                         r = r - 1;
%                     end
%                     count3 = count3-1;
%                 end
%             end
%             
%             rankusedr(tau1,tau2) = r;
%                        
%             e(abs(imag(e))>0.1) = 0;
%             e = real(e);
%                         
%             eig_tau = -1./ e; % 
%             eig_tau_original = eig_tau; 
%                         
%             [err1,ind1] = min(abs(eig_tau - tau1));
%             [err2,ind2] = min(abs(eig_tau - tau2));
%             
%             if abs(tau1-tau2) > 0.5
%                 if ind1 == ind2
%                     eig_tau(ind1) = inf; 
%                     if err1 > err2 
%                         [err1,ind1] = min(abs(eig_tau - tau1));
%                     else
%                         [err2,ind2] = min(abs(eig_tau - tau2));
%                     end
%                 end
%             end            
%          
%             tau1_calc = eig_tau_original(ind1);
%             tau2_calc = eig_tau_original(ind2);
%                        
%             tempmat(1,ii) = tau1_calc;
%             tempmat(2,ii) = tau2_calc;
%             
%         end
%         
%         tempmat(tempmat > max([tau1,tau2])*10) = 0; 
%         for jj = 1:2
%             row = tempmat(jj,:);
%             tempvec(jj) = sum(row) ./ sum(row~=0); 
%         end
%         
%         errel1_var_r(count1,count2) = abs(tempvec(1)-tau1)/tau1;
%         errel2_var_r(count1,count2) = abs(tempvec(2)-tau2)/tau2;
%         erreltot_var_r(count1,count2) = sqrt((abs(tempvec(1)-tau1)/tau1)^2 + (abs(tempvec(2)-tau2)/tau2)^2);
%         
%         fprintf('iterations completed varying r: %d / %d\n',(length(r2vals)*(count1-1))+count2,length(r1vals)*length(r2vals));
%         count2 = count2 + 1;
%         
%     end
%     
%     count1  = count1 + 1;
% 
% end
% 
% [R2,R1] = meshgrid(r2vals,r1vals);
% 
% R1_R2 = R1./R2;
% r1_r2 = reshape(R1_R2,1,[]);
% erreltot_var_r_vec = reshape(erreltot_var_r,1,[]);