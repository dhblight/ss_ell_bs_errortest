%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to test the reliability of ss_ell_bs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start timing so we can display progress 
tic

% seed random number generation based on time
rng('shuffle');

%% constants, variables, objects
numtrials = 20; % number of times we add noise and calculate the dmd 

nsph = 10; % number of small spheres to make big ellipsoid
r0 = 10e-9; % radius of small spheres
r1 = r0 * nsph^(1/3); % radius of final big spheres
ratio = 3; % ra/rb for the ellipse, ra > rb for prolate (needle-like) ellipsoid 
rb = r1 * ratio^(-1/3); % make it so that the volume of the ellipsoid and sphere are equal
ra = rb * ratio;

Mss = 180000; % molar mass of small spheres (g/mol)
Mbs = Mss * nsph; % " " " big spheres
Mell = Mss * nsph; % " " " ellipsoids

c0 = 1e-9; % initial concentration (mol/cm^3) 

ntau = 50; % there will be ntau^2 * numtrials * length(rankvals) iterations in total 
tau1vals = linspace(1,99,ntau); % relaxation time for small spheres -> ellipsoid
tau2vals = linspace(10,990,ntau); % relaxation time for ellipsoid -> big sphere

rankvals = 6; % ranks to use

n = 300; % number of spatial measurements
m = 50; % number of temporal measurements

sigma = 0.03; % multiplicative error

% temporary matrix used to calculate the means after all of the trials
tempmat = zeros(6,numtrials); 

% 3-dimensional objects,(x,y) plane: varying tau1 and tau2, z axis: varying rank
relerror = zeros(length(tau1vals),length(tau2vals),length(rankvals));
meanevals = zeros(length(tau1vals),length(tau2vals),length(rankvals));
maxevals = zeros(length(tau1vals),length(tau2vals),length(rankvals));
wfrac1 = zeros(length(tau1vals),length(tau2vals),length(rankvals));
wfrac2 = zeros(length(tau1vals),length(tau2vals),length(rankvals));
wfracconst = zeros(length(tau1vals),length(tau2vals),length(rankvals));

% the grid of co-ordinates for the data
qi = logspace(7.7,9.1761,n);
t = logspace(-1,3,m);
[qgrid,T] = meshgrid(qi,t);
qgrid = qgrid.';T = T.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% main program
count1 = 1;
for tau1_iter = tau1vals
    count2 = 1;
    for tau2_iter = tau2vals
        count3 = 1;
        
        % generate clean data
        
        % form factors
        Pss = sphere_p_q(qgrid,r0);
        Pell = elipse_p_q(qgrid,ra,rb);
        Pbs = sphere_p_q(qgrid,r1);
        
        % concentrations
        [Css,Cell,Cbs] = c_t_3species(c0,T,tau1_iter,tau2_iter,nsph);
        
        % intentisty
        Xclean = Mss^2 .* Pss .* Css + Mbs^2 .* Pbs .* Cbs + ...
            Mell^2 .* Pell .* Cell;
        
        for r = rankvals
            
            for kk = 1:numtrials
                
                % add noise
                G = sigma .* randn(n,m);
                X = Xclean.*(1+G);
                
                % DMD
                imode = 1; % don't project to POD modes
                [w,e,b,~,evals] = multioptdmd3(X,t,r,0.032);
                e(abs(imag(e))>0.1) = 0; % only care about real eigenvalues
                e = real(e);
                eig_tau = -1./ e; % tau values corresponding to the eigenvalues
                eig_tau_original = eig_tau; % eig_tau will be modified but we still want the values
                
                % find indices of eigenvalues and their errors
                [err1,ind1] = min(abs(eig_tau - tau1_iter));
                [err2,ind2] = min(abs(eig_tau - tau2_iter));
                
                % Can't use the same eigenvalue for both unless tau1 = tau2
                if tau1_iter ~= tau2_iter
                    if ind1 == ind2
                        eig_tau(ind1) = inf; % make sure this value isn't closest to tau so the next closest is used
                        if err1 > err2 % calculate the next closest value
                            [err1,ind1] = min(abs(eig_tau - tau1_iter)); % for tau1 if it was closer to tau2 originally
                        else
                            [err2,ind2] = min(abs(eig_tau - tau2_iter)); % or vice versa
                        end
                    end
                end
                
                % calculated tau value
                tau1_temp = eig_tau_original(ind1);
                tau2_temp = eig_tau_original(ind2);
                
                % weight fraction
                wfrac1temp = b(ind1)/sum(b);
                wfrac2temp = b(ind2)/sum(b);
                [~,indconst] = min(abs(e)); % index of constant mode 
                wfracconsttemp = b(indconst) / sum(b); % for comparison, constant mode normally has highest b value 
                
                % store info
                tempmat(1,kk) = tau1_temp;
                tempmat(2,kk) = tau2_temp;
                tempmat(3,kk) = evals;
                tempmat(4,kk) = wfrac1temp;
                tempmat(5,kk) = wfrac2temp;
                tempmat(6,kk) = wfracconsttemp;
                
            % display progress
            totiterations = length(tau1vals)*length(tau2vals)*length(rankvals)*numtrials;
            currentiteration = length(tau2vals)*length(rankvals)*numtrials*(count1-1)+...
                (length(rankvals)*numtrials*(count2-1)) + numtrials*(count3-1) + kk;
            timelapsed = toc;            
            [hours,mins,secs] = hourminsec(timelapsed);
            fprintf('iterations completed: %d / %d      progress: %.2f %%       time elapsed: %02d:%02d:%02d \n',...
                currentiteration,totiterations,currentiteration/totiterations*100,hours,mins,round(secs));            
           
            end
            
            % store info
            meanevals(count1,count2,count3) = mean(tempmat(3,:));
            maxevals(count1,count2,count3) = max(tempmat(3,:));
            
            wfrac1(count1,count2,count3) = mean(tempmat(4,:));
            wfrac2(count1,count2,count3) = mean(tempmat(5,:));
            wfracconst(count1,count2,count3) = mean(tempmat(6,:));
            
            tau1_calc = mean(tempmat(1,:),'omitnan');
            tau2_calc = mean(tempmat(2,:),'omitnan');
            errel1 = abs(tau1_calc-tau1_iter)/tau1_iter;
            errel2 = abs(tau2_calc-tau2_iter)/tau2_iter;
            erreltot = sqrt(errel1^2 + errel2^2);
            relerror(count1,count2,count3) = erreltot;                      
            
            count3 = count3+1;
        end
        count2 = count2+1;
    end
    count1 = count1 +1;
end

%% plot and/or save
%
%
%

% saves data as "prefix_dd-mmm-yyyy_HH-MM.mat"
formatout = 'dd-mmm-yyyy_HH-MM';
save(strcat('reliability_test','_',datestr(datetime('now'),formatout),'.mat'));
