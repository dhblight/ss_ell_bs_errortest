function[out] = ss_ell_bs_intensity(c0,nsph,r0,r1,ra,rb,qgrid,T,Mss,M,tau1,tau2)

t = T(1,:); % vector of times
dim = size(qgrid); 

% work out contribution of small spheres
Css = c0.*exp(-t./tau1);
Pss = sphere_p_q(qgrid,r0);
Iss = Mss^2 .* Pss .* Css;

% 3 dimensional object to store intensity due to ellipsoids formed at each
% t(i) e.g Ielltemp(7,15,3) = intensity due to ellipsoids formed at t =
% t(3) observerd at q = qi(7) and t = t(15)
Ielltemp = zeros(dim(1),dim(2),length(t));

% at each time calculate the concentration of big spheres that have been
% formed between t(i-1 )and t(i)  (or t0=0 and t(i) if i = 1) 
for ii = 1:length(t)      
    
    % re-initialize at the strt of each loop  
    Pell = zeros(dim(1),dim(2));
    
    % concentration of big spheres at time t(i)
    if ii == 1
        Cbs = (c0-Css(ii))/nsph; 
    else
        Cbs = (Css(ii-1)-Css(ii))/nsph;
    end
    
    % work out what the radii of the ellipsoids formed at time t(i) will
    % be for the rest of the experiment (uses r_t.m function), and create 
    % grids of those values
    ragrid = r_t(ra,r1,T(:,ii:end)-t(ii),tau2);
    rbgrid = r_t(rb,r1,T(:,ii:end)-t(ii),tau2);
    
    % work out the form factor of only the ellipsoids formed at t(i) for
    % all q and the remaining times
    Pell(:,ii:end) = elipse_p_q(qgrid(:,ii:end),ragrid,rbgrid);
    
    % calculate the intensity and store this as a layer in a 3-dimensional object
    Ielltemp(:,:,ii) = M^2 .* Pell .* Cbs;    
end

% the total intensity due to the ellipsoids is the sum of the intensities
% of ellipsoids formed at each t(i)
Iell = sum(Ielltemp,3);

% total intensity = intensity of small spheres + ellipsoids
out = Iell + Iss;
    
    

