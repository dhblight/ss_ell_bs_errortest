function[out] = sphere_p_q(q,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to caluculate the form factor of a sphere.
%
% Multiplied by 9 to get same results as elipse_p_q with ra = rb.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = 9 * ((sin(q.*r) - q.*r.*cos(q.*r)) ./ (q.*r).^3).^2;

end
