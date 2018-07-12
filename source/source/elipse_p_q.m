function[out] = elipse_p_q(q,Ra,Rb)

% function to caluculate the form factor of an ellipsoid

nu = Ra ./ Rb;

func1 = @(z) (3 .* ((sin(z) - z .* cos(z)) ./ z.^3)).^2;
func2 = @(x) func1(q .* Rb .* sqrt(1 + x.^2 .* (nu.^2 - 1)));

out = integral(func2, 0, 1,'ArrayValued',true);

end

