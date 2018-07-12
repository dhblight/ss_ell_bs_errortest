function[out] = rectpolydisperse_p_q(q,R,p)

% they don't match at boundary!!!

% p = sigmaR/R i.e. the root-mean-square deviation in radius divided by the
% average radius

cont = 3e-6;
width = sqrt(3) .* p .* R;
qr = q .* R;
qw = q .* width;

if qr < 0.1 % Guinier approx
    R2g = @(p,R) (3/5) .* R.^2 .* ((1 + 28.*p.^2 + 126.*p.^4 + 108.*p.^6 + 27.*p.^8) ./...
        (1 + 15.*p.^2 + 27.*p.^4 + (27./7).*p.^6));
    
    temp = cont^2 * 1e8 .* (4.*pi./3) .* R.^3 .* (1 + 15.*p.^2 + 27.*p.^4 +(27/7).*p.^6)...
        .* (1./(1 + 3.*p.^2)) .* exp((-1./3).*R2g(p,R).*q.^2);
    
    
else
    temp = -0.5*qw + qr*qr*qw + (qw^3)/3 - (5/2)*cos(2*qr)*sin(qw)*cos(qw)...
        + 0.5*qr*qr*cos(2*qr)*sin(2*qw) + 0.5*qw*qw*cos(2*qr)*sin(2*qw)...
        + qw*qr*sin(2*qr)*cos(2*qw) + 3*qw*(cos(qr)*cos(qw))^2 ...
        + 3*qw*(sin(qr)*sin(qw))^2 - 6*qr*cos(qr)*sin(qr)*cos(qw)*sin(qw);    
end

out = 8.*pi.*pi.*cont.*cont.*(1./(width.*q.^7)).* temp ;

end



