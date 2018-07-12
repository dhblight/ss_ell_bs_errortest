function[out1,out2,out3] = c_t_3species(c0,t,tau1,tau2,n)

lambda1 = 1/tau1;
lambda2 = 1/tau2;

out1 = c0.*exp(-t./tau1);
out2 = c0.*(lambda1./(lambda2-lambda1)).*(exp(-t./tau1)-exp(-t./tau2)).*(1./n);
out3 = c0./n - out1./n - out2;

end