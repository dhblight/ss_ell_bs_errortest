function[out1,out2] = c_t(c0,t,tau,fact)


out1 = c0*exp(-t./tau);
out2 = c0*(1-exp(-t./tau))*(1/fact);

end