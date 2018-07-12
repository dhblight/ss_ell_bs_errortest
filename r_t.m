function[out] = r_t(ri,r0,t,tau)

out = ri + (r0-ri).*(1-exp(-t./tau));

end