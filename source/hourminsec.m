function[h,m,s] = hourminsec(t)

% gives hours, mins, secs for hh:mm:ss format for a time t given in seconds

temp = t/3600;
h = floor(temp);
m = floor((temp-h)*60);
s = (temp-h-m/60)*3600;
end