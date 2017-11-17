function frac_delw = LongTermPlasticity(delt)
% Long Term Plasticity
a.ltp = 0.015;
tau.tlp = 13;
a.ltd = 0.021;
tau.ltd = 20;
if(delt>0)
    frac_delw = a.lpt*exp(-delt/tau.ltp);
else
    frac_delw = a.lpd*exp(delt/tau.ltd);
end
end