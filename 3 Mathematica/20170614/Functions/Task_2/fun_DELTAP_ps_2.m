% May 22, 2017
function DELTAP_ps...
    = fun_DELTAP_ps_2(tau_pms,tau_ps,tau_pps,g_p,mu)

temp1 = tau_ps-mu^2*g_p/tau_pps;
temp2 = 1+mu*g_p/tau_pps;

DELTAP_ps = temp2/(temp1);

return