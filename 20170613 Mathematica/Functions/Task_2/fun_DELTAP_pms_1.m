% May 22, 2017
function DELTAP_pms...
    = fun_DELTAP_pms_1(tau_pms,tau_ps,tau_pps,g_p,mu)

temp1 = 1+tau_ps/mu;
temp2 = tau_ps*tau_pms/mu;
temp3 = mu*(g_p*tau_pms/tau_pps+(1-g_p));

DELTAP_pms = temp1/(temp2-temp3);

return