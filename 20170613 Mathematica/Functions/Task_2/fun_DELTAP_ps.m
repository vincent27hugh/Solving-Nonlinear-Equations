% May 22, 2017
function DELTAP_ps...
    = fun_DELTAP_ps(epsilon,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu)

[~,DELTAP_ps,~]...
    = fun_DELTAP(epsilon,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);

return