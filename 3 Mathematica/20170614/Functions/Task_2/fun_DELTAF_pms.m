% May 23, 2017
% May 22, 2017
%%
function DELTAF_pms...
    = fun_DELTAF_pms(epsilon,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta)

[DELTAF_pms,~,~]...
    = fun_DELTAF(epsilon,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);

return