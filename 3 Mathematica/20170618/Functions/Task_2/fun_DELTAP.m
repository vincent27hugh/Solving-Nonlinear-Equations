% May 22, 2017
function [DELTAP_pms,DELTAP_ps,DELTAP_pps]...
    = fun_DELTAP(epsilon,epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,g_p,mu)

DELTAP_pms = NaN(size(epsilon));
DELTAP_ps = NaN(size(epsilon));
DELTAP_pps = NaN(size(epsilon));

for i = 1:numel(epsilon)
    if epsilon(i) >= epsilond_pms
        DELTAP_pms(i) = fun_DELTAP_pms_1(tau_pms,tau_ps,tau_pps,g_p,mu);
        
        DELTAP_ps(i) = (tau_pms*DELTAP_pms(i)-1)/mu;
        
        DELTAP_pps(i) = tau_pms*DELTAP_pms(i)/tau_pps;
    elseif epsilon(i) < epsilond_pms && ...
            epsilon(i) >= epsilond_ps
        DELTAP_pms(i) = 0;
        
        DELTAP_ps(i) = fun_DELTAP_ps_2(tau_pms,tau_ps,tau_pps,g_p,mu);
        
        DELTAP_pps(i) = (1+mu*DELTAP_ps(i))/tau_pps;
    elseif epsilon(i) < epsilond_ps && ...
            epsilon(i) >= epsilond_pps
        DELTAP_pms(i) = 0;
        
        DELTAP_ps(i) = 0;
        
        DELTAP_pps(i) = 1/tau_pps;
    else
        DELTAP_pms(i) = 0;
        
        DELTAP_ps(i) = 0;
        
        DELTAP_pps(i) = 0;
    end
end
return