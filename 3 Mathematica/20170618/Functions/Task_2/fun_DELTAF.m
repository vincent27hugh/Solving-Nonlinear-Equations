% May 22, 2017
function [DELTAF_pms,DELTAF_ps,DELTAF_pps]...
    = fun_DELTAF(epsilon,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta)

DELTAF_pms = NaN(size(epsilon));
DELTAF_ps = NaN(size(epsilon));
DELTAF_pps = NaN(size(epsilon));

for i = 1:numel(epsilon)
    [DELTAP_pms,DELTAP_ps,DELTAP_pps] = fun_DELTAP(epsilon(i),...
        epsilond_pms,epsilond_ps,epsilond_pps,...
        tau_pms,tau_ps,tau_pps,g_p,mu);
        
    if epsilon(i) >= epsilonc_pms
        DELTAF_ps(i) = fun_DELTAF_ps_1(r,lambda,mu);
        
        DELTAF_pms(i) = (1+mu*DELTAF_ps(i))/(r+lambda+mu);
        
        DELTAF_pps(i) = DELTAF_pms(i);
    elseif epsilon(i) < epsilonc_pms && ...
            epsilon(i) >= epsilonc_ps
        DELTAF_pms(i) = delta*DELTAP_pms;
        
        DELTAF_ps(i) = fun_DELTAF_ps_2(r,lambda,mu,g_p,delta,DELTAP_pms);
        
        DELTAF_pps(i) = (1+mu*DELTAF_ps(i))/(r+lambda+mu);
        
    elseif epsilon(i) < epsilonc_ps && ...
            epsilon(i) >= epsilonc_pps
        
        DELTAF_pms(i) = delta*DELTAP_pms;
        
        DELTAF_ps(i) = delta*DELTAP_ps;
        
        DELTAF_pps(i) = (1+mu*DELTAF_ps(i))/(r+lambda+mu);
        
    elseif epsilon(i) < epsilonc_pps

        DELTAF_pms(i) = delta*DELTAP_pms;
        
        DELTAF_ps(i) = delta*DELTAP_ps;
        
        DELTAF_pps(i) = delta*DELTAP_pps;
        
    end
end
return