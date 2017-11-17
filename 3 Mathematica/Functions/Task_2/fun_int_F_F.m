% symbolic Jun 13, 2017
% May 23, 2017
% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function intF = fun_int_F_F(epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta,typen,epsilon_u,pos)
% integration of (1-F(x))from a to b
% typen = 'i' 'ii' 'iii' 'iv'
syms x

if pos == 1
    DELTAF_pms=fun_DELTAF_pms(x,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);
    
    fun = DELTAF_pms.*(1-fun_F_x(x,typen,epsilon_u));

    intF=int(fun,epsilond_pms,epsilon_u);
    
elseif pos==2
    DELTAF_ps=fun_DELTAF_ps(x,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);
    
    fun = (DELTAF_ps.*(1-fun_F_x(x,typen,epsilon_u)));
    
    intF=int(fun,epsilond_ps,epsilon_u);
elseif pos ==3
    
    DELTAF_pps=fun_DELTAF_pps(x,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);
    
    fun = (DELTAF_pps.*(1-fun_F_x(x,typen,epsilon_u)));
    
    intF=int(fun,epsilond_pps,epsilon_u);
end

    

return