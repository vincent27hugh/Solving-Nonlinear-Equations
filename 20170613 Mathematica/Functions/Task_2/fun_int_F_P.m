% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function intF = fun_int_F_P(epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,g_p,mu,epsilon_u,typen,pos) 
% integration of (1-F(x))from a to b
% typen = 'i' 'ii' 'iii' 'iv'
syms x
if pos == 1
    DELTAP_pms=fun_DELTAP_pms(x,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    
    fun = DELTAP_pms.*(1-fun_F_x(x,typen,epsilon_u));

    intF=int(fun,epsilond_pms,epsilon_u);
elseif pos==2
    DELTAP_ps=fun_DELTAP_ps(x,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    
    fun = DELTAP_ps.*(1-fun_F_x(x,typen,epsilon_u));
    
    intF=int(fun,epsilond_ps,epsilon_u);
elseif pos ==3
    DELTAP_pps=fun_DELTAP_pps(x,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    
    fun = DELTAP_pps.*(1-fun_F_x(x,typen,epsilon_u));
    
    intF=int(fun,epsilond_pps,epsilon_u);
end

    

return