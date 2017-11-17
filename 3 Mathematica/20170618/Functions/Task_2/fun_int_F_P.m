% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function intF = fun_int_F_P(epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,g_p,mu,epsilon_u,typen,pos) 
% integration of (1-F(x))from a to b
% typen = 'i' 'ii' 'iii' 'iv'
if pos == 1
    DELTAP_pms=@(x)fun_DELTAP_pms(x,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    
    fun = @(x)(DELTAP_pms(x).*(1-fun_F_x(x,typen,epsilon_u)));

    intF=integral(fun,epsilond_pms,epsilon_u);
elseif pos==2
    DELTAP_ps=@(x)fun_DELTAP_ps(x,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    
    fun = @(x)(DELTAP_ps(x).*(1-fun_F_x(x,typen,epsilon_u)));
    
    intF=integral(fun,epsilond_ps,epsilon_u);
elseif pos ==3
    DELTAP_pps=@(x)fun_DELTAP_pps(x,epsilond_pms,epsilond_ps,...
    epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    
    fun = @(x)(DELTAP_pps(x).*(1-fun_F_x(x,typen,epsilon_u)));
    
    intF=integral(fun,epsilond_pps,epsilon_u);
end

    

return