% Task #20170302
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function Fn = fun_solveMar5_2(var,epsilon_d,epsilon_c,theta)
% var vetcor represent 
% [epsilon_dps;epsilon_cps;theta_ps]; 
% length(y)*3 variables
% Parameters; global 
global c_f c_p beta phi delta sigma lambda b r epsilon_u mu 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global step width typen
n=int16(2*width/step+1);
% length of the lattice
y = (-width:step:width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1
for i =1
    p = fun_p(y(i));
    theta_pps=theta(i+1);
    % theta_{p^+,sigma}
    q_theta = fun_q_theta(theta(i));
    % q(theta_{p,sigma})
    q_thetapps=fun_q_theta(theta_pps);
    %q(theta_{p^+,sigma})

    g_p=fun_gp(y(i));
    %
    g_pps=fun_gp(y(i+1));
    %g_{p^+,sigma}
    
    epsilon_cpms=epsilon_c(i);
    % epsilon_{p^-,sigma}^c
    epsilon_dpms=epsilon_d(i);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % epsilon_{p^-,sigma}^d
    epsilon_cpps=epsilon_c(i+1);
    % epsilon_{p^+,sigma}^c
    epsilon_dpps=epsilon_d(i+1);
    % epsilon_{p^+,sigma}^d
    
    intF_ucpm=fun_int_F(epsilon_cpms,epsilon_u,typen);%!!!!!!!!!!!!!!!!!!!!!!
    % int(epsilon_u,epsilon_{p^-,sigma}^c)
    
    intF_cpdpm=fun_int_F(epsilon_dpms,epsilon_c(i),typen);%!!!!!!!!!!!!!!!!!!!!!!
    % int(epsilon_{p,sigma}^c, epsilon_{p^-,sigma}^d)
   

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp111 = ((beta+phi*(1-beta))*c_f*(1-var(i)))/((1-beta)*(1-phi))...
        +beta*c_p*var(i)/(1-beta);
    temp112 = b-p+theta(i)*(temp111);
    
    part11 = temp112/sigma;

    temp121 = (1+(mu*g_p*delta)/(r+lambda+theta(i)*q_theta))*...
        intF_ucpm/(r+lambda+mu*(1-g_p));
    

    temp123 = intF_cpdpm/(r+lambda+theta(i)*q_theta);
    

    temp125 = phi*(1-fun_F_x(epsilon_c(i),typen))-r-lambda-mu;
    temp126 = delta*sigma*(epsilon_c(i)-epsilon_dpms)/(r+lambda+theta(i)*q_theta);

    temp127 = (1+(mu*g_p*delta)/(r+lambda+theta_pps*q_thetapps))*...
        (epsilon_c(i)-epsilon_cpps)/...
        (r+lambda+mu*(1-g_pps));
    temp128 = delta*(epsilon_cpps-epsilon_dpps)/...
        ((1-phi)*(r+lambda+theta_pps*q_thetapps));

    part12 = epsilon_c(i)+lambda*(temp121)+lambda*delta*(temp123)+...
        temp125*temp126/(sigma*(1-phi))+mu*g_p*(temp127+temp128);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Fn = [part11 - part12];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:n-1
    Fnm1=Fn;
    
    p = fun_p(y(i));
    
    theta_pps=theta(i+1);
    % theta_{p^+,sigma}
    q_theta = fun_q_theta(theta(i));
    % q(theta_{p,sigma})
    q_thetapps=fun_q_theta(theta_pps);
    %q(theta_{p^+,sigma})

    g_p=fun_gp(y(i));
    %
    g_pps=fun_gp(y(i+1));
    %g_{p^+,sigma}
    
    epsilon_cpms=epsilon_c(i-1);
    % epsilon_{p^-,sigma}^c
    epsilon_dpms=epsilon_d(i-1);
    % epsilon_{p^-,sigma}^d
    epsilon_cpps=epsilon_c(i+1);
    % epsilon_{p^+,sigma}^c
    epsilon_dpps=epsilon_d(i+1);
    % epsilon_{p^+,sigma}^d
    
    intF_ucpm=fun_int_F(epsilon_cpms,epsilon_u,typen) ;
    % int(epsilon_u,epsilon_{p^-,sigma}^c)
    intF_cppm=fun_int_F(epsilon_c(i),epsilon_cpms,typen) ;
    % int(epsilon_{p^-,sigma}^c, epsilon_{p,sigma}^c)
    intF_cpdpm=fun_int_F(epsilon_dpms,epsilon_c(i),typen) ;
    % int(epsilon_{p,sigma}^c, epsilon_{p^-,sigma}^d)
    intF_dpmp=fun_int_F(epsilon_d(i),epsilon_dpms,typen) ;
    % int(epsilon_{p^-,sigma}^d, epsilon_{p,sigma}^d)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp111 = ((beta+phi*(1-beta))*c_f*(1-var(i)))/((1-beta)*(1-phi))...
        +beta*c_p*var(i)/(1-beta);
    temp112 = b-p+theta(i)*(temp111);
    
    part11 = temp112/sigma;

    temp121 = (1+(mu*g_p*delta)/(r+lambda+theta(i)*q_theta))*...
        intF_ucpm/(r+lambda+mu*(1-g_p));
    temp122 = intF_cppm/(r+lambda);

    temp123 = intF_cpdpm/(r+lambda+theta(i)*q_theta);
    temp124 = intF_dpmp/(r+lambda+theta(i)*q_theta+mu*(1-g_p));

    temp125 = phi*(1-fun_F_x(epsilon_c(i),typen))-r-lambda-mu;
    temp126 = delta*sigma*(epsilon_c(i)-epsilon_dpms)/(r+lambda+theta(i)*q_theta)+...
        delta*sigma*(epsilon_dpms-epsilon_d(i))/(r+lambda+theta(i)*q_theta+mu*(1-g_p));

    temp127 = (1+(mu*g_p*delta)/(r+lambda+theta_pps*q_thetapps))*...
        (epsilon_c(i)-epsilon_cpps)/...
        (r+lambda+mu*(1-g_pps));
    temp128 = delta*(epsilon_cpps-epsilon_dpps)/...
        ((1-phi)*(r+lambda+theta_pps*q_thetapps));

    part12 = epsilon_c(i)+lambda*(temp121+temp122)+lambda*delta*(temp123+temp124)+...
        temp125*temp126/(sigma*(1-phi))+mu*g_p*(temp127+temp128);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fn = [Fnm1;part11 - part12];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i =n
for i = n
    Fnm1=Fn;
    
    p = fun_p(y(i));
    
    theta_pps=theta(i);%!!!!!!!!!!!!!!!
    % theta_{p^+,sigma}
    q_theta = fun_q_theta(theta(i));
    % q(theta_{p,sigma})
    q_thetapps=fun_q_theta(theta_pps);
    %q(theta_{p^+,sigma})

    g_p=fun_gp(y(i));
    %
    g_pps=fun_gp(y(i));
    %g_{p^+,sigma}

    intF_ucpm=fun_int_F(epsilon_c(i-1),epsilon_u,typen) ;
    % int(epsilon_u,epsilon_{p^-,sigma}^c)
    intF_cppm=fun_int_F(epsilon_c(i),epsilon_c(i-1),typen) ;
    % int(epsilon_{p^-,sigma}^c, epsilon_{p,sigma}^c)
    intF_cpdpm=fun_int_F(epsilon_d(i-1),epsilon_c(i),typen) ;
    % int(epsilon_{p,sigma}^c, epsilon_{p^-,sigma}^d)
    intF_dpmp=fun_int_F(epsilon_d(i),epsilon_d(i-1),typen) ;
    % int(epsilon_{p^-,sigma}^d, epsilon_{p,sigma}^d)
  

    epsilon_dpms=epsilon_d(i-1);
    % epsilon_{p^-,sigma}^d
    epsilon_cpps=epsilon_c(i);%~!!!!!!!!!!!!!!!!!!!!!
    % epsilon_{p^+,sigma}^c
    epsilon_dpps=epsilon_d(i);%~!!!!!!!!!!!!!!!!!!!!!!!!
    % epsilon_{p^+,sigma}^d
  

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   temp111 = ((beta+phi*(1-beta))*c_f*(1-var(i)))/((1-beta)*(1-phi))...
        +beta*c_p*var(i)/(1-beta);
    temp112 = b-p+theta(i)*(temp111);
    
    part11 = temp112/sigma;

    temp121 = (1+(mu*g_p*delta)/(r+lambda+theta(i)*q_theta))*...
        intF_ucpm/(r+lambda+mu*(1-g_p));
    temp122 = intF_cppm/(r+lambda);

    temp123 = intF_cpdpm/(r+lambda+theta(i)*q_theta);
    temp124 = intF_dpmp/(r+lambda+theta(i)*q_theta+mu*(1-g_p));

    temp125 = phi*(1-fun_F_x(epsilon_c(i),typen))-r-lambda-mu;
    temp126 = delta*sigma*(epsilon_c(i)-epsilon_dpms)/(r+lambda+theta(i)*q_theta)+...
        delta*sigma*(epsilon_dpms-epsilon_d(i))/(r+lambda+theta(i)*q_theta+mu*(1-g_p));

    
    temp128 = delta*(epsilon_cpps-epsilon_dpps)/...
        ((1-phi)*(r+lambda+theta_pps*q_thetapps));

    part12 = epsilon_c(i)+lambda*(temp121+temp122)+lambda*delta*(temp123+temp124)+...
        temp125*temp126/(sigma*(1-phi))+mu*g_p*(temp128);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fn = [Fnm1;part11 - part12];
end
return

