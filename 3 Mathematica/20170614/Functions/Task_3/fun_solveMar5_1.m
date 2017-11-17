% Task #20170302
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function Fn = fun_solveMar5_1(var)
% var vetcor represent 
% [epsilon_dps;epsilon_cps;theta_ps]; 
% length(y)*3 variables
% Parameters; global 
global c_f c_p beta phi delta sigma lambda b r epsilon_u mu 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global step width typen
n = int16(2*width/step+1);
% length of the lattice
y = (-width:step:width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1
for i =1
    p = fun_p(y(i));
    theta_pps=var(3*i+3);
    % theta_{p^+,sigma}
    q_theta = fun_q_theta(var(3*i));
    % q(theta_{p,sigma})
    q_thetapps=fun_q_theta(theta_pps);
    %q(theta_{p^+,sigma})

    g_p=fun_gp(y(i));
    %
    g_pps=fun_gp(y(i+1));
    %g_{p^+,sigma}
    
    epsilon_dpms=var(3*i-2);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % epsilon_{p^-,sigma}^d
   
    epsilon_dpps=var(3*i+1);
    % epsilon_{p^+,sigma}^d
    epsilon_cpms=var(3*i-1);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % epsilon_{p^-,sigma}^c
    
    intF_udpm=fun_int_F(epsilon_dpms,epsilon_u,typen);%!!!!!!!!!!!!!!!!!!!
    % int(epsilon_u, epsilon_{p^-,sigma}^d)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % var = [epsilon_d;epsilon_c;theta]
    part21 = (b-p)/sigma;

    temp221 =  intF_udpm/(r+lambda+var(3*i)*q_theta);
    
    temp223 = (var(3*i-2)-epsilon_dpps)/(r+lambda+theta_pps*q_thetapps+mu*(1-g_pps));

    part22 = var(3*i-2)+lambda*(temp221)+mu*g_p*temp223;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % var = [epsilon_d;epsilon_c;theta;alpha]
    temp311 = sigma*(epsilon_u-epsilon_cpms)/(r+lambda);
    
    temp313 = delta*sigma*(var(3*i-1)-var(3*i-2))/((1-phi)*(r+lambda+var(3*i)*q_theta));
    part31 = temp311+temp313;

    part32 = c_f/((1-beta)*(1-phi)*q_theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % var = [epsilon_d;epsilon_c;theta;alpha]
    temp411 = delta*sigma*(epsilon_u-epsilon_dpms)/...
        (r+lambda+var(3*i)*q_theta);
   
    part41 = temp411;

    part42 = c_p/((1-beta)*q_theta);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Fn = [part21 - part22;part31 - part32;part41-part42];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:n-1
    Fnm1=Fn;
    
    p = fun_p(y(i));
    
    theta_pps=var(3*i+3);
    % theta_{p^+,sigma}
    q_theta = fun_q_theta(var(3*i));
    % q(theta_{p,sigma})
    q_thetapps=fun_q_theta(theta_pps);
    %q(theta_{p^+,sigma})

    g_p=fun_gp(y(i));
    %
    g_pps=fun_gp(y(i+1));
    %g_{p^+,sigma}

    epsilon_dpms=var(3*i-5);
    % epsilon_{p^-,sigma}^d
    epsilon_dpps=var(3*i+1);
    % epsilon_{p^+,sigma}^d
    epsilon_cpms=var(3*i-4);
    % epsilon_{p^-,sigma}^c
    
    intF_dpmp=fun_int_F(var(3*i-2),epsilon_dpms,typen) ;
    % int(epsilon_{p^-,sigma}^d, epsilon_{p,sigma}^d)
    intF_udpm=fun_int_F(epsilon_dpms,epsilon_u,typen);
    % int(epsilon_u, epsilon_{p^-,sigma}^d)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % var = [epsilon_d;epsilon_c;theta]
    part21 = (b-p)/sigma;

    temp221 =  intF_udpm/(r+lambda+var(3*i)*q_theta);
    temp222 = intF_dpmp/(r+lambda+var(3*i)*q_theta+mu*(1-g_p));
    temp223 = (var(3*i-2)-epsilon_dpps)/(r+lambda+theta_pps*q_thetapps+mu*(1-g_pps));

    part22 = var(3*i-2)+lambda*(temp221+temp222)+mu*g_p*temp223;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % var = [epsilon_d;epsilon_c;theta;alpha]
    temp311 = sigma*(epsilon_u-epsilon_cpms)/(r+lambda);
    
    temp3121 = sigma+(mu*g_p*delta*sigma)/(r+lambda+var(3*i)*q_theta);
    temp3122 = epsilon_cpms-var(3*i-1);
    temp312 = (temp3121*temp3122)/(r+lambda+mu*(1-g_p));
    
    temp313 = delta*sigma*(var(3*i-1)-var(3*i-2))/((1-phi)*(r+lambda+var(3*i)*q_theta));
    part31 = temp311+temp312+temp313;

    part32 = c_f/((1-beta)*(1-phi)*q_theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % var = [epsilon_d;epsilon_c;theta;alpha]
    temp411 = delta*sigma*(epsilon_u-epsilon_dpms)/...
        (r+lambda+var(3*i)*q_theta);
    temp412 = delta*sigma*(epsilon_dpms-var(3*i-2))/...
        (r+lambda+var(3*i)*q_theta+mu*(1-g_p));
    part41 = temp411+temp412;

    part42 = c_p/((1-beta)*q_theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fn = [Fnm1;part21 - part22;part31 - part32;part41-part42];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i =n
for i = n
    Fnm1=Fn;
    
    p = fun_p(y(i));
    
    q_theta = fun_q_theta(var(3*i));
    % q(theta_{p,sigma})
    
    g_p=fun_gp(y(i));
    %

    intF_dpmp=fun_int_F(var(3*i-2),var(3*i-5),typen) ;
    % int(epsilon_{p^-,sigma}^d, epsilon_{p,sigma}^d)
    intF_udpm=fun_int_F(var(3*i-5),epsilon_u,typen);
    % int(epsilon_u, epsilon_{p^-,sigma}^d)

    epsilon_dpms=var(3*i-5);
    % epsilon_{p^-,sigma}^d
   
    epsilon_cpms=var(3*i-4);
    % epsilon_{p^-,sigma}^c

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % var = [epsilon_d;epsilon_c;theta]
    part21 = (b-p)/sigma;

    temp221 =  intF_udpm/(r+lambda+var(3*i)*q_theta);
    temp222 = intF_dpmp/(r+lambda+var(3*i)*q_theta+mu*(1-g_p));
   

    part22 = var(3*i-2)+lambda*(temp221+temp222);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % var = [epsilon_d;epsilon_c;theta;alpha]
    temp311 = sigma*(epsilon_u-epsilon_cpms)/(r+lambda);
    
    temp3121 = sigma+(mu*g_p*delta*sigma)/(r+lambda+var(3*i)*q_theta);
    temp3122 = epsilon_cpms-var(3*i-1);
    temp312 = (temp3121*temp3122)/(r+lambda+mu*(1-g_p));
    
    temp313 = delta*sigma*(var(3*i-1)-var(3*i-2))/((1-phi)*(r+lambda+var(3*i)*q_theta));
    part31 = temp311+temp312+temp313;

    part32 = c_f/((1-beta)*(1-phi)*q_theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % var = [epsilon_d;epsilon_c;theta;alpha]
    temp411 = delta*sigma*(epsilon_u-epsilon_dpms)/...
        (r+lambda+var(3*i)*q_theta);
    temp412 = delta*sigma*(epsilon_dpms-var(3*i-2))/...
        (r+lambda+var(3*i)*q_theta+mu*(1-g_p));
    part41 = temp411+temp412;

    part42 = c_p/((1-beta)*q_theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fn = [Fnm1;part21 - part22;part31 - part32;part41-part42];
end
return

