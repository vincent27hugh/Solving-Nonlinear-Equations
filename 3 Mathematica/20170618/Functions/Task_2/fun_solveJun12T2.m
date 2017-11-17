% Task1-2, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function Fn = fun_solveJun12T2(var,parameters)
% var vetcor represent 
% [epsilond_ps,SbarP_ps,theta_ps]; 

% length(y)*3 variables
% Parameters; global 
dbstop if error
% parameters = [c_f c_p beta phi delta sigma lambda b r epsilon_u mu step
% width typen]
A = parameters{1};
B1= parameters{2};
B2 = parameters{3};
c_f = parameters{4};
c_p = parameters{5};
beta = parameters{6};
phi = parameters{7};
delta = parameters{8};
sigma = parameters{9};
lambda = parameters{10};
b = parameters{11};
r = parameters{12};
epsilon_u = parameters{13};
mu = parameters{14};
step = parameters{15};
width = parameters{16};
pstar = parameters{17};
typen = parameters{18};

alpha_ps = parameters{19};
% alpha_{p,sigma}
casecode = parameters{20};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = int16(2*width/step+1);
% length of the lattice
y = (-width:step:width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch casecode
    case 'pstar' % 1
        pstar=var(13);
    case 'b' % 2
        b=var(13);
    case 'phi' % 3
        phi=var(13);
    case 'sigma' % 4
        sigma=var(13);
    case 'beta' % 5
        beta=var(13);
    case 'lambda' % 6
        lambda=var(13);
    case 'c_f' % 7
        c_f=var(13);
    case 'r' % 8
        r=var(13);
    case 'c_p' % 9
        c_p=var(13);
    case 'delta' % 10
        delta = var(13);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1
% p
pm = fun_p(y(1),pstar,width);
p = fun_p(y(2),pstar,width);
pp = fun_p(y(3),pstar,width);

% g_p
g_p=fun_gp(y(2),width);
% g_{p,sigma}

epsilond_pms=var(1);
% epsilon^d_{p^-,sigma}
epsilond_ps=var(3);
% epsilon^d_{p,sigma}
epsilond_pps=var(5);
% epsilon^d_{p^+,sigma}

theta_pms=var(2);
% theta_{p^-,sigma}
theta_ps=var(4);
% theta_{p,sigma}
theta_pps=var(6);
% theta_{p^+,sigma}

epsilonc_pms=var(7);
% epsilon^c_{p^-,sigma}
epsilonc_ps=var(12);
% epsilon^c_{p,sigma}
epsilonc_pps=var(17);
% epsilon^c_{p^+,sigma}

alpha_pms = var(8);
% alpha_{p^-,sigma}

alpha_pps = var(18);
% alpha_{p^+,sigma}

SP_pms_epc_pms = var(9);
% S^P_{p^-,sigma}(epsilon^c_{p^-,sigma})
SP_pms_epc_ps = var(10);
% S^P_{p^-,sigma}(epsilon^c_{p,sigma})
SP_pms_epc_pps = var(11);
% S^P_{p^-,sigma}(epsilon^c_{p^+,sigma})

SP_ps_epc_pms = var(14);
% S^P_{p,sigma}(epsilon^c_{p^-,sigma})
SP_ps_epc_ps = var(15);
% S^P_{p,sigma}(epsilon^c_{p,sigma})
SP_ps_epc_pps = var(16);
% S^P_{p,sigma}(epsilon^c_{p^+,sigma})

SP_pps_epc_pms = var(19);
% S^P_{p^+,sigma}(epsilon^c_{p^-,sigma})
SP_pps_epc_ps = var(20);
% S^P_{p^+,sigma}(epsilon^c_{p,sigma})
SP_pps_epc_pps = var(21);
% S^P_{p^+,sigma}(epsilon^c_{p^+,sigma})

SF_ps_epc_pms = var(22);
% S^F_{p,sigma}(epsilon^c_{p^-,sigma})
SF_pps_epc_pms = var(23);
% S^F_{p^+,sigma}(epsilon^c_{p^-,sigma})
SF_pps_epc_ps = var(24);
% S^F_{p^+,sigma}(epsilon^c_{p,sigma})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_theta_pms=fun_q_theta(theta_pms,...
A,B1,B2);
%q(theta_{p^-,sigma})
q_theta_ps = fun_q_theta(theta_ps,...
A,B1,B2);
% q(theta_{p,sigma})
q_theta_pps=fun_q_theta(theta_pps,...
A,B1,B2);
%q(theta_{p^+,sigma})

tau_pms = fun_tau(theta_pms,q_theta_pms,r,lambda,mu);
% tau_{p^-,sigma}
tau_ps = fun_tau(theta_ps,q_theta_ps,r,lambda,mu);
% tau_{p,sigma}
tau_pps = fun_tau(theta_pps,q_theta_pps,r,lambda,mu);
% tau_{p^+,sigma}

SbarP_pms=c_p/((1-beta)*q_theta_pms);
% Sbar^P_{p^-,sigma}
SbarP_ps=c_p/((1-beta)*q_theta_ps);
% Sbar^P_{p,sigma}
SbarP_pps=c_p/((1-beta)*q_theta_pps);
% Sbar^P_{p^+,sigma}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SbarF_pms = c_f/((1-beta)*(1-phi)*q_theta_pms);
% Sbar^F_{p^-,sigma}
SbarF_ps = c_f/((1-beta)*(1-phi)*q_theta_ps);
% Sbar^F_{p^-,sigma}
SbarF_pps = c_f/((1-beta)*(1-phi)*q_theta_pps);
% Sbar^F_{p^-,sigma}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intF_dp_u1=fun_int_F_P(epsilond_pms,epsilond_ps,epsilond_pps,...
tau_pms,tau_ps,tau_pps,g_p,mu,epsilon_u,typen,1);
% int(epsilon^d_{p,sigma},epsilon_u)

intF_dp_u2=fun_int_F_P(epsilond_pms,epsilond_ps,epsilond_pps,...
tau_pms,tau_ps,tau_pps,g_p,mu,epsilon_u,typen,2);
% int(epsilon^d_{p,sigma},epsilon_u)

intF_dp_u3=fun_int_F_P(epsilond_pms,epsilond_ps,epsilond_pps,...
tau_pms,tau_ps,tau_pps,g_p,mu,epsilon_u,typen,3);
% int(epsilon^d_{p,sigma},epsilon_u)

intF_dp_u4=fun_int_F_F(epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta,typen,epsilon_u,1);
% int(epsilon^d_{p,sigma},epsilon_u)

intF_dp_u5=fun_int_F_F(epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta,typen,epsilon_u,2);
% int(epsilon^d_{p,sigma},epsilon_u)

intF_dp_u6=fun_int_F_F(epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta,typen,epsilon_u,3);
% int(epsilon^d_{p,sigma},epsilon_u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELTAP_ps_dpms...
= fun_DELTAP_ps(epsilond_pms,epsilond_pms,epsilond_ps,epsilond_pps,...
tau_pms,tau_ps,tau_pps,g_p,mu);

temp111 = epsilond_pms;
temp112 = lambda*intF_dp_u1+mu*DELTAP_ps_dpms*(epsilond_pms-epsilond_ps);

part11 = temp111+temp112;

part12 = (b-pm)/sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part21 = tau_pps*SbarP_pps;

temp221 = delta*sigma*(epsilon_u-epsilond_pps);

temp222 = mu*SbarP_ps;

part22 = temp221+temp222;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
DELTAP_pps_dps...
= fun_DELTAP_pps(epsilond_ps,epsilond_pms,epsilond_ps,epsilond_pps,...
tau_pms,tau_ps,tau_pps,g_p,mu);

temp311 = epsilond_ps;
temp312 = lambda*intF_dp_u2+mu*g_p*DELTAP_pps_dps*(epsilond_ps-epsilond_pps);

part31 = temp311+temp312;

part32 = (b-p)/sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part41 = tau_ps*SbarP_ps;

temp421 = delta*sigma*(epsilon_u-epsilond_ps);
temp4221 = SbarP_pps-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps;
temp4222 = SbarP_pms;
temp422 = mu*(g_p*temp4221+(1-g_p)*temp4222);

part42 = temp421+temp422;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp511 = epsilond_pps;
temp512 = lambda*intF_dp_u3;

part51 = temp511+temp512;

part52 = (b-pp)/sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part61 = tau_pms*SbarP_pms;

temp621 = delta*sigma*(epsilon_u-epsilond_pms);
temp6221 = SbarP_ps;
temp6222 = -(1+mu*g_p/tau_pps)*delta*sigma*(epsilond_pms-epsilond_ps)...
    /(tau_ps-mu^2*g_p/tau_pps);
temp622 = mu*(temp6221+temp6222);

part62 = temp621+temp622;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELTAF_ps_cpms...
= fun_DELTAF_ps(epsilonc_pms,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);

F_epc_pms = fun_F_x (epsilonc_pms,typen,epsilon_u);

temp711 = (beta+phi*(1-beta))*c_f*(1-alpha_pms)/((1-beta)*(1-phi));
temp712 = beta*c_p*alpha_pms/(1-beta);

temp713 = lambda+(r+lambda*phi*F_epc_pms+mu)/(1-phi);

part71 = theta_pms*(temp711+temp712)+temp713*SP_pms_epc_pms;

temp721 = -b+pm;
temp722 = epsilonc_pms+lambda*intF_dp_u4;
temp723 = mu*DELTAF_ps_cpms*(epsilonc_pms-epsilonc_ps);

part72 = temp721 + sigma*(temp722+temp723);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELTAF_pps_cps...
= fun_DELTAF_pps(epsilonc_ps,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);

F_epc_ps = fun_F_x(epsilonc_ps,typen,epsilon_u);

temp811 = (beta+phi*(1-beta))*c_f*(1-alpha_ps)/((1-beta)*(1-phi));
temp812 = beta*c_p*alpha_ps/(1-beta);

temp813 = lambda+(r+lambda*phi*F_epc_ps+mu)/(1-phi);

part81 = theta_ps*(temp811+temp812)+temp813*SP_ps_epc_ps;

temp821 = -b+p;
temp822 = epsilonc_ps+lambda*intF_dp_u5;
temp823 = mu*g_p*DELTAF_pps_cps*(epsilonc_ps-epsilonc_pps);

part82 = temp821 + sigma*(temp822+temp823);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_epc_pps = fun_F_x(epsilonc_pps,typen,epsilon_u);

temp911 = (beta+phi*(1-beta))*c_f*(1-alpha_pps)/((1-beta)*(1-phi));
temp912 = beta*c_p*alpha_pps/(1-beta);

temp913 = lambda+(r+lambda*phi*F_epc_pps+mu)/(1-phi);

part91 = theta_pps*(temp911+temp912)+temp913*SP_pps_epc_pps;

temp921 = -b+pp;
temp922 = epsilonc_pps+lambda*intF_dp_u6;

part92 = temp921 + sigma*(temp922);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1001 = (r+lambda+mu)*SbarF_pps;

temp100201 = sigma*(epsilon_u-epsilonc_pps);
temp100202 = (r+lambda+mu)*SP_pps_epc_pps/(1-phi);
temp100203 = SbarF_ps-SP_ps_epc_pps;

part1002 = temp100201+temp100202+mu*temp100203;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1101 = (r+lambda+mu)*SbarF_ps;

temp110201 = sigma*(epsilon_u-epsilonc_ps);
temp110202 = (r+lambda+mu)*SP_ps_epc_ps/(1-phi);
temp110203 = g_p*(SbarF_pps-SF_pps_epc_ps)+(1-g_p)*(SbarF_pms-SP_pms_epc_ps);

part1102 = temp110201+temp110202+mu*temp110203;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1201 = (r+lambda+mu)*SbarF_pms;

temp120201 = sigma*(epsilon_u-epsilonc_pms);
temp120202 = (r+lambda+mu)*SP_pms_epc_pms/(1-phi);
temp120203 = SbarF_ps-SF_ps_epc_pms;

part1202 = temp120201+temp120202+mu*temp120203;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1301 = (r+lambda+mu)*SF_pps_epc_ps;

temp130201 = sigma*(epsilonc_ps-epsilonc_pps);
temp130202 = (r+lambda+mu)*SP_pps_epc_pps/(1-phi);
temp130203 = SP_ps_epc_ps/(1-phi)-SP_ps_epc_pps;

part1302 = temp130201+temp130202+mu*temp130203;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1401 = (r+lambda+mu)*SF_pps_epc_pms;

temp140201 = sigma*(epsilonc_pms-epsilonc_pps);
temp140202 = (r+lambda+mu)*SP_pps_epc_pps/(1-phi);
temp140203 = SF_ps_epc_pms-SP_ps_epc_pps;

part1402 = temp140201+temp140202+mu*temp140203;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1501 = (r+lambda+mu)*SF_ps_epc_pms;

temp150201 = sigma*(epsilonc_pms-epsilonc_ps);
temp150202 = (r+lambda+mu)*SP_ps_epc_ps/(1-phi);
temp150203 = g_p*(SF_pps_epc_pms - SF_pps_epc_ps)+...
    (1-g_p)*(SP_pms_epc_pms/(1-phi)-SP_pms_epc_ps);

part1502 = temp150201+temp150202+mu*temp150203;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1601 = tau_pps*SP_pps_epc_pps;

temp160201 = delta*sigma*(epsilonc_pps-epsilond_pps);
temp160202 = SP_ps_epc_pps;

part1602 = temp160201++mu*temp160202;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1701 = tau_ps*SP_ps_epc_ps;

temp170201 = delta*sigma*(epsilonc_ps-epsilond_ps);
temp170202 = g_p*(SP_pps_epc_ps-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps)+...
    (1-g_p)*(SP_pms_epc_ps);

part1702 = temp170201++mu*temp170202;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1801 = tau_pms*SP_pms_epc_pms;

temp180201 = delta*sigma*(epsilonc_pms-epsilond_pms);
temp180202 = SP_ps_epc_pms - ...
    (1+mu*g_p/tau_pps)*delta*sigma*...
    (epsilond_pms-epsilond_ps)/(tau_ps-mu^2*g_p/tau_pps);

part1802 = temp180201++mu*temp180202;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part1901 = tau_ps*SP_ps_epc_pps;

temp190201 = fun_indicator(epsilonc_pps,epsilond_ps);
temp190202 = delta*sigma*(epsilonc_pps-epsilond_ps);
temp190203 = g_p*(SP_pps_epc_pps-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps) + ...
    (1-g_p)*(SP_pms_epc_pps);

part1902 = temp190201*(temp190202+mu*temp190203);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part2001 = tau_pms*SP_pms_epc_pps;

temp200201 = fun_indicator(epsilonc_pps,epsilond_pms);
temp200202 = delta*sigma*(epsilonc_pps-epsilond_pms);
temp200203 = SP_ps_epc_pps-(1+mu*g_p/tau_pps)*delta*sigma*(epsilond_pms-epsilond_ps)/...
    (tau_ps-mu^2*g_p/tau_pps);

part2002 = temp200201*(temp200202+mu*temp200203);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part2101 = tau_pms*SP_pms_epc_ps;

temp210201 = fun_indicator(epsilonc_ps,epsilond_pms);
temp210202 = delta*sigma*(epsilonc_ps-epsilond_pms);
temp210203 = SP_ps_epc_ps-(1+mu*g_p/tau_pps)*delta*sigma*(epsilond_pms-epsilond_ps)/...
    (tau_ps-mu^2*g_p/tau_pps);

part2102 = temp210201*(temp210202+mu*temp210203);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part2201 = tau_pps*SP_pps_epc_ps;

temp220201 = delta*sigma*(epsilonc_ps-epsilond_pps);
temp220202 = SP_ps_epc_ps;

part2202 = temp220201++mu*temp220202;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part2301 = tau_pps*SP_pps_epc_pms;

temp230201 = delta*sigma*(epsilonc_pms-epsilond_pps);
temp230202 = SP_ps_epc_pms;

part2302 = temp230201++mu*temp230202;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part2401 = tau_ps*SP_ps_epc_pms;

temp240201 = delta*sigma*(epsilonc_pms-epsilond_ps);
temp240202 = g_p*(SP_pps_epc_pms-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps) + ...
    (1-g_p)*SP_pms_epc_pms;

part2402 = temp240201++mu*temp240202;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fn = [part11 - part12;part21 - part22;part31 - part32;...
    part41 - part42;part51 - part52;part61 - part62;...
    part71 - part72;part81 - part82;part91 - part92;...
    part1001 - part1002;part1101 - part1102;part1201 - part1202;...
    part1301 - part1302;part1401 - part1402;part1501 - part1502;...
    part1601 - part1602;part1701 - part1702;part1801 - part1802;...
    part1901 - part1902;part2001 - part2002;part2101 - part2102;...
    part2201 - part2202;part2301 - part2302;part2401 - part2402];

return

