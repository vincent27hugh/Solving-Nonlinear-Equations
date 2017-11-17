% Task1-2, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function Fn = fun_solveMay22T2P2(var,parameters)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilond = parameters{19};
theta = parameters{20};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = int16(2*width/step+1);
% length of the lattice
y = (-width:step:width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1
% p
pm = fun_p(y(1),pstar,width);
p = fun_p(y(2),pstar,width);
pp = fun_p(y(3),pstar,width);

% g_p
g_p=fun_gp(y(2),width);
% g_{p,sigma}

epsilond_pms=epsilond(1);
% epsilon^d_{p^-,sigma}
epsilond_ps=epsilond(2);
% epsilon^d_{p,sigma}
epsilond_pps=epsilond(3);
% epsilon^d_{p^+,sigma}

theta_pms=theta(1);
% theta_{p^-,sigma}
theta_ps=theta(2);
% theta_{p,sigma}
theta_pps=theta(3);
% theta_{p^+,sigma}

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilonc_pms=var(1);
% epsilon^c_{p^-,sigma}
epsilonc_ps=var(6);
% epsilon^c_{p,sigma}
epsilonc_pps=var(11);
% epsilon^c_{p^+,sigma}

alpha_pms = var(2);
% alpha_{p^-,sigma}
alpha_ps = var(7);
% alpha_{p,sigma}
alpha_pps = var(12);
% alpha_{p^+,sigma}

SP_pms_epc_pms = var(3);
% S^P_{p^-,sigma}(epsilon^c_{p^-,sigma})
SP_pms_epc_ps = var(4);
% S^P_{p^-,sigma}(epsilon^c_{p,sigma})
SP_pms_epc_pps = var(5);
% S^P_{p^-,sigma}(epsilon^c_{p^+,sigma})

SP_ps_epc_pms = var(8);
% S^P_{p,sigma}(epsilon^c_{p^-,sigma})
SP_ps_epc_ps = var(9);
% S^P_{p,sigma}(epsilon^c_{p,sigma})
SP_ps_epc_pps = var(10);
% S^P_{p,sigma}(epsilon^c_{p^+,sigma})

SP_pps_epc_pms = var(13);
% S^P_{p^+,sigma}(epsilon^c_{p^-,sigma})
SP_pps_epc_ps = var(14);
% S^P_{p^+,sigma}(epsilon^c_{p,sigma})
SP_pps_epc_pps = var(15);
% S^P_{p^+,sigma}(epsilon^c_{p^+,sigma})

SF_ps_epc_pms = var(16);
% S^F_{p,sigma}(epsilon^c_{p^-,sigma})
SF_pps_epc_pms = var(17);
% S^F_{p^+,sigma}(epsilon^c_{p^-,sigma})
SF_pps_epc_ps = var(18);
% S^F_{p^+,sigma}(epsilon^c_{p,sigma})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SbarF_pms = c_f/((1-beta)*(1-phi)*q_theta_pms);
% Sbar^F_{p^-,sigma}
SbarF_ps = c_f/((1-beta)*(1-phi)*q_theta_ps);
% Sbar^F_{p^-,sigma}
SbarF_pps = c_f/((1-beta)*(1-phi)*q_theta_pps);
% Sbar^F_{p^-,sigma}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intF_dp_u1=fun_int_F_F(epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta,typen,epsilon_u,1);
% int(epsilon^d_{p,sigma},epsilon_u)

intF_dp_u2=fun_int_F_F(epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta,typen,epsilon_u,2);
% int(epsilon^d_{p,sigma},epsilon_u)

intF_dp_u3=fun_int_F_F(epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta,typen,epsilon_u,3);
% int(epsilon^d_{p,sigma},epsilon_u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_DELTAP='off';
if strcmp(plot_DELTAP,'on')
    f1 = figure;
    t=linspace(-100,epsilon_u,100);

    s1 = subplot(2,2,1);
    hold on
    t_DELTAP_pms= fun_DELTAP_pms(t,epsilond_pms,epsilond_ps,...
        epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    plot(t,t_DELTAP_pms,'LineWidth',2);
    plot(epsilond_pms,0,'or',epsilond_ps,0,'or',epsilond_pps,0,'or');
    xlabel('$\epsilon$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    ylabel('$\Delta_{p^{-},\sigma}(\epsilon)$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    xlim([min(t),max(t)]);
    hold off

    s2 = subplot(2,2,2);
    hold on
    t_DELTAP_ps= fun_DELTAP_ps(t,epsilond_pms,epsilond_ps,...
        epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    plot(t,t_DELTAP_ps,'LineWidth',2);
    plot(epsilond_pms,0,'or',epsilond_ps,0,'or',epsilond_pps,0,'or');
    xlabel('$\epsilon$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    ylabel('$\Delta_{p,\sigma}(\epsilon)$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    hold off
    xlim([min(t),max(t)]);

    s3 = subplot(2,2,3);
    hold on
    t_DELTAP_pps= fun_DELTAP_pps(t,epsilond_pms,epsilond_ps,...
        epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    plot(t,t_DELTAP_pps,'LineWidth',2);
    plot(epsilond_pms,0,'or',epsilond_ps,0,'or',epsilond_pps,0,'or');
    xlabel('$\epsilon$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    ylabel('$\Delta_{p^{+},\sigma}(\epsilon)$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    hold off
    xlim([min(t),max(t)]);

    grid(s1,'on');grid(s2,'on');grid(s3,'on');
    set(gcf,'Position',get(0,'Screensize'));
    pause(0.5);
    close(f1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_int='off';
if strcmp(plot_int,'on')
    f1 = figure;
    t=linspace(-100,epsilon_u,100);
    temp =1-fun_F_x(t,typen,epsilon_u);
    
    s1 = subplot(2,2,1);
    hold on
    t_DELTAP_pms= fun_DELTAP_pms(t,epsilond_pms,epsilond_ps,...
        epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    
    plot(t,t_DELTAP_pms.*temp,'LineWidth',2);
    plot(epsilond_pms,0,'or',epsilond_ps,0,'or',epsilond_pps,0,'or');
    xlabel('$\epsilon$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    ylabel('$\Delta_{p^{-},\sigma}(\epsilon)*(1-F(\epsilon)$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    xlim([min(t),max(t)]);
    hold off

    s2 = subplot(2,2,2);
    hold on
    t_DELTAP_ps= fun_DELTAP_ps(t,epsilond_pms,epsilond_ps,...
        epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    plot(t,t_DELTAP_ps.*temp,'LineWidth',2);
    plot(epsilond_pms,0,'or',epsilond_ps,0,'or',epsilond_pps,0,'or');
    xlabel('$\epsilon$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    ylabel('$\Delta_{p,\sigma}(\epsilon)*(1-F(\epsilon)$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    hold off
    xlim([min(t),max(t)]);

    s3 = subplot(2,2,3);
    hold on
    t_DELTAP_pps= fun_DELTAP_pps(t,epsilond_pms,epsilond_ps,...
        epsilond_pps,tau_pms,tau_ps,tau_pps,g_p,mu);
    plot(t,t_DELTAP_pps.*temp,'LineWidth',2);
    plot(epsilond_pms,0,'or',epsilond_ps,0,'or',epsilond_pps,0,'or');
    xlabel('$\epsilon$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    ylabel('$\Delta_{p^{+},\sigma}(\epsilon)*(1-F(\epsilon)$','interpreter',...
        'Latex','FontSize',18,'FontWeight','bold');
    hold off
    xlim([min(t),max(t)]);

    grid(s1,'on');grid(s2,'on');grid(s3,'on');
    set(gcf,'Position',get(0,'Screensize'));
    pause(0.5);
    close(f1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELTAF_ps_cpms...
= fun_DELTAF_ps(epsilonc_pms,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);

F_epc_pms = fun_F_x (epsilonc_pms,typen,epsilon_u);

temp111 = (beta+phi*(1-beta))*c_f*(1-alpha_pms)/((1-beta)*(1-phi));
temp112 = beta*c_p*alpha_pms/(1-beta);

temp113 = lambda+(r+lambda*phi*F_epc_pms+mu)/(1-phi);

part11 = theta_pms*(temp111+temp112)+temp113*SP_pms_epc_pms;

temp121 = -b+pm;
temp122 = epsilonc_pms+lambda*intF_dp_u1;
temp123 = mu*DELTAF_ps_cpms*(epsilonc_pms-epsilonc_ps);

part12 = temp121 + sigma*(temp122+temp123);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELTAF_pps_cps...
= fun_DELTAF_pps(epsilonc_ps,epsilonc_pms,epsilonc_ps,epsilonc_pps,...
    epsilond_pms,epsilond_ps,epsilond_pps,...
    tau_pms,tau_ps,tau_pps,...
    r,lambda,g_p,mu,delta);

F_epc_ps = fun_F_x(epsilonc_ps,typen,epsilon_u);

temp211 = (beta+phi*(1-beta))*c_f*(1-alpha_ps)/((1-beta)*(1-phi));
temp212 = beta*c_p*alpha_ps/(1-beta);

temp213 = lambda+(r+lambda*phi*F_epc_ps+mu)/(1-phi);

part21 = theta_ps*(temp211+temp212)+temp213*SP_ps_epc_ps;

temp221 = -b+p;
temp222 = epsilonc_ps+lambda*intF_dp_u2;
temp223 = mu*g_p*DELTAF_pps_cps*(epsilonc_ps-epsilonc_pps);

part22 = temp221 + sigma*(temp222+temp223);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_epc_pps = fun_F_x(epsilonc_pps,typen,epsilon_u);

temp311 = (beta+phi*(1-beta))*c_f*(1-alpha_pps)/((1-beta)*(1-phi));
temp312 = beta*c_p*alpha_pps/(1-beta);

temp313 = lambda+(r+lambda*phi*F_epc_pps+mu)/(1-phi);

part31 = theta_pps*(temp311+temp312)+temp313*SP_pps_epc_pps;

temp321 = -b+pp;
temp322 = epsilonc_pps+lambda*intF_dp_u3;

part32 = temp321 + sigma*(temp322);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part41 = (r+lambda+mu)*SbarF_pps;

temp421 = sigma*(epsilon_u-epsilonc_pps);
temp422 = (r+lambda+mu)*SP_pps_epc_pps/(1-phi);
temp423 = SbarF_ps-SP_ps_epc_pps;

part42 = temp421+temp422+mu*temp423;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part51 = (r+lambda+mu)*SbarF_ps;

temp521 = sigma*(epsilon_u-epsilonc_ps);
temp522 = (r+lambda+mu)*SP_ps_epc_ps/(1-phi);
temp523 = g_p*(SbarF_pps-SF_pps_epc_ps)+(1-g_p)*(SbarF_pms-SP_pms_epc_ps);

part52 = temp521+temp522+mu*temp523;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part61 = (r+lambda+mu)*SbarF_pms;

temp621 = sigma*(epsilon_u-epsilonc_pms);
temp622 = (r+lambda+mu)*SP_pms_epc_pms/(1-phi);
temp623 = SbarF_ps-SF_ps_epc_pms;

part62 = temp621+temp622+mu*temp623;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part71 = (r+lambda+mu)*SF_pps_epc_ps;

temp721 = sigma*(epsilonc_ps-epsilonc_pps);
temp722 = (r+lambda+mu)*SP_pps_epc_pps/(1-phi);
temp723 = SP_ps_epc_ps/(1-phi)-SP_ps_epc_pps;

part72 = temp721+temp722+mu*temp723;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part81 = (r+lambda+mu)*SF_pps_epc_pms;

temp821 = sigma*(epsilonc_pms-epsilonc_pps);
temp822 = (r+lambda+mu)*SP_pps_epc_pps/(1-phi);
temp823 = SF_ps_epc_pms-SP_ps_epc_pps;

part82 = temp821+temp822+mu*temp823;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part91 = (r+lambda+mu)*SF_ps_epc_pms;

temp921 = sigma*(epsilonc_pms-epsilonc_ps);
temp922 = (r+lambda+mu)*SP_ps_epc_ps/(1-phi);
temp923 = g_p*(SF_pps_epc_pms - SF_pps_epc_ps)+...
    (1-g_p)*(SP_pms_epc_pms/(1-phi)-SP_pms_epc_ps);

part92 = temp921+temp922+mu*temp923;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part101 = tau_pps*SP_pps_epc_pps;

temp1021 = delta*sigma*(epsilonc_pps-epsilond_pps);
temp1022 = SP_ps_epc_pps;

part102 = temp1021++mu*temp1022;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part111 = tau_ps*SP_ps_epc_ps;

temp1121 = delta*sigma*(epsilonc_ps-epsilond_ps);
temp1122 = g_p*(SP_pps_epc_ps-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps)+...
    (1-g_p)*(SP_pms_epc_ps);

part112 = temp1121++mu*temp1122;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part121 = tau_pms*SP_pms_epc_pms;

temp1221 = delta*sigma*(epsilonc_pms-epsilond_pms);
temp1222 = SP_ps_epc_pms - ...
    (1+mu*g_p/tau_pps)*delta*sigma*...
    (epsilond_pms-epsilond_ps)/(tau_ps-mu^2*g_p/tau_pps);

part122 = temp1221++mu*temp1222;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part131 = tau_ps*SP_ps_epc_pps;

temp1321 = fun_indicator(epsilonc_pps,epsilond_ps);
temp1322 = delta*sigma*(epsilonc_pps-epsilond_ps);
temp1323 = g_p*(SP_pps_epc_pps-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps) + ...
    (1-g_p)*(SP_pms_epc_pps);

part132 = temp1321*(temp1322+mu*temp1323);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part141 = tau_pms*SP_pms_epc_pps;

temp1421 = fun_indicator(epsilonc_pps,epsilond_pms);
temp1422 = delta*sigma*(epsilonc_pps-epsilond_pms);
temp1423 = SP_ps_epc_pps-(1+mu*g_p/tau_pps)*delta*sigma*(epsilond_pms-epsilond_ps)/...
    (tau_ps-mu^2*g_p/tau_pps);

part142 = temp1421*(temp1422+mu*temp1423);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part151 = tau_pms*SP_pms_epc_ps;

temp1521 = fun_indicator(epsilonc_ps,epsilond_pms);
temp1522 = delta*sigma*(epsilonc_ps-epsilond_pms);
temp1523 = SP_ps_epc_ps-(1+mu*g_p/tau_pps)*delta*sigma*(epsilond_pms-epsilond_ps)/...
    (tau_ps-mu^2*g_p/tau_pps);

part152 = temp1521*(temp1522+mu*temp1523);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part161 = tau_pps*SP_pps_epc_ps;

temp1621 = delta*sigma*(epsilonc_ps-epsilond_pps);
temp1622 = SP_ps_epc_ps;

part162 = temp1621++mu*temp1622;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part171 = tau_pps*SP_pps_epc_pms;

temp1721 = delta*sigma*(epsilonc_pms-epsilond_pps);
temp1722 = SP_ps_epc_pms;

part172 = temp1721++mu*temp1722;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part181 = tau_ps*SP_ps_epc_pms;

temp1821 = delta*sigma*(epsilonc_pms-epsilond_ps);
temp1822 = g_p*(SP_pps_epc_pms-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps) + ...
    (1-g_p)*SP_pms_epc_pms;

part182 = temp1821++mu*temp1822;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fn = [part11 - part12;part21 - part22;part31 - part32;...
    part41 - part42;part51 - part52;part61 - part62;...
    part71 - part72;part81 - part82;part91 - part92;...
    part101 - part102;part111 - part112;part121 - part122;...
    part131 - part132;part141 - part142;part151 - part152;...
    part161 - part162;part171 - part172;part181 - part182];

return

