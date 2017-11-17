% Task1-2, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function Fn = fun_solveMay22T2P1(var,parameters)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = int16(2*width/step+1);
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

temp411 = epsilond_ps;
temp412 = lambda*intF_dp_u2+mu*g_p*DELTAP_pps_dps*(epsilond_ps-epsilond_pps);

part41 = temp411+temp412;

part42 = (b-p)/sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part51 = tau_ps*SbarP_ps;

temp521 = delta*sigma*(epsilon_u-epsilond_ps);
temp5221 = SbarP_pps-delta*sigma*(epsilond_ps-epsilond_pps)/tau_pps;
temp5222 = SbarP_pms;
temp522 = mu*(g_p*temp5221+(1-g_p)*temp5222);

part52 = temp521+temp522;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp711 = epsilond_pps;
temp712 = lambda*intF_dp_u3;

part71 = temp711+temp712;

part72 = (b-pp)/sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part81 = tau_pms*SbarP_pms;

temp821 = delta*sigma*(epsilon_u-epsilond_pms);
temp8221 = SbarP_ps;
temp8222 = -(1+mu*g_p/tau_pps)*delta*sigma*(epsilond_pms-epsilond_ps)...
    /(tau_ps-mu^2*g_p/tau_pps);
temp822 = mu*(temp8221+temp8222);

part82 = temp821+temp822;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fn = [part11 - part12;part41 - part42;part71 - part72;...
    part21 - part22;part51 - part52;part81 - part82];

return

