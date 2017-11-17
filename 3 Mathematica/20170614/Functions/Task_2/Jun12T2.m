% Task1-2, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilond,epsilonc,theta,alpha,...
    SbarP,SbarF,SP_pms_epc,SP_ps_epc,SP_pps_epc,...
    SF_epc,ceq,exitflag,sol] = Jun12T2(parameters)
dbstop if error
%{ 
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;step;...
                width;pstar;typen};
%}
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
casecode = parameters{20};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = int16(2*width/step+1);
% length of the lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;step;...
                width;pstar;typen;alpha_ps;casecode};
% F1 = @(var)fun_solveJun12T2P1(var,para1);
% var = [epsilon_d;epsilon_c;theta]
% theta should be larger than 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch casecode
    case 'pstar' % 1
        sol0=pstar;
    case 'b' % 2
        sol0=b;
    case 'phi' % 3
        sol0=phi;
    case 'sigma' % 4
        sol0=sigma;
    case 'beta' % 5
        sol0=beta;
    case 'lambda' % 6
        sol0=lambda;
    case 'c_f' % 7
        sol0=c_f;
    case 'r' % 8
        sol0=r;
    case 'c_p' % 9
        sol0=c_p;
    case 'delta' % 10
        sol0=delta;
end
%%%%%%%%%%%%%%%%%%
t11 = linspace(-1,-10,n);
t21 = linspace(1,5,n);
var01=reshape([t11;t21],[],1);

t12 = linspace(-1,-10,n);
t22 = linspace(0.5,0.5,n);
t32 = linspace(5,5,n);
t42 = linspace(5,5,n);
t52 = linspace(5,5,n);
t12_t52=reshape([t12;t22;t32;t42;t52],[],1);
t62=linspace(5,5,n);
var02 = [t12_t52;t62'];

var02(7)=sol0;

var0=[var01;var02];

options = optimoptions(@fmincon,'Algorithm','interior-point' ,...
     'Display','iter',...
     'PlotFcn',@optimplotx,...
     'MaxFunctionEvaluations',Inf,...
     'MaxIter',250,...
     'ConstraintTolerance',1e-12,...
     'OptimalityTolerance',1e-12,...
     'StepTolerance',1e-12);
%   'interior-point' (default)
%   'trust-region-reflective'
%   A gradient to be supplied in the objective function
%   'sqp'
%   'active-set'
Afmin=[];bfmin=[];Aeqfmin=[];beqfmin=[];
%%%%
t1=-Inf(1,n);
t3=-Inf(1,n);
lbfmin1 = reshape([t1;t3],[],1);

%%%%
t1=-Inf(1,n);
t2=-Inf(1,n);
t3=-Inf(1,n);
t4=-Inf(1,n);
t5=-Inf(1,n);
t1_t5 = reshape([t1;t2;t3;t4;t5],[],1);
t6 = -Inf(1,n);
lbfmin2 = [t1_t5;t6'];

lbfmin=[lbfmin1;lbfmin2];
%%%%
t1=repmat(epsilon_u,1,n);
% for epsilon_d
t3=Inf(1,n);
% for theta
ubfmin1 = reshape([t1;t3],[],1);
%%%%
t1=repmat(epsilon_u,1,n);
% for epsilon_c
t2=Inf(1,n);
% for alpha
t3=Inf(1,n);
t4=Inf(1,n);
t5=Inf(1,n);
t1_t5 = reshape([t1;t2;t3;t4;t5],[],1);
t6 = Inf(1,n);
ubfmin2 = [t1_t5;t6'];

ubfmin=[ubfmin1;ubfmin2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,...
    Aeqfmin,beqfmin,lbfmin,ubfmin,...
    @(x)fun_solveJun12T2_nonlcon(x,para),options);
% Find minimum of constrained nonlinear multivariable function
[~,ceq]=fun_solveJun12T2_nonlcon(fs,para);

epsilond=fs(1:2:2*n);
theta=fs(2:2:2*n);
%%%
q_theta_pms=fun_q_theta(theta(1),...
A,B1,B2);
%q(theta_{p^-,sigma})
q_theta_ps = fun_q_theta(theta(2),...
A,B1,B2);
% q(theta_{p,sigma})
q_theta_pps=fun_q_theta(theta(3),...
A,B1,B2);
%q(theta_{p^+,sigma})

SbarP_pms=c_p/((1-beta)*q_theta_pms);
% Sbar^P_{p^-,sigma}
SbarP_ps=c_p/((1-beta)*q_theta_ps);
% Sbar^P_{p,sigma}
SbarP_pps=c_p/((1-beta)*q_theta_pps);
% Sbar^P_{p^+,sigma}

SbarP = [SbarP_pms;SbarP_ps;SbarP_pps];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_theta_pms=fun_q_theta(theta(1),...
A,B1,B2);
%q(theta_{p^-,sigma})
q_theta_ps = fun_q_theta(theta(2),...
A,B1,B2);
% q(theta_{p,sigma})
q_theta_pps=fun_q_theta(theta(3),...
A,B1,B2);
%q(theta_{p^+,sigma})

SbarF_pms = c_f/((1-beta)*(1-phi)*q_theta_pms);
% Sbar^F_{p^-,sigma}
SbarF_ps = c_f/((1-beta)*(1-phi)*q_theta_ps);
% Sbar^F_{p^-,sigma}
SbarF_pps = c_f/((1-beta)*(1-phi)*q_theta_pps);
% Sbar^F_{p^-,sigma}

SbarF = [SbarF_pms,SbarF_ps,SbarF_pps];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilonc=fs(7:5:17);
alpha=fs(8:10:18);

sol = fs(13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP_pms_epc_pms=fs(9);
SP_pms_epc_ps=fs(10);
SP_pms_epc_pps=fs(11);

SP_pms_epc = [SP_pms_epc_pms,SP_pms_epc_ps,SP_pms_epc_pps];

SP_ps_epc_pms=fs(14);
SP_ps_epc_ps=fs(15);
SP_ps_epc_pps=fs(16);

SP_ps_epc = [SP_ps_epc_pms,SP_ps_epc_ps,SP_ps_epc_pps];

SP_pps_epc_pms=fs(19);
SP_pps_epc_ps=fs(20);
SP_pps_epc_pps=fs(21);

SP_pps_epc = [SP_pps_epc_pms,SP_pps_epc_ps,SP_pps_epc_pps];

SF_ps_epc_pms=fs(22);
SF_pps_epc_pms=fs(23);
SF_pps_epc_ps=fs(24);

SF_epc = [SF_ps_epc_pms,SF_pps_epc_pms,SF_pps_epc_ps];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return