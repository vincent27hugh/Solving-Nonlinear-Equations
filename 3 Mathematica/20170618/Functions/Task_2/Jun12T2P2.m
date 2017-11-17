% Task1-2, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilonc,alpha,SbarF,SP_pms_epc,SP_ps_epc,SP_pps_epc,...
    SF_epc,ceqP2,exitflag] = Jun12T2P2(parameters)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilond = parameters{19};
theta = parameters{20};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = int16(2*width/step+1);
% length of the lattice
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;step;...
                width;pstar;typen;epsilond;theta};
F = @(var)fun_solveMay22T2P2(var,para);
% var = [epsilon_d;epsilon_c;theta]
% theta should be larger than 0
t1 = linspace(-1,-10,n);
t2 = linspace(0.5,0.5,n);
t3 = linspace(5,5,n);
t4 = linspace(5,5,n);
t5 = linspace(5,5,n);
t1_t5=reshape([t1;t2;t3;t4;t5],[],1);
t6=linspace(5,5,n);
var0 = [t1_t5;t6'];
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','iter',...
    'PlotFcn',@optimplotx,...
    'MaxFunctionEvaluations',Inf,...
     'MaxIteration',1500,...
     'FunctionTolerance',1e-12,...
     'StepTolerance',1e-12);
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt'
% Option to display output
% LM Algo is better

[fs,~,exitflag] = fsolve(F,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]
ceqP2=fun_solveMay22T2P2(fs,para);
%{
t1 = linspace(-1600,-1800,n);
t2 = linspace(50,50,n);
t3 = linspace(50,50,n);
t4 = linspace(50,50,n);
t5 = linspace(50,50,n);
t1_t5=reshape([t1;t2;t3;t4;t5],[],1);
t6=linspace(50,50,n);
var0 = [t1_t5;t6'];
options = optimoptions(@fmincon,'Algorithm','interior-point' ,...
    'Display','iter',...
     'PlotFcn',@optimplotx,...
     'MaxFunctionEvaluations',Inf,...
     'MaxIter',500,...
     'ConstraintTolerance',1e-12,...
     'OptimalityTolerance',1e-12,...
     'StepTolerance',1e-12);
%{
     'interior-point' (default)
'trust-region-reflective'
 % A gradient to be supplied in the objective function
'sqp'
'active-set'
     %}
Afmin=[];bfmin=[];Aeqfmin=[];beqfmin=[];
%%%%
t1=-Inf(1,n);
t2=-Inf(1,n);
t3=-Inf(1,n);
t4=-Inf(1,n);
t5=-Inf(1,n);
t1_t5 = reshape([t1;t2;t3;t4;t5],[],1);
t6 = -Inf(1,n);
lbfmin = [t1_t5;t6'];
%%%%
t1=repmat(epsilon_u,1,n);
t2=Inf(1,n);
t3=Inf(1,n);
t4=Inf(1,n);
t5=Inf(1,n);
t1_t5 = reshape([t1;t2;t3;t4;t5],[],1);
t6 = Inf(1,n);
ubfmin = [t1_t5;t6'];
%%%%
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;step;...
                width;pstar;typen;epsilond;theta};
[fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,...
    Aeqfmin,beqfmin,lbfmin,ubfmin,...
    @(x)fun_solveMay22T2P2_nonlcon(x,para),options);
% Find minimum of constrained nonlinear multivariable function
[~,ceqP2]=fun_solveMay22T2P2_nonlcon(fs,para);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
epsilonc=fs(1:5:5*n);
alpha=fs(2:5:5*n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP_pms_epc_pms=fs(3);
SP_pms_epc_ps=fs(4);
SP_pms_epc_pps=fs(5);

SP_pms_epc = [SP_pms_epc_pms,SP_pms_epc_ps,SP_pms_epc_pps];

SP_ps_epc_pms=fs(8);
SP_ps_epc_ps=fs(9);
SP_ps_epc_pps=fs(10);

SP_ps_epc = [SP_ps_epc_pms,SP_ps_epc_ps,SP_ps_epc_pps];

SP_pps_epc_pms=fs(13);
SP_pps_epc_ps=fs(14);
SP_pps_epc_pps=fs(15);

SP_pps_epc = [SP_pps_epc_pms,SP_pps_epc_ps,SP_pps_epc_pps];

SF_ps_epc_pms=fs(16);
SF_pps_epc_pms=fs(17);
SF_pps_epc_ps=fs(18);

SF_epc = [SF_ps_epc_pms,SF_pps_epc_pms,SF_pps_epc_ps];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
F2 = @(var)fun_solveMar5_2(var,epsilon_d,epsilon_c,theta);
% var = [alpha]
t = 0.5;
var0=repmat(t,n,1);
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','off'); 
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt'
% Option to display output

[fs2,~] = fsolve(F2,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]
alpha=fs2;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return