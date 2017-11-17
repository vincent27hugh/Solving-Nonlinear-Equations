% Task1-2, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilond,theta,ceqP1,SbarP,exitflag] = Jun12T2P1(parameters)
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

n = int16(2*width/step+1);
% length of the lattice
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;step;...
                width;pstar;typen};
F = @(var)fun_solveMay22T2P1(var,para);
% var = [epsilon_d;epsilon_c;theta]
% theta should be larger than 0
t1 = linspace(-1,-10,n);
t2 = linspace(1,5,n);
var0=reshape([t1;t2],[],1);
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','iter'); 
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt'
% Option to display output
% LM Algo is better

[fs,~,exitflag] = fsolve(F,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]
ceqP1=fun_solveMay22T2P1(fs,para);
%{
t1 = linspace(-1600,-1800,n);
t2 = linspace(50,50,n);
var0=reshape([t1;t2],[],1);
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
t3=-Inf(1,n);
lbfmin = reshape([t1;t3],[],1);
%%%%
t1=repmat(epsilon_u,1,n);
t3=Inf(1,n);
ubfmin = reshape([t1;t3],[],1);
%%%%
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;step;...
                width;pstar;typen};
[fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,...
    Aeqfmin,beqfmin,lbfmin,ubfmin,...
    @(x)fun_solveMay22T2P1_nonlcon(x,para),options);
% Find minimum of constrained nonlinear multivariable function
[~,ceqP1]=fun_solveMay22T2P1_nonlcon(fs,para);
%}
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