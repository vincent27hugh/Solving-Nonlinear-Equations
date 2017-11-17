% Symbolic Task1, Jun13, 2017
% June 12, 2017
% June 1,2017
% May 22-24, 2017
% 20170411 PM
% Mar18,2017
% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path0 = pwd;
add_path=strcat(path0,'\Functions\Task_1\');
addpath(add_path);

path_fig = strcat(path0,'\Figures\Task_1\');
path_data = strcat(path0,'\Data\Task_1\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = sym('A');
B1= sym('B1');
B2 = sym('B2');
cf = sym('cf');
cp = sym('cp');
beta = sym('beta');
phi = sym('phi');
delta = sym('delta');
sigma = sym('sigma');
lambda = sym('lambda');
b = sym('b');
r = sym('r');
epsilonu = sym('epsilonu');
mu = sym('mu');
pstar = sym('pstar');

typen = 'III';

%%%%%%%%%%%%%%%%%%%%%%%%
epsilond=sym('epsilond');
epsilonc=sym('epsilonc');
theta=sym('theta');
alpha=sym('alpha');

p=pstar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qtheta = fun_q_theta(theta,A,B1,B2);

intFdu = fun_int_F(epsilond,epsilonu,typen,epsilonu);
intFcu = fun_int_F(epsilonc,epsilonu,typen,epsilonu);
%{
%test
test = fun_int_F(epsilond,epsilonu,typen,epsilonu);
test=subs(test,epsilond,-10);
test=subs(test,epsilonu,1);
vpa(test,5)
% ans = 10.001
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part11 = epsilond+lambda*intFdu/(r+lambda+theta*qtheta);

part12 = (b-p)/sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp211 = delta*(r+lambda-phi*(1-fun_F_x(epsilonc,typen,epsilonu)))...
    *(epsilonc-epsilond);
temp212 = r+lambda+theta*qtheta;
temp213 = (lambda/(r+lambda)-lambda*delta/temp212)*intFcu;
part21 = epsilonc-temp211/((1-phi)*temp212)+temp213;

temp2211 = ((beta+phi*(1-beta))*cf*(1-alpha))/((1-beta)*(1-phi));
temp2212 = (beta*cp*alpha)/(1-beta);
temp221 = theta*(temp2211 + temp2212);
part22 = delta*epsilond+((1-delta)*(b-p)+temp221)/delta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part31 = 1/qtheta;

temp321 = (1-beta)*(1-phi)/cf;
temp322 = sigma*(epsilonu-epsilonc)/(r+lambda)+...
    delta*sigma*(epsilonc-epsilond)/((1-phi)*(r+lambda+theta*qtheta));
part32 = temp321*temp322;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part41 = 1/qtheta;

temp421 = (1-beta)*delta*sigma*(epsilonu-epsilond);
temp422 = cp*(r+lambda+theta*qtheta);
part42 = temp421/temp422;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = [part11-part12,part21-part22,part31 - part32,part41-part42];
% function vector
%{
%   test
Ftest=F;
Ftest=subs(Ftest,...
    [A,B1,B2],...
    [1.355,18,25]);
Ftest=subs(Ftest,mu,0.08);
Ftest=subs(Ftest,epsilonu,1);
Ftest=subs(Ftest,[cf,cp,beta,phi,delta,sigma,lambda,b,r,pstar],...
    [0.213,0.113,0.72,0.1,0.5,0.0436,0.1,0.2,0.012,1]);
% 10 parameters
Ftest=subs(Ftest,epsilond,-5);
Ftest=subs(Ftest,epsilonc,-1);
Ftest=subs(Ftest,theta,1);
Ftest=subs(Ftest,alpha,0.5);

vpa(Ftest,5)
%   ans =
 
%    13.69
%   2.3107
% -0.26126
%  0.51708
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var=[epsilond,epsilonc,theta,alpha,...
    pstar,b,phi,sigma,beta,...
    lambda,cf,r,cp,delta];
% [1,14]

Jacfvar = jacobian(F,var);
% jacobian matrix
% [4,14]
% [dF1_dx1,...dF1_dxN,
%  
%  dFN_dx1,.....dFN_dxN]

%{
Ftest=Jacfvar;

Ftest=subs(Ftest,...
    [A,B1,B2],...
    [1.355,18,25]);
Ftest=subs(Ftest,mu,0.08);
Ftest=subs(Ftest,epsilonu,1);
Ftest=subs(Ftest,[cf,cp,beta,phi,delta,sigma,lambda,b,r,pstar],...
    [0.213,0.113,0.72,0.1,0.5,0.0436,0.1,0.2,0.012,1]);
% 10 parameters
Ftest=subs(Ftest,epsilond,-5);
Ftest=subs(Ftest,epsilonc,-1);
Ftest=subs(Ftest,theta,1);
Ftest=subs(Ftest,alpha,0.5);
vpa(Ftest,5)
%{
ans =
 
[  0.93218,       0, -0.088335,       0, 22.936, -22.936,       0, -420.84,       0,   3.1827,       0, -0.23283,       0,        0]
[ -0.49145, 0.20694,  -0.90403, 0.68333,    1.0,    -1.0, 0.37728,       0,   -4.46, -0.77137, -2.9683,  -10.449, -2.5714,   3.5004]
[ 0.019535, 0.44103,   0.55157,       0,      0,       0,  1.0235, -22.919,  3.5688,   8.2776,  4.6914,   8.2776,       0, -0.15628]
[ 0.036822,       0,    0.5885,       0,      0,       0,       0, -5.0672, 0.78904,   0.1506,       0,   0.1506,  1.9551, -0.44186]
 
%}
%}

%{
matlabFunction(Jacfvar,...
    'vars',{A,B1,B2,mu,epsilonu,...
    cf,cp,beta,phi,delta,sigma,lambda,b,r,pstar,...
    epsilond,epsilonc,theta,alpha},...
    'file',strcat('Jacobian_Jun13_',typen));
% 19 input arguments
%}
%{
% test
Jacfvar = Jacobian_Jun13_III(1.355,18,25,0.08,1,...
    0.213,0.113,0.72,0.1,0.5,0.0436,0.1,0.2,0.012,1,...
    -5,-1,1,0.5)
vpa(Jacfvar,5)
%{
ans =
 
[  0.93218,       0, -0.088335,       0, 22.936, -22.936,       0, -420.84,       0,   3.1827,       0, -0.23283,       0,        0]
[ -0.49145, 0.20694,  -0.90403, 0.68333,    1.0,    -1.0, 0.37728,       0,   -4.46, -0.77137, -2.9683,  -10.449, -2.5714,   3.5004]
[ 0.019535, 0.44103,   0.55157,       0,      0,       0,  1.0235, -22.919,  3.5688,   8.2776,  4.6914,   8.2776,       0, -0.15628]
[ 0.036822,       0,    0.5885,       0,      0,       0,       0, -5.0672, 0.78904,   0.1506,       0,   0.1506,  1.9551, -0.44186]
 
%}
%}