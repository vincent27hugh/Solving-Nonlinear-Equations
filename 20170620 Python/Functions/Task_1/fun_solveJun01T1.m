% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function F = fun_solveJun01T1(var,parameters)
% var vetcor represent [epsilon_d;epsilon_c;theta]; 3 variables 
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
pstar = parameters{15};
typen = parameters{16};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_d=var(1);
epsilon_c=var(2);
theta=var(3);
alpha=var(4);
p=pstar;

q_theta = fun_q_theta(theta,A,B1,B2);

int_Fdu = fun_int_F(epsilon_d,epsilon_u,typen,epsilon_u);
int_Fcu = fun_int_F(epsilon_c,epsilon_u,typen,epsilon_u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part11 = epsilon_d+lambda*int_Fdu/(r+lambda+theta*q_theta);

part12 = (b-p)/sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp211 = delta*(r+lambda-phi*(1-fun_F_x(epsilon_c,typen,epsilon_u)))...
    *(epsilon_c-epsilon_d);
temp212 = r+lambda+theta*q_theta;
temp213 = (lambda/(r+lambda)-lambda*delta/temp212)*int_Fcu;
part21 = epsilon_c-temp211/((1-phi)*temp212)+temp213;

temp2211 = ((beta+phi*(1-beta))*c_f*(1-alpha))/((1-beta)*(1-phi));
temp2212 = (beta*c_p*alpha)/(1-beta);
temp221 = theta*(temp2211 + temp2212);
part22 = delta*epsilon_d+((1-delta)*(b-p)+temp221)/delta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part31 = 1/q_theta;

temp321 = (1-beta)*(1-phi)/c_f;
temp322 = sigma*(epsilon_u-epsilon_c)/(r+lambda)+...
    delta*sigma*(epsilon_c-epsilon_d)/((1-phi)*(r+lambda+theta*q_theta));
part32 = temp321*temp322;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part41 = 1/q_theta;

temp421 = (1-beta)*delta*sigma*(epsilon_u-epsilon_d);
temp422 = c_p*(r+lambda+theta*q_theta);
part42 = temp421/temp422;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = [part11-part12;part21-part22;part31 - part32;part41-part42];

return

