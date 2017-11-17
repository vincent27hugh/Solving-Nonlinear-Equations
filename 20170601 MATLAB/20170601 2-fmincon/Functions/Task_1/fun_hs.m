% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_s = fun_hs(epsilon_c,epsilon_d,theta,alpha,u,lambda,typen,...
    A,B1,B2,epsilon_u)
% parameter

q_theta=fun_q_theta(theta,A,B1,B2);

F_epsilonc = fun_F_x(epsilon_c,typen,epsilon_u);
F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

temp1 = (1-u)*(1-alpha)*lambda*(F_epsilonc-F_epsilond);
temp2 = lambda*F_epsilonc+(1-alpha)*theta*q_theta;

h_s =temp1 / temp2;

return