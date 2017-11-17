% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function u = fun_u(epsilon_d,theta,lambda,typen,...
    A,B1,B2,epsilon_u)
% parameter

q_theta=fun_q_theta(theta,A,B1,B2);

F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

u = lambda*F_epsilond/(lambda*F_epsilond+theta*q_theta);

return