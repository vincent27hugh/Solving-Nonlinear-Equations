% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_sm = fun_hsm(epsilon_d,h_s,theta,lambda,typen,....
    A,B1,B2,epsilon_u)
% parameter

q_theta = fun_q_theta(theta,A,B1,B2);

F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

h_sm = h_s*(lambda*F_epsilond+theta*q_theta);

return