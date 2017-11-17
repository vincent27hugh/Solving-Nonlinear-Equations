% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_nm = fun_hnm(h_n,theta,alpha,epsilon_d,lambda,typen,...
    A,B1,B2,epsilon_u)
% parameter

q_theta=fun_q_theta(theta,A,B1,B2);

F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

h_nm = h_n*(lambda*F_epsilond+(1-alpha)*theta*q_theta);

return