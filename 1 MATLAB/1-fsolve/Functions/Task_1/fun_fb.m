% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function fb = fun_fb(alpha,theta,A,B1,B2)
% parameter

q_theta=fun_q_theta(theta,A,B1,B2);

fb = (1-alpha)*theta*q_theta;
return