% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function fp = fun_fp(u,h_s,h_n,alpha,theta,A,B1,B2)
% parameter

q_theta=fun_q_theta(theta,A,B1,B2);

fp = (u+h_s+h_n)*(1-alpha)*theta*q_theta;

return