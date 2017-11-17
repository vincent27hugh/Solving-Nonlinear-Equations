% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_np = fun_hnp(h_s,theta,alpha,u,A,B1,B2)

q_theta=fun_q_theta(theta,A,B1,B2);

h_np = (u+h_s)*alpha*theta*q_theta;

return