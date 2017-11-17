% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function um = fun_um(theta,u,A,B1,B2)

q_theta = fun_q_theta(theta,A,B1,B2);

um = (u)*theta*q_theta;

return