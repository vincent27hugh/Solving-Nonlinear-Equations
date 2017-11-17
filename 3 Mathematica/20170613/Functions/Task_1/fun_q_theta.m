% Jun 13 2017 symbolic
% Task #20160308
% Related to May16
% edited in Oct11,2016
function q_theta = fun_q_theta(theta,A,B1,B2)
temp1 = theta.^(1/B2);
temp2 = (temp1).^(-B1);
q_theta = A*temp2;
return