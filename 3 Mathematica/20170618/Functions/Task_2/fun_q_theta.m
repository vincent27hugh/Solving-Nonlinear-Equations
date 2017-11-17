% Task #20160308
% Related to May16
% edited in Oct11,2016
function q_theta = fun_q_theta(theta,...
    A,B1,B2)

temp1 = nthroot(theta,B2);
temp2 = (temp1).^(B1);
temp3 = 1./temp2;
q_theta = A.*temp3;

return