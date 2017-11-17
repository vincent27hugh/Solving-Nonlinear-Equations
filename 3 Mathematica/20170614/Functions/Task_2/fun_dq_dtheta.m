% Task #20160308
% Related to May16
% edited in Oct11,2016
function dq_dtheta = fun_dq_dtheta(theta,...
    A,B1,B2)
if isa(theta,'sym')
    syms q_theta
    
    temp2 = (theta)^(B1/B2+1);
    temp3 = 1/temp2;
    dq_dtheta = (-A*B1/B2)*temp3;
else
    temp1 = nthroot(theta,B2);
    temp2 = (temp1)^(B1+B2);
    temp3 = 1/temp2;
    dq_dtheta = (-A*B1/B2)*temp3;
end
return