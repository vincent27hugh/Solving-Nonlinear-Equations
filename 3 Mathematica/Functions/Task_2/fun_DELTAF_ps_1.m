% mAY 23, 2017
% May 22, 2017
function DELTAF_ps = fun_DELTAF_ps_1(r,lambda,mu)

temp1 = r+lambda+mu;

temp2 = 1+mu/temp1;
temp3 = temp1-mu^2/temp1;

DELTAF_ps = temp2/temp3;

return