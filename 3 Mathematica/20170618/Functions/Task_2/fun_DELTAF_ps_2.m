% mAY 23, 2017
% May 22, 2017
function DELTAF_ps = fun_DELTAF_ps_2(r,lambda,mu,g_p,delta,DELTAP_pms)

temp1 = r+lambda+mu;

temp2 = 1+mu*g_p/temp1+mu*(1-g_p)*delta*DELTAP_pms;
temp3 = temp1-mu^2*g_p/temp1;

DELTAF_ps = temp2/temp3;

return