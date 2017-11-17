% May 22, 2017
function tau = fun_tau(theta,q_theta,r,lambda,mu)

tau = r+lambda+theta*q_theta+mu;

return