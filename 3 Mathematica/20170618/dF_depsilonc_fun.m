function dF_depsilonc=dF_depsilonc_fun(A,B1,B2,mu,epsilon_u,...
    c_f,c_p,beta,phi,delta,sigma,lambda,b,r,pstar,...
    epsilon_d,epsilon_c,theta,alpha)

%    This function was generated by ToMatlab Package from Mathematica
%    18-Jun-2017

B=B1./B2;

dF_depsilonc=[0,1+(-0.479179E0).*exp(1).^((-1/2).*log(2).^(-1) ...
  .*(0.346574E0+log(1+(-1).*epsilon_c)).^2).*delta.*(1+(-1).* ...
  epsilon_c).^(-1).*(epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^(1+(-1).*B)+lambda) ...
  .^(-1).*(1+(-1).*phi).^(-1).*phi+(-1).*delta.*(r+A.*theta.^(1+ ...
  (-1).*B)+lambda).^(-1).*(1+(-1).*phi).^(-1).*(r+lambda+(-1).* ...
  phi.*(0.5E0+0.5E0.*erf((2.*log(2)).^(-1/2).*( ...
  0.346574E0+log(1+(-1).*epsilon_c)))))+(lambda.*(r+lambda).^(-1)+( ...
  -1).*delta.*lambda.*(r+A.*theta.^(1+(-1).*B)+lambda).^(-1)).*(( ...
  -0.5E0)+(-0.479179E0).*exp(1).^((-1).*(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1).*epsilon_c)).^2)+ ...
  0.479179E0.*exp(1).^((-1).*(0.294353E0+( ...
  -0.849322E0).*log(0.1E1+(-0.1E1).*epsilon_c)).^2).*( ...
  0.1E1+(-0.1E1).*epsilon_c).^(-1)+(-0.5E0).*erf( ...
  0.294353E0+0.849322E0.*log(0.1E1+(-0.1E1).*epsilon_c))),( ...
  -1).*c_f.^(-1).*(1+(-1).*beta).*((-1).*(r+lambda).^(-1).*sigma+ ...
  delta.*(r+A.*theta.^(1+(-1).*B)+lambda).^(-1).*sigma.*(1+(-1).*phi) ...
  .^(-1)).*(1+(-1).*phi),0];

return