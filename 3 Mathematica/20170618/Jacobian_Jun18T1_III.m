function Jac=Jacobian_Jun18T1_III(A,B1,B2,mu,epsilon_u,...
    c_f,c_p,beta,phi,delta,sigma,lambda,b,r,pstar,...
    epsilon_d,epsilon_c,theta,alpha)
%    This function was generated by ToMatlab Package from Mathematica
%    18-Jun-2017

B=B1./B2;

Jac=[1+lambda.*(r+A.*theta.^(1+(-1).*B)+lambda) ...
  .^(-1).*((-0.5E0)+( ...
  -0.479179E0).*exp(1).^((-1).*( ...
  0.294353E0+0.849322E0.*log( ...
  0.1E1+(-0.1E1).*epsilon_d)).^2)+ ...
  0.479179E0.*exp(1).^((-1).*( ...
  0.294353E0+(-0.849322E0).*log( ...
  0.1E1+(-0.1E1).*epsilon_d)).^2).*( ...
  0.1E1+(-0.1E1).*epsilon_d).^(-1)+( ...
  -0.5E0).*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_d))),0,(-1).*A.*(1+(-1).*B) ...
  .*theta.^((-1).*B).*lambda.*(r+A.*theta.^( ...
  1+(-1).*B)+lambda).^(-2).*(0.5E0.*( ...
  0.1E1+(-0.1E1).*epsilon_d).^0.1E1+( ...
  -0.5E0).*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1+0.5E0.*erf(0.294353E0+ ...
  (-0.849322E0).*log(0.1E1+( ...
  -0.1E1).*epsilon_d))+0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_d).^0.1E1.*erf( ...
  0.294353E0+0.849322E0.*log( ...
  0.1E1+(-0.1E1).*epsilon_d))+(-0.5E0) ...
  .*erf(0.294353E0+(-0.849322E0) ...
  .*log(0.1E1+(-0.1E1).*epsilon_u))+( ...
  -0.5E0).*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u))),0,sigma.^(-1),(-1).*sigma.^( ...
  -1),0,(b+(-1).*pstar).*sigma.^(-2) ...
  ,0,(-1).*lambda.*(r+A.*theta.^(1+(-1).* ...
  B)+lambda).^(-2).*(0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_d).^0.1E1+(-0.5E0).* ...
  (0.1E1+(-0.1E1).*epsilon_u).^0.1E1+ ...
  0.5E0.*erf(0.294353E0+( ...
  -0.849322E0).*log(0.1E1+( ...
  -0.1E1).*epsilon_d))+0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_d).^0.1E1.*erf( ...
  0.294353E0+0.849322E0.*log( ...
  0.1E1+(-0.1E1).*epsilon_d))+(-0.5E0) ...
  .*erf(0.294353E0+(-0.849322E0) ...
  .*log(0.1E1+(-0.1E1).*epsilon_u))+( ...
  -0.5E0).*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u)))+(r+A.*theta.^(1+(-1).*B)+ ...
  lambda).^(-1).*(0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_d).^0.1E1+(-0.5E0).* ...
  (0.1E1+(-0.1E1).*epsilon_u).^0.1E1+ ...
  0.5E0.*erf(0.294353E0+( ...
  -0.849322E0).*log(0.1E1+( ...
  -0.1E1).*epsilon_d))+0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_d).^0.1E1.*erf( ...
  0.294353E0+0.849322E0.*log( ...
  0.1E1+(-0.1E1).*epsilon_d))+(-0.5E0) ...
  .*erf(0.294353E0+(-0.849322E0) ...
  .*log(0.1E1+(-0.1E1).*epsilon_u))+( ...
  -0.5E0).*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u))),0,(-1).*lambda.*(r+A.*theta.^( ...
  1+(-1).*B)+lambda).^(-2).*(0.5E0.*( ...
  0.1E1+(-0.1E1).*epsilon_d).^0.1E1+( ...
  -0.5E0).*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1+0.5E0.*erf(0.294353E0+ ...
  (-0.849322E0).*log(0.1E1+( ...
  -0.1E1).*epsilon_d))+0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_d).^0.1E1.*erf( ...
  0.294353E0+0.849322E0.*log( ...
  0.1E1+(-0.1E1).*epsilon_d))+(-0.5E0) ...
  .*erf(0.294353E0+(-0.849322E0) ...
  .*log(0.1E1+(-0.1E1).*epsilon_u))+( ...
  -0.5E0).*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u))),0,0;(-1).*delta+delta.*(r+A.* ...
  theta.^(1+(-1).*B)+lambda).^(-1).*(1+( ...
  -1).*phi).^(-1).*(r+lambda+(-1).*phi.*( ...
  0.5E0+0.5E0.*erf((2.*log(2)) ...
  .^(-1/2).*(0.346574E0+log(1+( ...
  -1).*epsilon_c))))),1+(-0.479179E0).* ...
  exp(1).^((-1/2).*log(2).^(-1) ...
  .*(0.346574E0+log(1+(-1).*epsilon_c)) ...
  .^2).*delta.*(1+(-1).*epsilon_c).^(-1).*( ...
  epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^(1+(-1) ...
  .*B)+lambda).^(-1).*(1+(-1).*phi).^( ...
  -1).*phi+(-1).*delta.*(r+A.*theta.^(1+( ...
  -1).*B)+lambda).^(-1).*(1+(-1).*phi) ...
  .^(-1).*(r+lambda+(-1).*phi.*(0.5E0+ ...
  0.5E0.*erf((2.*log(2)).^(-1/2) ...
  .*(0.346574E0+log(1+(-1).*epsilon_c)) ...
  )))+(lambda.*(r+lambda).^(-1)+(-1).*delta.* ...
  lambda.*(r+A.*theta.^(1+(-1).*B)+lambda).^( ...
  -1)).*((-0.5E0)+(-0.479179E0) ...
  .*exp(1).^((-1).*(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_c)).^2)+0.479179E0.*exp(1) ...
  .^((-1).*(0.294353E0+( ...
  -0.849322E0).*log(0.1E1+( ...
  -0.1E1).*epsilon_c)).^2).*(0.1E1+( ...
  -0.1E1).*epsilon_c).^(-1)+(-0.5E0).* ...
  erf(0.294353E0+0.849322E0.* ...
  log(0.1E1+(-0.1E1).*epsilon_c))),(-1) ...
  .*delta.^(-1).*(c_p.*alpha.*(1+(-1).*beta) ...
  .^(-1).*beta+c_f.*(1+(-1).*alpha).*(1+ ...
  (-1).*beta).^(-1).*(1+(-1).*phi).^( ...
  -1).*(beta+(1+(-1).*beta).*phi))+A.*( ...
  1+(-1).*B).*delta.*(epsilon_c+(-1).*epsilon_d).* ...
  theta.^((-1).*B).*(r+A.*theta.^(1+(-1) ...
  .*B)+lambda).^(-2).*(1+(-1).*phi).^( ...
  -1).*(r+lambda+(-1).*phi.*(0.5E0+ ...
  0.5E0.*erf((2.*log(2)).^(-1/2) ...
  .*(0.346574E0+log(1+(-1).*epsilon_c)) ...
  )))+A.*(1+(-1).*B).*delta.*theta.^(( ...
  -1).*B).*lambda.*(r+A.*theta.^(1+(-1).* ...
  B)+lambda).^(-2).*(0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_c).^0.1E1+(-0.5E0).* ...
  (0.1E1+(-0.1E1).*epsilon_u).^0.1E1+ ...
  0.5E0.*erf(0.294353E0+( ...
  -0.849322E0).*log(0.1E1+( ...
  -0.1E1).*epsilon_c))+0.5E0.*(0.1E1+( ...
  -0.1E1).*epsilon_c).^0.1E1.*erf( ...
  0.294353E0+0.849322E0.*log( ...
  0.1E1+(-0.1E1).*epsilon_c))+(-0.5E0) ...
  .*erf(0.294353E0+(-0.849322E0) ...
  .*log(0.1E1+(-0.1E1).*epsilon_u))+( ...
  -0.5E0).*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u))),(-1).*delta.^(-1).*theta.*( ...
  c_p.*(1+(-1).*beta).^(-1).*beta+(-1) ...
  .*c_f.*(1+(-1).*beta).^(-1).*(1+( ...
  -1).*phi).^(-1).*(beta+(1+(-1).*beta) ...
  .*phi)),(-1).*((-1)+delta).*delta.^(-1), ...
  (-1).*(1+(-1).*delta).*delta.^(-1),( ...
  -1).*delta.^(-1).*theta.*(c_f.*(1+(-1) ...
  .*alpha).*(1+(-1).*phi).^(-1)+c_f.*( ...
  1+(-1).*alpha).*(1+(-1).*beta).^(-1) ...
  .*(1+(-1).*phi).^(-2).*(beta+(1+( ...
  -1).*beta).*phi))+(-1).*delta.*(epsilon_c+(-1) ...
  .*epsilon_d).*(r+A.*theta.^(1+(-1).*B)+lambda) ...
  .^(-1).*(1+(-1).*phi).^(-2).*(r+ ...
  lambda+(-1).*phi.*(0.5E0+0.5E0.*erf(( ...
  2.*log(2)).^(-1/2).*( ...
  0.346574E0+log(1+(-1).*epsilon_c))))) ...
  +delta.*(epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^( ...
  1+(-1).*B)+lambda).^(-1).*(1+(-1).* ...
  phi).^(-1).*(0.5E0+0.5E0.*erf(( ...
  2.*log(2)).^(-1/2).*( ...
  0.346574E0+log(1+(-1).*epsilon_c)))), ...
  0,(-1).*delta.^(-1).*theta.*(c_f.*(1+( ...
  -1).*alpha).*(1+(-1).*beta).^(-1)+ ...
  c_p.*alpha.*(1+(-1).*beta).^(-1)+c_p.* ...
  alpha.*(1+(-1).*beta).^(-2).*beta+c_f.*( ...
  1+(-1).*alpha).*(1+(-1).*beta).^(-2) ...
  .*(1+(-1).*phi).^(-1).*(beta+(1+( ...
  -1).*beta).*phi)),(-1).*delta.*(epsilon_c+(-1) ...
  .*epsilon_d).*(r+A.*theta.^(1+(-1).*B)+lambda) ...
  .^(-1).*(1+(-1).*phi).^(-1)+delta.*( ...
  epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^(1+(-1) ...
  .*B)+lambda).^(-2).*(1+(-1).*phi).^( ...
  -1).*(r+lambda+(-1).*phi.*(0.5E0+ ...
  0.5E0.*erf((2.*log(2)).^(-1/2) ...
  .*(0.346574E0+log(1+(-1).*epsilon_c)) ...
  )))+((-1).*lambda.*(r+lambda).^(-2)+(r+ ...
  lambda).^(-1)+delta.*lambda.*(r+A.*theta.^(1+( ...
  -1).*B)+lambda).^(-2)+(-1).*delta.*(r+ ...
  A.*theta.^(1+(-1).*B)+lambda).^(-1)).*( ...
  0.5E0.*(0.1E1+(-0.1E1).*epsilon_c) ...
  .^0.1E1+(-0.5E0).*(0.1E1+( ...
  -0.1E1).*epsilon_u).^0.1E1+0.5E0.* ...
  erf(0.294353E0+(-0.849322E0).* ...
  log(0.1E1+(-0.1E1).*epsilon_c))+ ...
  0.5E0.*(0.1E1+(-0.1E1).*epsilon_c) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_c))+(-0.5E0).*erf( ...
  0.294353E0+(-0.849322E0).*log( ...
  0.1E1+(-0.1E1).*epsilon_u))+(-0.5E0) ...
  .*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u))),(-1).*(1+(-1).*alpha).*(1+ ...
  (-1).*beta).^(-1).*delta.^(-1).*theta.*( ...
  1+(-1).*phi).^(-1).*(beta+(1+(-1).* ...
  beta).*phi),(-1).*delta.*(epsilon_c+(-1).*epsilon_d) ...
  .*(r+A.*theta.^(1+(-1).*B)+lambda).^( ...
  -1).*(1+(-1).*phi).^(-1)+delta.*(epsilon_c+ ...
  (-1).*epsilon_d).*(r+A.*theta.^(1+(-1).* ...
  B)+lambda).^(-2).*(1+(-1).*phi).^(-1) ...
  .*(r+lambda+(-1).*phi.*(0.5E0+0.5E0.* ...
  erf((2.*log(2)).^(-1/2).*( ...
  0.346574E0+log(1+(-1).*epsilon_c))))) ...
  +((-1).*lambda.*(r+lambda).^(-2)+delta.*lambda.*( ...
  r+A.*theta.^(1+(-1).*B)+lambda).^(-2)) ...
  .*(0.5E0.*(0.1E1+(-0.1E1).*epsilon_c) ...
  .^0.1E1+(-0.5E0).*(0.1E1+( ...
  -0.1E1).*epsilon_u).^0.1E1+0.5E0.* ...
  erf(0.294353E0+(-0.849322E0).* ...
  log(0.1E1+(-0.1E1).*epsilon_c))+ ...
  0.5E0.*(0.1E1+(-0.1E1).*epsilon_c) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_c))+(-0.5E0).*erf( ...
  0.294353E0+(-0.849322E0).*log( ...
  0.1E1+(-0.1E1).*epsilon_u))+(-0.5E0) ...
  .*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u))),(-1).*alpha.*(1+(-1).*beta) ...
  .^(-1).*beta.*delta.^(-1).*theta,(-1).*(( ...
  -1).*b+pstar).*delta.^(-1)+(-1).* ...
  epsilon_d+delta.^(-2).*((b+(-1).*pstar).* ...
  (1+(-1).*delta)+theta.*(c_p.*alpha.*(1+(-1) ...
  .*beta).^(-1).*beta+c_f.*(1+(-1).*alpha) ...
  .*(1+(-1).*beta).^(-1).*(1+(-1).* ...
  phi).^(-1).*(beta+(1+(-1).*beta).*phi))) ...
  +(-1).*(epsilon_c+(-1).*epsilon_d).*(r+A.* ...
  theta.^(1+(-1).*B)+lambda).^(-1).*(1+( ...
  -1).*phi).^(-1).*(r+lambda+(-1).*phi.*( ...
  0.5E0+0.5E0.*erf((2.*log(2)) ...
  .^(-1/2).*(0.346574E0+log(1+( ...
  -1).*epsilon_c)))))+(-1).*lambda.*(r+A.* ...
  theta.^(1+(-1).*B)+lambda).^(-1).*( ...
  0.5E0.*(0.1E1+(-0.1E1).*epsilon_c) ...
  .^0.1E1+(-0.5E0).*(0.1E1+( ...
  -0.1E1).*epsilon_u).^0.1E1+0.5E0.* ...
  erf(0.294353E0+(-0.849322E0).* ...
  log(0.1E1+(-0.1E1).*epsilon_c))+ ...
  0.5E0.*(0.1E1+(-0.1E1).*epsilon_c) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_c))+(-0.5E0).*erf( ...
  0.294353E0+(-0.849322E0).*log( ...
  0.1E1+(-0.1E1).*epsilon_u))+(-0.5E0) ...
  .*(0.1E1+(-0.1E1).*epsilon_u) ...
  .^0.1E1.*erf(0.294353E0+ ...
  0.849322E0.*log(0.1E1+(-0.1E1) ...
  .*epsilon_u)));c_f.^(-1).*(1+(-1).*beta) ...
  .*delta.*(r+A.*theta.^(1+(-1).*B)+lambda) ...
  .^(-1).*sigma,(-1).*c_f.^(-1).*(1+( ...
  -1).*beta).*((-1).*(r+lambda).^(-1).* ...
  sigma+delta.*(r+A.*theta.^(1+(-1).*B)+lambda) ...
  .^(-1).*sigma.*(1+(-1).*phi).^(-1)) ...
  .*(1+(-1).*phi),A.^(-1).*B.*theta.^( ...
  (-1)+B)+A.*(1+(-1).*B).*c_f.^( ...
  -1).*(1+(-1).*beta).*delta.*(epsilon_c+(-1) ...
  .*epsilon_d).*theta.^((-1).*B).*(r+A.* ...
  theta.^(1+(-1).*B)+lambda).^(-2).*sigma,0, ...
  0,0,c_f.^(-1).*(1+(-1).*beta).*((( ...
  -1).*epsilon_c+epsilon_u).*(r+lambda).^(-1).*sigma+ ...
  delta.*(epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^(1+ ...
  (-1).*B)+lambda).^(-1).*sigma.*(1+(-1) ...
  .*phi).^(-1))+(-1).*c_f.^(-1).*( ...
  1+(-1).*beta).*delta.*(epsilon_c+(-1).*epsilon_d).* ...
  (r+A.*theta.^(1+(-1).*B)+lambda).^(-1) ...
  .*sigma.*(1+(-1).*phi).^(-1),(-1).* ...
  c_f.^(-1).*(1+(-1).*beta).*(((-1) ...
  .*epsilon_c+epsilon_u).*(r+lambda).^(-1)+delta.*(epsilon_c+( ...
  -1).*epsilon_d).*(r+A.*theta.^(1+(-1).*B) ...
  +lambda).^(-1).*(1+(-1).*phi).^(-1)) ...
  .*(1+(-1).*phi),c_f.^(-1).*(((-1) ...
  .*epsilon_c+epsilon_u).*(r+lambda).^(-1).*sigma+delta.*( ...
  epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^(1+(-1) ...
  .*B)+lambda).^(-1).*sigma.*(1+(-1).*phi) ...
  .^(-1)).*(1+(-1).*phi),(-1).* ...
  c_f.^(-1).*(1+(-1).*beta).*((-1).* ...
  ((-1).*epsilon_c+epsilon_u).*(r+lambda).^(-2).*sigma+ ...
  (-1).*delta.*(epsilon_c+(-1).*epsilon_d).*(r+A.* ...
  theta.^(1+(-1).*B)+lambda).^(-2).*sigma.*( ...
  1+(-1).*phi).^(-1)).*(1+(-1).*phi) ...
  ,c_f.^(-2).*(1+(-1).*beta).*(((-1) ...
  .*epsilon_c+epsilon_u).*(r+lambda).^(-1).*sigma+delta.*( ...
  epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^(1+(-1) ...
  .*B)+lambda).^(-1).*sigma.*(1+(-1).*phi) ...
  .^(-1)).*(1+(-1).*phi),(-1).* ...
  c_f.^(-1).*(1+(-1).*beta).*((-1).* ...
  ((-1).*epsilon_c+epsilon_u).*(r+lambda).^(-2).*sigma+ ...
  (-1).*delta.*(epsilon_c+(-1).*epsilon_d).*(r+A.* ...
  theta.^(1+(-1).*B)+lambda).^(-2).*sigma.*( ...
  1+(-1).*phi).^(-1)).*(1+(-1).*phi) ...
  ,0,(-1).*c_f.^(-1).*(1+(-1).*beta) ...
  .*(epsilon_c+(-1).*epsilon_d).*(r+A.*theta.^(1+( ...
  -1).*B)+lambda).^(-1).*sigma;c_p.^(-1).* ...
  (1+(-1).*beta).*delta.*(r+A.*theta.^(1+( ...
  -1).*B)+lambda).^(-1).*sigma,0,A.^(-1) ...
  .*B.*theta.^((-1)+B)+A.*(1+(-1).* ...
  B).*c_p.^(-1).*(1+(-1).*beta).*delta.* ...
  ((-1).*epsilon_d+epsilon_u).*theta.^((-1).*B).*( ...
  r+A.*theta.^(1+(-1).*B)+lambda).^(-2).* ...
  sigma,0,0,0,0,(-1).*c_p.^(-1).*(1+( ...
  -1).*beta).*delta.*((-1).*epsilon_d+epsilon_u).*(r+ ...
  A.*theta.^(1+(-1).*B)+lambda).^(-1), ...
  c_p.^(-1).*delta.*((-1).*epsilon_d+epsilon_u).*( ...
  r+A.*theta.^(1+(-1).*B)+lambda).^(-1).* ...
  sigma,c_p.^(-1).*(1+(-1).*beta).*delta.*(( ...
  -1).*epsilon_d+epsilon_u).*(r+A.*theta.^(1+(-1) ...
  .*B)+lambda).^(-2).*sigma,0,c_p.^(-1).*( ...
  1+(-1).*beta).*delta.*((-1).*epsilon_d+epsilon_u).* ...
  (r+A.*theta.^(1+(-1).*B)+lambda).^(-2) ...
  .*sigma,c_p.^(-2).*(1+(-1).*beta).*delta.* ...
  ((-1).*epsilon_d+epsilon_u).*(r+A.*theta.^(1+( ...
  -1).*B)+lambda).^(-1).*sigma,(-1).* ...
  c_p.^(-1).*(1+(-1).*beta).*((-1).* ...
  epsilon_d+epsilon_u).*(r+A.*theta.^(1+(-1).*B)+ ...
  lambda).^(-1).*sigma];

return
