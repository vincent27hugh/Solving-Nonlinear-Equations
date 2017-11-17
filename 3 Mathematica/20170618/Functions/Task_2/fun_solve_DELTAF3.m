% May 22, 2017
function Fn = fun_solve_DELTAF3(var,para)
% var = (DELTAF_ps, DELTAF_pps)
r = para{1};
lambda = para{2};

mu = para{4};
DELTAP_ps = para{5};

part31 = (r+lambda+mu)*var;
part32 = 1+mu*DELTAP_ps;

Fn = [part31-part32];
return