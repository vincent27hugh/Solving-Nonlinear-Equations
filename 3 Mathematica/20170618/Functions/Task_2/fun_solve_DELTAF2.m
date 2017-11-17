% May 22, 2017
function Fn = fun_solve_DELTAF2(var,para)
% var = (DELTAF_ps, DELTAF_pps)
r = para{1};
lambda = para{2};
g_p = para{3};
mu = para{4};
DELTAP_pms = para{5};

part21 = (r+lambda+mu)*var(1);
part22 = 1+mu*(g_p*var(2)+(1-g_p)*DELTAP_pms);

part31 = (r+lambda+mu)*var(2);
part32 = 1+mu*var(1);

Fn = [part21-part22;part31-part32];
return