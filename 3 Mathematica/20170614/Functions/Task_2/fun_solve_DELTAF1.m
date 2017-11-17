% May 22, 2017
function Fn = fun_solve_DELTAF1(var,para)
% var = (DELTAF_pms, DELTAF_ps, DELTAF_pps)
r = para{1};
lambda = para{2};
g_p = para{3};
mu = para{4};

part11 = (r+lambda+mu)*var(1);
part12 = 1+mu*var(2);

part21 = (r+lambda+mu)*var(2);
part22 = 1+mu*(g_p*var(3)+(1-g_p)*var(1));

part31 = (r+lambda+mu)*var(3);
part32 = 1+mu*var(2);

Fn = [part11-part12;part21-part22;part31-part32];
return