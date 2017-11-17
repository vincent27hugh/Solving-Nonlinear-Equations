% May 22, 2017
function Fn = fun_solve_DELTAP1(var,para)
% var = (DELTAP_pms, DELTAP_ps, DELTAP_pps)
tau_pms = para{1};
tau_ps = para{2};
tau_pps = para{3};

g_p = para{4};
mu = para{5};

part11 = tau_pms*var(1);
part12 = 1+mu*var(2);

part21 = tau_ps*var(2);
part22 = 1+mu*(g_p*var(3)+(1-g_p)*var(1));

part31 = tau_pps*var(3);
part32 = 1+mu*var(2);

Fn = [part11-part12;part21-part22;part31-part32];
return