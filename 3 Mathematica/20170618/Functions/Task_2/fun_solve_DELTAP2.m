% May 22, 2017
function Fn = fun_solve_DELTAP2(var,para)
% var = (DELTAP_ps, DELTAP_pps)

tau_ps = para{2};
tau_pps = para{3};

g_p = para{4};
mu = para{5};

part21 = tau_ps*var(1);
part22 = 1+mu*g_p*var(2);

part31 = tau_pps*var(2);
part32 = 1+mu*var(1);

Fn = [part21-part22;part31-part32];
return