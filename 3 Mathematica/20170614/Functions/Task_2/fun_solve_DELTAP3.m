% May 22, 2017
function Fn = fun_solve_DELTAP3(var,para)
% var = (DELTAP_pps)

tau_pps = para{3};


part31 = tau_pps*var;
part32 = 1;

Fn = [part31-part32];
return