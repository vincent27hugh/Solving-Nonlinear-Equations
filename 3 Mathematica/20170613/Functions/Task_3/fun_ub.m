% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function ub = fun_ub(epsilon_d,lambda,typen,epsilon_u)

F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

ub = lambda*(F_epsilond);

return