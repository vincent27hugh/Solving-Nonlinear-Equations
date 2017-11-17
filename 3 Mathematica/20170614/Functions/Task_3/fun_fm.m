% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function fm = fun_fm(f,epsilon_c,lambda,typen,epsilon_u)

F_epsilonc = fun_F_x(epsilon_c,typen,epsilon_u);

fm = f*lambda*F_epsilonc;

return