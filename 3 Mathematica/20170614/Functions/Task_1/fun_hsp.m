% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_sp = fun_hsp(epsilon_c,epsilon_d,h_s,h_n,u,lambda,typen)
% parameter

F_epsilonc = fun_F_x(epsilon_c,typen,epsilon_u);
F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

h_sp = (1-u-h_s-h_n)*lambda*(F_epsilonc-F_epsilond);

return