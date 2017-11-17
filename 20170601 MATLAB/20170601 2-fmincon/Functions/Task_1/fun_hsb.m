% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_sb = fun_hsb(epsilon_c,epsilon_d,lambda,typen,...
    epsilon_u)
% parameter 
F_epsilonc = fun_F_x(epsilon_c,typen,epsilon_u);
F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

h_sb = lambda*(F_epsilonc-F_epsilond);

return