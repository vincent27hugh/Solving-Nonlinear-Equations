%Mar18,2017
% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [c,ceq] = fun_solveMar4T1_nonlcon(var,para)
c=[]; % no nonlinear inequality
ceq = fun_solveMar4T1(var,para); 
% the fsolve objective is fmincon constrain
return

