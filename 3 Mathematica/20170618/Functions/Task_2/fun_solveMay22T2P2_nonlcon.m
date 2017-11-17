% Task 2, May22, 2017
% Task1-2, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [c,ceq] = fun_solveMay22T2P2_nonlcon(var,para)
c=[]; % no nonlinear inequality
ceq = fun_solveMay22T2P2(var,para); 
% the fsolve objective is fmincon constrain
return

