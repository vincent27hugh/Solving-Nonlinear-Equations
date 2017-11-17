% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function up = fun_up(epsilon_d,u,lambda,typen,epsilon_u)

F_epsilond = fun_F_x(epsilon_d,typen,epsilon_u);

up = (1-u)*lambda*F_epsilond;

return