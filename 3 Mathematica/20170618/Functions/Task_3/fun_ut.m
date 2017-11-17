% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function u = fun_ut(utm1,epsilon_dt,thetat,...
    lambda,typen,A,B1,B2,epsilon_u)
% parameter

q_theta_t=fun_q_theta(thetat,A,B1,B2);

F_epsilond_t = fun_F_x(epsilon_dt,typen,epsilon_u);

u = utm1+(1-utm1)*lambda*F_epsilond_t-utm1*thetat*q_theta_t;

return