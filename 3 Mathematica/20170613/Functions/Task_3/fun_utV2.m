% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function u = fun_utV2(utm1,epsilon_dt,thetat,epsilon_dtm1,epsilon_ctm1,...
    hs_tm1,hn_tm1,lambda,typen,A,B1,B2,epsilon_u)
% parameter

q_theta_t=fun_q_theta(thetat,A,B1,B2);

F_epsilond_t = fun_F_x(epsilon_dt,typen,epsilon_u);
F_epsilond_tm1 = fun_F_x(epsilon_dtm1,typen,epsilon_u);
F_epsilonc_tm1 = fun_F_x(epsilon_ctm1,typen,epsilon_u);

temp = (F_epsilonc_tm1-F_epsilond_tm1)*(hs_tm1+hn_tm1)/(1-F_epsilond_tm1);
temp2 = (F_epsilond_t-F_epsilonc_tm1)*(1-utm1-hstm1-hntm1)/(1-F_epsilonc_tm1);

u = utm1+(1-utm1)*lambda*F_epsilond_t-utm1*thetat*q_theta_t+temp+temp2;

return