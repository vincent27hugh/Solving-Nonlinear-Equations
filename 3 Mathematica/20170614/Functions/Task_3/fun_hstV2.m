% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_s = fun_hstV2(hstm1,hntm1,epsilon_ctm1,...
    epsilon_dt,thetat,epsilon_ct,utm1,lambda,typen,A,B1,B2,epsilon_u)

q_theta_t=fun_q_theta(thetat,A,B1,B2);

F_epsilonc_t = fun_F_x(epsilon_ct,typen,epsilon_u);
F_epsilond_t = fun_F_x(epsilon_dt,typen,epsilon_u);
F_epsilonc_tm1 = fun_F_x(epsilon_ctm1,typen,epsilon_u);

temp1 = hstm1;
temp2 = (1-utm1-hstm1-hntm1)*lambda*(F_epsilonc_t-F_epsilond_t);
temp3 = hstm1*(lambda*F_epsilond_t+thetat*q_theta_t);
temp4 = (F_epsilonc_t-F_epsilond_t)*(1-utm1-hstm1-hntm1)/...
    (1-F_epsilonc_tm1);

h_s =temp1 + temp2 - temp3+temp4;

return