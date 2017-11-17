% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function h_n = fun_hnt(hntm1,hstm1,epsilon_dt,...
    thetat,alphat,utm1,lambda,typen,A,B1,B2,epsilon_u)

q_theta_t=fun_q_theta(thetat,A,B1,B2);

F_epsilond_t = fun_F_x(epsilon_dt,typen,epsilon_u);

temp1 = hntm1;
temp2 = (utm1+hstm1)*alphat*thetat*q_theta_t;
temp3 = hntm1*(lambda*F_epsilond_t+(1-alphat)*thetat*q_theta_t);

h_n =temp1 + temp2 -temp3;

return