% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function intF = fun_int_F(a,b,typen,epsilon_u) 
% integration of (1-F(x))from a to b
% typen = 'i' 'ii' 'iii' 'iv'

temp = integral(@(x)(fun_F_x(x,typen,epsilon_u)),a,b);

intF = (b-a)-temp;


return