% Jun 13, symbolic
% Updated in Mar18,2 017
% Task #20170302
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function F_x =fun_F_x(x,typen,epsilon_u)

if strcmp(typen,'O')
    temp1 = (pi/4)/sin(pi/4)-x;
    temp2 = (pi/2)/sin(pi/2)-((pi/4)/sin(pi/4))^2;
    temp3 = 1+(temp1/sqrt(temp2))^(-4);
    F_x = 1-(temp3)^(-1);
elseif strcmp(typen,'I')
    temp = (x+epsilon_u)/(2*epsilon_u);
   F_x = heaviside(x+epsilon_u)*(temp)+(1-temp)*heaviside(x-epsilon_u);
elseif strcmp(typen,'II')
    temp=sqrt(1/pi)-x*sqrt(.5-1/pi);
    F_x =  1-erf(temp);
elseif strcmp(typen,'III')
    F_x=.5-.5*erf((log(-x+1)+.5*log(2))/sqrt(2*log(2)));
elseif strcmp(typen,'IV')
    F_x =((3-2*x)/sqrt(3))^(-3);      
end
return