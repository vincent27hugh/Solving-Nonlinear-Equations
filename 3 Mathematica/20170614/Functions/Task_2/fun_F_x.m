% Updated in Mar18,2 017
% Task #20170302
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function F_x =fun_F_x (x,typen,epsilon_u)
if size(x,1)>1&&size(x,2)>1
    disp('Warning:input of fun_F_x must be a vector!');
    return;
end
F_x = NaN(size(x));

if strcmp(typen,'O')
    for i = 1:length(x)
        if x(i)<=epsilon_u
            temp1 = (pi/4)/sin(pi/4)-x(i);
            temp2 = (pi/2)/sin(pi/2)-((pi/4)/sin(pi/4))^2;
            temp3 = 1+(temp1/sqrt(temp2))^(-4);
            F_x(i) = 1-(temp3)^(-1);
        end
    end
elseif strcmp(typen,'I')
    for i = 1:length(x)
       if x(i)<-epsilon_u
           F_x(i)=0;
       elseif x(i)>epsilon_u
           F_x(i)=1;
       else
           F_x(i)=(x(i)+epsilon_u)/(2*epsilon_u);
       end
    end
elseif strcmp(typen,'II')
    for i = 1:length(x)
        if x(i)<=epsilon_u
            temp=sqrt(1/pi)-x(i)*sqrt(.5-1/pi);
            F_x(i) = 1-erf(temp);
        end
    end
elseif strcmp(typen,'III')
    for i = 1:length(x)
        if x(i)<epsilon_u
            F_x(i)=.5-.5*erf((log(-x(i)+1)+.5*log(2))/sqrt(2*log(2)));
        end
    end
elseif strcmp(typen,'IV')
    for i = 1:length(x)
        if x(i)<=epsilon_u
            F_x(i) = ((3-2*x(i))/sqrt(3))^(-3);
        end
    end
end

return