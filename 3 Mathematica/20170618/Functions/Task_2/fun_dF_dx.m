% Updated in Mar18,2 017
% Task #20170302
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function dF_dx =fun_dF_dx (x,typen,epsilon_u)
if size(x,1)>1&&size(x,2)>1
    disp('Warning:input of fun_dF_dx must be a vector!');
    return;
end
dF_dx = NaN(size(x));

if strcmp(typen,'O')
    for i = 1:length(x)
        if x(i)<=epsilon_u
            dF_dx(i) = -6.89/...
                ((1/(1.72*x(i) - 1.91)^4 + 1.0)^2*(1.72*x - 1.91)^5);
        end
    end
elseif strcmp(typen,'I')
    for i = 1:length(x)
        dF_dx(i) = (0.5*heaviside(epsilon_u + x(i)))/epsilon_u -...
        1.0*dirac(1.0*x(i) - epsilon_u)*((0.5*(epsilon_u + x(i)))/epsilon_u - 1.0) - ...
        (0.5*heaviside(1.0*x(i) - 1.0*epsilon_u))/epsilon_u +...
        (0.5*dirac(epsilon_u + x(i))*(epsilon_u + x(i)))/epsilon_u;
    end
elseif strcmp(typen,'II')
    for i = 1:length(x)
        if x(i)<=epsilon_u 
            dF_dx(i) = 0.481*exp(-1.0*(0.426*x(i) - 0.564)^2);
        end
    end
elseif strcmp(typen,'III')
    for i = 1:length(x)
        if x(i)<epsilon_u
            dF_dx(i)=-(0.479*exp(-1.0*(0.849*log(1.0 - 1.0*x(i)) + 0.294)^2))...
                /(x(i) - 1.0);
        end
    end
elseif strcmp(typen,'IV')
    for i = 1:length(x)
        if x(i)<=epsilon_u
            dF_dx(i)=31.2/(2.0*x(i) - 3.0)^4;
        end
    end
end

return