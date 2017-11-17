function ind = fun_indicator(x1,x2,type)
% ind = 1, if x1>=x2
% ind = 0, if x1<x2
%
if nargin==1
    fprintf('Number of input should be: %d\n',2);
    return
elseif nargin==2
    ind = heaviside(x1-x2);
    if ind == 0.5
        ind = 1;
    end
elseif nargin==3
    ind = heaviside(x1-x2);
    if ind == 0.5
        if strcmp(type,'NaN')
            ind = NaN;
        elseif strcmp(type,'zero')
            ind = 0;
        elseif strcmp(type,'one')
            ind = 1;
        end
    end
end
return