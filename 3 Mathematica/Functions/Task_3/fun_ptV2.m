% Task #20170302
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [pt,ynum] = fun_ptV2(prob,width,step,T,Tb,p_1,ynum_1,p_Tp1,ynum_Tp1,pstar,mu)

% yt is coloumn vector
pt=NaN(Tb,1);
ynum=NaN(Tb,1);

pt(1)=p_1;
ynum(1)=ynum_1;

for i = 2:T
    %prob is a single uniformly distributed random number in the interval (0,1).
    pt(i)=p_1;
    ynum(i)=ynum_1;
end
pt(T+1)=p_Tp1;
ynum(T+1)=ynum_Tp1;
for i = T+2:Tb
    %prob is a single uniformly distributed random number in the interval (0,1).
    g_p = fun_gp(log(pt(i-1)/pstar),width);
    if prob(i)< mu*g_p
        pt(i) = exp(log(pt(i-1))+step);
        ynum(i) = ynum(i-1)+1;
    elseif prob(i)>=mu*g_p&&prob(i)<mu
        pt(i) = exp(log(pt(i-1))-step);
        ynum(i) = ynum(i-1)-1;
    elseif prob(i)>=mu
        pt(i)=pt(i-1);
        ynum(i) = ynum(i-1);
    end
end
return