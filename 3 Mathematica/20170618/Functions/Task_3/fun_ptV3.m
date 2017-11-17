% Task #20170302
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [pt,ynum] = fun_ptV3(T1,T2,Tb,p_1,ynum_1,...
    p_T1p1,ynum_T1p1,p_T2p1,ynum_T2p1)

% yt is coloumn vector
pt=NaN(Tb,1);
ynum=NaN(Tb,1);

pt(1)=p_1;
ynum(1)=ynum_1;

for i = 2:T1
    %prob is a single uniformly distributed random number in the interval (0,1).
    pt(i)=p_1;
    ynum(i)=ynum_1;
end
pt(T1+1)=p_T1p1;
ynum(T1+1)=ynum_T1p1;
for i = T1+2:T2
    %prob is a single uniformly distributed random number in the interval (0,1). 
    pt(i)=p_T1p1;
    ynum(i) = ynum_T1p1;
end
pt(T2+1)=p_T2p1;
ynum(T2+1)=ynum_T2p1;
for i = T2+2:Tb
    %prob is a single uniformly distributed random number in the interval (0,1). 
    pt(i)=p_T2p1;
    ynum(i) = ynum_T2p1;
end
return