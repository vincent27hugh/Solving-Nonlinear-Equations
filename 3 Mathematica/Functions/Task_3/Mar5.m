% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilon_d,epsilon_c,theta,alpha] = Mar5
global typen epsilon_u width step
n = int16(2*width/step+1);
% length of the lattice

F = @(var)fun_solveMar5_1(var);
% var = [epsilon_d;epsilon_c;theta]
% theta should be larger than 0
% var = [epsilon_d;epsilon_c;theta]
% theta should be larger than 0
t1 = linspace(-5,-500,n);
t2 = linspace(-5,-500,n);
t3 = linspace(0.5,5,n);
t4 = [t1;t2;t3];
var0=reshape(t4,[],1);
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','iter',...
    'PlotFcn',@optimplotfirstorderopt,...
    'StepTolerance',1e-10); 
% for fsolve function
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt',0.05
% Option to display output
% LM Algo is better
% 'FunctionTolerance',1e-16

[fs,~] = fsolve(F,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_d=fs(1:3:3*n);
epsilon_c=fs(2:3:3*n);
theta=fs(3:3:3*n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F2 = @(var)fun_solveMar5_2(var,epsilon_d,epsilon_c,theta);
% var = [alpha]
t = 0.5;
var0=repmat(t,n,1);
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','off'); 
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt'
% Option to display output

[fs2,~] = fsolve(F2,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]

alpha=fs2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n
    if ~strcmp(typen,'I')
        if epsilon_d(i)>=epsilon_u
            epsilon_d(i)=epsilon_u;
        end
        if epsilon_c(i)>=epsilon_u
            epsilon_c(i)=epsilon_u;
        end
    end
    if theta(i)<0
        theta(i)=0;
    end
    if alpha(i)<0
        alpha(i)=0;
    elseif alpha(i)>1
        alpha(i)=1;
    end
end
return