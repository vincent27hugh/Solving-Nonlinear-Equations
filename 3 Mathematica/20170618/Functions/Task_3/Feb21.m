% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilon_d,epsilon_c,theta,alpha] =Feb21
global step width 
n = 2*width/step+1;
% length of the lattice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = @(var)fun_solve33(var);
% var = [epsilon_d;epsilon_c;theta]
% theta should be larger than 0
t = [-0.5;-0.5;0.5;0.5];
var0=repmat(t,n,1);
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','off'); 
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt'
% Option to display output

[fs,~] = fsolve(F,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]

epsilon_d=fs(1:4:4*n);
epsilon_c=fs(2:4:4*n);
theta=fs(3:4:4*n);
alpha=fs(4:4:4*n);
return