% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parameters=myParameter_Jun01V6(typen,epsilon_u)
%% set parameters
sprintf('Initialization parameters for type %s of function F(x)\n',typen)

% Parameter22PM for Type_II
A=1.355;
B1=18;
B2=25;
B=B1/B2;
%
r=0.012;
%case(1) r=(0.001:0.001:0.05)
c_p=0.113;
%case(2) c_p=(0.001:0.001:2)
beta=0.72;
%case(3) (0.1:0.01:0.7)
phi=0.1;
%case(4) (0.01:0.01:0.5)
delta=0.5;
%case(5) (0.1:0.01:0.8)
sigma=0.0436;
%case(6) (1:0.5:5)
lambda=0.1;
%case(7) (0.1:0.05:0.5)
pstar=1;
%case(8) (1:0.1:5)
b=0.2;
%case(9) (0.1:0.01:2)
c_f=0.213;
%case(10) (0.001:0.001:2)
mu=0.08;

Names={'A';'B';'r';'c^p';'beta';'phi';'delta';'sigma';'lambda';...
    'pstar';'b';'c^f';'mu';'epsilon_u'};
Values=[A;B;r;c_p;beta;phi;delta;sigma;lambda;pstar;b;c_f;mu;epsilon_u];
Tab=table(Values,'RowNames',Names);
disp(Tab);

parameters = ...
    {A;B;r;c_p;beta;phi;delta;sigma;lambda;pstar;b;c_f;mu;B1;B2};
return