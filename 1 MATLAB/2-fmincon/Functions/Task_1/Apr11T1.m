% Mar18,2017
% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilon_d,epsilon_c,theta,alpha,ceq,exitflag] = Apr11T1(parameters)
dbstop if error
%{ 
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;step;...
                width;pstar;typen};
%}
A = parameters{1};
B1= parameters{2};
B2 = parameters{3};
c_f = parameters{4};
c_p = parameters{5};
beta = parameters{6};
phi = parameters{7};
delta = parameters{8};
sigma = parameters{9};
lambda = parameters{10};
b = parameters{11};
r = parameters{12};
epsilon_u = parameters{13};
mu = parameters{14};
pstar = parameters{15};
typen = parameters{16};
taskcode = parameters{17};
str_s = parameters{18};
str_para = parameters{19};
path_fig = parameters{20};
variable = parameters{21};
num = parameters{22};
%%%%%%
exitflag=-10;
v_epd=linspace(-1e3,epsilon_u-1e-3,5);
v_theta=linspace(1e-1,1e3,5);
verror=NaN(numel(v_epd),numel(v_theta));
i=1;
while exitflag~=1 && i<=length(v_epd)
    j=1;
    while exitflag~=1&&j<=length(v_theta)
        var0 = [v_epd(i);v_epd(i);v_theta(j);0.5];
        % epd;epc;theta;alpha
        options = optimoptions(@fmincon,'Algorithm','interior-point',...
            'Display','iter',...
             'MaxFunctionEvaluations',Inf,...
             'MaxIteration',300,...
             'ConstraintTolerance',1e-9);
        %{
             'interior-point' (default)
        'trust-region-reflective'
         % A gradient to be supplied in the objective function
        'sqp'
        'active-set'
             
         'PlotFcn',@optimplotx,...
             %}
        Afmin=[];bfmin=[];Aeqfmin=[];beqfmin=[];
        lbfmin = [-Inf;-Inf;0;0];
        if strcmp(typen,'I')
            ubfmin = [epsilon_u;epsilon_u;Inf;1];
        else
            ubfmin = [epsilon_u;epsilon_u;Inf;1];
        end
        %%%%
        para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;...
                        pstar;typen};
        [fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,Aeqfmin,beqfmin,lbfmin,ubfmin,...
            @(x)fun_solveMar4T1_nonlcon(x,para),options);

        [~,ceq]=fun_solveMar4T1_nonlcon(fs,para);

        verror(i,j)=sqrt(sum(ceq.^2));
        
        j=j+1;
    end
    i=i+1;
end

[M,I]=min(verror(:));
[Ir,Ic]=ind2sub(size(verror),I);
% Subscripts from linear index
    
plot_error='on';
if strcmp(plot_error,'on')
    f1=figure;
    hold on 
    surf(v_theta,v_epd,verror);
    view([0,90]);
    colorbar;
    
    scatter3(v_theta(Ic),v_epd(Ir),M,'filled','r');
    xlabel('$\theta$','interpreter','Latex','FontSize',18,'FontWeight','bold');
    ylabel('$\epsilon^d$','interpreter','Latex','FontSize',18,'FontWeight','bold');
    zlabel('$error$','interpreter','Latex','FontSize',18,'FontWeight','bold');
    hold off
    
    str_title = strcat('$',taskcode,',',typen,',',...
        str_s,'=',num2str(variable(num)),',',str_para,'$');
    title(str_title,'interpreter','Latex');
    
    legend('Error','Minimum','Location','northeastoutside')
    %%%%%%%
    set(gcf,'Position',get(0,'Screensize'));%Maximize figure
    str_fig2 = strcat(path_fig,taskcode,'_',typen,'_',str_s,'_',str_para,'_fimconTest_',num2str(num));
    saveas(gcf,str_fig2,'tif');
    savefig(gcf,str_fig2);
    
    close(f1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var0 = [v_epd(Ir);v_epd(Ir);v_theta(Ic);0.5];
% epd;epc;theta;alpha
options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'Display','iter',...
     'PlotFcn',@optimplotx,...
     'MaxFunctionEvaluations',Inf,...
     'MaxIteration',300,...
     'ConstraintTolerance',1e-9);
%{
     'interior-point' (default)
'trust-region-reflective'
 % A gradient to be supplied in the objective function
'sqp'
'active-set'
     %}
Afmin=[];bfmin=[];Aeqfmin=[];beqfmin=[];
lbfmin = [-Inf;-Inf;0;0];
if strcmp(typen,'I')
    ubfmin = [epsilon_u;epsilon_u;Inf;1];
else
    ubfmin = [epsilon_u;epsilon_u;Inf;1];
end
%%%%
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;...
                pstar;typen};
[fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,Aeqfmin,beqfmin,lbfmin,ubfmin,...
    @(x)fun_solveMar4T1_nonlcon(x,para),options);
%%%%%
epsilon_d=fs(1);
epsilon_c=fs(2);
theta=fs(3);
alpha=fs(4);
% size(ceq) = [2*n,1]
%{
F = @(var)fun_solveMar4_1(var);
% var = [epsilon_d;epsilon_c;theta]
% theta should be larger than 0
var0 = [-100;-10;5];
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','iter',...
    'PlotFcn',@optimplotx); 
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt'
% Option to display output

% @optimplotfval plots the function value.
% @optimplotfirstorderopt plots the first-order optimality measure.
% @optimplotx plots the current point.

[fs,~] = fsolve(F,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]

epsilon_d=fs(1);
epsilon_c=fs(2);
theta=fs(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F2 = @(var)fun_solveMar4_2(var,epsilon_d,epsilon_c,theta);
% var = [alpha]

var0 = 0.5;
% initial condition
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','iter'); 
% 'trust-region-dogleg' (default)
% 'trust-region-reflective'
% 'levenberg-marquardt'
% Option to display output

[fs2,~] = fsolve(F2,var0,options);
% fs is solution of [epsilon_d;epsilon_c;theta;alpha]

alpha=fs2;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
if alpha<0
    alpha=0;
elseif alpha>1
    alpha=1;
end
%}
return