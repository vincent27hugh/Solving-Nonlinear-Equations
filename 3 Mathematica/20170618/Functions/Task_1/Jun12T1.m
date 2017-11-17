% Mar18,2017
% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilon_d,epsilon_c,theta,sol,ceq,exitflag] = Jun12T1(parameters)
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

alpha = parameters{23};
casecode=parameters{24};
%%%%%%
%{
%%%%%%%%%%%%%%%%%%
% OPTIONS of fsolve
options = optimoptions('fsolve','algorithm','levenberg-marquardt'...
    ,'Display','iter',...
    'PlotFcn',@optimplotx,...
    'MaxFunctionEvaluations',Inf,...
     'MaxIteration',1500,...
     'FunctionTolerance',1e-18,...
     'StepTolerance',1e-18);
%       'trust-region-dogleg' (default)
%       'trust-region-reflective'
%       'levenberg-marquardt'
%       Option to display output
%       @optimplotfval plots the function value.
%       @optimplotfirstorderopt plots the first-order optimality measure.
%       @optimplotx plots the current point.
%}
%%%%%%%%%%%%%%%%%%
% OPTIONS of fmincon
options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'Display','iter',...
    'PlotFcn',@optimplotx,...
     'MaxFunctionEvaluations',Inf,...
     'MaxIteration',500,...
     'ConstraintTolerance',1e-9);
%       'interior-point' (default)
%       'trust-region-reflective'
%       A gradient to be supplied in the objective function
%       'sqp'
%       'active-set'    
%%%%%%%%%%%%%%%%
% For fmoncon
Afmin=[];bfmin=[];Aeqfmin=[];beqfmin=[];
lbfmin = [-Inf;-Inf;0;-Inf];
if strcmp(typen,'I')
    ubfmin = [epsilon_u;epsilon_u;Inf;Inf];
else
    ubfmin = [epsilon_u;epsilon_u;Inf;Inf];
end
%%%%
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;...
                pstar;typen;casecode;alpha};
%{
% fsolve
F = @(var)fun_solveJun12T1(var,para);
% var = [epsilon_d;epsilon_c;theta]
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch casecode
    case 'pstar' % 1
        sol0=pstar;
    case 'b' % 2
        sol0=b;
    case 'phi' % 3
        sol0=phi;
    case 'sigma' % 4
        sol0=sigma;
    case 'beta' % 5
        sol0=beta;
    case 'lambda' % 6
        sol0=lambda;
    case 'c_f' % 7
        sol0=c_f;
    case 'r' % 8
        sol0=r;
    case 'c_p' % 9
        sol0=c_p;
    case 'delta' % 10
        sol0=delta;
end
%%%%%%%%%%%%%%%%%%
v_epd=[-5;-1;-7];
v_epc=[-1;-5;-7];
v_theta=[1;5;7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_var0_test='off';
%%%%%
if strcmp(str_var0_test,'on')
    exitflag=-10;
    
    verror=NaN(numel(v_epd),numel(v_theta),numel(v_epc));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=1;
    % for fsolve
    %while exitflag<=0 && i<=length(v_epd)
    % for fmincon
    while exitflag~=1 && i<=length(v_epd)
        j=1;
        % fsolve
        %while exitflag<=0&&j<=length(v_theta)
        % fmincon
        while exitflag~=1 && j<=length(v_theta)
            k=1;
            % fsolve
            %while exitflag<=0&&k<=length(v_epc)
            while exitflag~=1 && k<=length(v_epc)
                var0 = [v_epd(i);v_epc(k);v_theta(j);sol0];
                % epd;epc;theta;alpha
                %{
                [fs,~,exitflag] = fsolve(F,var0,options);
                % fs is solution of 
                ceq=fun_solveJun12T1(fs,para);
                %}
                 [fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,...
                    Aeqfmin,beqfmin,lbfmin,ubfmin,...
                    @(x)fun_solveJun12T1_nonlcon(x,para),options);

                [~,ceq]=fun_solveJun12T1_nonlcon(fs,para);

                verror(i,j,k)=sqrt(sum(ceq.^2));

                k=k+1;
            end

            j=j+1;
        end
        i=i+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [MinValue,Index]=min(verror(:));
    [Ix,Iy,Iz]=ind2sub(size(verror),Index);
    % Subscripts from linear index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    plot_error='off';
    if strcmp(plot_error,'on')
        f1=figure;
        hold on 
        surf(v_theta,v_epd,verror);
        view([0,90]);
        colorbar;

        scatter3(v_theta(Iy),v_epd(Ix),MinValue,'filled','r');
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
        str_fig2 = strcat(path_fig,taskcode,'_',typen,'_',str_s,'_',...
            str_para,'_fimconTest_',num2str(num));
        saveas(gcf,str_fig2,'tif');
        savefig(gcf,str_fig2);

        close(f1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    var0 = [v_epd(Ix);v_epc(Iz);v_theta(Iy);sol0];
    % epd;epc;theta;alpha
else
    var0=[v_epd(1);v_epc(1);v_theta(1);sol0];
end

% Fmincon
[fs,~,exitflag] = fmincon(@(var)1,var0,Afmin,bfmin,Aeqfmin,beqfmin,lbfmin,ubfmin,...
    @(x)fun_solveJun12T1_nonlcon(x,para),options);
[~,ceq]=fun_solveJun12T1_nonlcon(fs,para);

%{
% fsolve
[fs,~,exitflag] = fsolve(F,var0,options);
% fs is solution of 
ceq=fun_solveJun12T1(fs,para);
%}
%%%%%
epsilon_d=fs(1);
epsilon_c=fs(2);
theta=fs(3);
sol=fs(4);
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