% Mar18,2017
% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
function [epsilon_d,epsilon_c,theta,alpha,...
    ceq,flag,str_sol_alg] = Jun01T1(parameters)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;b;r;epsilon_u;mu;...
                pstar;typen};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver='fmincon';%%%%%%%%%%%
% define the solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
if strcmp(solver,'fsolve')
    % Solve system of nonlinear equations
    options = optimoptions('fsolve');
    %   Option to display output
    options.Algorithm = 'trust-region-dogleg';
    str_algo = 'TRD';
    %   'trust-region-dogleg' (default)
    %   'trust-region-reflective'
    %   'levenberg-marquardt'
    options.MaxFunctionEvaluations = Inf;
    options.MaxIterations = 500;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = @(var)fun_solveJun01T1(var,para);
    % var = [epsilon_d;epsilon_c;theta]
elseif strcmp(solver,'fmincon')
    % Find minimum of constrained nonlinear multivariable functions
    options = optimoptions(@fmincon);
    options.Algorithm = 'interior-point';
    str_algo = 'IP';
    %   'interior-point' (default)
    %   'trust-region-reflective'
    %   A gradient to be supplied in the objective function
    %   'sqp'
    %   'active-set'
    options.MaxFunctionEvaluations = Inf;
    options.MaxIterations = 500;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Afmin=[];bfmin=[];Aeqfmin=[];beqfmin=[];
    lbfmin = [-Inf;-Inf;-Inf;-Inf];
    if strcmp(typen,'I')
        ubfmin = [epsilon_u;epsilon_u;Inf;Inf];
    else
        ubfmin = [epsilon_u;epsilon_u;Inf;Inf];
    end
elseif strcmp(solver,'lsqnonlin')
    % Solve nonlinear least-squares (nonlinear data fitting problems)
    options = optimoptions(@lsqnonlin);
    options.Algorithm = 'trust-region-reflective';
    str_algo = 'TRR';
    %   'trust-region-reflective' (default)
    %   A gradient to be supplied in the objective function
    %   'levenberg-marquardt'
    options.MaxFunctionEvaluations = Inf;
    options.MaxIterations = 500;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = @(var)fun_solveJun01T1(var,para);
    % var = [epsilon_d;epsilon_c;theta]
    lblsq = [-Inf;-Inf;-Inf;-Inf];
    if strcmp(typen,'I')
        ublsq = [epsilon_u;epsilon_u;Inf;Inf];
    else
        ublsq = [epsilon_u;epsilon_u;Inf;Inf];
    end
end

options.Display = 'iter';
options.PlotFcn = @optimplotx;
%   @optimplotfval plots the function value.
%   @optimplotfirstorderopt plots the first-order optimality measure.
%   @optimplotx plots the current point.

options.FunctionTolerance = 1e-18;
options.StepTolerance = 1e-18;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial guess
v_epd=[-5;-1;-7];
v_epc=[-1;-5;-7];
v_theta=[1;5;7];

str_var0_test='off';
%%%%%
if strcmp(str_var0_test,'on')
    exitflag=-10;
    verror=NaN(numel(v_epd),numel(v_theta),numel(v_epc));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=1;
    if strcmp(solver,'fsolve')
        while exitflag<=0 && i<=length(v_epd)
            j=1;
            while exitflag<=0&&j<=length(v_theta)
                k=1;
                while exitflag<=0&&k<=length(v_epc)
                    var0 = [v_epd(i);v_epc(k);v_theta(j);0.5];
                    % epd;epc;theta;alpha

                    [fs,~,exitflag] = fsolve(F,var0,options);
                    % fs is solution of 
                    ceq=fun_solveJun01T1(fs,para);
              
                    verror(i,j,k)=sqrt(sum(ceq.^2));
                    k=k+1;
                end
                j=j+1;
            end
            i=i+1;
        end
    elseif strcmp(solver,'fmincon')||...
            strcmp(solver,'lsqnonlin')
        while exitflag~=1 && i<=length(v_epd)
            j=1;
            while exitflag~=1&&j<=length(v_theta)
                k=1;
                while exitflag~=1&&k<=length(v_epc)
                    var0 = [v_epd(i);v_epc(k);v_theta(j);0.5];
                    % epd;epc;theta;alpha

                    [fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,...
                        Aeqfmin,beqfmin,lbfmin,ubfmin,...
                        @(x)fun_solveJun01T1_nonlcon(x,para),options);

                    [~,ceq]=fun_solveJun01T1_nonlcon(fs,para);
                    
                    verror(i,j,k)=sqrt(sum(ceq.^2));
                    k=k+1;
                end
                j=j+1;
            end
            i=i+1;
        end
    elseif strcmp(solver,'lsqnonlin')
        while exitflag~=1 && i<=length(v_epd)
            j=1;
            while exitflag~=1&&j<=length(v_theta)
                k=1;
                while exitflag~=1&&k<=length(v_epc)
                    var0 = [v_epd(i);v_epc(k);v_theta(j);0.5];
                    % epd;epc;theta;alpha

                    [fs,~,~,exitflag]=lsqnonlin(F,var0,lblsq,ublsq,options);
                    ceq=fun_solveJun01T1(fs,para);
                    
                    verror(i,j,k)=sqrt(sum(ceq.^2));
                    k=k+1;
                end
                j=j+1;
            end
            i=i+1;
        end
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
    var0 = [v_epd(Ix);v_epc(Iz);v_theta(Iy);0.5];
    % epd;epc;theta;alpha
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    var0=[v_epd(1);v_epc(1);v_theta(1);0.5];
end

if strcmp(solver,'fsolve')
    [fs,~,exitflag] = fsolve(F,var0,options);
    % fs is solution of 
    ceq=fun_solveJun01T1(fs,para);
elseif strcmp(solver,'fmincon')
    [fs,~,exitflag] = fmincon(@(var)0,var0,Afmin,bfmin,Aeqfmin,beqfmin,lbfmin,ubfmin,...
        @(x)fun_solveJun01T1_nonlcon(x,para),options);
    [~,ceq]=fun_solveJun01T1_nonlcon(fs,para);
elseif strcmp(solver,'lsqnonlin')
    [fs,~,~,exitflag]=lsqnonlin(F,var0,lblsq,ublsq,options);
    ceq=fun_solveJun01T1(fs,para);
end
%%%%%
epsilon_d=fs(1);
epsilon_c=fs(2);
theta=fs(3);
alpha=fs(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(solver,'fsolve')
    flag= exitflag>0;
elseif strcmp(solver,'fmincon')||...
        strcmp(solver,'lsqnonlin')
    flag= exitflag==1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_sol_alg=strcat(solver,'-',str_algo);
return