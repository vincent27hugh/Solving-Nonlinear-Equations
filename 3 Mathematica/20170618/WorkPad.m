close all

taskcode='Jun12T1';
loadcode='Jun01T1';
num_para='V27';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sym_solver='Jacobian_Jun13T1_III';
sym_solver='symbolic_Jun13T1_III';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_start =cputime;
% elasped CPU time when starting the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path0 = pwd;
if strcmp(sym_solver,'Jacobian_Jun13T1_III')
    add_path=strcat(path0,'\Functions\Task_1\');
    addpath(add_path);
    
    % mu = 0.08
    % step =0.1
    path_fig = strcat(path0,'\Figures\Task_1\');
    path_data = strcat(path0,'\Data\Task_1\');
elseif strcmp(sym_solver,'symbolic_Jun13T1_III')
    add_path_sym=strcat(path0,'\Functions\Task_1_sym\');
    addpath(add_path_sym);
    
    % mu = 0.08
    % step =0.1
    path_fig = strcat(path0,'\Figures\Task_1_sym\');
    path_data = strcat(path0,'\Data\Task_1_sym\');
end


str_para0=strcat('paraJun01',num_para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fsolve algorithm is 'interior-point'
str_algo='fmincon-IP';
str_para=strcat(str_para0,'-',str_algo);
%% 
cell_typen = {'I';'II';'III';'IV';'O'};% 
cell_case = {'pstar';'b';'phi';'sigma';'beta';'lambda';'c_f';'r';'c_p';'delta'};
% 10 cases

v_tt = 3;% 1:length(cell_typen);
v_cc = 1;%(1:10);

for tt = v_tt
    % F(x) type: 'O' 'I' 'II' 'III' 'IV'
    typen=cell_typen{tt};
    for c = v_cc
        close all
        % case code: 'p';'b';'phi';'sigma';'beta';'lambda';'c_f';'r';'c_p';'delta'
        casecode=cell_case{c};
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % epsilon_u is bar{epsilon}
        % change epsilon_u with different type of function of F(x)
        switch typen
            case 'I'
                epsilon_u=sqrt(3);
            case 'II'
                epsilon_u=sqrt(2)/(sqrt(pi-2));
            case 'III'
                epsilon_u=1;
            case 'IV'
                epsilon_u=(3-sqrt(3))/2;
            case 'O'
                epsilon_u=(pi/4)/(sin(pi/4));
        end
        str_func=strcat('myParameter_Jun01',num_para);
        func_handle=str2func(str_func);
        parameters=func_handle(typen,epsilon_u);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % parameters=[A;B;r;c_p;beta;phi;delta;sigma;lambda;pstar;b;c_f;mu;B1;B2];
        A=parameters{1};
        B=parameters{2};
        r=parameters{3};
        c_p=parameters{4};
        beta=parameters{5};
        phi=parameters{6};
        delta=parameters{7};
        sigma=parameters{8};
        lambda=parameters{9};
        pstar=parameters{10};
        b=parameters{11};
        c_f=parameters{12};
        mu=parameters{13};
        B1=parameters{14};
        B2=parameters{15};
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

epsilon_d=-5;epsilon_c=-1;theta=1;alpha=0.5;

Jacfvar = Jacobian_Jun13T1_III(A,B1,B2,mu,epsilon_u,...
    c_f,c_p,beta,phi,delta,sigma,lambda,b,r,pstar,...
    epsilon_d,epsilon_c,theta,alpha);

vpa(Jacfvar(:,2),5)
% dF/depsilond

dF_depsilonc=dF_depsilonc_fun(A,B1,B2,mu,epsilon_u,...
    c_f,c_p,beta,phi,delta,sigma,lambda,b,r,pstar,...
    epsilon_d,epsilon_c,theta,alpha);

vpa(dF_depsilonc',5)