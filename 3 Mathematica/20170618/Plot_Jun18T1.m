% function main_Jun01_Task1(num_para)
% June 12, 2017
% June 1,2017
% May 22-24, 2017
% 20170411 PM
% Mar18,2017
% Task #20170221
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all

taskcode='Jun12T1';
loadcode='Jun01T1';
num_para='V27';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sym_solver='Jacobian_Jun13T1_III';
% sym_solver='symbolic_Jun13T1_III';
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
v_cc = (1:10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % use the parameters of Mar4
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
        % case
        % change parameters with different cases
        % vt is column vector
        switch casecode
            case 'pstar'% 1
                vt=(0.18:0.01:2.5);
                str_s = 'pstar';
                str_var = 'p^*';
            case 'b' % 2
                vt=(0.01:0.01:1);
                str_s = 'b';
                str_var = 'b';
            case 'phi' % 3
                % vt=(0.05:0.01:0.5);
                vt=(0.3:0.01:0.95);
                str_s = 'phi';
                str_var = '\phi';
            case 'sigma' % 4
                vt=(0.01:0.01:1.2);
                str_s = 'sigma';
                str_var = '\sigma';
            case 'beta' % 5
                % vt=(0.1:0.05:0.7);
                vt=(0.4:0.01:0.99);
                str_s = 'beta';
                str_var = '\beta';
            case 'lambda' % 6
                vt=(0.02:0.01:0.63);
                str_s = 'lambda';
                str_var = '\lambda';
            case 'c_f' % 7
                % vt=(0.01:0.01:0.5);
                vt=(0.12:0.01:1.3);
                str_s = 'c_F';
                str_var = 'c^F';
            case 'r' % 8
                vt=(0.005:0.001:0.15);
                str_s = 'r';
                str_var = 'r';
            case 'c_p' % 9
                % vt=(0.01:0.005:0.2);
                vt=(0.01:0.001:0.2);
                str_s = 'c_P';
                str_var = 'c^P';
            case 'delta' % 10
                vt=(0.21:0.01:0.7);
                str_s = 'delta';
                str_var = '\delta';
        end
        
        str_fig=strcat(taskcode,'_',typen,'_',str_s);
        str_load=strcat(loadcode,'_',typen,'_',str_s);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(strcat(path_data,str_load,'_',str_para,'_ep_d.mat'),'vepsilon_d');
        load(strcat(path_data,str_load,'_',str_para,'_ep_c.mat'),'vepsilon_c');
        load(strcat(path_data,str_load,'_',str_para,'_theta.mat'),'vtheta');
        load(strcat(path_data,str_load,'_',str_para,'_alpha.mat'),'valpha');
            
        load(strcat(path_data,str_load,'_',str_para,'_ceq.mat'),'vceq');
        load(strcat(path_data,str_load,'_',str_para,'_exitflag.mat'),'vexitflag');
        
        load(strcat(path_data,str_load,'_',str_para,'_u.mat'),'vu');
        load(strcat(path_data,str_load,'_',str_para,'_hs.mat'),'vhs');
        load(strcat(path_data,str_load,'_',str_para,'_hn.mat'),'vhn');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        invA_B = NaN(4,numel(vt));
        % dF/dx
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:numel(valpha)
            %%%%%%%%%%%%%%%%%
            Dis1 = sprintf('We are coputing the case %s \nAnd function F(x) type %s...'...
                ,casecode,typen);
            Dis2 = sprintf('The parameters are:\n');
            Names={'epsilon_u';'A';'B';'r';'c^p';'beta';'phi';'delta';'sigma';'lambda';...
                'pstar';'b';'c^f'};
            Values=[epsilon_u;A;B;r;c_p;beta;phi;delta;sigma;lambda;pstar;b;c_f];
            Tab=table(Values,'RowNames',Names);
            fprintf(Dis1);
            fprintf(Dis2);
            disp(Tab);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            switch casecode
                case 'pstar' % 1
                    pstar=vt(i);
                case 'b' % 2
                    b=vt(i);
                case 'phi' % 3
                    phi=vt(i);
                case 'sigma' % 4
                    sigma=vt(i);
                case 'beta' % 5
                    beta=vt(i);
                case 'lambda' % 6
                    lambda=vt(i);
                case 'c_f' % 7
                    c_f=vt(i);
                case 'r' % 8
                    r=vt(i);
                case 'c_p' % 9
                    c_p=vt(i);
                case 'delta' % 10
                    delta = vt(i);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(sym_solver,'Jacobian_Jun13T1_III')
                Jacfvar = Jacobian_Jun13T1_III(A,B1,B2,mu,epsilon_u,...
                    c_f,c_p,beta,phi,delta,sigma,lambda,b,r,pstar,...
                    vepsilon_d(i),vepsilon_c(i),vtheta(i),valpha(i));
            elseif strcmp(sym_solver,'symbolic_Jun13T1_III')
                Jacfvar = symbolic_Jun13T1_III(A,B1,B2,mu,epsilon_u,...
                    c_f,c_p,beta,phi,delta,sigma,lambda,b,r,pstar,...
                    vepsilon_d(i),vepsilon_c(i),vtheta(i),valpha(i));
            end
        %{
        F = [part11-part12,part21-part22,part31 - part32,part41-part42];

        var=[epsilond,epsilonc,theta,alpha,...
            pstar,b,phi,sigma,beta,...
            lambda,cf,r,cp,delta];
        % [1,14]

        Jacfvar = jacobian(F,var);
        % jacobian matrix
        % [4,14]
        % [dF1_dx1,...dF1_dxN,
        %  
        %  dFN_dx1,.....dFN_dxN]
        %}
            MatB=Jacfvar(:,1:4);
            % [4,4] 
            MatA=Jacfvar(:,5:14);
            % [4,10]
            
            MatC = MatB\MatA; 
            % MatC2 = inv(MatB)*MatA;
            % [4,10]
            switch casecode
                case 'pstar' % 1
                    invA_B(:,i)=MatC(:,1);
                case 'b' % 2
                    invA_B(:,i)=MatC(:,2);
                case 'phi' % 3
                    invA_B(:,i)=MatC(:,3);
                case 'sigma' % 4
                    invA_B(:,i)=MatC(:,4);
                case 'beta' % 5
                    invA_B(:,i)=MatC(:,5);
                case 'lambda' % 6
                    invA_B(:,i)=MatC(:,6);
                case 'c_f' % 7
                    invA_B(:,i)=MatC(:,7);
                case 'r' % 8
                    invA_B(:,i)=MatC(:,8);
                case 'c_p' % 9
                    invA_B(:,i)=MatC(:,9);
                case 'delta' % 10
                    invA_B(:,i)=MatC(:,10);
            end
        end
        save(strcat(path_data,str_fig,'_',str_para,'_invA_B.mat'),'invA_B','-v7.3');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alpha_flag=fun_alphaflag(valpha);
        theta_flag=vtheta>0;
        flag = vexitflag==1&alpha_flag&theta_flag;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure;
        %
        s1=subplot(2,2,1);
        plot(valpha(flag),invA_B(1,flag),'LineWidth',2);
        xlabel('$\alpha$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        %ylabel(strcat('$\frac{\partial F_1}{\partial {',str_var,'}}$'),...
        %    'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel(strcat('$B^{-1}A(',str_var,')$'),...
            'interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(valpha(flag))
            xlim([min(valpha(flag)),max(valpha(flag))]);
        else
            xlim([min(valpha),max(valpha)]);
        end
        str_title = {strcat('$',taskcode,',',typen,',','$'),...
            strcat('$',str_para0,',','$'),...
            strcat('$',str_algo,'$')};
        title(str_title,'interpreter','Latex');
        %
        s2=subplot(2,2,2);
        plot(valpha(flag),invA_B(2,flag),'LineWidth',2);
        ylabel(strcat('$B^{-1}A(',str_var,')$'),...
            'interpreter','Latex','FontSize',18,'FontWeight','bold');
        xlabel('$\alpha$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(valpha(flag))
            xlim([min(valpha(flag)),max(valpha(flag))]);
        else
            xlim([min(valpha),max(valpha)]);
        end
        str_title = {strcat('$',taskcode,',',typen,',','$'),...
            strcat('$',str_para0,',','$'),...
            strcat('$',str_algo,'$')};
        title(str_title,'interpreter','Latex');
        %
        s3=subplot(2,2,3);
        plot(valpha(flag),invA_B(3,flag),'LineWidth',2);
        xlabel('$\alpha$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel(strcat('$B^{-1}A(',str_var,')$'),...
            'interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(valpha(flag))
            xlim([min(valpha(flag)),max(valpha(flag))]);
        else
            xlim([min(valpha),max(valpha)]);
        end
        str_title = {strcat('$',taskcode,',',typen,',','$'),...
            strcat('$',str_para0,',','$'),...
            strcat('$',str_algo,'$')};
        title(str_title,'interpreter','Latex');
        %
        s4=subplot(2,2,4);
        plot(valpha(flag),invA_B(4,flag),'LineWidth',2);
        xlabel('$\alpha$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel(strcat('$B^{-1}A(',str_var,')$'),...
            'interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(valpha(flag))
            xlim([min(valpha(flag)),max(valpha(flag))]);
        else
            xlim([min(valpha),max(valpha)]);
        end
        str_title = {strcat('$',taskcode,',',typen,',','$'),...
            strcat('$',str_para0,',','$'),...
            strcat('$',str_algo,'$')};
        title(str_title,'interpreter','Latex');
        %
        grid(s1,'on');grid(s2,'on');grid(s3,'on');
        grid(s4,'on');
        %
        set(gcf,'Position',get(0,'Screensize'));%Maximize figure
        str_fig2 = strcat(path_fig,str_fig,'_',str_para,'_invA_B_flag');
        saveas(gcf,str_fig2,'tif');
        savefig(gcf,str_fig2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
    end
end

close all
time_passby=cputime-time_start;
DCPU = sprintf...
    ('\n Total CPU Time is :%12.3f \n'...
    , time_passby);
disp(DCPU);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Send email to notify the fufillment
to = 'huweihugh@gmail.com';
% the efrom address you want to send message to
% you can send to multiple email address, use ; to seperate
subject = strcat('Program_',mfilename,'_done') ;

txt1 = sprintf('%s is done and total cputime is: %12.3f seconds! \n'...
    ,mfilename,time_passby);

txt2=[];
for ii=1:numel(v_tt)
    txt2 = strcat(txt2,'types-',cell_typen(v_tt(ii)),';');
end
txt2r = sprintf('\n The types of function F(x) include:\n %s \n',char(txt2));
txt3=[];
for ii = 1:numel(v_cc)
    txt3 = strcat(txt3,'case-',cell_case(v_cc(ii)),';');
end
txt3r = sprintf('\n The cases of variables include:\n %s \n',char(txt3));

txt4 = sprintf('\n This is automatically sending mail, please DO NOT reply!\n');

message = {txt1;txt2r;txt3r;txt4};

attachment_path = 'D:\GN Luo\20170522\Figures\Task_1\';
filename = 'Apr11_I_b_paraApr11-IP_solution.tif';
attachment = strcat(attachment_path,filename);
% 'folder/attach1.html','attach2.doc'

mysendemail(to,subject,message);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beep 
pause(1)
mydialog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(sym_solver,'Jacobian_Jun13T1_III')
    rmpath(add_path);
elseif strcmp(sym_solver,'symbolic_Jun13T1_III')
    rmpath(add_path_sym);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return