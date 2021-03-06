% function main_Jun01_Task1(num_para)
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

taskcode='Jun01T1';
% taskcode='May22T4_1';
num_para='V21';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_start =cputime;
% elasped CPU time when starting the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path0 = pwd;
add_path=strcat(path0,'\Functions\Task_1\');
addpath(add_path);

if strcmp(taskcode,'Jun01T1')
    % mu = 0.08
    % step =0.1
    path_fig = strcat(path0,'\Figures\Task_1\');
    path_data = strcat(path0,'\Data\Task_1\');
    str_para0=strcat('paraJun01',num_para);
elseif strcmp(taskcode,'May22T4_1')
    % mu = 0.012
    % step = 0.15
    path_fig = strcat(path0,'\Figures\Task_4\',...
        'Task_1\');
    path_data = strcat(path0,'\Data\Task_4\',...
        'Task_1\');
    str_para0='paraMay22';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fsolve algorithm is 'interior-point'
str_algo='LM';
str_para=strcat(str_para0,'-',str_algo);
%% 
cell_typen = {'I';'II';'III';'IV';'O'};% 
cell_case = {'pstar';'b';'phi';'sigma';'beta';'lambda';'c_f';'r';'c_p';'delta'};
% 10 cases

v_tt = 3;% 1:length(cell_typen);
v_cc = 7;% [1;2;3;4;5;6];
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
        if strcmp(taskcode,'Jun01T1')
            % use the parameters of Mar4
            str_func=strcat('myParameter_Jun01',num_para);
            func_handle=str2func(str_func);
            parameters=func_handle(typen,epsilon_u);
            
        elseif strcmp(taskcode,'May22T4_1')
            % use the parameters of May22
            parameters=myParameter_May22T2(typen,epsilon_u);
        end
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
                vt=(0.5:0.01:2);
                str_s = 'pstar';
                str_var = '$p^*$';
            case 'b' % 2
                vt=(0.1:0.01:0.5);
                str_s = 'b';
                str_var = '$b$';
            case 'phi' % 3
                vt=(0.05:0.01:0.5);
                str_s = 'phi';
                str_var = '$\phi$';
            case 'sigma' % 4
                vt=(0.01:0.01:0.5);
                str_s = 'sigma';
                str_var = '$\sigma$';
            case 'beta' % 5
                vt=(0.1:0.05:0.7);
                str_s = 'beta';
                str_var = '$\beta$';
            case 'lambda' % 6
                vt=(0.01:0.01:0.5);
                str_s = 'lambda';
                str_var = '$\lambda$';
            case 'c_f' % 7
                % vt=(0.01:0.01:0.5);
                vt=(0.12:0.01:0.5);
                str_s = 'c_F';
                str_var = '$c^F$';
            case 'r' % 8
                vt=(0.005:0.001:0.03);
                str_s = 'r';
                str_var = '$r$';
            case 'c_p' % 9
                vt=(0.01:0.005:0.2);
                str_s = 'c_P';
                str_var = '$c^P$';
            case 'delta' % 10
                vt=(0.1:0.01:0.7);
                str_s = 'delta';
                str_var = '$\delta$';
        end
        str_fig=strcat(taskcode,'_',typen,'_',str_s);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vepsilon_d=NaN(size(vt));
        vepsilon_c=NaN(size(vt));
        vtheta=NaN(size(vt));
        valpha=NaN(size(vt));
        vexitflag=NaN(size(vt));
        
        vceq=NaN(4,numel(vt));
        
        vu=NaN(size(vt));
        vhs=NaN(size(vt));
        vhn=NaN(size(vt));
        for i = 1:numel(vt)
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
            %%%%%%%%%%% (i)
            % the most important part
            para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;...
                b;r;epsilon_u;mu;pstar;typen;taskcode;...
                str_s;str_para;path_fig;vt;i};
            [vepsilon_d(i),vepsilon_c(i),vtheta(i),valpha(i),vceq(:,i),...
                vexitflag(i)]=Jun01T1(para);
            % solve the equations and get the solution
            

            vu(i) = fun_u(vepsilon_d(i),vtheta(i),lambda,typen,...
                A,B1,B2,epsilon_u);
            vhs(i)=fun_hs(vepsilon_c(i),vepsilon_d(i),vtheta(i),valpha(i),...
                vu(i),lambda,typen,A,B1,B2,epsilon_u);
            vhn(i)=fun_hn(vepsilon_c(i),vtheta(i),valpha(i),vu(i),...
                lambda,typen,A,B1,B2,epsilon_u);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save variables
        save(strcat(path_data,str_fig,'_',str_para,'_ep_d.mat'),'vepsilon_d','-v7.3');
        save(strcat(path_data,str_fig,'_',str_para,'_ep_c_1s.mat'),'vepsilon_c','-v7.3');
        save(strcat(path_data,str_fig,'_',str_para,'_the_1s.mat'),'vtheta','-v7.3');
        save(strcat(path_data,str_fig,'_',str_para,'_alp_1s.mat'),'valpha','-v7.3');
            
        save(strcat(path_data,str_fig,'_',str_para,'_ceq_1s.mat'),'vceq','-v7.3');
        
        save(strcat(path_data,str_fig,'_',str_para,'_u.mat'),'vu','-v7.3');
        save(strcat(path_data,str_fig,'_',str_para,'_hs.mat'),'vhs','-v7.3');
        save(strcat(path_data,str_fig,'_',str_para,'_hn.mat'),'vhn','-v7.3');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alpha_flag=fun_alphaflag(valpha);
        theta_flag=vtheta>0;
        flag = vexitflag>=0&alpha_flag&theta_flag;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure;
        %
        s1=subplot(2,3,1);
        q_theta=fun_q_theta(vtheta,A,B1,B2);
        tqt=vtheta.*q_theta;
        plot(vt(flag),tqt(flag),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$\theta q(\theta)$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s2=subplot(2,3,2);
        plot(vt(flag),valpha(flag),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$\alpha$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s3=subplot(2,3,3);
        vF_epd=fun_F_x(vepsilon_d,typen,epsilon_u);
        vF_epc=fun_F_x(vepsilon_c,typen,epsilon_u);
        hold on 
        plot(vt(flag),vF_epd(flag),'LineWidth',2);
        plot(vt(flag),vF_epc(flag),'LineWidth',2);
        hold off
        leg = legend('$F(\epsilon^d)$','$F(\epsilon^c)$');
        set(leg,'interpreter','Latex');
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$F(\epsilon)$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s4=subplot(2,3,4);
        plot(vt(flag),vu(flag),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$u$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s5=subplot(2,3,5);
        plot(vt(flag),vhs(flag),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$h_s$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s6=subplot(2,3,6);
        plot(vt(flag),vhn(flag),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$h_n$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        grid(s1,'on');grid(s2,'on');grid(s3,'on');
        grid(s4,'on');grid(s5,'on');grid(s6,'on');
        %
        set(gcf,'Position',get(0,'Screensize'));%Maximize figure
        str_fig2 = strcat(path_fig,str_fig,'_',str_para,'_sol_flag');
        saveas(gcf,str_fig2,'tif');
        savefig(gcf,str_fig2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure;
        %
        s1=subplot(2,3,1);
        q_theta=fun_q_theta(vtheta,A,B1,B2);
        tqt=vtheta.*q_theta;
        plot(vt(vexitflag>=0),tqt(vexitflag>=0),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$\theta q(\theta)$','interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        if ~isempty(vt(vexitflag>=0))
            xlim([min(vt(vexitflag>=0)),max(vt(vexitflag>=0))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s2=subplot(2,3,2);
        plot(vt(vexitflag>=0),valpha(vexitflag>=0),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$\alpha$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(vexitflag>=0))
            xlim([min(vt(vexitflag>=0)),max(vt(vexitflag>=0))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s3=subplot(2,3,3);
        vF_epd=fun_F_x(vepsilon_d,typen,epsilon_u);
        vF_epc=fun_F_x(vepsilon_c,typen,epsilon_u);
        hold on 
        plot(vt(vexitflag>=0),vF_epd(vexitflag>=0),'LineWidth',2);
        plot(vt(vexitflag>=0),vF_epc(vexitflag>=0),'LineWidth',2);
        hold off
        leg = legend('$F(\epsilon^d)$','$F(\epsilon^c)$');
        set(leg,'interpreter','Latex');
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$F(\epsilon)$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(vexitflag>=0))
            xlim([min(vt(vexitflag>=0)),max(vt(vexitflag>=0))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s4=subplot(2,3,4);
        plot(vt(vexitflag>=0),vu(vexitflag>=0),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$u$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(vexitflag>=0))
            xlim([min(vt(vexitflag>=0)),max(vt(vexitflag>=0))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s5=subplot(2,3,5);
        plot(vt(vexitflag>=0),vhs(vexitflag>=0),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$h_s$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(vexitflag>=0))
            xlim([min(vt(vexitflag>=0)),max(vt(vexitflag>=0))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s6=subplot(2,3,6);
        plot(vt(vexitflag>=0),vhn(vexitflag>=0),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$h_n$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(vexitflag>=0))
            xlim([min(vt(vexitflag>=0)),max(vt(vexitflag>=0))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        grid(s1,'on');grid(s2,'on');grid(s3,'on');
        grid(s4,'on');grid(s5,'on');grid(s6,'on');
        %
        set(gcf,'Position',get(0,'Screensize'));%Maximize figure
        str_fig2 = strcat(path_fig,str_fig,'_',str_para,'_solution');
        saveas(gcf,str_fig2,'tif');
        savefig(gcf,str_fig2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ERRORS
        figure;
        nFn=(1:4);
        surf(vt,nFn,abs(vceq));
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('Fn');zlabel('Error');
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        set(gcf,'Position',get(0,'Screensize'));
        str_fig2 = strcat(path_fig,str_fig,'_',str_para,'_ceq');
        saveas(gcf,str_fig2,'tif');
        savefig(gcf,str_fig2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLAGS
        figure;
        s1 = subplot(2,2,1);
        plot(vt,vexitflag>=0,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$exit flag$','interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s2 = subplot(2,2,2);
        plot(vt,alpha_flag,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$flag:\alpha\subseteq(0,1)$','interpreter',...
            'Latex','FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s3 = subplot(2,2,3);
        plot(vt,theta_flag,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$flag:\theta\subseteq(0,\infty)$','interpreter',...
            'Latex','FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        s4 = subplot(2,2,4);
        plot(vt,flag,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$flag:\exists solution $','interpreter',...
            'Latex','FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = strcat('$',taskcode,',',typen,',',str_para,'$');
        title(str_title,'interpreter','Latex');
        %
        grid(s1,'on');grid(s2,'on');grid(s3,'on');grid(s4,'on');
        %
        set(gcf,'Position',get(0,'Screensize'));%Maximize figure
        str_fig2 = strcat(path_fig,str_fig,'_',str_para,'_exitflag');
        saveas(gcf,str_fig2,'tif');
        savefig(gcf,str_fig2);
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
%{
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
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
beep 
pause(1)
mydialog
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath(add_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return