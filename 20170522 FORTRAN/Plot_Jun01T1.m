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

taskcode='Jun01T1FOR';
loadcode='Jun01T1';
num_para='V27';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_start =cputime;
% elasped CPU time when starting the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path0 = pwd;

path_fig = strcat(path0,'\Figures\Task_1\');
path_data = strcat(path0,'\Task_1\Task_1\');
str_para0=strcat('paraJun01_',num_para);
str_algo='NEQNF-PH';%Powell Hybrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
cell_typen = {'I';'II';'III';'IV';'O'};% 
cell_case = {'pstar';'b';'phi';'sigma';'beta';'lambda';'c_f';'r';'c_p';'delta'};
% 10 cases

v_tt = 3;% 1:length(cell_typen);
v_cc = (9:10);
for tt = v_tt
    % F(x) type: 'O' 'I' 'II' 'III' 'IV'
    typen=cell_typen{tt};
    for c = v_cc
        close all
        % case code: 'p';'b';'phi';'sigma';'beta';'lambda';'c_f';'r';'c_p';'delta'
        casecode=cell_case{c};
      
        %%
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % case
        % change parameters with different cases
        % vt is column vector
        switch casecode
            case 'pstar'% 1
                vt=(0.18:0.01:2.5);
                str_s = 'pstar';
                str_var = '$p^*$';
            case 'b' % 2
                vt=(0.01:0.01:1);
                str_s = 'b';
                str_var = '$b$';
            case 'phi' % 3
                vt=(0.3:0.01:0.95);
                str_s = 'phi';
                str_var = '$\phi$';
            case 'sigma' % 4
                vt=(0.01:0.01:1.2);
                str_s = 'sigma';
                str_var = '$\sigma$';
            case 'beta' % 5
                vt=(0.4:0.01:0.99);
                str_s = 'beta';
                str_var = '$\beta$';
            case 'lambda' % 6
                vt=(0.01:0.01:1.2);
                str_s = 'lambda';
                str_var = '$\lambda$';
            case 'c_f' % 7
                vt=(0.12:0.01:1.3);
                str_s = 'c_F';
                str_var = '$c^F$';
            case 'r' % 8
                vt=(0.005:0.001:0.15);
                str_s = 'r';
                str_var = '$r$';
            case 'c_p' % 9
                vt=(0.01:0.001:0.2);
                str_s = 'c_P';
                str_var = '$c^P$';
            case 'delta' % 10
                vt=(0.01:0.01:0.7);
                str_s = 'delta';
                str_var = '$\delta$';
        end
        str_fig=strcat(taskcode,'_',typen,'_',str_s);
        str_load=strcat(loadcode,'_',typen,'_',str_s);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_ep_d.txt'));
        dat=textscan(fid,'%f');
        vepsilon_d=dat{1};
        
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_ep_c.txt'));
        dat=textscan(fid,'%f');
        vepsilon_c=dat{1};
        
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_theta.txt'));
        dat=textscan(fid,'%f');
        vtheta=dat{1};
        
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_alpha.txt'));
        dat=textscan(fid,'%f');
        valpha=dat{1};
        
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_exitflag.txt'));
        dat=textscan(fid,'%s');
        vexitflag=strcmp(dat{1},'T');
        
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_u.txt'));
        dat=textscan(fid,'%f');
        vu=dat{1};
        
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_hs.txt'));
        dat=textscan(fid,'%f');
        vhs=dat{1};
        
        fid=fopen(strcat(path_data,str_load,'_',str_para0,'_hn.txt'));
        dat=textscan(fid,'%f');
        vhn=dat{1};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alpha_flag=(valpha>0)&(valpha<1);
        theta_flag=vtheta>0;
        flag = (vexitflag==1)&alpha_flag&theta_flag;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure;
        %
        s1=subplot(2,3,1);
        %q_theta=fun_q_theta(vtheta,A,B1,B2);
        %tqt=vtheta.*q_theta;
        plot(vt(flag),vtheta(flag),'LineWidth',2);
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$\theta$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = {strcat('$',taskcode,',',typen,'$'),...
            strcat('$',str_para0,'$'),...
            strcat('$',str_algo,'$')};
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
        str_title = {strcat('$',taskcode,',',typen,'$'),...
            strcat('$',str_para0,'$'),...
            strcat('$',str_algo,'$')};
        title(str_title,'interpreter','Latex');
        %
        s3=subplot(2,3,3);
        %vF_epd=fun_F_x(vepsilon_d,typen,epsilon_u);
        %vF_epc=fun_F_x(vepsilon_c,typen,epsilon_u);
        hold on 
        plot(vt(flag),vepsilon_d(flag),'LineWidth',2);
        plot(vt(flag),vepsilon_c(flag),'LineWidth',2);
        hold off
        leg = legend('$\epsilon^d$','$\epsilon^c$','Location','northeastoutside');
        set(leg,'interpreter','Latex');
        xlabel(str_var,'interpreter','Latex','FontSize',18,'FontWeight','bold');
        ylabel('$\epsilon$','interpreter','Latex','FontSize',18,'FontWeight','bold');
        if ~isempty(vt(flag))
            xlim([min(vt(flag)),max(vt(flag))]);
        else
            xlim([min(vt),max(vt)]);
        end
        str_title = {strcat('$',taskcode,',',typen,'$'),...
            strcat('$',str_para0,'$'),...
            strcat('$',str_algo,'$')};
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
        str_title = {strcat('$',taskcode,',',typen,'$'),...
            strcat('$',str_para0,'$'),...
            strcat('$',str_algo,'$')};
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
        str_title = {strcat('$',taskcode,',',typen,'$'),...
            strcat('$',str_para0,'$'),...
            strcat('$',str_algo,'$')};
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
        str_title = {strcat('$',taskcode,',',typen,'$'),...
            strcat('$',str_para0,'$'),...
            strcat('$',str_algo,'$')};
        title(str_title,'interpreter','Latex');
        
        %
        grid(s1,'on');grid(s2,'on');grid(s3,'on');
        grid(s4,'on');grid(s5,'on');grid(s6,'on');
        
        %
        set(gcf,'Position',get(0,'Screensize'));%Maximize figure
        str_fig2 = strcat(path_fig,str_fig,'_',str_para0,'_sol_flag');
        saveas(gcf,str_fig2,'tif');
        savefig(gcf,str_fig2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLAGS
        figure;
        s1 = subplot(2,2,1);
        plot(vt(1:numel(vexitflag)),vexitflag,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$exit flag$','interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = {strcat('$',taskcode,',',typen,'$'),...
            strcat('$',str_para0,'$')};
        title(str_title,'interpreter','Latex');
        %
        s2 = subplot(2,2,2);
        plot(vt(1:numel(alpha_flag)),alpha_flag,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$flag:\alpha\subseteq(0,1)$','interpreter',...
            'Latex','FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = strcat('$',taskcode,',',typen,',',str_para0,'$');
        title(str_title,'interpreter','Latex');
        %
        s3 = subplot(2,2,3);
        plot(vt(1:numel(theta_flag)),theta_flag,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$flag:\theta\subseteq(0,\infty)$','interpreter',...
            'Latex','FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = strcat('$',taskcode,',',typen,',',str_para0,'$');
        title(str_title,'interpreter','Latex');
        %
        s4 = subplot(2,2,4);
        plot(vt(1:numel(flag)),flag,'LineWidth',2);
        xlabel(str_var,'interpreter','Latex',...
            'FontSize',18,'FontWeight','bold');
        ylabel('$flag:\exists solution $','interpreter',...
            'Latex','FontSize',18,'FontWeight','bold');
        xlim([min(vt),max(vt)]);
        str_title = strcat('$',taskcode,',',typen,',',str_para0,'$');
        title(str_title,'interpreter','Latex');
        %
        grid(s1,'on');grid(s2,'on');grid(s3,'on');grid(s4,'on');
        %
        set(gcf,'Position',get(0,'Screensize'));%Maximize figure
        str_fig2 = strcat(path_fig,str_fig,'_',str_para0,'_exitflag');
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
%rmpath(add_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return