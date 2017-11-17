% June 12, 2017
% June 01, 2017
% Task 2, May 22, 2017
% Task3 April 01, 2017
% Task2,1 sigma, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taskcode = 'Jun12T2';
% taskcode = 'May22T4_2P1';
num_para='V6';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path0 = pwd;
add_path=strcat(path0,'\Functions\Task_2\');
addpath(add_path);


% mu = 0.08
% step =0.1
path_fig = strcat(path0,'\Figures\Task_2\');
path_data = strcat(path0,'\Data\Task_2\');
str_para0=strcat('paraJun01',num_para);

% fsolve algorithm is 'levenberg-marquardt'
str_algo='IP';
str_para=strcat(str_para0,'-',str_algo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_start =cputime;
% elasped CPU time when starting the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(strcat(path_data,taskcode,'_celltype.mat')...
    ,'cell_typen');
load(strcat(path_data,taskcode,'_cellcase.mat')...
    ,'cell_case');
load(strcat(path_data,taskcode,'_y.mat')...
    ,'y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_Fcn1=2*numel(y);
Num_Fcn2=18;
Num_Fcn=Num_Fcn1+Num_Fcn2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_tt = 3;%1:length(cell_typen) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_cc = 1:length(cell_case)-1;
v_cc2 = length(cell_case);%1:length(cell_case)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = v_tt 
    % F(x) type: 'O' 'I' 'II' 'III' 'IV'
    typen=cell_typen{tt};
    
    str_fig=strcat(taskcode,'_',typen);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocation  
    load(strcat(path_data,str_fig,'_',str_para,'_epd_pms.mat')...
        ,'cellepsilond_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_epd_ps.mat')...
        ,'cellepsilond_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_epd_pps.mat')...
        ,'cellepsilond_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_SbarP_pms.mat')...
        ,'cellSbarP_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_SbarP_ps.mat')...
        ,'cellSbarP_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_SbarP_pps.mat')...
        ,'cellSbarP_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_theta_pms.mat')...
        ,'celltheta_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_theta_ps.mat')...
        ,'celltheta_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_theta_pps.mat')...
        ,'celltheta_pps');

    load(strcat(path_data,str_fig,'_',str_para,'_ceq.mat')...
        ,'cellceq');
    load(strcat(path_data,str_fig,'_',str_para,'_exitflag.mat')...
        ,'cellexitflag');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(strcat(path_data,str_fig,'_',str_para,'_epc_pms.mat')...
        ,'cellepsilonc_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_epc_ps.mat')...
        ,'cellepsilonc_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_epc_pps.mat')...
        ,'cellepsilonc_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_SbarF_pms.mat')...
        ,'cellSbarF_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_SbarF_ps.mat')...
        ,'cellSbarF_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_SbarF_pps.mat')...
        ,'cellSbarF_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_alpha_pms.mat')...
        ,'cellalpha_pms');
    
    load(strcat(path_data,str_fig,'_',str_para,'_alpha_pps.mat')...
        ,'cellalpha_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_SP_pms_epc_pms.mat')...
        ,'cellSP_pms_epc_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_SP_pms_epc_ps.mat')...
        ,'cellSP_pms_epc_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_SP_pms_epc_pps.mat')...
        ,'cellSP_pms_epc_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_SP_ps_epc_pms.mat')...
        ,'cellSP_ps_epc_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_SP_ps_epc_ps.mat')...
        ,'cellSP_ps_epc_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_SP_ps_epc_pps.mat')...
        ,'cellSP_ps_epc_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_SP_pps_epc_pms.mat')...
        ,'cellSP_pps_epc_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_SP_pps_epc_ps.mat')...
        ,'cellSP_pps_epc_ps');
    load(strcat(path_data,str_fig,'_',str_para,'_SP_pps_epc_pps.mat')...
        ,'cellSP_pps_epc_pps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_SF_ps_epc_pms.mat')...
        ,'cellSF_ps_epc_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_SF_pps_epc_pms.mat')...
        ,'cellSF_pps_epc_pms');
    load(strcat(path_data,str_fig,'_',str_para,'_SF_pps_epc_ps.mat')...
        ,'cellSF_pps_epc_ps');
    
    load(strcat(path_data,str_fig,'_',str_para,'_var.mat')...
        ,'cellvar');
    
    for c = v_cc
        casecode=cell_case{c};
        for c2 = v_cc2
            close all
            % case code: '1' '2' '3' '4' '5' '6' '7' '8' '9' 
            casecode2=cell_case{c2};
            if ~strcmp(casecode,casecode2)
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % epsilon_u is bar{epsilon}
                % change epsilon_u with different type of function of F(x)
                
                vepsilond_pms = cellepsilond_pms{c,c2};
                vepsilond_ps = cellepsilond_ps{c,c2};
                vepsilond_pps = cellepsilond_pps{c,c2};
                
                vSbarP_pms = cellSbarP_pms{c,c2};
                vSbarP_ps = cellSbarP_ps{c,c2};
                vSbarP_pps = cellSbarP_pps{c,c2};
                
                vtheta_pms = celltheta_pms{c,c2};
                vtheta_ps = celltheta_ps{c,c2};
                vtheta_pps = celltheta_pps{c,c2};

                vceq = cellceq{c,c2};
                vexitflag = cellexitflag{c,c2};
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                vepsilonc_pms = cellepsilonc_pms{c,c2};
                vepsilonc_ps = cellepsilonc_ps{c,c2};
                vepsilonc_pps = cellepsilonc_pps{c,c2};
                
                vSbarF_pms = cellSbarF_pms{c,c2};
                vSbarF_ps = cellSbarF_ps{c,c2};
                vSbarF_pps = cellSbarF_pps{c,c2};
                
                valpha_pms = cellalpha_pms{c,c2};
                
                valpha_pps = cellalpha_pps{c,c2};
                
                vSP_pms_epc_pms=cellSP_pms_epc_pms{c,c2};
                vSP_pms_epc_ps=cellSP_pms_epc_ps{c,c2};
                vSP_pms_epc_pps=cellSP_pms_epc_pps{c,c2};

                vSP_ps_epc_pms=cellSP_ps_epc_pms{c,c2};
                vSP_ps_epc_ps=cellSP_ps_epc_ps{c,c2};
                vSP_ps_epc_pps=cellSP_ps_epc_pps{c,c2};
                
                vSP_pps_epc_pms=cellSP_pps_epc_pms{c,c2};
                vSP_pps_epc_ps=cellSP_pps_epc_ps{c,c2};
                vSP_pps_epc_pps=cellSP_pps_epc_pps{c,c2};
                
                vSF_ps_epc_pms=cellSF_ps_epc_pms{c,c2};
                vSF_pps_epc_pms=cellSF_pps_epc_pms{c,c2};
                vSF_pps_epc_ps=cellSF_pps_epc_ps{c,c2};

                vt = cellvar{c,c2};
                %%%%%%%%%%%%%%%%%%%%
                if strcmp(typen,'I')
                    epsilon_u=sqrt(3);
                elseif strcmp(typen,'II')
                    epsilon_u=sqrt(2)/(sqrt(pi-2));
                elseif strcmp(typen,'III')
                    epsilon_u=1;
                elseif strcmp(typen,'IV')
                    epsilon_u=(3-sqrt(3))/2;
                elseif strcmp(typen,'O')
                    epsilon_u=(pi/4)/(sin(pi/4));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % use the parameters of Mar4
                str_func=strcat('myParameter_Jun01',num_para);
                func_handle=str2func(str_func);
                parameters=func_handle(typen,epsilon_u);
                
                % parameters=[A;B;r;c_p;beta;phi;delta;sigma;lambda;
                % pstar;b;c_f;mu;B1;B2];
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
                
                Orig_mu=mu;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % case
                % change parameters with different cases
                % vt is column vector
                if strcmp(casecode,'pstar') % 1
                    %vt=(0.5:0.01:2);
                    str_s = 'pstar';
                    str_var = '$p^*$';
                elseif strcmp(casecode,'b') % 2
                    %vt=(0.1:0.01:0.5);
                    %vt=(0.1:0.01:0.4);
                    str_s = 'b';
                    str_var = '$b$';                        
                elseif strcmp(casecode,'phi') % 3
                    %vt=(0.3:0.01:0.8);
                    str_s = 'phi';
                    str_var = '$\phi$';
                elseif strcmp(casecode,'sigma') % 4
                    %vt=(0.01:0.01:0.5);
                    str_s = 'sigma';
                    str_var = '$\sigma$';
                elseif strcmp(casecode,'beta') % 5
                    %vt=(0.4:0.01:0.95);
                    str_s = 'beta';
                    str_var = '$\beta$';
                elseif strcmp(casecode,'lambda') % 6
                    %vt=(0.01:0.01:0.5);
                    str_s = 'lambda';
                    str_var = '$\lambda$';
                elseif strcmp(casecode,'c_f') % 7
                    % vt=(0.01:0.01:0.5);
                    %vt=(0.12:0.01:0.5);
                    str_s = 'c_F';
                    str_var = '$c^F$';
                elseif strcmp(casecode,'r') % 8
                    %vt=(0.005:0.001:0.03);
                    str_s = 'r';
                    str_var = '$r$';
                elseif strcmp(casecode,'c_p') % 9
                    %vt=(0.01:0.001:0.2);
                    str_s = 'c_P';
                    str_var = '$c^P$';
                elseif strcmp(casecode,'delta') % 10
                    %vt=(0.1:0.01:0.7);
                    str_s = 'delta';
                    str_var = '$\delta$';
                elseif strcmp(casecode,'original') % original case
                    %vt=Orig_mu;
                    str_s = 'OrigPara';
                    str_var = '';
                end

                if strcmp(casecode2,'pstar') % 1
                    vt2=[0.6;0.8;1;1.2];
                    str_s2 = 'pstar';
                    str_var2 = '$p^*$';
                elseif strcmp(casecode2,'b') % 2
                    vt2=[0.1;0.2;0.3;0.4];
                    str_s2 = 'b';
                    str_var2 = '$b$';                        
                elseif strcmp(casecode2,'phi') % 3
                    vt2=[0.02;0.05;0.1;0.15];
                    str_s2 = 'phi';
                    str_var2 = '$\phi$';
                elseif strcmp(casecode2,'sigma') % 4
                    vt2=[0.001;0.0048;0.01;0.02];
                    str_s2 = 'sigma';
                    str_var2 = '$\sigma$';
                elseif strcmp(casecode2,'beta') % 5
                    vt2=[0.32;0.52;0.72;0.82];
                    str_s2 = 'beta';
                    str_var2 = '$\beta$';
                elseif strcmp(casecode2,'lambda') % 6
                    vt2=[0.05;0.1;0.2;0.3];
                    str_s2 = 'lambda';
                    str_var2 = '$\lambda$';
                elseif strcmp(casecode2,'c_f') % 7
                    % vt=[0.09;0.15;0.213;0.3];
                    vt2=[0.15;0.213;0.3];
                    str_s2 = 'c_F';
                    str_var2 = '$c^F$';
                elseif strcmp(casecode2,'r') % 8
                    vt2=[0.006;0.008;0.01;0.012];
                    str_s2 = 'r';
                    str_var2 = '$r$';
                elseif strcmp(casecode2,'c_p') % 9
                    vt2=[0.05;0.1;0.15;0.2];
                    str_s2 = 'c_P';
                    str_var2 = '$c^P$';
                elseif strcmp(casecode2,'delta') % 10
                    vt2=[0.3;0.4;0.5;0.6];
                    str_s2 = 'delta';
                    str_var2 = '$\delta$';
                elseif strcmp(casecode2,'original') % original case
                    vt2=Orig_mu;
                    str_s2 = 'OrigPara';
                    str_var2 = '$OrigPara$';
                end

                str_figt=strcat(path_fig,taskcode,'_',typen,'_',str_s,...
                    '_',str_s2);
                
                valpha_ps=repmat((0.01:0.01:0.99),numel(vt2),1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% plot
                alpha_pms_flag = fun_alphaflag(valpha_pms);
                
                alpha_pps_flag = fun_alphaflag(valpha_ps);
                
                theta_pms_flag = vtheta_pms>0;
                theta_ps_flag = vtheta_ps>0;
                theta_pps_flag = vtheta_pps>0;
                
                flag1 = vexitflag==1;
                
                flag = alpha_pms_flag&alpha_pps_flag&...
                    theta_pms_flag&theta_ps_flag&theta_pps_flag&...
                    flag1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                figure;
                %%%%%%%%%%%%%%%%%%%%%%
                s1=subplot(3,3,1);
                hold on 
                str_leg=cell(size(vt2));
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    vttheta_pps = vtheta_pps(kk,:);
                    plot(vtalpha_ps(flag1(kk,:)),vttheta_pps(flag1(kk,:)),'LineWidth',2);
                   
                    str_leg{kk}=strcat('$',str_s2,'=',...
                            num2str(vt2(kk)),'$'); 
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                
                if v_cc2~=length(cell_case)
                    leg=legend(str_leg,'Location','northeastoutside');
                    set(leg,'interpreter','Latex');
                end
                
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\theta_{p^{+},\sigma}$',...
                    'interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s2=subplot(3,3,2);
                hold on 
                str_leg=cell(size(vt2));
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    vtalpha_pps = valpha_pps(kk,:);
                    plot(vtalpha_ps(flag1(kk,:)),vtalpha_pps(flag1(kk,:)),'LineWidth',2);

                    str_leg{kk}=strcat('$',str_s2,'=',...
                        num2str(vt2(kk)),'$');
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                
                if v_cc2~=length(cell_case)
                    leg=legend(str_leg,'Location','northeastoutside');
                    set(leg,'interpreter','Latex');
                end
                
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\alpha_{p^{+},\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s3=subplot(3,3,3);
                hold on 
                str_leg=cell(numel(vt2)*2,1);
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    
                    vtepsilond_pps = vepsilond_pps(kk,:);
                    vtepsilonc_pps = vepsilonc_pps(kk,:);
                    
                    %vF_epd_pps=fun_F_x (vtepsilond_pps,typen,epsilon_u);
                    %vF_epc_pps=fun_F_x (vtepsilonc_pps,typen,epsilon_u);
                    
                    plot(vtalpha_ps(flag1(kk,:)),vtepsilond_pps(flag1(kk,:)),'LineWidth',2);
                    
                    if v_cc2==length(cell_case)
                        str_leg{2*kk-1}=strcat('$\epsilon^d$');
                    else
                        str_leg{2*kk-1}=strcat('$\epsilon^d,',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                    
                    
                    plot(vtalpha_ps(flag1(kk,:)),vtepsilonc_pps(flag1(kk,:)),'LineWidth',2);
                    
                    if v_cc2==length(cell_case)
                        str_leg{2*kk}=strcat('$\epsilon^c$');
                    else
                        str_leg{2*kk}=strcat('$\epsilon^c,',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                
                leg=legend(str_leg,'Location','northeastoutside');
                set(leg,'interpreter','Latex');
                
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\epsilon_{p^{+},\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s4=subplot(3,3,4);
                hold on 
                str_leg=cell(size(vt2));
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    vttheta_ps = vtheta_ps(kk,:);
                    plot(vtalpha_ps(flag1(kk,:)),vttheta_ps(flag1(kk,:)),'LineWidth',2);
                    
                    str_leg{kk}=strcat('$',str_s2,...
                        '=',num2str(vt2(kk)),'$');
                    
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                
                if v_cc2~=length(cell_case)
                    leg=legend(str_leg,'Location','northeastoutside');
                    set(leg,'interpreter','Latex');
                end
                
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\theta_{p,\sigma}$',...
                    'interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s5=subplot(3,3,5);
                hold on 
                str_leg=cell(size(vt2));
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    plot(vtalpha_ps(flag1(kk,:)),vt(flag1(kk,:)),'LineWidth',2);
                    
                    str_leg{kk}=strcat('$',str_s2,...
                        '=',num2str(vt2(kk)),'$');
                    
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                
                if v_cc2~=length(cell_case)
                    leg=legend(str_leg,'Location','northeastoutside');
                    set(leg,'interpreter','Latex');
                end
                
                ylabel(str_var,'interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s6=subplot(3,3,6);
                hold on 
                str_leg=cell(numel(vt2)*2,1);
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    
                    vtepsilond_ps = vepsilond_ps(kk,:);
                    vtepsilonc_ps = vepsilonc_ps(kk,:);
                    
                    %vF_epd_ps=fun_F_x (vtepsilond_ps,typen,epsilon_u);
                    %vF_epc_ps=fun_F_x (vtepsilonc_ps,typen,epsilon_u);
                    
                    plot(vtalpha_ps(flag1(kk,:)),vtepsilond_ps(flag1(kk,:)),'LineWidth',2);
                    
                    if v_cc2==length(cell_case)
                        str_leg{2*kk-1}=strcat('$\epsilon^d$');
                    else
                        str_leg{2*kk-1}=strcat('$\epsilon^d,',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                    
                    plot(vtalpha_ps(flag1(kk,:)),vtepsilonc_ps(flag1(kk,:)),'LineWidth',2);
                    
                    if v_cc2==length(cell_case)
                        str_leg{2*kk}=strcat('$\epsilon^c$');
                    else
                        str_leg{2*kk}=strcat('$\epsilon^c,',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                leg=legend(str_leg,'Location','northeastoutside');
                set(leg,'interpreter','Latex');
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\epsilon_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s7=subplot(3,3,7);
                hold on 
                str_leg=cell(size(vt2));
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    vttheta_pms = vtheta_pms(kk,:);
                    plot(vtalpha_ps(flag1(kk,:)),vttheta_pms(flag1(kk,:)),'LineWidth',2);
                    if v_cc2==length(cell_case)
                        str_leg{kk}=strcat(str_var2);
                    else
                        str_leg{kk}=strcat('$',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                
                if v_cc2~=length(cell_case)
                    leg=legend(str_leg,'Location','northeastoutside');
                    set(leg,'interpreter','Latex');
                end
                
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\theta_{p^{-},\sigma}$',...
                    'interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s8=subplot(3,3,8);
                hold on 
                str_leg=cell(size(vt2));
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    vtalpha_pms = valpha_pms(kk,:);
                    plot(vtalpha_ps(flag1(kk,:)),vtalpha_pms(flag1(kk,:)),'LineWidth',2);
                    if v_cc2==length(cell_case)
                        str_leg{kk}=strcat(str_var2);
                    else
                        str_leg{kk}=strcat('$',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                
                if v_cc2~=length(cell_case)
                    leg=legend(str_leg,'Location','northeastoutside');
                    set(leg,'interpreter','Latex');
                end
                
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\alpha_{p^{-},\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %%%%%%%%%%%%%%%%%%%%%%
                s9=subplot(3,3,9);
                hold on 
                str_leg=cell(numel(vt2)*2,1);
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    vtepsilond_pms = vepsilond_pms(kk,:);
                    vtepsilonc_pms = vepsilonc_pms(kk,:);
                    
                    %vF_epd_pms=fun_F_x (vtepsilond_pms,typen,epsilon_u);
                    %vF_epc_pms=fun_F_x (vtepsilonc_pms,typen,epsilon_u);
                    
                    plot(vtalpha_ps(flag1(kk,:)),vtepsilond_pms(flag1(kk,:)),'LineWidth',2);
                    
                    if v_cc2==length(cell_case)
                        str_leg{2*kk-1}=strcat('$\epsilon^d$');
                    else
                        str_leg{2*kk-1}=strcat('$\epsilon^d,',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                    
                    plot(vtalpha_ps(flag1(kk,:)),vtepsilonc_pms(flag1(kk,:)),'LineWidth',2);
                    
                    if v_cc2==length(cell_case)
                        str_leg{2*kk}=strcat('$\epsilon^c$');
                    else
                        str_leg{2*kk}=strcat('$\epsilon^c,',str_s2,...
                            '=',num2str(vt2(kk)),'$');
                    end
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                leg=legend(str_leg,'Location','northeastoutside');
                set(leg,'interpreter','Latex');
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$\epsilon_{p^{-},\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                %
                grid(s1,'on');grid(s2,'on');grid(s3,'on');
                grid(s4,'on');grid(s5,'on');grid(s6,'on');
                grid(s7,'on');grid(s8,'on');grid(s9,'on');
                %
                set(gcf,'Position',get(0,'Screensize'));
                str_fig2 = strcat(str_figt,'_',str_para,'_solution');
                saveas(gcf,str_fig2,'tif');
                savefig(gcf,str_fig2);            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% plot of EXITFLAG
                figure;
                %
                hold on 
                str_leg=cell(numel(vt2),1);
                for kk = 1:numel(vt2)
                    vtalpha_ps = valpha_ps(kk,:);
                    plot(vtalpha_ps,flag(kk,:),'LineWidth',2);
                    str_leg{kk}=strcat('$Flag,',...
                        str_s2,'=',num2str(vt2(kk)),'$'); 
                end
                hold off
                xlim([min(vtalpha_ps),max(vtalpha_ps)]);
                leg=legend(str_leg,'Location','northeastoutside');
                set(leg,'interpreter','Latex');
                xlabel('$\alpha_{p,\sigma}$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                ylabel('$Flag$','interpreter',...
                    'Latex','FontSize',18,'FontWeight','bold');
                
                str_title = {strcat('$',taskcode,',',typen,'$'),...
                    strcat('$',str_para,'$')};
                title(str_title,'interpreter','Latex');
                
                grid on 
                
                descr={'Exit flag:';...
                    '1:';...
                    'Local minimum found';...
                    'that satisfies the';...
                    'constraints.';...
                    '0:';...
                    'No feasible point';...
                    'was found'};
                dim = [0,0.95,0,0];
                annotation('textbox',dim,'string',descr,...
                    'FontSize',8,'FontWeight','bold',...
                    'FitBoxToText','on');
                
                set(gcf,'Position',get(0,'Screensize'));
                str_fig2 = strcat(str_figt,'_',str_para,'_exitflag');
                saveas(gcf,str_fig2,'tif');
                savefig(gcf,str_fig2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                N_h = 1;
                nFn =(1:Num_Fcn);
                for kk = 1:numel(vt2)
                    figure;
                    
                    surf(vt,nFn,abs(vceq(Num_Fcn*kk-(Num_Fcn-1):...
                        Num_Fcn*kk,:)));
                    xlabel(str_var,'interpreter',...
                        'Latex','FontSize',18,'FontWeight','bold');
                    ylabel('$No. of Function$','interpreter',...
                        'Latex','FontSize',18,'FontWeight','bold');
                    zlabel('$Error$','interpreter',...
                        'Latex','FontSize',18,'FontWeight','bold');
                    title(strcat('Error of constraint functions-',...
                        taskcode,',',typen,...
                        '-',str_s2,'_',num2str(vt2(kk))));
                    %
                    set(gcf,'Position',get(0,'Screensize'));
                    str_fig2 = strcat(str_figt,'_',str_para,'_ceq',...
                        str_s2,'_eq_',num2str(kk));
                    saveas(gcf,str_fig2,'tif');
                    savefig(gcf,str_fig2);
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end

close all
time_passby=cputime-time_start;
DCPU = sprintf('\n !!!ATTENTION!!!\n Total CPU Time is :%12.3f \n !!!ATTENTION!!!\n', time_passby);
disp(DCPU);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Send email to notify the fufillment
to = 'huweihugh@gmail.com';
% the efrom address you want to send message to
% you can send to multiple email address, use ; to seperate
subject = strcat('Program_',taskcode,'_done') ;

txt = sprintf('%s is done and total cputime is: %12.3f seconds!'...
    ,taskcode,time_passby);

message = {txt,...
    '',...
    'This is automatically sending mail,',...
    'please DO NOT reply!'};

attachment_path = 'D:\GN Luo\20170522\Figures\Task_1\';
filename = 'Apr11_I_b_paraApr11-IP_solution.tif';
attachment = strcat(attachment_path,filename);
% {'folder/attach1.html','attach2.doc'}

mysendemail(to,subject,message);

pause(1)
beep
% delete(gcp);
mydialog;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath(add_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%