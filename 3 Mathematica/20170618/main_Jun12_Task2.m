% June 12, 2017
% June 01, 2017
% Task 2, May 22, 2017
% Task3 April 01, 2017
% Task2,1 sigma, Mar18, 2017
% Task2, Mar4, 2017
% Related to May16/Oct17,2016
% edited in Feb21,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Parallel Computing Tollbox
% Using 'Large Scale' in fmincon 
% Specify Gradient of nonlinear constraint in fmincon
%% Main part
clear
close all

taskcode = 'Jun12T2';
% taskcode = 'Jun01T4_2P1';
num_para='V6';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parpool('local',3);
% gcp('nocreate')
% start a parallel pool of 4 workers using the local profile
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
cell_typen = {'I';'II';'III';'IV';'O'};
% {0;1;2;3;4};
% cell for the code of type of function F(x)
save(strcat(path_data,taskcode,'_celltype.mat')...
    ,'cell_typen','-v7.3');

cell_case = {'pstar';'b';'phi';'sigma';'beta';'lambda';'c_f';...
    'r';'c_p';'delta';'original'};
% cell_case = {'pstar';'b'};
% 10 cases
save(strcat(path_data,taskcode,'_cellcase.mat')...
    ,'cell_case','-v7.3');

%%%%%%%%%%%%%%%%%%%
N_h=1;%%%%%%%%%%%%
% the number of points in the right half of x>0
%%%%%%%%%%%%%%%%%%%

width=0.1;

% the width of the half chain that x>0 
step=width/N_h;
% the step size
y=(-width:step:width);

% y is the vartiable of the location of each points on the chain
save(strcat(path_data,taskcode,'_y.mat')...
    ,'y','-v7.3');
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
    cellepsilond_pms=cell(length(cell_case),length(cell_case));
    cellepsilond_ps=cell(length(cell_case),length(cell_case));
    cellepsilond_pps=cell(length(cell_case),length(cell_case));
    
    cellSbarP_pms=cell(length(cell_case),length(cell_case));
    cellSbarP_ps=cell(length(cell_case),length(cell_case));
    cellSbarP_pps=cell(length(cell_case),length(cell_case));
    
    celltheta_pms=cell(length(cell_case),length(cell_case));
    celltheta_ps=cell(length(cell_case),length(cell_case));
    celltheta_pps=cell(length(cell_case),length(cell_case));
    
    cellceq=cell(length(cell_case),length(cell_case));
    cellexitflag = cell(length(cell_case),length(cell_case));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellepsilonc_pms=cell(length(cell_case),length(cell_case));
    cellepsilonc_ps=cell(length(cell_case),length(cell_case));
    cellepsilonc_pps=cell(length(cell_case),length(cell_case));
    
    cellalpha_pms=cell(length(cell_case),length(cell_case));
    
    cellalpha_pps=cell(length(cell_case),length(cell_case));
    
    cellSbarF_pms=cell(length(cell_case),length(cell_case));
    cellSbarF_ps=cell(length(cell_case),length(cell_case));
    cellSbarF_pps=cell(length(cell_case),length(cell_case));
    
    cellSP_pms_epc_pms=cell(length(cell_case),length(cell_case));
    cellSP_pms_epc_ps=cell(length(cell_case),length(cell_case));
    cellSP_pms_epc_pps=cell(length(cell_case),length(cell_case));

    cellSP_ps_epc_pms=cell(length(cell_case),length(cell_case));
    cellSP_ps_epc_ps=cell(length(cell_case),length(cell_case));
    cellSP_ps_epc_pps=cell(length(cell_case),length(cell_case));
    
    cellSP_pps_epc_pms=cell(length(cell_case),length(cell_case));
    cellSP_pps_epc_ps=cell(length(cell_case),length(cell_case));
    cellSP_pps_epc_pps=cell(length(cell_case),length(cell_case));
    
    cellSF_ps_epc_pms=cell(length(cell_case),length(cell_case));
    cellSF_pps_epc_pms=cell(length(cell_case),length(cell_case));
    cellSF_pps_epc_ps=cell(length(cell_case),length(cell_case));
    
    cellvar=cell(length(cell_case),length(cell_case));
    for c = v_cc
    % parfor c=v_cc
    casecode=cell_case{c};
        for c2 = v_cc2
            % case code: '1' '2' '3' '4' '5' '6' '7' '8' '9' 
            casecode2=cell_case{c2};
            if ~strcmp(casecode,casecode2)
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % epsilon_u is bar{epsilon}
                % change epsilon_u with different type of function of F(x)
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
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% 
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
                    str_var2 = '';
                end
                
                %valpha_ps=NaN(numel(vt2),numel(vt));
                valpha_ps=repmat((0.01:0.01:0.99),numel(vt2),1);
                % preallocation
                vepsilond_pms=NaN(numel(vt2),size(valpha_ps,2));
                vepsilond_ps=NaN(numel(vt2),size(valpha_ps,2));
                vepsilond_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                vSbarP_pms=NaN(numel(vt2),size(valpha_ps,2));
                vSbarP_ps=NaN(numel(vt2),size(valpha_ps,2));
                vSbarP_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                vtheta_pms=NaN(numel(vt2),size(valpha_ps,2));
                vtheta_ps=NaN(numel(vt2),size(valpha_ps,2));
                vtheta_pps=NaN(numel(vt2),size(valpha_ps,2));

                
                % preallocation
                vepsilonc_pms=NaN(numel(vt2),size(valpha_ps,2));
                vepsilonc_ps=NaN(numel(vt2),size(valpha_ps,2));
                vepsilonc_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                vSbarF_pms=NaN(numel(vt2),size(valpha_ps,2));
                vSbarF_ps=NaN(numel(vt2),size(valpha_ps,2));
                vSbarF_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                valpha_pms=NaN(numel(vt2),size(valpha_ps,2));
                
                valpha_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                vSP_pms_epc_pms=NaN(numel(vt2),size(valpha_ps,2));
                vSP_pms_epc_ps=NaN(numel(vt2),size(valpha_ps,2));
                vSP_pms_epc_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                vSP_ps_epc_pms=NaN(numel(vt2),size(valpha_ps,2));
                vSP_ps_epc_ps=NaN(numel(vt2),size(valpha_ps,2));
                vSP_ps_epc_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                vSP_pps_epc_pms=NaN(numel(vt2),size(valpha_ps,2));
                vSP_pps_epc_ps=NaN(numel(vt2),size(valpha_ps,2));
                vSP_pps_epc_pps=NaN(numel(vt2),size(valpha_ps,2));
                
                vSF_ps_epc_pms=NaN(numel(vt2),size(valpha_ps,2));
                vSF_pps_epc_pms=NaN(numel(vt2),size(valpha_ps,2));
                vSF_pps_epc_ps=NaN(numel(vt2),size(valpha_ps,2));

                vceq = NaN(Num_Fcn*numel(vt2),size(valpha_ps,2));
                vexitflag = NaN(numel(vt2),size(valpha_ps,2));
                
                vt = NaN(numel(vt2),size(valpha_ps,2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % main computation
                for i = 1:size(valpha_ps,2)

                    %{
                    if strcmp(casecode,'pstar') % 1
                        pstar=vt(i);
                    elseif strcmp(casecode,'b') % 2
                        b=vt(i);
                    elseif strcmp(casecode,'phi') % 3
                        phi=vt(i);
                    elseif strcmp(casecode,'sigma') % 4
                        sigma=vt(i);
                    elseif strcmp(casecode,'beta') % 5
                        beta=vt(i);
                    elseif strcmp(casecode,'lambda') % 6
                        lambda=vt(i);
                    elseif strcmp(casecode,'c_f') % 7
                        c_f=vt(i);
                    elseif strcmp(casecode,'r') % 8
                        r=vt(i);
                    elseif strcmp(casecode,'c_p') % 9
                        c_p=vt(i);
                    elseif strcmp(casecode,'delta') % 10
                        delta=vt(i);
                    elseif strcmp(casecode,'original') % 11
                        mu=vt(i);
                    end
                    %}
                    for j = 1:numel(vt2)
                        
                        if strcmp(casecode2,'pstar') % 1
                            pstar=vt2(j);
                        elseif strcmp(casecode2,'b') % 2
                            b=vt2(j);
                        elseif strcmp(casecode2,'phi') % 3
                            phi=vt2(j);
                        elseif strcmp(casecode2,'sigma') % 4
                            sigma=vt2(j);
                        elseif strcmp(casecode2,'beta') % 5
                            beta=vt2(j);
                        elseif strcmp(casecode2,'lambda') % 6
                            lambda=vt2(j);
                        elseif strcmp(casecode2,'c_f') % 7
                            c_f=vt2(j);
                        elseif strcmp(casecode2,'r') % 8
                            r=vt2(j);
                        elseif strcmp(casecode2,'c_p') % 9
                            c_p=vt2(j);
                        elseif strcmp(casecode2,'delta') % 10
                            delta=vt2(j);
                        elseif strcmp(casecode2,'original') % 11
                            mu=vt2(j);
                        end
                        
                        alpha_ps = valpha_ps(j,i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        tic
                        % Mar5 is the most important part
                        para ={A;B1;B2;c_f;c_p;beta;phi;delta;sigma;lambda;...
                            b;r;epsilon_u;mu;step;...
                            width;pstar;typen;alpha_ps;casecode};

                        [epsilon_d,epsilon_c,theta,alpha,...
                            SbarP,SbarF,SP_pms_epc,SP_ps_epc,SP_pps_epc,...
                            SF_epc,ceq,exitflag,sol] = Jun12T2(para);
                        toc      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        vepsilond_pms(j,i)=epsilon_d(1);
                        vepsilond_ps(j,i)=epsilon_d(2);
                        vepsilond_pps(j,i)=epsilon_d(3);

                        vtheta_pms(j,i)=theta(1);
                        vtheta_ps(j,i)=theta(2);
                        vtheta_pps(j,i)=theta(3);
                        
                        vceq((Num_Fcn*j-(Num_Fcn-1):Num_Fcn*j),i)=ceq;
                        vexitflag(j,i)=exitflag;
                        
                        vSbarP_pms(j,i)=SbarP(1);
                        vSbarP_ps(j,i)=SbarP(2);
                        vSbarP_pps(j,i)=SbarP(3);
                        %%%%%%%%%%%%
                        vepsilonc_pms(j,i)=epsilon_c(1);
                        vepsilonc_ps(j,i)=epsilon_c(2);
                        vepsilonc_pps(j,i)=epsilon_c(3);

                        valpha_pms(j,i)=alpha(1);
                        
                        valpha_pps(j,i)=alpha(2);
                        
                        vSbarF_pms(j,i)=SbarF(1);
                        vSbarF_ps(j,i)=SbarF(2);
                        vSbarF_pps(j,i)=SbarF(3);
                        
                        vSP_pms_epc_pms(j,i)=SP_pms_epc(1);
                        vSP_pms_epc_ps(j,i)=SP_pms_epc(2);
                        vSP_pms_epc_pps(j,i)=SP_pms_epc(3);
                        
                        vSP_ps_epc_pms(j,i)=SP_ps_epc(1);
                        vSP_ps_epc_ps(j,i)=SP_ps_epc(2);
                        vSP_ps_epc_pps(j,i)=SP_ps_epc(3);
                        
                        vSP_pps_epc_pms(j,i)=SP_pps_epc(1);
                        vSP_pps_epc_ps(j,i)=SP_pps_epc(2);
                        vSP_pps_epc_pps(j,i)=SP_pps_epc(3);
                        
                        vSF_ps_epc_pms(j,i)=SF_epc(1);
                        vSF_pps_epc_pms(j,i)=SF_epc(2);
                        vSF_pps_epc_ps(j,i)=SF_epc(3);
                        
                        vt(j,i)=sol;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                end
                cellepsilond_pms{c,c2}=vepsilond_pms;
                cellepsilond_ps{c,c2}=vepsilond_ps;
                cellepsilond_pps{c,c2}=vepsilond_pps;
                
                cellSbarP_pms{c,c2}=vSbarP_pms;
                cellSbarP_ps{c,c2}=vSbarP_ps;
                cellSbarP_pps{c,c2}=vSbarP_pps;
                
                celltheta_pms{c,c2}=vtheta_pms;
                celltheta_ps{c,c2}=vtheta_ps;
                celltheta_pps{c,c2}=vtheta_pps;
                
                cellceq{c,c2}=vceq;
                %
                cellexitflag{c,c2} = vexitflag;
                %%%%%%%%%%%
                cellepsilonc_pms{c,c2}=vepsilonc_pms;
                cellepsilonc_ps{c,c2}=vepsilonc_ps;
                cellepsilonc_pps{c,c2}=vepsilonc_pps;
                
                cellalpha_pms{c,c2}=valpha_pms;
                
                cellalpha_pps{c,c2}=valpha_pps;
                
                cellSbarF_pms{c,c2}=vSbarF_pms;
                cellSbarF_ps{c,c2}=vSbarF_ps;
                cellSbarF_pps{c,c2}=vSbarF_pps;
                
                cellSP_pms_epc_pms{c,c2} = vSP_pms_epc_pms;
                cellSP_pms_epc_ps{c,c2} = vSP_pms_epc_ps;
                cellSP_pms_epc_pps{c,c2} = vSP_pms_epc_pps;
                
                cellSP_ps_epc_pms{c,c2}=vSP_ps_epc_pms;
                cellSP_ps_epc_ps{c,c2}=vSP_ps_epc_ps;
                cellSP_ps_epc_pps{c,c2}=vSP_ps_epc_pps;
                
                cellSP_pps_epc_pms{c,c2}=vSP_pps_epc_pms;
                cellSP_pps_epc_ps{c,c2}=vSP_pps_epc_ps;
                cellSP_pps_epc_pps{c,c2}=vSP_pps_epc_pps;
                
                cellSF_ps_epc_pms{c,c2}=vSF_ps_epc_pms;
                cellSF_pps_epc_pms{c,c2}=vSF_pps_epc_pms;
                cellSF_pps_epc_ps{c,c2}=vSF_pps_epc_ps;
                
                cellvar{c,c2}=vt;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end 
        end   
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save variables
    save(strcat(path_data,str_fig,'_',str_para,'_epd_pms.mat')...
        ,'cellepsilond_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_epd_ps.mat')...
        ,'cellepsilond_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_epd_pps.mat')...
        ,'cellepsilond_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_SbarP_pms.mat')...
        ,'cellSbarP_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SbarP_ps.mat')...
        ,'cellSbarP_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SbarP_pps.mat')...
        ,'cellSbarP_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_theta_pms.mat')...
        ,'celltheta_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_theta_ps.mat')...
        ,'celltheta_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_theta_pps.mat')...
        ,'celltheta_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_ceq.mat')...
            ,'cellceq','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_exitflag.mat')...
        ,'cellexitflag','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_epc_pms.mat')...
        ,'cellepsilonc_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_epc_ps.mat')...
        ,'cellepsilonc_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_epc_pps.mat')...
        ,'cellepsilonc_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_alpha_pms.mat')...
        ,'cellalpha_pms','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_alpha_pps.mat')...
        ,'cellalpha_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_SbarF_pms.mat')...
        ,'cellSbarF_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SbarF_ps.mat')...
        ,'cellSbarF_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SbarF_pps.mat')...
        ,'cellSbarF_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_SP_pms_epc_pms.mat')...
        ,'cellSP_pms_epc_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SP_pms_epc_ps.mat')...
        ,'cellSP_pms_epc_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SP_pms_epc_pps.mat')...
        ,'cellSP_pms_epc_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_SP_ps_epc_pms.mat')...
        ,'cellSP_ps_epc_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SP_ps_epc_ps.mat')...
        ,'cellSP_ps_epc_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SP_ps_epc_pps.mat')...
        ,'cellSP_ps_epc_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_SP_pps_epc_pms.mat')...
        ,'cellSP_pps_epc_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SP_pps_epc_ps.mat')...
        ,'cellSP_pps_epc_ps','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SP_pps_epc_pps.mat')...
        ,'cellSP_pps_epc_pps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_SF_ps_epc_pms.mat')...
        ,'cellSF_ps_epc_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SF_pps_epc_pms.mat')...
        ,'cellSF_pps_epc_pms','-v7.3');
    save(strcat(path_data,str_fig,'_',str_para,'_SF_pps_epc_ps.mat')...
        ,'cellSF_pps_epc_ps','-v7.3');
    
    save(strcat(path_data,str_fig,'_',str_para,'_var.mat')...
        ,'cellvar','-v7.3');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
run('Plot_Jun12T2.m');
%%%%%%
rmpath(add_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%