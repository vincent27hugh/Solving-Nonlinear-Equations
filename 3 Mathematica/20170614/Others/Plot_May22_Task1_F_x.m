clear
close all

global epsilon_u

path0=pwd;
addpath(strcat(path0,'\Functions\'));
str_var='y';
cell_typen = {'I';'II';'III';'IV';'O'};% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = 1:length(cell_typen)
    typen=cell_typen{tt};
    switch typen
        case 'I'
            epsilon_u=sqrt(3);
            y1 = (-20:0.1:epsilon_u);
            F1 = fun_F_x(y1,typen);
            intF1=NaN(size(y1));
            for i = 1:length(y1)
                intF1(i) = fun_int_F(y1(i),epsilon_u,typen);
            end
        case 'II'
            epsilon_u=sqrt(2)/(sqrt(pi-2));
            y2 = (-20:0.1:epsilon_u);
            F2 = fun_F_x(y2,typen);
            intF2=NaN(size(y2));
            for i = 1:length(y2)
                intF2(i) = fun_int_F(y2(i),epsilon_u,typen);
            end
        case 'III'
            epsilon_u=1;
            y3 = (-20:0.1:epsilon_u-0.1);
            F3 = fun_F_x(y3,typen);
            intF3=NaN(size(y3));
            for i = 1:length(y3)
                intF3(i) = fun_int_F(y3(i),epsilon_u,typen);
            end
        case 'IV'
            epsilon_u=(3-sqrt(3))/2;
            y4 = (-20:0.1:epsilon_u);
            F4 = fun_F_x(y4,typen);
            intF4=NaN(size(y4));
            for i = 1:length(y4)
                intF4(i) = fun_int_F(y4(i),epsilon_u,typen);
            end
        case 'O'
            epsilon_u=(pi/4)/(sin(pi/4));
            y0 = (-20:0.1:epsilon_u);
            F0 = fun_F_x(y0,typen);
            intF0=NaN(size(y0));
            for i = 1:length(y0)
                intF0(i) = fun_int_F(y0(i),epsilon_u,typen);
            end
    end
end

path1 = strcat(path0,'\Figures\');
%%
figure;
%
subplot(2,3,1);
plot(y0,F0,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y0),max(y0)]);
ylabel('F(x)- Type O','FontSize',12,'FontWeight','bold');
%
subplot(2,3,2);
plot(y1,F1,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y1),max(y1)]);
ylabel('F(x)- Type I','FontSize',12,'FontWeight','bold');
%
subplot(2,3,3);
plot(y2,F2,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y2),max(y2)]);
ylabel('F(x)- Type II','FontSize',12,'FontWeight','bold');
%
subplot(2,3,4);
plot(y3,F3,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y3),max(y3)]);
ylabel('F(x)- Type III','FontSize',12,'FontWeight','bold');
%
subplot(2,3,5);
plot(y4,F4,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y4),max(y4)]);
ylabel('F(x)- Type IV','FontSize',12,'FontWeight','bold');
%
set(gcf,'Position',get(0,'Screensize'));
str_fig2 = strcat(path1,'Mar3','_','Plot_of_F_x');
saveas(gcf,str_fig2,'tif');
savefig(gcf,str_fig2);

%%
figure;
%
subplot(2,3,1);

plot(y0,intF0,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y0),max(y0)]);
ylabel('integ(1-F(x))- Type O','FontSize',12,'FontWeight','bold');
%
subplot(2,3,2);
plot(y1,intF1,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y1),max(y1)]);
ylabel('integ(1-F(x))- Type I','FontSize',12,'FontWeight','bold');
%
subplot(2,3,3);
plot(y2,intF2,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y2),max(y2)]);
ylabel('integ(1-F(x))- Type II','FontSize',12,'FontWeight','bold');
%
subplot(2,3,4);
plot(y3,intF3,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y3),max(y3)]);
ylabel('integ(1-F(x))- Type III','FontSize',12,'FontWeight','bold');
%
subplot(2,3,5);
plot(y4,intF4,'LineWidth',2);
xlabel(str_var,'FontSize',12,'FontWeight','bold');
xlim([min(y4),max(y4)]);
ylabel('integ(1-F(x))- Type IV','FontSize',12,'FontWeight','bold');
%
set(gcf,'Position',get(0,'Screensize'));
str_fig2 = strcat(path1,'Mar3','_','Plot_of_integ(1-F_x)');
saveas(gcf,str_fig2,'tif');
savefig(gcf,str_fig2);

rmpath(strcat(path0,'\Functions\'));