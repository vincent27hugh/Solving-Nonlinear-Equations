%{
var=[epsilond,epsilonc,theta,alpha,...
    pstar,b,phi,sigma,beta,...
    lambda,cf,r,cp,delta];
% [1,14]
%}
for j=1:14
    if j==1
        str_txt='dF_depsilond';
    elseif j==2
        str_txt='dF_depsilonc';
    elseif j==3
        str_txt='dF_dtheta';
    elseif j==4
        str_txt='dF_dalpha';
    elseif j==5
        str_txt='dF_dpstar';
    elseif j==6
        str_txt='dF_db';
    elseif j==7
        str_txt='dF_dphi';
    elseif j==8
        str_txt='dF_dsigma';
    elseif j==9
        str_txt='dF_dbeta';
    elseif j==10
        str_txt='dF_dlambda';
    elseif j==11
        str_txt='dF_dcf';
    elseif j==12
        str_txt='dF_dr';
    elseif j==13
        str_txt='dF_dcp';
    elseif j==14
        str_txt='dF_ddelta';
    end
    fid=fopen(strcat(str_txt,'.txt'),'wt');
    for i = 1:4
        fprintf(fid,'dF(%d)/dx: %s ; \n ',i,char(vpa(Jacfvar(i,j),5)));
    end
    fclose(fid);
end