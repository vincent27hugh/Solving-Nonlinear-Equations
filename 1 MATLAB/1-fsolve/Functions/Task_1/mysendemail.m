% Send efrom example

%{
sendfrom('user@otherdomain.com',...
         'Test subject','Test message',...
         {'folder/attach1.html','attach2.doc'});
%}
%%
function mysendemail(to,subject,message,attachment)
dbstop if error

% define these variable appropriately
from = '2846043323@qq.com';
% your efrom address to send out efrom
% password = 'jhioitqmkwetdhff';
password = 'ebkgzgcrrpgedgfe';
% your efrom password
% If your from is gfrom, you need to set the "less secure apps"
% if your from is qqfrom, you need to activate POPS/SMTP Service and get
% the Authorization code, which you need to verift by phone

setpref('Internet','SMTP_Server','smtp.qq.com');

setpref('Internet','E_mail',from);

setpref('Internet','SMTP_Username',from);
setpref('Internet','SMTP_Password',password);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class',...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

setpref('Internet','E_mail_Charset','SJIS');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(to)
    to = strsplit(to,';');
end
% for multiple recipients
if nargin <=3 
    attachment = {};
end

sendmail(to,subject,message,attachment);

return
