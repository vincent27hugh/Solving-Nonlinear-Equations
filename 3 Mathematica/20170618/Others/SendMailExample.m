to = {'huweihugh@gmail.com'};
% the efrom address you want to send message to
% you can send to multiple email address, use ; to seperate
subject = 'Test3';

txt = 'Successful!';

message = {txt,...
    '',...
    'This is automatically sending mail,',...
    'please DO NOT reply!'};

attachment_path = 'D:\GN Luo\20170522\Figures\Task_1\';
filename = 'Apr11_I_b_paraApr11-IP_solution.tif';
attachment = strcat(attachment_path,filename);
% {'folder/attach1.html','attach2.doc'}

mysendemail(to,subject,message);