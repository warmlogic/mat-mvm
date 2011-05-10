function send_gmail(subject,mail_message,recipients,attachments)
% SEND_GMAIL sends a mail message through Gmail
%
% function send_gmail(subject,mail_message,recipients,attachments)
% input arguments:
%
% subject:      A string - subject of the mail
%
% mail_message: A cell array of strings - Message body. One cell per line.
%
% recipents:    A cell array of strings - Email addresses of recipients
%               e.g., {'mail1@mail.com', 'mail2@mail.com'}.
%               Default: {'matt.mollison@gmail.com'}
%
% attachments:  A cell array of strings - path to file to attach, one cell
%                                         per attachment. Use [] for none.
%                                         Default: []
%
% NB: You must modify send_gmail.m to inlcude your email address and
% password and then do an encrpyted save using pcode: 'pcode send_gmail'.
% The pcoded file takes precedence over the .m file. Be sure to remove your
% password from the .m file after pcoding.
%
% Adapted from: http://www.amirwatad.com/blog/archives/2009/01/31/sending-emails-with-matlab/
% and: http://www.mathworks.com/support/solutions/en/data/1-3PRRDV/
%
% See also: SENDMAIL

if nargin < 4
  attachments = [];
  if nargin < 3
    recipients = {'matt.mollison@gmail.com'};
  end
end    

% Define these:
address = 'matt.mollison@gmail.com'; %Your GMail email address
password = ''; %Your GMail password

if isempty(password)
  error('Modify %s to inlcude your password and save it encrpyted to your MATLAB directory using ''pcode send_gmail''.',mfilename);
  %password = input('\n\nType pass wd and press ''enter'': ','s');
end

% Don't touch unless you need to change the Email supplier (currently Gmail)
setpref('Internet','E_mail',address);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',address);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Send the email
sendmail(recipients,subject,mail_message,attachments);

end
