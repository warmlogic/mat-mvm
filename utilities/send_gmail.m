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
%               Default: the hardcoded from email address
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
% Adapted from:
% http://www.amirwatad.com/blog/archives/2009/01/31/sending-emails-with-matlab/
% http://www.mathworks.com/support/solutions/en/data/1-3PRRDV/
%
% See also: SENDMAIL

%% You must define these:
address = 'matt.mollison@gmail.com'; % your Gmail address
password = ''; % your Gmail password

%% Send the email

% make sure the password is set
if isempty(password)
  error('Modify %s to inlcude your password and save it encrpyted using ''pcode send_gmail''.',mfilename);
end

% set some defaults
if nargin < 4
  attachments = [];
  if nargin < 3
    recipients = {address};
  end
end    

% make sure recipients is a cell
if ischar(recipients)
  recipients = {recipients};
end

% tell us who we're sending to
fprintf('Sending mail to:%s...',sprintf(repmat(' ''%s''',1,length(recipients)),recipients{:}));

% don't touch this unless you need to change the email supplier
setpref('Internet','E_mail',address);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',address);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% send the email
sendmail(recipients,subject,mail_message,attachments);

fprintf('Done.\n');

end
