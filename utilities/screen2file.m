function screen2file(filename,files)
%SCREEN2FILE Generate an image file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
%
%   SCREEN2FILE(filename,files) saves the current figure to the file "filename".
%
%   files is a struct containing the fields figPrintFormat and figPrintRes.
%
%    Sean P. McCarthy
%    Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved
%
% modified by MVM

if nargin < 1
  error('Not enough input arguments!')
end

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;

set(gcf,'PaperUnits','inches',...
  'PaperPosition',newpos)

print(gcf,sprintf('-d%s',files.figPrintFormat),sprintf('-r%d',files.figPrintRes),filename);
drawnow

set(gcf,'Units',oldscreenunits,...
  'PaperUnits',oldpaperunits,...
  'PaperPosition',oldpaperpos);

