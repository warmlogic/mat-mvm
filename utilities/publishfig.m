function publishfig(handle,cleartitle,ticFontSize,labelFontSize,figFontName,legendFontSize)
%PUBLISHFIG - Prepared figure for publication.
% 
% Prepares a figure for publication or presentation by making all
% the fonts larger.  The default is to work on the current figure
% and to clear the title away.
%
% FUNCTION:
%   publishfig(handle,cleartitle,ticFontSize,labelFontSize,figFontName,legendFontSize)
%
% INPUT ARGS:
%   handle        = gcf (default) % handle to figure to modify
%   cleartitle    = 1 (default) % whether to clear the title or not
%   ticFontSize   = 20 (default) % size of tic mark font
%   labelFontSize = 24 (default) % size of label font
%   legendFontSize= 16 (default) % size of legend font
%   figFontName   = 'Helvetica' (default) % some journals require a
%                   specific font (e.g., 'Arial', 'Courier', 'Times',
%                   'FixedWidth') % see LISTFONTS
%
% See also: LISTFONTS

if ~exist('handle','var') || isempty(handle)
  handle = gcf;
end
if ~exist('cleartitle','var') || isempty(cleartitle)
  cleartitle = 1;
end
if ~exist('ticFontSize','var') || isempty(ticFontSize)
  ticFontSize = 20;
end
if ~exist('labelFontSize','var') || isempty(labelFontSize)
  labelFontSize = 24;
end
if ~exist('figFontName','var') || isempty(figFontName)
  figFontName = 'Helvetica';
end
if ~exist('legendFontSize','var') || isempty(legendFontSize)
  legendFontSize = 16;
end

% remove the title
if cleartitle
  title('');
end

% set the tics on all children
h = get(handle,'Children');
for j=1:length(h)
  if isprop(h(j),'FontSize') && isprop(h(j),'FontWeight') && isprop(h(j),'FontName')
    set(h(j),'FontSize',ticFontSize,'FontWeight','Bold','FontName',figFontName);
  end
end

% set the title and labels
h = get(handle,'Children');
for j = 1:length(h)
  if isprop(h(j),'Title');
    set(get(h(j),'Title'),'FontSize',labelFontSize,'FontName',figFontName);
  end
  
  if isprop(h(j),'xlabel');
    set(get(h(j),'xlabel'),'FontSize',labelFontSize,'FontWeight','Bold','FontName',figFontName);
  end
  
  if isprop(h(j),'ylabel');
    set(get(h(j),'ylabel'),'FontSize',labelFontSize,'FontWeight','Bold','FontName',figFontName);
  end
end


% h = get(get(handle,'Children'),'Title');
% if length(h) > 1 
%   for j=1:length(h)
%     if isprop(h{j},'FontSize')
%       set(h{j},'FontSize',24);
%     end
%   end
% else
  
%   if isprop(h,'FontSize')
%     set(h,'FontSize',24);
%   end
% end


% % set the xlabel
% h = get(get(handle,'Children'),'xlabel');
% if length(h) > 1 
%   for j=1:length(h)
%     set(h{j},'FontSize',24,'FontWeight','Bold');
%   end
% else
%   set(h,'FontSize',24,'FontWeight','Bold');
% end

% % set the ylabel
% h = get(get(handle,'Children'),'ylabel');
% if length(h) > 1 
%   for j=1:length(h)
%     set(h{j},'FontSize',24,'FontWeight','Bold');
%   end
% else
%   set(h,'FontSize',24,'FontWeight','Bold');
% end

% redo the legend
h_legend = legend;
set(h_legend,'FontSize',legendFontSize);
% legend

return

% fix the bounding box
% hndl2=axes('position',[0,0,1,1],...
%        'units','normalized',...
%        'UserData','timestamp',...
%        'visible','off');


% save the figure
%saveas(handle,figs(i).name,'fig')

% export to deps
%[fpath,fname,fext] = fileparts(figs(i).name);
%print('-depsc2',[fpath fname '.eps'])

