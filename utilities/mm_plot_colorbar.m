function mm_plot_colorbar(minmax,orient,label,filename,cmap,figPrintFormat,figPrintRes,figFontName)
%MM_PLOT_COLORBAR - Plot a colorbar and print it to file
%
% mm_plot_colorbar(minmax,orient,label,filename,cmap,figPrintFormat,figPrintRes)
%
% Input:
%  minmax         = [-1.5 1.5] (required)
%  orient         = 'vert' (default); or 'horiz' % see COLORBAR
%  label          = 'Voltage (\muV)' (default)
%  filename       = 'colorbar' (pwd; default); or e.g.,'~/Desktop/colorbar'
%                   % no extension necessary
%  cmap           = the colormap to use; hot(64) (default); or, e.g., 'jet',
%                   or jet(64), 'hot', hot(128), etc.
%                   % see COLORMAP
%  figPrintFormat = 'epsc2' (default) % don't include '-d'; see PRINT
%  figPrintRes    = figure resolution; 150 (default) % see PRINT
%  figFontName    = 'Helvetica' (default); some journals require a specific
%                   font (e.g., 'Arial', 'Courier', 'Times', 'FixedWidth')
%                   % see LISTFONTS
%
% This function will automatically add the minmax to the filename, to one
% decimal place. However, it will replace the decimal with a 'p' so there
% are no extra periods in the file name because some programs can't handle
% that (e.g., LaTeX)
%
% convert this to an eps with eps2eps (installed with MacTeX) to remove the
% whitespace. see inline code for a bash script to do this, plus conversion
% to pdf.
%
% See also: COLORBAR, COLORMAP, PRINT, LISTFONTS

% Bash script for conversion: epstopdf_batch.sh
%
% #!/bin/bash
% 
% if [ $# -le 0 ]; then
%     echo
%     echo "Usage: $(basename $0) file1.eps [file2.eps ...]"
%     echo
%     echo "  This script allows batch conversion of eps to pdf."
%     echo
%     exit 1
% fi
% 
% # set the filelist
% files=$*
% 
% # loop and add each file
% for i in ${files[*]} ; do
%   echo "$i -->> ${i%.eps}.pdf"
%   eval "eps2eps $i temp_$i ; epstopdf temp_$i"
%   eval "mv temp_${i%.eps}.pdf ${i%.eps}.pdf"
%   eval "rm temp_$i"
% done

if ~exist('minmax','var') || isempty(minmax)
  error('Must set minmax variable');
end

if ~exist('orient','var') || isempty(orient)
  orient = 'vert';
end
if ~strcmp(orient,'vert') && ~strcmp(orient,'horiz')
  error('You set ''orient'' to ''%s''. Must set it to ''vert'' or ''horiz''',orient);
end
if ~exist('label','var') || isempty(label)
  label = 'Voltage (\muV)';
end
if ~exist('filename','var') || isempty(filename)
  filename = 'colorbar';
end
if ~exist('cmap','var') || isempty(cmap)
  cmap = hot(64);
end
if ~exist('figPrintFormat','var') || isempty(figPrintFormat)
  figPrintFormat = 'epsc2';
end
if ~exist('figPrintRes','var') || isempty(figPrintRes)
  figPrintRes = 150;
end
if ~exist('figFontName','var') || isempty(figFontName)
  figFontName = 'Helvetica';
end

% set the colormap
colormap(cmap);
% plot the dummy data
imagesc(0,0,0,minmax);

% turn on the colorbar
h = colorbar(orient);
if strcmp(orient,'vert')
  set(get(h,'YLabel'),'string',label);
elseif strcmp(orient,'horiz')
  set(get(h,'XLabel'),'string',label);
end

% make the fonts bigger
publishfig(gcf,1,[],[],figFontName);

% turn off the main plot
cla
axis off

% add the limits to the file name
minmax_str{1} = strrep(sprintf('%.1f',minmax(1)),'.','p');
minmax_str{2} = strrep(sprintf('%.1f',minmax(2)),'.','p');
[pathstr,name,ext] = fileparts(filename);
filename = fullfile(pathstr,sprintf('%s_%s_%s_%s',name,orient,minmax_str{1},minmax_str{2}));
% put the extension back if there was one
if ~isempty(ext)
  filename = sprintf('%s%s',filename,ext);
end

print(gcf,sprintf('-d%s',figPrintFormat),sprintf('-r%d',figPrintRes),filename);

end

