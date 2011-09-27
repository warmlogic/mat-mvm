function mm_plot_colorbar(minmax,label,filename,orient,cmap,figPrintFormat)
%MM_PLOT_COLORBAR - Plot a colorbar and print it to file
%
% Input:
%  minmax = [-1.5 1.5];
%  label = 'Voltage (\muV)';
%  orient = 'vert'; % or 'horiz'
%  cmap = 'jet'; % or 'hot'
%  filename = '~/Desktop/colormap'; % no extension necessary
%  figPrintFormat = 'epsc2'; % see PRINT
%
% This function will automatically add the minmax to the filename, to one
% decimal place. However, it will replace the decimal with a 'p' so there
% are no extra periods in the file name because some programs can't handle
% that (e.g., LaTeX)
%
% convert this to an eps with eps2eps (installed with MacTeX) to remove the
% whitespace. see inline code for a bash script to do this, plus conversion
% to pdf.

% Bash script for conversion: epstopdf_batch.sh
%
% #!/bin/bash
% 
% #files=(`ls *.png`)
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
% # set the base command
% CMD="epstopdf "
% 
% # loop and add each file
% for i in ${files[*]} ; do
%   echo "$i -->> ${i%.eps}.pdf"
%   eval "eps2eps $i temp_$i ; epstopdf temp_$i"
%   eval "mv temp_${i%.eps}.pdf ${i%.eps}.pdf"
%   eval "rm temp_$i"
%   #eval $CMD $i
% done

if nargin < 6
  figPrintFormat = 'epsc2';
  if nargin < 5
    cmap = 'jet';
    if nargin < 4
      orient = 'vert';
    end
  end
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
publishfig(gca,1);

% turn off the main plot
cla
axis off

% add the limits to the file name
minmax_str{1} = strrep(sprintf('%.1f',minmax(1)),'.','p');
minmax_str{2} = strrep(sprintf('%.1f',minmax(2)),'.','p');
[pathstr,name,ext] = fileparts(filename);
filename = fullfile(pathstr,sprintf('%s_%s_%s',name,minmax_str{1},minmax_str{2}));
% put the extension back if there was one
if ~isempty(ext)
  filename = sprintf('%s%s',filename,ext);
end

saveas(gcf,filename,figPrintFormat);

end

