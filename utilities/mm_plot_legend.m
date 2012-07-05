function mm_plot_legend(label,linestyle,linewidth,filename,figPrintFormat,figPrintRes,figFontName)
%MM_PLOT_LEGEND - Plot a legend and print it to file
%
% mm_plot_legend(label,linestyle,linewidth,filename,figPrintFormat,figPrintRes)
%
% Input:
%  label          = {'a','b','c'} (required) % see LEGEND
%  linestyle      = {'r-','b--','k:','g-','c--','m:','y-'} (default) % see PLOT
%  linewidth      = 1.0 (default) % see PLOT
%  filename       = 'legend' (default); or, e.g., '~/Desktop/legend';
%                   % no extension necessary
%  figPrintFormat = 'epsc2' (default) % don't include '-d'; see PRINT
%  figPrintRes    = figure resolution; 150 (default) % see PRINT
%  figFontName    = 'Helvetica' (default); some journals require a specific
%                   font (e.g., 'Arial', 'Courier', 'Times', 'FixedWidth')
%                   % see LISTFONTS
%
% This function will automatically add the labels to the filename.
%
% convert this to an eps with eps2eps (installed with MacTeX) to remove the
% whitespace. see inline code for a bash script to do this, plus conversion
% to pdf.
%
% See also: LEGEND, PLOT, PRINT, LISTFONTS

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

if ~exist('label','var') || isempty(label)
  error('Must set label variable');
end
if ~exist('linestyle','var') || isempty(linestyle)
  linestyle = {'r-','b--','k:','g-','c--','m:','y-'};
  %error('Must set label linestyle');
end

if ~exist('linewidth','var') || isempty(linewidth)
  linewidth = 1.0;
end
if ~exist('filename','var') || isempty(filename)
  filename = 'legend';
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

if length(label) > length(linestyle)
  error('The number of labels and linestyles does not match.');
end

% initialize to store the labels for the file name
filename_labels = cell(size(label));

figure
for i = 1:length(label)
  plot([0,0],linestyle{i},'LineWidth',linewidth);
  hold on
  
  % get the label without punctuation or whitespace
  filename_labels{i} = cell2mat(regexp(label{i}, '\w+','match'));
end
hold off

% turn on the legend
legend(label);

% make the fonts bigger
publishfig(gcf,1,[],[],figFontName);

% turn off the main plot
cla
axis off

% add the labels to the file name
[pathstr,name,ext] = fileparts(filename);
% make the filename labels into a string
filename_labels = sprintf(repmat('_%s',1,length(filename_labels)),filename_labels{:});
filename = fullfile(pathstr,sprintf('%s%s',name,filename_labels));
% put the extension back if there was one
if ~isempty(ext)
  filename = sprintf('%s%s',filename,ext);
end

print(gcf,sprintf('-d%s',figPrintFormat),sprintf('-r%d',figPrintRes),filename);

end
