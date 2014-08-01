function [data,badChan_str] = mm_ft_artifact_repairChan(data,badChan_str,elecfile,contData,ylim,blocksize)

if nargin < 6
  if strcmp(contData,'no')
    blocksize = [];
  elseif strcmp(contData,'yes')
    blocksize = 30;
  end
  if nargin < 5
    ylim = [];
  end
end

cfgChannelRepair = [];
cfgChannelRepair.elecfile = elecfile;
cfgChannelRepair.continuous = contData;
cfgChannelRepair.blocksize = blocksize;

% viewmode?
repair_viewmode = [];
while isempty(repair_viewmode) || (repair_viewmode ~= 0 && repair_viewmode ~= 1)
  repair_viewmode = input('\nDo you want to plot in butterfly (1) or vertical (0) mode? (1 or 0, then press ''return''):\n\n');
end
if repair_viewmode
  cfgChannelRepair.viewmode = 'butterfly';
  if ~isempty(ylim)
    cfgChannelRepair.ylim = ylim;
  else
    cfgChannelRepair.ylim = [-150 150];
  end
else
  cfgChannelRepair.viewmode = 'vertical';
  if ~isempty(ylim)
    cfgChannelRepair.ylim = ylim;
  else
    cfgChannelRepair.ylim = [-10 10];
  end
end

% subset?
repair_chanNum = [];
while isempty(repair_chanNum) || (repair_chanNum ~= 0 && repair_chanNum ~= 1)
  repair_chanNum = input('\nDo you want to plot all channels (1) or a particular subset of channels (0)? (1 or 0, then press ''return''):\n\n');
end
if repair_chanNum
  cfgChannelRepair.channel = 'all';
else
  channel = ft_channelselection('gui', data.label);
  cfgChannelRepair.channel = channel;
  fprintf('(Manually close any empty figure windows.)\n');
end

if strcmp(cfgChannelRepair.viewmode,'butterfly')
  fprintf('\nUse the ''i'' key and mouse to identify channels in the data browser. Note any consistently bad channels.\n');
end
fprintf('\nUse the ''q'' key to quit the data browser when finished. Then channel selection will begin.\n');
cfgChannelRepair = ft_databrowser(cfgChannelRepair, data);
% bug when calling rejectartifact right after databrowser, pause first
pause(1);

rejArt_repair_really = [];
while isempty(rejArt_repair_really) || (rejArt_repair_really ~= 0 && rejArt_repair_really ~= 1)
  rejArt_repair_really = input('\nWere there channels to repair? (1 or 0, then press ''return''):\n\n');
end

if rejArt_repair_really
  badchannel = ft_channelselection('gui', data.label);
  fprintf('(Manually close any empty figure windows.)\n');
  badChan_str = cat(1,badChan_str,badchannel);
  cfgChannelRepair.channel = 'all';
  cfgChannelRepair.badchannel = badchannel;
  cfgChannelRepair.method = 'spline';
  fprintf('Repairing channels%s using method=''%s''...\n',sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}),cfgChannelRepair.method);
  data = ft_channelrepair(cfgChannelRepair, data);
end
