function [dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg,dirs,files,prefix)
%MM_FT_SETSAVEDIRS Sets the saving directories and creates them if they
%don't exist. Also sets up figure options.
%
% [dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg,dirs,files,prefix)
%
% prefix           = gets put at the beginning of the directory within the
%                    dir named for the event values. Use, e.g., 'tla' for
%                    time-lock analysis, 'pow' for a power analysis, 'tfr'
%                    for time-freq representation, 'conn' for a
%                    connectivity analysis.
%
% files.saveFigs   = 1 (default) or 0
% files.figFileExt = 'png' (default), 'eps', 'pdf', 'jpg'
%

if nargin < 6
  prefix = [];
end

if ~isfield(dirs,'dataDir')
  error('Must specift dirs.dataDir.');
end

if ~exist(fullfile(dirs.dataroot,dirs.dataDir),'dir')
  error('%s does not exist; dirs.dataDir must be a real directory (fullfile(dirs.dataroot,dirs.dataDir))',fullfile(dirs.dataroot,dirs.dataDir));
end

if ~isfield(dirs,'saveDirStem')
  dirs.saveDirStem = dirs.dataDir;
end

%% name of the folder to save the FT data in; get all the event names

if ~isfield(exper,'eventValuesExtra')
  evStr = sprintf(repmat('%s_',1,length(exper.eventValues)),exper.eventValues{:});
else
  if ~isfield(exper.eventValuesExtra,'newValue')
    evStr = sprintf(repmat('%s_',1,length(exper.eventValues)),exper.eventValues{:});
  elseif isfield(exper.eventValuesExtra,'newValue') && (isfield(exper.eventValuesExtra,'onlyKeepExtras') && exper.eventValuesExtra.onlyKeepExtras == 1)
    evStr = sort(cat(2,exper.eventValuesExtra.newValue{:}));
    evStr = sprintf(repmat('%s_',1,length(evStr)),evStr{:});
  elseif isfield(exper.eventValuesExtra,'newValue') && (~isfield(exper.eventValuesExtra,'onlyKeepExtras') || (isfield(exper.eventValuesExtra,'onlyKeepExtras') && exper.eventValuesExtra.onlyKeepExtras == 0))
    evStr = sort(cat(2,exper.eventValues,cat(2,exper.eventValuesExtra.newValue{:})));
    evStr = sprintf(repmat('%s_',1,length(evStr)),evStr{:});
  end
end
% remove the underscore
evStr = evStr(1:end-1);

% denote whether the trial counts are being equated
evStrDir = sprintf('%s_eq%d',evStr,exper.equateTrials);

%% set the directory name, dependent upon the analysis type

if strcmp(ana.ftFxn,'ft_timelockanalysis')
  if isempty(prefix)
    prefix = 'tla';
  end
%   if ~strcmp(ana.ftype,'tla')
%     warning([mfilename,':incompatiblePrefixName'],'The ana.ftype ''%s'' might not be appropriate for the ftFxn %s',ana.ftype,ana.ftFxn);
%   end
    
  if isfield(ana,'otherFxn') && ~isempty(ana.otherFxn{1})
    add_prefix = [];
    for i = 1:length(ana.cfg_other)
      add_prefix = cat(2,sprintf('%s_%s',add_prefix,ana.cfg_other{i}.ftype));
    end
    prefix = [prefix,add_prefix];
  end
  
  infoDir = sprintf('%s_%d_%d',prefix,exper.prepost(1)*1000,exper.prepost(2)*1000);
  
elseif strcmp(ana.ftFxn,'ft_freqanalysis')
  if isempty(prefix)
    prefix = 'tfr';
  end
%   if ~strcmp(ana.ftype,'tfr') && ~strcmp(ana.ftype,'coh') && ~strcmp(ana.ftype,'plv') && ~strcmp(ana.ftype,'fourier')
%     warning([mfilename,':incompatiblePrefixName'],'The ana.ftype ''%s'' might not be appropriate for the ftFxn %s',ana.ftype,ana.ftFxn);
%   end
  
  if isfield(ana,'otherFxn') && ~isempty(ana.otherFxn{1})
    add_prefix = [];
    for i = 1:length(ana.cfg_other)
      add_prefix = cat(2,sprintf('%s_%s',add_prefix,ana.cfg_other{i}.ftype));
    end
    prefix = [prefix,add_prefix];
  end
  
  % extra string in filename for identification
  if strcmp(cfg.method,'wavelet')
    extraInfo = sprintf('_w%d',cfg.width);
  elseif strcmp(cfg.method,'mtmconvol') || strcmp(cfg.method,'mtmfft')
    extraInfo = sprintf('_%s',cfg.taper);
  else
    extraInfo = '';
  end
  
  if isfield(cfg,'foi')
    infoDir = sprintf('%s_%s%s_%s_%d_%d_%d_%d',prefix,cfg.method,extraInfo,cfg.output,cfg.toi(1)*1000,cfg.toi(end)*1000,round(cfg.foi(1)),round(cfg.foi(end)));
  elseif isfield(cfg,'foilim')
    infoDir = sprintf('%s_%s%s_%s_%d_%d_%d_%d',prefix,cfg.method,extraInfo,cfg.output,cfg.toi(1)*1000,cfg.toi(end)*1000,round(cfg.foilim(1)),round(cfg.foilim(end)));
  end
else
  error('mm_ft_setSaveDirs:unknownFtFxn','Unsure of how to name the saving directory for %s',ana.ftFxn);
end

% set the save directory
dirs.saveDirRaw = fullfile(dirs.dataroot,dirs.saveDirStem,'ft_data',evStrDir,'ft_raw');
dirs.saveDirProc = fullfile(dirs.dataroot,dirs.saveDirStem,'ft_data',evStrDir,infoDir);

% denote that these data will be averaged
if isfield(cfg,'keeptrials') && strcmp(cfg.keeptrials,'no')
  dirs.saveDirProc = [dirs.saveDirProc,'_avg'];
end

% directory in which to save the raw data
if ~exist(dirs.saveDirRaw,'dir')
  mkdir(dirs.saveDirRaw)
end

% directory in which to save the processed data
if ~exist(dirs.saveDirProc,'dir')
  mkdir(dirs.saveDirProc)
end

% directory in which to save figures
dirs.saveDirFigs = fullfile(dirs.saveDirProc,'figs');
if ~exist(dirs.saveDirFigs,'dir')
  mkdir(dirs.saveDirFigs)
end

%% set some options for figures

% do we want to save them?
if ~isfield(files,'saveFigs')
  files.saveFigs = 1;
end

% set the default type to png
if ~isfield(files,'figFileExt')
  files.figFileExt = 'png';
end

if strcmp(files.figFileExt,'eps')
  files.figPrintFormat = '-depsc2';
elseif strcmp(files.figFileExt,'pdf')
  files.figPrintFormat = '-dpdf';
elseif strcmp(files.figFileExt,'png')
  files.figPrintFormat = '-dpng';
elseif strcmp(files.figFileExt,'jpg')
  files.figPrintFormat = '-djpg90';
end

end
