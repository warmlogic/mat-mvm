function [dirs,files] = mm_ft_setSaveDirs_multiSes(exper,ana,cfg,dirs,files,prefix,allSesEvInDir)
%MM_FT_SETSAVEDIRS_MULTISES Sets the saving directories and creates them if they
%don't exist. Also sets up figure options.
%
% [dirs,files] = mm_ft_setSaveDirs_multiSes(exper,ana,cfg,dirs,files,prefix)
%
% prefix           = gets put at the beginning of the directory within the
%                    dir named for the event values. Use, e.g., 'tla' for
%                    time-lock analysis, 'pow' for a power analysis, 'tfr'
%                    for time-freq representation, 'conn' for a
%                    connectivity analysis.
%
% allSesEvInDir    = true/false; whether to include the names of all events
%                    from all sessions in the directory name (default:
%                    false). false means directory will start with 'data'
%
%
% DEFAULTS (if you don't set these):
%  files.saveFigs   = 1 (default), or 0 if you don't want to print figures
%  files.figPrintFormat = 'epsc2' (default), or 'tiff', 'png', 'pdf', 'jpeg90'
%                         (do not include '-d') % see PRINT
%  files.figPrintRes    = 150 (default); some journals require higher DPI
%                         resolution % see PRINT
%  files.figFontName    = 'Helvetica' (default); some journals require a
%                         specific font (e.g., 'Arial', 'Courier', 'Times',
%                         'FixedWidth') % see LISTFONTS
%
% See also: PRINT, LISTFONTS

if ~exist('allSesEvInDir','var') || isempty(allSesEvInDir)
  allSesEvInDir = false;
end

if ~exist('prefix','var')
  prefix = [];
end

if ~isfield(dirs,'dataDir')
  error('Must specify dirs.dataDir.');
end

if ~exist(fullfile(dirs.dataroot,dirs.dataDir),'dir')
  error('%s does not exist; dirs.dataDir must be a real directory, located at fullfile(dirs.dataroot,dirs.dataDir).',fullfile(dirs.dataroot,dirs.dataDir));
end

if ~isfield(dirs,'saveDirStem')
  dirs.saveDirStem = dirs.dataDir;
end

% if ~isfield(exper,'equateTrials')
%   exper.equateTrials = 0;
% end

if ~isfield(ana,'artifact') || (isfield(ana,'artifact') && ~isfield(ana.artifact,'type'))
  ana.artifact.type = {'none'};
elseif isfield(ana,'artifact') && isfield(ana.artifact,'type') && ischar(ana.artifact.type)
  ana.artifact.type = {ana.artifact.type};
end

%% make the field names legal

% some of the many characters that cannot be used in a struct field, to be
% replaced with an underscore
illegalStructFieldChars = {' ','-','+','=','!','@','#','$','%','^','&','*','?','~','(',')','{','}','[',']'};
replaceIllegalCharWith = '_';

% appent this character in front of a field that starts with a number
appendInFrontOfNum = 'a';

for ses = 1:length(exper.sessions)
  for i = 1:length(exper.eventValues{ses})
    % look for illegal characters
    for ic = 1:length(illegalStructFieldChars)
      if ~isempty(strfind(exper.eventValues{ses}{i},illegalStructFieldChars{ic}))
        exper.eventValues{ses}{i} = strrep(exper.eventValues{ses}{i},illegalStructFieldChars{ic},replaceIllegalCharWith);
      end
    end
    % cannot start a struct field with a number
    if isstrprop(exper.eventValues{ses}{i}(1),'digit')
      exper.eventValues{ses}{i} = [appendInFrontOfNum,exper.eventValues{ses}{i}];
    end
  end
end

%% name of the folder to save the FT data in

if allSesEvInDir
  % get all the event names
  evFlat = sort(cellflat(exper.eventValues));
  evStr = sprintf(repmat('%s_',1,length(evFlat)),evFlat{:});
  % remove the trailing underscore
  evStr = evStr(1:end-1);
else
  evStr = 'data';
end

% denote the artifact type
if ~isempty(ana.artifact.type)
  if length(ana.artifact.type) > 1
    artStr = ana.artifact.type{1};
    for i = 2:length(ana.artifact.type)
      artStr = cat(2,artStr,'_',ana.artifact.type{i});
    end
  else
    artStr = ana.artifact.type{1};
  end
  if isfield(ana.artifact,'continuousICA')
    if ana.artifact.continuousICA
      artStr = cat(2,'continuousICA_',artStr);
    end
  end
  artStr = cat(2,'_art_',artStr);
else
  artStr = '';
end

% % denote whether the trial counts are being equated
% evStrDir = sprintf('%s_eq%d%s',evStr,exper.equateTrials,artStr);
evStrDir = sprintf('%s%s',evStr,artStr);

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
  
  %infoDir = sprintf('%s_%d_%d',prefix,exper.prepost(1)*1000,exper.prepost(2)*1000);
  infoDir = sprintf('%s',prefix);
  
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
    %infoDir = sprintf('%s_%s%s_%s_%d_%d_%d_%d',prefix,cfg.method,extraInfo,cfg.output,round(cfg.toi(1)*1000),round(cfg.toi(end)*1000),round(cfg.foi(1)),round(cfg.foi(end)));
    infoDir = sprintf('%s_%s%s_%s_%d_%d',prefix,cfg.method,extraInfo,cfg.output,round(cfg.foi(1)),round(cfg.foi(end)));
  elseif isfield(cfg,'foilim')
    %infoDir = sprintf('%s_%s%s_%s_%d_%d_%d_%d',prefix,cfg.method,extraInfo,cfg.output,round(cfg.toi(1)*1000),round(cfg.toi(end)*1000),round(cfg.foilim(1)),round(cfg.foilim(end)));
    infoDir = sprintf('%s_%s%s_%s_%d_%d',prefix,cfg.method,extraInfo,cfg.output,round(cfg.foilim(1)),round(cfg.foilim(end)));
  end
else
  error('%s:unknownFtFxn','Unsure of how to name the saving directory for %s',mfilename,ana.ftFxn);
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

%% set some default options for figures

% do we want to save them?
if ~isfield(files,'saveFigs')
  files.saveFigs = 1;
end

% set default figure file format to color encapsulated PS (vector graphics)
if ~isfield(files,'figPrintFormat')
  files.figPrintFormat = 'epsc2';
  %files.figPrintFormat = 'png';
else
  if strcmp(files.figPrintFormat(1:2),'-d')
    files.figPrintFormat = files.figPrintFormat(3:end);
  end
end

% default resolution for printed figures (DPI)
if ~isfield(files,'figPrintRes')
  files.figPrintRes = 150; % 150 is the Matlab default
end

% default font for figures
if ~isfield(files,'figFontName')
  files.figFontName = 'Helvetica'; % Helvetica is the Matlab default
end

end
