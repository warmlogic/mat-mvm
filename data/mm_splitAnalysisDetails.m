function mm_splitAnalysisDetails(adFile,badEvPerValue,rm_eq)

%function mm_splitAnalysisDetails(adFile,badEvPerValue,rm_eq)
%
% adFile        = path to and name of original analysisDetails.mat
%
% badEvPerValue = number of possible events per event value, listed in
%                 alphabetical order because these event values got sorted
%                 when processing using the old method

if ~exist('rm_eq','var')
  rm_eq = false;
end

if ~exist('badEvPerVal','var')
  badEvPerValue = [];
end

ad = load(adFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
if isfield(ad,'exper')
  exper = ad.exper;
else
  error('exper variable not found');
end
if isfield(ad,'ana')
  ana = ad.ana; %#ok<NASGU>
else
  error('ana variable not found');
end
if isfield(ad,'dirs')
  dirs = ad.dirs;
  if rm_eq
    if ~isempty(strfind(dirs.saveDirRaw,'_eq0'))
      dirs.saveDirRaw = strrep(dirs.saveDirRaw,'_eq0','');
    end
    if ~isempty(strfind(dirs.saveDirRaw,'_eq1'))
      dirs.saveDirRaw = strrep(dirs.saveDirRaw,'_eq1','');
    end
    
    if ~isempty(strfind(dirs.saveDirProc,'_eq0'))
      dirs.saveDirProc = strrep(dirs.saveDirProc,'_eq0','');
    end
    if ~isempty(strfind(dirs.saveDirProc,'_eq1'))
      dirs.saveDirProc = strrep(dirs.saveDirProc,'_eq1','');
    end
    
    if ~isempty(strfind(dirs.saveDirFigs,'_eq0'))
      dirs.saveDirFigs = strrep(dirs.saveDirFigs,'_eq0','');
    end
    if ~isempty(strfind(dirs.saveDirFigs,'_eq1'))
      dirs.saveDirFigs = strrep(dirs.saveDirFigs,'_eq1','');
    end
  end
else
  error('dirs variable not found');
end
if isfield(ad,'files')
  files = ad.files; %#ok<NASGU>
else
  error('files variable not found');
end
if isfield(ad,'cfg_proc')
  cfg_proc = ad.cfg_proc; %#ok<NASGU>
else
  error('cfg_proc variable not found');
end
if isfield(ad,'cfg_pp')
  cfg_pp = ad.cfg_pp; %#ok<NASGU>
else
  error('cfg_pp variable not found');
end

exper_orig = exper;

if ~isempty(badEvPerValue)
  evCount = 1;
end

for sub = 1:length(exper_orig.subjects)
  for ses = 1:length(exper_orig.sessions)
    fprintf('%s %s...\n',exper_orig.subjects{sub},exper_orig.sessions{ses});
    
    %     % turn the session name into a string for easier printing
    %     if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
    %       sesStr = exper.sessions{ses}{1};
    %       for i = 2:length(exper.sessions{ses})
    %         sesStr = cat(2,sesStr,'_',exper.sessions{ses}{i});
    %       end
    %     elseif ~iscell(exper.sessions{ses})
    %       sesStr = exper.sessions{ses};
    %     elseif iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1
    %       sesStr = cell2mat(exper.sessions{ses});
    %     end
    
    exper.subjects = exper_orig.subjects(sub);
    exper.sessions = {exper_orig.sessions(ses)};
    exper.sesStr = exper_orig.sesStr(ses);
    
    exper.eventValues = {exper_orig.eventValues};
    
    prepost = exper_orig.prepost(1,:);
    if size(exper_orig.prepost,1) > 1
      for i = 2:length(exper_orig.eventValues)
        prepost = cat(1,prepost,exper_orig.prepost(i,:));
      end
    end
    exper.prepost = {prepost};
    
    exper = rmfield(exper,'badChan');
    exper = rmfield(exper,'nTrials');
    exper = rmfield(exper,'badEv');
    
    exper.badChan.(exper_orig.sesStr{ses}) = exper_orig.badChan(sub,ses);
    for evVal = 1:length(exper_orig.eventValues)
      exper.nTrials.(exper_orig.sesStr{ses}).(exper_orig.eventValues{evVal}) = exper_orig.nTrials.(exper_orig.eventValues{evVal})(sub,ses);
    end
    
    if ~isempty(exper_orig.badEv{sub,ses})
      manArtFile = fullfile(dirs.saveDirRaw,exper_orig.subjects{sub},exper_orig.sesStr{ses},sprintf('%s_%s_manArtFT.mat',exper_orig.subjects{sub},exper_orig.sesStr{ses}));
      if exist(manArtFile,'file') && isempty(badEvPerValue)
        evNumColNum = 4;
        maf = load(manArtFile);
        for evVal = 1:length(exper_orig.eventValues)
          exper.badEv.(exper_orig.sesStr{ses}).(exper_orig.eventValues{evVal}) = {exper_orig.badEv{sub,ses}(maf.cfg_manArt.trl(:,evNumColNum) == evVal)};
        end
      elseif ~exist(manArtFile,'file') && ~isempty(badEvPerValue)
        % predefined number of events
        if sum(badEvPerValue) == length(exper_orig.badEv{sub,ses})
          for evVal = 1:length(exper_orig.eventValues)
            exper.badEv.(exper_orig.sesStr{ses}).(exper_orig.eventValues{evVal}) = {exper_orig.badEv{sub,ses}(evCount:badEvPerValue(evVal))};
            evCount = evCound + badEvPerValue(evVal);
          end
        else
          warning('Sum across number of events per value does not equal the number of bad events');
          keyboard
          %exper.badEv.(exper_orig.sesStr{ses}).(exper_orig.eventValues{evVal}) = {[]};
        end
      else
        warning('Unable to set badEv for individual event values, simply including entire badEv vector.');
        exper.badEv.(exper_orig.sesStr{ses}).(exper_orig.eventValues{evVal}) = {[]};
        exper.badEv.(exper_orig.sesStr{ses}).badEv = exper_orig.badEv(sub,ses);
      end
    else
      exper.badEv.(exper_orig.sesStr{ses}).(exper_orig.eventValues{evVal}) = {[]};
    end
    
    raw_dir = fullfile(dirs.saveDirRaw,exper_orig.subjects{sub},exper_orig.sesStr{ses});
    proc_dir = fullfile(dirs.saveDirProc,exper_orig.subjects{sub},exper_orig.sesStr{ses});
    
    if ~exist(raw_dir,'dir')
      mkdir(raw_dir);
    end
    if ~exist(proc_dir,'dir')
      mkdir(proc_dir);
    end
    
    raw_subDetailsFile = fullfile(raw_dir,'subjectDetails.mat');
    proc_subDetailsFile = fullfile(proc_dir,'subjectDetails.mat');
    
    save(raw_subDetailsFile,'exper','ana','dirs','files','cfg_pp');
    save(proc_subDetailsFile,'exper','ana','dirs','files','cfg_pp','cfg_proc');
    
  end
end

end