%% run FieldTrip's automatic artifact detection on the data

% log for info
resultsDir = fullfile(dirs.localDir,dirs.dataDir,'ftAuto_results',subject);
if ~exist(resultsDir, 'dir')
  resultsDir = mkdir(fullfile(dirs.localDir,dirs.dataDir,'ftAuto_results'),subject);
end
resultsFileLoc = fullfile(dirs.localDir,dirs.dataDir,'ftAuto_results',subject);
resultsFile = fopen(fullfile(resultsFileLoc,'ftAuto_results.txt'), 'w');

if rejArt_ftAuto
  if ~isfield(ana.artifact,'trlpadding')
    ana.artifact.trlpadding = 0;
  end
  if ~isfield(ana.artifact,'artpadding')
    ana.artifact.artpadding = 0.1;
  end
  if ~isfield(ana.artifact,'fltpadding')
    ana.artifact.fltpadding = 0;
  end
  if ~exist('badChan_str','var')
    badChan_str = {};
  end
  
  if ~isfield(ana.artifact,'basic_art') || isempty(ana.artifact.basic_art)
    ana.artifact.basic_art = true;
  end
  if ~isfield(ana.artifact,'basic_art_z')
    ana.artifact.basic_art_z = 30;
  end
  
  if ~isfield(ana.artifact,'jump_art') || isempty(ana.artifact.jump_art)
    ana.artifact.jump_art = true;
  end
  if ~isfield(ana.artifact,'jump_art_z')
    ana.artifact.jump_art_z = 50;
  end
  
  if ~isfield(ana.artifact,'eog_art') || isempty(ana.artifact.eog_art)
    ana.artifact.eog_art = true;
  end
  if ~isfield(ana.artifact,'eog_art_z')
    ana.artifact.eog_art_z = 3.5;
  end
  
  if ~isfield(ana.artifact,'thresh') || isempty(ana.artifact.thresh)
    ana.artifact.thresh = true;
  end
  if ~isfield(ana.artifact,'threshmin')
    ana.artifact.threshmin = -100;
  end
  if ~isfield(ana.artifact,'threshmax')
    ana.artifact.threshmax = 100;
  end
  if ~isfield(ana.artifact,'threshrange')
    ana.artifact.threshrange = 150;
  end
  
  % get the trial definition for automated FT artifact rejection
  trl = ft_findcfg(data.cfg,'trl');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for threshold artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ana.artifact.thresh
    if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
      cfg = [];
      cfg.continuous = 'no';
      % get the trial definition for automated FT artifact rejection
      cfg.trl = trl;
      
      % % don't exclude eye channels because we want to reject any blinks
      % % that ICA didn't catch
      % cfg.artfctdef.threshold.channel = {'all'};
      % exclStr = '';
      
      % % exclude eye channels - assumes we're using EGI's HCGSN
      % cfg.artfctdef.threshold.channel = {'all', '-E25', '-E8', '-E127', '-E126', '-E128', '-E125'};
      % exclStr = ' (excludes eye channels)';
      
      % % exclude eye channels and all the channels around the periphery - assumes we're using EGI's HCGSN
      % cfg.artfctdef.threshold.channel = {'all', '-E1', '-E8', '-E14', '-E17', '-E21', '-E25', '-E32', '-E38', '-E43', '-E44', '-E48', '-E49', '-E56', '-E57', '-E63', '-E64', '-E68', '-E69', '-E73', '-E74', '-E81', '-E82', '-E88', '-E89', '-E94', '-E95', '-E99', '-E100', '-E107', '-E113', '-E114', '-E119', '-E120', '-E121', '-E125', '-E126', '-E127', '-E128'};
      % exclStr = ' (excludes eye channels and peripheral channels)';
      
      % % exclude eye channels and neighbors - assumes we're using EGI's HCGSN
      % cfg.artfctdef.threshold.channel = {'all', ...
      %   '-E48', '-E128', '-E127', '-E126', '-E125', '-E119', ...
      %   '-E43', '-E32', '-E25', '-E21', '-E17', '-E14', '-E8', '-E1', '-E120', ...
      %   '-E26', '-E22', '-E15', '-E9', '-E2', ...
      %   '-E23', '-E18', '-E16', '-E10', '-E3', ...
      %   '-E19', '-E11', '-E4'};
      % exclStr = ' (excludes eye channels and neighbors)';
      
      % exclude eye channels and neighbors and all peripheral channels - assumes we're using EGI's HCGSN
      cfg.artfctdef.threshold.channel = {'all', '-E1', '-E8', '-E14', '-E17', '-E21', '-E25', '-E32', '-E38', '-E43', '-E44', '-E48', '-E49', '-E56', '-E57', '-E63', '-E64', '-E68', '-E69', '-E73', '-E74', '-E81', '-E82', '-E88', '-E89', '-E94', '-E95', '-E99', '-E100', '-E107', '-E113', '-E114', '-E119', '-E120', '-E121', '-E125', '-E126', '-E127', '-E128', ...
        '-E26', '-E22', '-E15', '-E9', '-E2', ...
        '-E23', '-E18', '-E16', '-E10', '-E3', ...
        '-E19', '-E11', '-E4'};
      exclStr = ' (excludes eye channels and neighbors and peripheral channels)';
      
      cfg.artfctdef.threshold.bpfilter = 'yes';
      cfg.artfctdef.threshold.bpfreq = [0.3 30];
      cfg.artfctdef.threshold.bpfiltord = 4;
      
      cfg.artfctdef.threshold.min = ana.artifact.threshmin;
      cfg.artfctdef.threshold.max = ana.artifact.threshmax;
      cfg.artfctdef.threshold.range = ana.artifact.threshrange;
      
      fprintf('\nUsing EGI HydroCel GSN...\nChecking for voltages above %.1f uV and below %.1f uV or out of peak-to-peak range %.1f uV%s...\n',cfg.artfctdef.threshold.max,cfg.artfctdef.threshold.min,cfg.artfctdef.threshold.range,exclStr);
      
      % auto mark zvalue artifacts
      [cfg, artifact_thresh] = ft_artifact_threshold(cfg, data);
    else
      warning('Not using EGI HydroCel GSN 128/129 electrode file! Threshold artifacts are not being assessed!!');
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for basic zvalue artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ana.artifact.basic_art
    cfg = [];
    cfg.continuous = 'no';
    % get the trial definition for automated FT artifact rejection
    cfg.trl = trl;
    
    cfg.artfctdef.zvalue.channel = 'all';
    cfg.artfctdef.zvalue.cutoff = ana.artifact.basic_art_z;
    cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
    cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
    cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
    %cfg.artfctdef.zvalue.fltpadding = 0;
    
    % interactive artifact viewer
    %cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
    
    fprintf('Checking for (basic) zvalue artifacts at z=%.1f...\n',cfg.artfctdef.zvalue.cutoff);
    
    % auto mark some artifacts
    [cfg, artifact_basic] = ft_artifact_zvalue(cfg, data);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for jump artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ana.artifact.jump_art
    cfg = [];
    cfg.continuous = 'no';
    % get the trial definition for automated FT artifact rejection
    cfg.trl = trl;
    
    % cutoff and padding
    % select a set of channels on which to run the artifact detection
    cfg.artfctdef.zvalue.channel = 'all';
    cfg.artfctdef.zvalue.cutoff = ana.artifact.jump_art_z;
    cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
    cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
    cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
    %cfg.artfctdef.zvalue.fltpadding = 0;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative = 'yes';
    cfg.artfctdef.zvalue.medianfilter = 'yes';
    cfg.artfctdef.zvalue.medianfiltord = 9;
    cfg.artfctdef.zvalue.absdiff = 'yes';
    
    % interactive artifact viewer
    %cfg.artfctdef.zvalue.interactive = ft_autoCheckArt_interactive;
    
    fprintf('\nChecking for jump artifacts at z=%.1f...\n',cfg.artfctdef.zvalue.cutoff);
    
    % auto mark jump artifacts
    [cfg, artifact_jump] = ft_artifact_zvalue(cfg, data);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for EOG artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ana.artifact.eog_art
    if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
      cfg = [];
      cfg.trl = trl;
      cfg.continuous = 'no';
      
      % cutoff and padding
      % select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
      %cfg.artfctdef.zvalue.channel = 'all';
      cfg.artfctdef.zvalue.channel = {'E127','E126','E128','E125','E8','E25'};
      cfg.artfctdef.zvalue.cutoff      = ana.artifact.eog_art_z;
      cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
      if strcmp(cfg.continuous,'yes')
        cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      elseif strcmp(cfg.continuous,'no')
        cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
      end
      cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
      
      % algorithmic parameters
      cfg.artfctdef.zvalue.bpfilter   = 'yes';
      cfg.artfctdef.zvalue.bpfilttype = 'but';
      cfg.artfctdef.zvalue.bpfreq     = [1 15];
      cfg.artfctdef.zvalue.bpfiltord  = 4;
      cfg.artfctdef.zvalue.hilbert    = 'yes';
      
      % interactive artifact viewer
      %cfg.artfctdef.zvalue.interactive = 'yes';
      
      fprintf('\nChecking for EOG artifacts at z=%.1f...\n',cfg.artfctdef.zvalue.cutoff);
      [cfg, artifact_EOG] = ft_artifact_zvalue(cfg,data);
    else
      warning('Not using EGI HydroCel GSN 128/129 electrode file! EOG artifacts are not being assessed!!');
    end
  end
  
%% figures & log for results

% log information

% total segments marked bad b/c of blinks
totalBlinks = sum(logical(sum(artifact_EOG,2)));
fprintf(resultsFile,'\nBlink artifacts (automatically rejected): %d\n', totalBlinks);

% total segments marked bad b/c of eye movement
totalEyeMove = sum(logical(sum(foundEyeMove,2)));
fprintf(resultsFile,'\nEye movement artifacts (automatically rejected): %d\n', totalEyeMove);

% total segments marked bad b/c of too many bad channels
totalTooManyBadChan = sum(logical(sum(foundTooManyBadChan,2)));
fprintf(resultsFile,'\nToo many bad channels (automatically rejected): %d\n', totalTooManyBadChan);

% total segments marked bad
totalbadevents = sum(logical(sum(foundArt,2)));
fprintf(resultsFile,'\nTotal number of segments with artifacts: %d\n',totalbadevents);

% # of bad channels (repaired & rejected?)
totalbadchans = sum(foundArt(:));
fprintf(resultsFile,'\nTotal number of bad channels across all events: %d\n',totalbadchans);

% channels bad for the entire recording
if ~isempty(fullyRepairChan_str)
    for i = 1:length(fullyRepairChan_str)
        fprintf(resultsFile,'\nChannels repaired across entire recording: %d\n',fullyRepairChan_str{i});
    end
else
    fprintf(resultsFile,'\nChannels repaired across entire recording: 0\n');
end

% prints a figure for each trial in foundArt, with subplot for each
% artifact type labeled by artifact type and trial
for tr = 1:nTrial
    if any(foundArt(tr,:))
        h = figure;
        figName = num2str(tr);
        hold on
        subplot(3,1,[1,1]);
        plot(data.trial{tr},'LineWidth',1)
        ylim([-100 100]);
        title(sprintf('foundArt_%s', num2str(tr)));
        if any(foundBlink(tr,:))
            subplot(3,1,[2,2]);
            hold on
            plot(data.trial{tr}(eog_upper_left,:),'m','LineWidth',1)
            plot(data.trial{tr}(eog_lower_left,:),'k','LineWidth',1)
            plot(data.trial{tr}(eog_horiz_left,:),'b','LineWidth',2)
            plot(data.trial{tr}(eog_horiz_right,:),'r','LineWidth',2)
            hold off
            ylim([-100 100]);
            title(sprintf('foundBlink_%s', num2str(tr)));
        end
        if any(foundEyeMove(tr,:))
            subplot(3,1,[3,3]);
            hold on
            plot(data.trial{tr}(eog_upper_left,:),'m','LineWidth',1)
            plot(data.trial{tr}(eog_lower_left,:),'k','LineWidth',1)
            plot(data.trial{tr}(eog_horiz_left,:),'b','LineWidth',2)
            plot(data.trial{tr}(eog_horiz_right,:),'r','LineWidth',2)
            hold off
            ylim([-100 100]);
            title(sprintf('foundEyeMove_%s', num2str(tr)));
        end
        hold off
        saveas(h,fullfile(resultsFileLoc, figName),'png');
        close(h);
    end
end

%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % reject the automatically defined artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  if ana.artifact.basic_art
    cfg.artfctdef.basic.artifact = artifact_basic;
  end
  if ana.artifact.jump_art
    cfg.artfctdef.jump.artifact = artifact_jump;
  end
  %cfg.artfctdef.muscle.artifact = artifact_muscle;
  if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
    if ana.artifact.eog_art
      cfg.artfctdef.eog.artifact = artifact_EOG;
    end
    if ana.artifact.thresh
      cfg.artfctdef.threshold.artifact = artifact_thresh;
    end
  end
  
  % initialize to store whether there was an artifact for each trial
  if ~exist('badEv','var') || isempty(badEv)
    combineArtLists = false;
    %badEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
    badEv = zeros(size(data.sampleinfo,1), 1);
  else
    combineArtLists = true;
  end
  autoEv = zeros(size(data.sampleinfo,1), 1);
  
  % find out what kind of artifacts we're dealing with
  theseArt = {};
  if isfield(cfg,'artfctdef')
    fn = fieldnames(cfg.artfctdef);
    for i = 1:length(fn)
      if isstruct(cfg.artfctdef.(fn{i})) && isfield(cfg.artfctdef.(fn{i}),'artifact') && ~isempty(cfg.artfctdef.(fn{i}).artifact)
        theseArt = cat(2,theseArt,fn{i});
      end
    end
  end
  % find out which samples were marked as artifacts
  if ~isempty(theseArt)
    artSamp = single(zeros(max(data.sampleinfo(:)),1));
    for i = 1:length(theseArt)
      for j = 1:size(cfg.artfctdef.(theseArt{i}).artifact,1)
        % mark that it was a particular type of artifact
        artSamp(cfg.artfctdef.(theseArt{i}).artifact(j,1):cfg.artfctdef.(theseArt{i}).artifact(j,2)) = find(ismember(theseArt,theseArt{i}));
      end
    end
    % save a list of trials with artifact status
    for k = 1:size(data.sampleinfo,1)
      if any(artSamp(data.sampleinfo(k,1):data.sampleinfo(k,2)) > 0)
        autoEv(k,1) = 1;
      end
    end
  end
  
  if combineArtLists
    % put the new artifacts into the old list
    rCount = 0;
    for i = 1:size(badEv,1)
      %if badEv(i,2) == 0
      if badEv(i,1) == 0
        rCount = rCount + 1;
        %if autoEv(rCount,2) == 1
        if autoEv(rCount) == 1
          %badEv(i,2) = 1;
          badEv(i,1) = 1;
        end
      end
    end
    if ~isempty(theseArt)
      if ~exist('artfctdef','var')
        artfctdef = cfg_manArt.artfctdef;
      else
        for i = 1:length(theseArt)
          if isfield(artfctdef,theseArt{i})
            artfctdef.(theseArt{i}).artifact = cat(1,artfctdef.(theseArt{i}).artifact,cfg.artfctdef.(theseArt{i}).artifact);
          else
            artfctdef.(theseArt{i}).artifact = cfg.artfctdef.(theseArt{i}).artifact;
          end
        end
      end
    end
  else
    badEv = autoEv;
    artfctdef = cfg.artfctdef;
  end
  
  cfg.artfctdef.reject = ana.artifact.reject;
  data = ft_rejectartifact(cfg,data);
end