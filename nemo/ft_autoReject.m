%% run FieldTrip's automatic artifact detection on the data

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
  
  % get the trial definition for automated FT artifact rejection
  trl = ft_findcfg(data.cfg,'trl');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for jump artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection
  cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.cutoff = 20;
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
  cfg.artfctdef.zvalue.interactive = 'yes';
  
  [cfg, artifact_jump] = ft_artifact_zvalue(cfg,data);
  
  %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   % look for muscle artifacts
  %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %   cfg = [];
  %   cfg.trl = trl;
  %   cfg.continuous = 'no';
  %
  %   % cutoff and padding
  %   % select a set of channels on which to run the artifact detection
  %   cfg.artfctdef.zvalue.channel = 'all';
  %   cfg.artfctdef.zvalue.cutoff      = 40;
  %   cfg.artfctdef.zvalue.trlpadding = ana.artifact.trlpadding;
  %   if strcmp(cfg.continuous,'yes')
  %     cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
  %   elseif strcmp(cfg.continuous,'no')
  %     cfg.artfctdef.zvalue.artpadding = ana.artifact.artpadding;
  %   end
  %   cfg.artfctdef.zvalue.fltpadding = ana.artifact.fltpadding;
  %
  %   % algorithmic parameters
  %   cfg.artfctdef.zvalue.bpfilter    = 'yes';
  %   if (sampleRate / 2) < 140
  %     cfg.artfctdef.zvalue.bpfreq      = [110 ((sampleRate / 2) - 1)];
  %   else
  %     cfg.artfctdef.zvalue.bpfreq      = [110 140];
  %   end
  %   cfg.artfctdef.zvalue.bpfiltord   = 6;
  %   cfg.artfctdef.zvalue.bpfilttype  = 'but';
  %   cfg.artfctdef.zvalue.hilbert     = 'yes';
  %   cfg.artfctdef.zvalue.boxcar      = 0.2;
  %
  %   % interactive artifact viewer
  %   %cfg.artfctdef.zvalue.interactive = 'yes';
  %
  %   [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,data);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % look for EOG artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.trl = trl;
  cfg.continuous = 'no';
  
  % cutoff and padding
  % select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
  %cfg.artfctdef.zvalue.channel = 'all';
  cfg.artfctdef.zvalue.channel = {'E127','E126','E128','E125'};
  cfg.artfctdef.zvalue.cutoff      = 4;
  cfg.artfctdef.zvalue.trlpadding = ana.artifact.artpadding;
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
  
  [cfg, artifact_EOG] = ft_artifact_zvalue(cfg,data);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % reject the automatically defined artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.artfctdef.jump.artifact = artifact_jump;
  %cfg.artfctdef.muscle.artifact = artifact_muscle;
  cfg.artfctdef.eog.artifact = artifact_EOG;
  
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
  fn = fieldnames(cfg.artfctdef);
  theseArt = {};
  for i = 1:length(fn)
    if isstruct(cfg.artfctdef.(fn{i})) && isfield(cfg.artfctdef.(fn{i}),'artifact') && ~isempty(cfg.artfctdef.(fn{i}).artifact)
      theseArt = cat(2,theseArt,fn{i});
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

if ~exist('badEv','var') || isempty(badEv)
  badEv = zeros(size(data.sampleinfo,1), 1);
end

if ~exist('artfctdef','var')
  artfctdef = [];
end