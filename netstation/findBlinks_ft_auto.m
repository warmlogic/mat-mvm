%%Define the trials w/artifacts

% get the trial definition for automated FT artifact rejection
  trl = ft_findcfg(data.cfg,'trl');
  badEv = zeros(size(data.sampleinfo,1), 1);


  cfg = [];
  cfg.artfctdef.threshold = 70;
  cfg.artfctdef.params = [.5, .5, .975, .025];
  fast_start = 0;
  slow_start = mean(data(1:10));
  
  for i = 1:length(data)
       % update the running averages
        if i > 1
            fast(i) = a*fast(i-1) + b*(data(i)-slow(i-1));
            slow(i) = c*slow(i-1) + d*data(i);
        else
            fast(i) = a*fast_start + b*(data(i)-slow_start);
            slow(i) = c*slow_start + d*data(i);
        end
        if abs(fast(i)) >= cfg.artfctdef.threshold
            badEv(i) = 
         
  end
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % reject the automatically defined artifacts
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %cfg.artfctdef.jump.artifact = artifact_jump;
  %cfg.artfctdef.muscle.artifact = artifact_muscle;
  %cfg.artfctdef.eog.artifact = artifact_EOG;
  
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