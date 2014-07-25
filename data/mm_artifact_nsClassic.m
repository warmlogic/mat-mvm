function [data,badChan_str,badEv,artfctdefEv,artfctdefSamp] = mm_artifact_nsClassic(data,exper,ana,elecfile,badChan_str,badEv,artfctdefEv,artfctdefSamp)

% function [data,badChan_str,badEv,artfctdefEv,artfctdefSamp] = mm_artifact_nsClassic(data,ana,elecfile,badChan_str,badEv,artfctdefEv,artfctdefSamp)
%
% Simulate Net Station's Artifact Detection Classic
%
% Set up in wrapper script:
%
% ana.artifact.type = {'nsClassic'};
%
% % time (in seconds) to check; use [-Inf Inf] to check the entire trial
% ana.artifact.checkArtSec = [-Inf Inf];
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% To detect bad channels:
%
% % repair the channel in all trials if ch is bad in >20% of segments
% ana.artifact.repairChan_percentBadTrials = 20;
%
% Will attempt to repair channels using spline interpolation on individual
% trials if fewer than ana.artifact.rejectTrial_nBadChan channels are bad
% and no bad channels are neighbors (if ana.artifact.allowBadNeighborChan =
% false)
%
% Does not repair eye channels; trials with blinks are always rejected
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% To detect bad segments:
%
% % fast average, reject the trial if over 100 uV
% ana.artifact.fast_threshold = 100;
%
% % differential average, reject the trial if over 50 uV
% ana.artifact.diff_threshold = 50;
%
% % reject the trial if 10+ bad channels
% ana.artifact.rejectTrial_nBadChan = 10;
%
% % reject the trial if neighboring channels are bad
% ana.artifact.allowBadNeighborChan = false;
%
% Again, trials with blinks are always rejected.
%

%% set up

%fast avg threshold
if ~isfield(ana.artifact,'fast_threshold')
  ana.artifact.fast_threshold = 100;
end
%differential avg threshold
if ~isfield(ana.artifact,'diff_threshold')
  ana.artifact.diff_threshold = 50;
end

if ~isfield(ana.artifact,'blink_threshold')
  ana.artifact.blink_threshold = 70;
end

%values for running avg
params = [.5, .5, .975, .025];

nTrial = length(data.trial);
nChan = size(data.trial{1},1);

if ~isfield(ana.artifact,'checkArtSec')
  ana.artifact.checkArtSec = [-Inf Inf];
end
tbeg = nearest(data.time{1}, ana.artifact.checkArtSec(1));
tend = nearest(data.time{1}, ana.artifact.checkArtSec(2));
nSamp = length(tbeg:tend);
sampNums = tbeg:tend;

% reject a trial  if it has more than X bad channels
if ~isfield(ana.artifact,'rejectTrial_nBadChan')
  ana.artifact.rejectTrial_nBadChan = 10;
end

% allow bad neighboring channels
if ~isfield(ana.artifact,'allowBadNeighborChan')
  ana.artifact.allowBadNeighborChan = false;
end

% fully repair the channel if it is bad in more than X% of trials (do not
% include eye channels by default)
if ~isfield(ana.artifact,'repairChan_percentBadTrials')
  ana.artifact.repairChan_percentBadTrials = 20;
end
nTrialThresh = fix(nTrial * (ana.artifact.repairChan_percentBadTrials / 100));

if ~isfield(ana.artifact,'doNotRepairEyes')
  ana.artifact.doNotRepairEyes = true;
end

if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
  % IMPORTANT: assumes we're using EGI's HCGSN
  
  eyeChan = {'E25', 'E8', 'E127', 'E126', 'E128', 'E125', 'E17', 'E21', 'E14'};
  
  eyeAndNeighbChan = {...
    'E48', 'E128', 'E127', 'E126', 'E125', 'E119', ...
    'E43', 'E32', 'E25', 'E21', 'E17', 'E14', 'E8', 'E1', 'E120', ...
    'E26', 'E22', 'E15', 'E9', 'E2', ...
    'E23', 'E18', 'E16', 'E10', 'E3', ...
    'E19', 'E11', 'E4'};
  
  %   % exclude eye channels and neighbors and all periphery channels
  %   eyeAndNeighbAndPeriphChan = {...
  %     'E1', 'E8', 'E14', 'E17', 'E21', 'E25', 'E32', 'E38', 'E43', 'E44', 'E48', 'E49', 'E56', 'E57', 'E63', 'E64', 'E68', 'E69', 'E73', 'E74', 'E81', 'E82', 'E88', 'E89', 'E94', 'E95', 'E99', 'E100', 'E107', 'E113', 'E114', 'E119', 'E120', 'E121', 'E125', 'E126', 'E127', 'E128', ...
  %     'E26', 'E22', 'E15', 'E9', 'E2', ...
  %     'E23', 'E18', 'E16', 'E10', 'E3', ...
  %     'E19', 'E11', 'E4'};
else
  warning('Must set up eye channels for your electrode layout!');
  keyboard
end

%% Detect bad channels

% denote channels that are chronically bad; must be repaired on all trials
fullyRepairChan = false(nChan,1);

% set up matrices to store whether artifacts were found; we always reject
% blinks and trials with too many bad channels
foundBlink = false(nTrial,1);
foundEyeMove = false(nTrial,1);
foundTooManyBadChan = false(nTrial,1);

% whether we reject the trial due to these depends on the number of bad
% channels
foundThresh = false(nTrial,nChan);
foundDead = false(nTrial,nChan);

% params
a = params(1);
b = params(2);
c = params(3);
d = params(4);

if ~ana.artifact.allowBadNeighborChan
  % check on neighbors
  cfg_nb = [];
  % cfg_nb.method = 'triangulation';
  cfg_nb.method = 'distance';
  cfg_nb.neighbourdist = 3.5;
  cfg_nb.elec = ana.elec;
  %if strcmp(cfg_ft.avgoverchan,'no')
  cfg_nb.neighbours = ft_prepare_neighbours(cfg_nb);
  
  hasBadNeighbor = false(nTrial,nChan);
end

fprintf('Checking for artifacts...');
for tr = 1:nTrial
  % init the two running averages
  fast = zeros(nChan,nSamp);
  slow = zeros(nChan,nSamp);
  
  % see if any channels had zero variance on this trial
  if ana.flatRef
    % assumes reference channel is at the end
    foundDead(tr,:) = cat(1,var(data.trial{tr}(~ismember(data.label,exper.refChan),:),0,2) == 0,false);
  else
    foundDead(tr,:) = var(data.trial{tr},0,2) == 0;
  end
  
  fast_start = 0;
  slow_start = mean(data.trial{tr}(:,tbeg:tbeg+10),2);
  
  for i = 1:nSamp
    % update the running averages
    if i > 1
      fast(:,i) = a*fast(:,i-1) + b*(data.trial{tr}(:,sampNums(i))-slow(:,i-1));
      slow(:,i) = c*slow(:,i-1) + d*data.trial{tr}(:,sampNums(i));
    else
      fast(:,i) = a*fast_start + b*(data.trial{tr}(:,sampNums(i))-slow_start);
      slow(:,i) = c*slow_start + d*data.trial{tr}(:,sampNums(i));
    end
  end
  
  % check for any threshold violations at every sample
  % by subtracting min mV from max mV and comparing the difference to
  % threshold
  badChanFast = any(abs(fast) >= ana.artifact.fast_threshold,2);
  
  % DIFFERENTIAL AVG is difference between a slower & faster moving avg
  badChanDiff = any((abs(fast)-abs(slow)) >= ana.artifact.diff_threshold,2);
  
  foundThresh(tr,:) = logical(badChanFast + badChanDiff);
  
  if ~ana.artifact.allowBadNeighborChan
    if any(foundThresh(tr,:))
      % check on neighboring channels
      
      % get all the channels that were bad
      badChans = data.label(foundThresh(tr,:));
      
      for ch = 1:length(badChans)
        % find the index of this bad channel
        chanIdx = ismember(cfg_nb.elec.label,badChans{ch});
        % see if any of the neighboring channels were bad
        if any(ismember(badChans,cfg_nb.neighbours(chanIdx).neighblabel))
          hasBadNeighbor(tr,ismember(data.label,badChans{ch})) = true;
        end
      end
    end
  end
  
end
fprintf('Done.\n');

%% set up which trials to reject and which channels to fix

foundArt = logical(foundThresh + foundDead);

% repair channel on all trials if it was bad on more than a given percent
fullyRepairChan(sum(foundArt,1) > nTrialThresh) = true;

% were any eye channels bad?
if ~fullyRepairChan(25)
  eog_upper_left = 25;
elseif ~fullyRepairChan(21)
  eog_upper_left = 21;
elseif ~fullyRepairChan(32)
  eog_upper_left = 32;
else
  eog_upper_left = [];
end
if ~fullyRepairChan(8)
  eog_upper_right = 8;
elseif ~fullyRepairChan(14)
  eog_upper_right = 14;
elseif ~fullyRepairChan(1)
  eog_upper_right = 1;
else
  eog_upper_right = [];
end
if ~fullyRepairChan(127)
  eog_lower_left = 127;
else
  eog_lower_left = [];
end
if ~fullyRepairChan(126)
  eog_lower_right = 126;
else
  eog_lower_right = [];
end
if ~fullyRepairChan(125)
  eog_horiz_left = 125;
else
  eog_horiz_left = [];
end
if ~fullyRepairChan(128)
  eog_horiz_right = 128;
else
  eog_horiz_right = [];
end

for tr = 1:nTrial
  %if any(foundArt(tr,ismember(data.label,eyeAndNeighbChan)))
  if any(foundArt(tr,:))
    fast = zeros(nChan,nSamp);
    slow = zeros(nChan,nSamp);
    
    fast_start = 0;
    slow_start = mean(data.trial{tr}(:,tbeg:tbeg+10),2);
    
    for i = 1:nSamp
      % update the running averages
      if i > 1
        fast(:,i) = a*fast(:,i-1) + b*(data.trial{tr}(:,sampNums(i))-slow(:,i-1));
        slow(:,i) = c*slow(:,i-1) + d*data.trial{tr}(:,sampNums(i));
      else
        fast(:,i) = a*fast_start + b*(data.trial{tr}(:,sampNums(i))-slow_start);
        slow(:,i) = c*slow_start + d*data.trial{tr}(:,sampNums(i));
      end
    end
    
    if ~isempty(eog_lower_left) && ~isempty(eog_upper_left)
      if any(abs(fast(eog_lower_left,:)) + abs(fast(eog_upper_left,:)) > 2*ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif ~isempty(eog_lower_right) && ~isempty(eog_upper_right)
      if any(abs(fast(eog_lower_right,:)) + abs(fast(eog_upper_right,:)) > 2*ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif ~isempty(eog_lower_left) && isempty(eog_upper_left)
      if any(abs(fast(eog_lower_left,:)) > ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif isempty(eog_lower_left) && ~isempty(eog_upper_left)
      if any(abs(fast(eog_upper_left,:)) > ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif ~isempty(eog_lower_right) && isempty(eog_upper_right)
      if any(abs(fast(eog_lower_right,:)) > ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif isempty(eog_lower_right) && ~isempty(eog_upper_right)
      if any(abs(fast(eog_upper_right,:)) > ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif ~isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
      if any(abs(fast(eog_horiz_left,:)) + abs(fast(eog_horiz_right,:)) > 2*ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif ~isempty(eog_horiz_left) && isempty(eog_horiz_right)
      if any(abs(fast(eog_horiz_left,:)) > ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
      if any(abs(fast(eog_horiz_right,:)) > ana.artifact.blink_threshold)
        foundBlink(tr) = true;
      end
    elseif isempty(eog_lower_left) && isempty(eog_upper_left) && isempty(eog_lower_right) && isempty(eog_upper_right)
      warning('Cannot check for eyeblinks because all eye channels are bad!');
    end
    
    if ~foundBlink(tr)
      if ~isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
        if any(abs(fast(eog_horiz_left,:)) + abs(fast(eog_horiz_right,:)) > 2*ana.artifact.blink_threshold)
          foundEyeMove(tr) = true;
        end
      elseif ~isempty(eog_horiz_left) && isempty(eog_horiz_right)
        if any(abs(fast(eog_horiz_left,:)) > 1*ana.artifact.blink_threshold)
          foundEyeMove(tr) = true;
        end
      elseif isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
        if any(abs(fast(eog_horiz_right,:)) > 1*ana.artifact.blink_threshold)
          foundEyeMove(tr) = true;
        end
      elseif isempty(eog_horiz_left) && isempty(eog_horiz_right)
        warning('Cannot check for eye movements because both horizontal eye channels are bad!');
      end
    end
    
  end
end

fprintf('**********************************************************************');
totalbadevents = sum(logical(sum(foundArt,2)));
fprintf('\nTotal number of bad events: %d\n',totalbadevents);
totalbadchans = sum(foundArt(:));
fprintf('\nTotal number of bad channels across all events: %d\n',totalbadchans);
fprintf('**********************************************************************\n');

if ana.artifact.doNotRepairEyes
  fullyRepairChan(ismember(data.label,eyeChan)) = false;
end
% turn it into a cell array
fullyRepairChan_str = data.label(fullyRepairChan);
% don't repair the reference channel
if ismember(exper.refChan,fullyRepairChan_str)
  fullyRepairChan_str = fullyRepairChan_str(~ismember(fullyRepairChan_str,exper.refChan));
end

% reject events with ana.artifact.rejectTrial_nBadChan or more bad channels
if ana.artifact.rejectTrial_nBadChan > 0
  foundTooManyBadChan(sum(foundArt,2) >= ana.artifact.rejectTrial_nBadChan) = true;
end

% see if there are bad neighbors, excluding eye channels
if ~ana.artifact.allowBadNeighborChan
  foundBadNeighborChan = logical(sum(hasBadNeighbor(:,~ismember(data.label,eyeAndNeighbChan)),2));
end

%% determine whether to repair channels on each trial

cfgChannelRepair = [];
cfgChannelRepair.channel = 'all';
cfgChannelRepair.method = 'spline';
cfgChannelRepair.elecfile = elecfile;

for tr = 1:nTrial
  trials = false(nTrial,1);
  trials(tr) = true;
  cfgChannelRepair.trials = trials;
  
  if any(foundArt(tr,~fullyRepairChan)) && ~foundBlink(tr) && ~foundEyeMove(tr) && ~foundTooManyBadChan(tr)
    if ana.artifact.allowBadNeighborChan || (~ana.artifact.allowBadNeighborChan && ~foundBadNeighborChan(tr))
      theseBadChan = data.label(foundArt(tr,:));
      theseBadChan = theseBadChan(~ismember(theseBadChan,fullyRepairChan_str));
      
      if ~isempty(theseBadChan)
        cfgChannelRepair.badchannel = theseBadChan;
        fprintf('\nTrial %d: using method=''%s'' to repair %d channels:%s\n',tr,cfgChannelRepair.method,length(cfgChannelRepair.badchannel),sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}));
        
        repaired_data = ft_channelrepair(cfgChannelRepair, data);
        data.trial{tr} = repaired_data.trial{1};
      end
    end
  end
end

%% do the rejection

if ~isfield(ana.artifact,'reject')
  ana.artifact.reject = 'complete';
end

cfg = [];
cfg.artfctdef.reject = ana.artifact.reject;

cfg.artfctdef.blink.artifact = data.sampleinfo(foundBlink,:);
cfg.artfctdef.eyemove.artifact = data.sampleinfo(foundEyeMove,:);
cfg.artfctdef.manybadchan.artifact = data.sampleinfo(foundTooManyBadChan,:);
if ~ana.artifact.allowBadNeighborChan
  cfg.artfctdef.badneighborchan.artifact = data.sampleinfo(foundBadNeighborChan,:);
end

% initialize to store whether there was an artifact for each trial
if isempty(badEv)
  combineArtLists = false;
  %badEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
  badEv = false(size(data.sampleinfo,1), 1);
else
  combineArtLists = true;
end

% find out what kind of artifacts we're dealing with
fn = fieldnames(cfg.artfctdef);
theseArt = {};
for i = 1:length(fn)
  if isstruct(cfg.artfctdef.(fn{i})) && isfield(cfg.artfctdef.(fn{i}),'artifact')
    theseArt = cat(2,theseArt,fn{i});
    if ~ismember(fn{i},artfctdefEv.types)
      artfctdefEv.types = cat(2,artfctdefEv.types,fn{i});
    else
      warning('''%s'' is already in the artifact types, if you continue you will overwrite the previous artifact information for this type.',fn{i});
      keyboard
    end
    artfctdefEv.(fn{i}) = false(size(badEv));
  end
end

% store artifacts for the only events we are checking
foundArtEv = false(size(data.sampleinfo,1),length(theseArt));

% find out which samples were marked as artifacts
if ~isempty(theseArt)
  artSamp = false(max(data.sampleinfo(:)),length(theseArt));
  for i = 1:length(theseArt)
    for j = 1:size(cfg.artfctdef.(theseArt{i}).artifact,1)
      % mark that it was a particular type of artifact
      artSamp(cfg.artfctdef.(theseArt{i}).artifact(j,1):cfg.artfctdef.(theseArt{i}).artifact(j,2),ismember(theseArt,theseArt{i})) = true;
    end
  end
  % save a list of trials with artifact status
  for k = 1:size(data.sampleinfo,1)
    if any(any(artSamp(data.sampleinfo(k,1):data.sampleinfo(k,2),:),1),2)
      foundArtEv(k,any(artSamp(data.sampleinfo(k,1):data.sampleinfo(k,2),:),1)) = true;
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
      if any(foundArtEv(rCount,:))
        %badEv(i,2) = 1;
        badEv(i,1) = 1;
        setTheseArt = theseArt(foundArtEv(rCount,:));
        for a = 1:length(setTheseArt)
          artfctdefEv.(setTheseArt{a})(i) = true;
        end
      end
    end
  end
  if ~isempty(theseArt)
    if isempty(artfctdefSamp)
      artfctdefSamp = cfg.artfctdef;
    else
      for i = 1:length(theseArt)
        if isfield(artfctdefSamp,theseArt{i})
          artfctdefSamp.(theseArt{i}).artifact = cat(1,artfctdefSamp.(theseArt{i}).artifact,cfg.artfctdef.(theseArt{i}).artifact);
        else
          artfctdefSamp.(theseArt{i}).artifact = cfg.artfctdef.(theseArt{i}).artifact;
        end
      end
    end
  end
else
  badEv = logical(sum(foundArtEv,2));
  artfctdefSamp = cfg.artfctdef;
  for a = 1:length(theseArt)
    artfctdefEv.(theseArt{a}) = foundArtEv(:,a);
  end
end

data = ft_rejectartifact(cfg, data);

%% do the channel repair on all remaining trials

if ~isempty(fullyRepairChan_str)
  cfgChannelRepair = [];
  cfgChannelRepair.channel = 'all';
  cfgChannelRepair.badchannel = fullyRepairChan_str;
  cfgChannelRepair.method = 'spline';
  cfgChannelRepair.elecfile = elecfile;
  fprintf('Repairing channels%s using method=''%s''...\n',sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}),cfgChannelRepair.method);
  data = ft_channelrepair(cfgChannelRepair, data);
  fprintf('Done.\n')
  
  badChan_str = cat(1,badChan_str,fullyRepairChan_str);
end

end
