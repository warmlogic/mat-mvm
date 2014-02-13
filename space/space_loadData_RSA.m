% load data for SPACE experiment; find pairs of stimulus presentations

verbose = false;

% set this to where you store the data
dataroot = pwd;

% list the conditions - NB: recognition misses (RgM) are not split by rc/fo
conditions = {'img_RgH_rc_spac_p1', 'img_RgH_rc_spac_p2', 'img_RgH_rc_mass_p1', 'img_RgH_rc_mass_p2', 'img_RgH_fo_spac_p1', 'img_RgH_fo_spac_p2', 'img_RgH_fo_mass_p1', 'img_RgH_fo_mass_p2', ...
  'img_RgM_spac_p1', 'img_RgM_spac_p2', 'img_RgM_mass_p1', 'img_RgM_mass_p2', ...
  'word_RgH_rc_spac_p1', 'word_RgH_rc_spac_p2', 'word_RgH_rc_mass_p1', 'word_RgH_rc_mass_p2', 'word_RgH_fo_spac_p1', 'word_RgH_fo_spac_p2', 'word_RgH_fo_mass_p1', 'word_RgH_fo_mass_p2', ...
  'word_RgM_spac_p1', 'word_RgM_spac_p2', 'word_RgM_mass_p1', 'word_RgM_mass_p2'};

% list the conditions without presentation number
cond_noP = {'img_RgH_rc_spac', 'img_RgH_rc_mass', 'img_RgH_fo_spac', 'img_RgH_fo_mass', ...
  'word_RgM_spac', 'word_RgM_mass', ...
  'word_RgH_rc_spac', 'word_RgH_rc_mass', 'word_RgH_fo_spac', 'word_RgH_fo_mass', ...
  'word_RgM_spac', 'word_RgM_mass'};

% the field where the EEG data is stored
param = 'trial';

% initialize
data = struct;

% load in the data and put it in the data struct
for i = 1:length(conditions)
  load(fullfile(dataroot,conditions{i}));
  data.(conditions{i}) = timelock;
end

% column numbers in trialinfo for selecting trials

% differentiate between multiple study phases
cols.phaseCount = 4;

% % trial number - do not need
% cols.trialNum = 5;

% unique number assigned to this stimulus
cols.stimNum = 6;
% image category number (1 = face, 2 = house)
cols.categNum = 7;

% % word-image pair number, only unique within a phase - do not need
% cols.pairNum = 13;

% find pairs of P1 and P2 for individual stimuli
for c = 1:length(cond_noP)
  p1_str = sprintf('%s_p1',cond_noP{c});
  p2_str = sprintf('%s_p2',cond_noP{c});
  
  for i = 1:size(data.(p1_str).(param),1)
    % get info about the P1 stimulus
    p1_trlInd = i;
    p1_phaseCount = data.(p1_str).trialinfo(p1_trlInd,cols.phaseCount);
    p1_stimNum = data.(p1_str).trialinfo(p1_trlInd,cols.stimNum);
    p1_categNum = data.(p1_str).trialinfo(p1_trlInd,cols.categNum);
    
    % see if this stimulus has a P2 (if not, was rejected for artifact)
    p2_trlInd = find(...
      data.(p2_str).trialinfo(:,cols.phaseCount) == p1_phaseCount & ...
      data.(p2_str).trialinfo(:,cols.stimNum) == p1_stimNum & ...
      data.(p2_str).trialinfo(:,cols.categNum) == p1_categNum);
    
    if ~isempty(p2_trlInd)
      if length(p2_trlInd) == 1
        % we found a match of this stimulus
        
        % get the data
        p1_data = squeeze(data.(p1_str).(param)(p1_trlInd,:,:));
        p2_data = squeeze(data.(p2_str).(param)(p2_trlInd,:,:));
        
        % after this point, compare similarity of P1 and P2...
        
      else
        error('Found more than one P2 match for P1 index %d!',i);
      end
    else
      if verbose
        fprintf('No match found for P1 index %d.\n',i);
      end
    end
  end
end

