function [data,fullyRepairChan_str,badEv,artfctdef] = mm_artifact_nsClassic_wZtran(subject,dirs,data,ana,elecfile,badChan_str,badEv)

% function [data,fullyRepairChan_str,badEv,artfctdef] = mm_artifact_nsClassic(subject,dirs,data,ana,elecfile,badChan_str,badEv)
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
% log for info
resultsDir = fullfile(dirs.localDir,dirs.dataDir,'nsClassic_wztran_results',subject);
if ~exist(resultsDir, 'dir')
    mkdir(fullfile(dirs.localDir,dirs.dataDir,'nsClassic_wztran_results'),subject);
end
resultsFileLoc = fullfile(dirs.localDir,dirs.dataDir,'nsClassic_wztran_results',subject);
resultsFile = fopen(fullfile(resultsFileLoc,'nsClassic_wztran_results.txt'), 'w');

% z basic artifact threshold
if ~isfield(ana.artifact,'basic_art_z')
    ana.artifact.basic_art_z = 30;
end

% z eog artifact threshold for blinks
if ~isfield(ana.artifact,'eog_art_z')
    ana.artifact.eog_art_z = 3.5;
end

% z jump artifact threshold for eye muscle movement
if ~isfield(ana.artifact,'jump_art_z')
    ana.artifact.jump_art_z = 50;
end

% %fast avg threshold
% if ~isfield(ana.artifact,'fast_threshold')
%   ana.artifact.fast_threshold = 100;
% end
% %differential avg threshold
% if ~isfield(ana.artifact,'diff_threshold')
%   ana.artifact.diff_threshold = 50;
% end
%
% if ~isfield(ana.artifact,'blink_threshold')
%   ana.artifact.blink_threshold = 70;
% end

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

% reject a trial  if it has more than X bad channels
if ~isfield(ana.artifact,'rejectTrial_nBadChan')
    ana.artifact.rejectTrial_nBadChan = 10;
end

% allow bad neighboring channels
if ~isfield(ana.artifact,'allowBadNeighborChan')
    ana.artifact.allowBadNeighborChan = true;
end

% fully repair the channel if it is bad in more than X% of trials (do not
% include eye channels by default)
if ~isfield(ana.artifact,'repairChan_percentBadTrials')
    ana.artifact.repairChan_percentBadTrials = 20;
end
nTrialThresh = fix(nTrial * (ana.artifact.repairChan_percentBadTrials / 100));

if ~isfield(ana.artifact,'doNotRepairEyes')
    ana.artifact.doNotRepairEyes = false;
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
    cfg_nb.neighbourdist = 3;
    cfg_nb.elec = ana.elec;
    %if strcmp(cfg_ft.avgoverchan,'no')
    cfg_nb.neighbours = ft_prepare_neighbours(cfg_nb);
    
    hasBadNeighbor = false(nTrial,nChan);
end

fprintf('Checking for artifacts...');
% for tr = 1:nTrial
%   % init the two running averages
%   fast = zeros(nChan,nSamp);
%   slow = zeros(nChan,nSamp);
%
%   % see if any channels had zero variance on this trial
%   foundDead(tr,:) = var(data.trial{tr},0,2) == 0;
%
%   fast_start = 0;
%   slow_start = mean(data.trial{tr}(:,tbeg:tbeg+10),2);
%
%   for i = tbeg:tend
%     % update the running averages
%     if i > 1
%       fast(:,i) = a*fast(:,i-1) + b*(data.trial{tr}(:,i)-slow(:,i-1));
%       slow(:,i) = c*slow(:,i-1) + d*data.trial{tr}(:,i);
%     else
%       fast(:,i) = a*fast_start + b*(data.trial{tr}(:,i)-slow_start);
%       slow(:,i) = c*slow_start + d*data.trial{tr}(:,i);
%     end
%   end
%
%   % check for any threshold violations at every sample
%   % by subtracting min mV from max mV and comparing the difference to
%   % threshold
%   badChanFast = any(abs(fast) >= ana.artifact.fast_threshold,2);
%
%   % DIFFERENTIAL AVG is difference between a slower & faster moving avg
%   badChanDiff = any((abs(fast)-abs(slow)) >= ana.artifact.diff_threshold,2);
%
%   foundThresh(tr,:) = logical(badChanFast + badChanDiff);
%
%   if ~ana.artifact.allowBadNeighborChan
%     if any(foundThresh(tr,:))
%       % check on neighboring channels
%
%       % get all the channels that were bad
%       badChans = data.label(foundThresh(tr,:));
%
%       for ch = 1:length(badChans)
%         % find the index of this bad channel
%         chanIdx = ismember(cfg_nb.elec.label,badChans{ch});
%         % see if any of the neighboring channels were bad
%         if any(ismember(badChans,cfg_nb.neighbours(chanIdx).neighblabel))
%           hasBadNeighbor(tr,ismember(data.label,badChans{ch})) = true;
%         end
%       end
%     end
%   end
%
% end
% fprintf('Done.\n');

for tr = 1:nTrial
    
    % init the two running averages
    fast = zeros(nChan,nSamp);
    slow = zeros(nChan,nSamp);
    
    % init difference
    datadiff = zeros(nChan,nSamp);
    % init std dev
    datastddev = zeros(nChan,nSamp);
    % init z score
    %zdata = zeros(nChan,nSamp);
    
    % see if any channels had zero variance on this trial
    foundDead(tr,:) = var(data.trial{tr},0,2) == 0;
    
    fast_start = 0;
    slow_start = mean(data.trial{tr}(:,tbeg:tbeg+10),2);
    
    for i = tbeg:tend
        % update the running average
        if i > 1
            fast(:,i) = a*fast(:,i-1) + b*(data.trial{tr}(:,i)-slow(:,i-1));
            slow(:,i) = c*slow(:,i-1) + d*data.trial{tr}(:,i);
            
            % update difference of data-fast for that sample
            datadiff(:,i) = data.trial{tr}(:,i)-fast(:,i-1);
            
            % update standard deviation
            % datastddev(:,i) = sqrt(i/(i-1)*sum(((data.trial{tr}(:,i)-fast(:,i-1)).^2),2));
            datastddev(:,i) = std(data.trial{tr}(:,tbeg:tbeg+10),0,2);
            
        else
            fast(:,i) = a*fast_start + b*(data.trial{tr}(:,i)-slow_start);
            slow(:,i) = c*slow_start + d*data.trial{tr}(:,i);
            
            % difference of data-fast for that sample
            datadiff(:,i) = data.trial{tr}(:,i)-fast(:,i);

            % standard deviation
            datastddev(:,i) = std(data.trial{tr}(:,tbeg:tbeg+10),0,2);
             
        end %if  
    end %i
    
    % smooth datadiff & stddev before calculating zdata
    datadiff = ft_preproc_smooth(datadiff, tend)*tend;
    datastddev = ft_preproc_smooth(datastddev, tend)*tend;
    
    % z score matrix
    zdata = datadiff./datastddev;
    
    
    % check for any threshold violations at every sample
    % compare z scores of data to threshold
    foundThresh(tr,:) = any(abs(zdata) >= ana.artifact.basic_art_z,2);
    
    
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
    
end %tr
fprintf('Done.\n');

% %% from ft_artifact_zvalue
%
% % z score collapsing across all channels per trial
%
% % read the data and apply preprocessing options
% sumvalue = zeros(nChan, nTrial);
% sumsqared = zeros(nChan, nTrial);
% numsample = zeros(nChan, nTrial);
% fprintf('searching trials');
%
% for tr = 1:nTrial
%     fprintf('.');
%
%     % store per trial the sum and the sum-of-squares
%     sumvalue(:,tr) = sum(data.trial{tr},2);
%     sumsqared(:,tr) = sum(data.trial{tr}.^2,2);
%     numsample(:,tr) = size(data.trial{tr},2);
%
% end % for tr
%
% % smooth the data (instead of zscores by sample - better or worse?)
% sumvalue = ft_preproc_smooth(sumvalue, nSamp)*nSamp;
% sumsqared = ft_preproc_smooth(sumsqared, nSamp)*nSamp;
% numsample = ft_preproc_smooth(numsample, nSamp)*nSamp;
%
%
% % compute the average and the standard deviation
% datavg = sumvalue./numsample;
% datstd = sqrt(sumsqared./numsample - (sumvalue./numsample).^2);
%
% zmax = cell(1, nTrial);
% zsum = cell(1, nTrial);
% zindx = cell(1, nTrial);
%
% % create a vector that indexes the trials, or is all 1, in order
% % to a per trial z-scoring, or use a static std and mean
% indvec = 1:nTrial;
%
% for tr = 1:nTrial
%
%     % initialize some matrices
%     zmax{tr}  = -inf + zeros(1,size(data.trial{tr},2));
%     zsum{tr}  = zeros(1,size(data.trial{tr},2));
%     zindx{tr} = zeros(1,size(data.trial{tr},2));
%
%     nsmp          = size(data.trial{tr},2);
%     zdata         = (data.trial{tr} - datavg(:,indvec(tr)*ones(1,nsmp)))./datstd(:,indvec(tr)*ones(1,nsmp));  % convert the filtered data to z-values
%     zsum{tr}   = nansum(zdata,1);                   % accumulate the z-values over channels
%     [zmax{tr},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
%     zindx{tr}      = data.label(ind);                % also remember the channel number that has the largest z-value
%
% end % for tr
%
%
% for tr = 1:nTrial
%     % z score of all channels per trial
%     zsum{tr} = zsum{tr} ./ sqrt(nChan);
%
%     % compare z-transformed data to threshold
%     foundThresh = any(zsum(tr) >= ana.artifact.basic_art_z,2);
%     % see if any channels had zero variance on this trial
%     foundDead(tr,:) = var(data.trial{tr},0,2) == 0;
%
%     if ~ana.artifact.allowBadNeighborChan
%         if any(foundThresh(tr,:))
%             % check on neighboring channels
%
%             % get all the channels that were bad
%             badChans = data.label(foundThresh(tr,:));
%
%             for ch = 1:length(badChans)
%                 % find the index of this bad channel
%                 chanIdx = ismember(cfg_nb.elec.label,badChans{ch});
%                 % see if any of the neighboring channels were bad
%                 if any(ismember(badChans,cfg_nb.neighbours(chanIdx).neighblabel))
%                     hasBadNeighbor(tr,ismember(data.label,badChans{ch})) = true;
%                 end
%             end
%         end
%     end
%
% end
%
% fprintf('Done.\n');

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
        
        % init difference
        datadiff = zeros(nChan,nSamp);
        % init std dev
        datastddev = zeros(nChan,nSamp);
        % init z score
        %zdata = zeros(nChan,nSamp);
        
        fast_start = 0;
        slow_start = mean(data.trial{tr}(:,tbeg:tbeg+10),2);
        
        for i = tbeg:tend
            % update the running average
            if i > 1
                fast(:,i) = a*fast(:,i-1) + b*(data.trial{tr}(:,i)-slow(:,i-1));
                slow(:,i) = c*slow(:,i-1) + d*data.trial{tr}(:,i);
                
                % update difference of data-fast for that sample
                datadiff(:,i) = data.trial{tr}(:,i)-fast(:,i-1);
                
                % update standard deviation
                % datastddev(:,i) = sqrt(i/(i-1)*sum(((data.trial{tr}(:,i)-fast(:,i-1)).^2),2));
                datastddev(:,i) = std(data.trial{tr}(:,tbeg:tbeg+10),0,2);
                
            else
                fast(:,i) = a*fast_start + b*(data.trial{tr}(:,i)-slow_start);
                slow(:,i) = c*slow_start + d*data.trial{tr}(:,i);
                
                % difference of data-fast for that sample
                datadiff(:,i) = data.trial{tr}(:,i)-fast(:,i);
                
                % standard deviation
                datastddev(:,i) = std(data.trial{tr}(:,tbeg:tbeg+10),0,2);
                
            end %if
        end %i
        
        % smooth datadiff & stddev before calculating zdata
        datadiff = ft_preproc_smooth(datadiff, tend)*tend;
        datastddev = ft_preproc_smooth(datastddev, tend)*tend;
        
        % z score matrix
        zdata = datadiff./datastddev;
        
        if ~isempty(eog_lower_left) && ~isempty(eog_upper_left)
            if any(abs(zdata(eog_lower_left,:)) + abs(zdata(eog_upper_left,:)) > 2*ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif ~isempty(eog_lower_right) && ~isempty(eog_upper_right)
            if any(abs(zdata(eog_lower_right,:)) + abs(zdata(eog_upper_right,:)) > 2*ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif ~isempty(eog_lower_left) && isempty(eog_upper_left)
            if any(abs(zdata(eog_lower_left,:)) > ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif isempty(eog_lower_left) && ~isempty(eog_upper_left)
            if any(abs(zdata(eog_upper_left,:)) > ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif ~isempty(eog_lower_right) && isempty(eog_upper_right)
            if any(abs(zdata(eog_lower_right,:)) > ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif isempty(eog_lower_right) && ~isempty(eog_upper_right)
            if any(abs(zdata(eog_upper_right,:)) > ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif ~isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
            if any(abs(zdata(eog_horiz_left,:)) + abs(zdata(eog_horiz_right,:)) > 2*ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif ~isempty(eog_horiz_left) && isempty(eog_horiz_right)
            if any(abs(zdata(eog_horiz_left,:)) > ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
            if any(abs(zdata(eog_horiz_right,:)) > ana.artifact.eog_art_z)
                foundBlink(tr) = true;
            end
        elseif isempty(eog_lower_left) && isempty(eog_upper_left) && isempty(eog_lower_right) && isempty(eog_upper_right)
            warning('Cannot check for eyeblinks because all eye channels are bad!');
        end
        
        if ~foundBlink(tr)
            if ~isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
                if any(abs(zdata(eog_horiz_left,:)) + abs(zdata(eog_horiz_right,:)) > 2*ana.artifact.jump_art_z)
                    % keyboard
                    foundEyeMove(tr) = true;
                end
            elseif ~isempty(eog_horiz_left) && isempty(eog_horiz_right)
                if any(abs(zdata(eog_horiz_left,:)) > 1*ana.artifact.jump_art_z)
                    foundEyeMove(tr) = true;
                end
            elseif isempty(eog_horiz_left) && ~isempty(eog_horiz_right)
                if any(abs(zdata(eog_horiz_right,:)) > 1*ana.artifact.jump_art_z)
                    foundEyeMove(tr) = true;
                end
            elseif isempty(eog_horiz_left) && isempty(eog_horiz_right)
                warning('Cannot check for eye movements because both horizontal eye channels are bad!');
            end
        end
        
        % more plots for checking zdata
        if any(foundArt(tr,:))
            h = figure;
            figName = sprintf('zdata_%s', num2str(tr));
            hold on
            subplot(3,1,[1,1]);
            plot(zdata(:,:),'LineWidth',1)
            ylim([-10 10]);
            title(sprintf('foundArt_%s', num2str(tr)));
            if any(foundBlink(tr,:))
                subplot(3,1,[2,2]);
                hold on
                plot(zdata(eog_upper_left,:),'m','LineWidth',1)
                plot(zdata(eog_lower_left,:),'k','LineWidth',1)
                plot(zdata(eog_horiz_left,:),'b','LineWidth',2)
                plot(zdata(eog_horiz_right,:),'r','LineWidth',2)
                hold off
                ylim([-10 10]);
                title(sprintf('foundBlink_%s', num2str(tr)));
            end
            if any(foundEyeMove(tr,:))
                subplot(3,1,[3,3]);
                hold on
                plot(zdata(eog_upper_left,:),'m','LineWidth',1)
                plot(zdata(eog_lower_left,:),'k','LineWidth',1)
                plot(zdata(eog_horiz_left,:),'b','LineWidth',2)
                plot(zdata(eog_horiz_right,:),'r','LineWidth',2)
                hold off
                ylim([-10 10]);
                title(sprintf('foundEyeMove_%s', num2str(tr)));
            end
            hold off
            saveas(h,fullfile(resultsFileLoc, figName),'png');
            close(h);
        end
        
        
    end
end




if ana.artifact.doNotRepairEyes
    fullyRepairChan(ismember(data.label,eyeChan)) = false;
end
% turn it into a cell array
fullyRepairChan_str = cat(1,badChan_str,data.label(fullyRepairChan));

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
                fprintf(resultsFile,'\nTrial %d: Repairing %d channels:%s\n',tr,length(cfgChannelRepair.badchannel),sprintf(repmat(' %s',1,length(cfgChannelRepair.badchannel)),cfgChannelRepair.badchannel{:}));
                
                repaired_data = ft_channelrepair(cfgChannelRepair, data);
                data.trial{tr} = repaired_data.trial{1};
            end
        end
    end
end

%% figures & log for results

% log information

% channels repaired by segment - prints during repair
fprintf(resultsFile,'\nEnd of channels repaired by trial.\n\n');

% total segments marked bad b/c of blinks
totalBlinks = sum(logical(sum(foundBlink,2)));
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

%% do the rejection

if ~isfield(ana.artifact,'reject')
    ana.artifact.reject = 'complete';
end

cfg = [];
cfg.artfctdef.reject = ana.artifact.reject;

cfg.artfctdef.blink.artifact = data.sampleinfo(foundBlink,:);
cfg.artfctdef.eyemove.artifact = data.sampleinfo(foundEyeMove,:);
cfg.artfctdef.manybadchan.artifact = data.sampleinfo(foundTooManyBadChan,:);

% totalBlinks = sum(logical(sum(foundBlink,2)));
% fprintf(resultsFile,'Blink artifacts: %d\n', totalBlinks);
% totalEyeMove = sum(foundEyeMove,:);
% fprintf(resultsFile,'Eye movement artifacts: %d\n', totalEyeMove);
% totalTooManyBadChan = sum(foundTooManyBadChan,:);
% fprintf(resultsFile,'Too many bad channels: %d\n', totalTooManyBadChan);

if ~ana.artifact.allowBadNeighborChan
    cfg.artfctdef.badneighborchan.artifact = data.sampleinfo(foundBadNeighborChan,:);
    fprintf(resultsFile,'Blink artifacts: %d\n',data.sampleinfo(foundBadNeighborChan,:));
end

% initialize to store whether there was an artifact for each trial
if ~exist('badEv','var') || isempty(badEv)
    combineArtLists = false;
    %badEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
    badEv = zeros(size(data.sampleinfo,1), 1);
else
    combineArtLists = true;
end
manualEv = zeros(size(data.sampleinfo,1), 1);

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
            manualEv(k,1) = 1;
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
            %if manualEv(rCount,2) == 1
            if manualEv(rCount) == 1
                %badEv(i,2) = 1;
                badEv(i,1) = 1;
            end
        end
    end
    if ~isempty(theseArt)
        if ~exist('artfctdef','var')
            artfctdef = cfg.artfctdef;
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
    badEv = manualEv;
    artfctdef = cfg.artfctdef;
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
    fprintf(resultsFile,'Done.\n');
    fclose(resultsFile);
    fprintf('Done.\n')
end

end
