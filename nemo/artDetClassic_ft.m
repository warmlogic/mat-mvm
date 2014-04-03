%% Artifact Detection Classic in FT?? %%
%
% To detect bad channels:
%   By segment:
%       fast average, bad if over 100 mV
%       differential average, bad if over 50 mV
%       zero variance across channel (dead)
%   Whole recording:
%       if ch is bad in >20% of segments
%
% Reject bad channels w/rmBadChan
%
% To detect bad segments:
%   10+ bad channels
%   ft_artifact_threshold w/thresholds of 70 mV
%
% Reject bad segments w/ft_rejectartifact

%% Detect bad channels

fast_threshold = 100; %fast avg threshold
diff_threshold = 50; %differential avg threshold
params = [.5, .5, .975, .025]; %values for running avg

% Fast average:
% init the two running averages
fast = zeros(1,length(data));
slow = zeros(1,length(data));
var = zeros(1,length(data));
%ind = logical(zeros(1,length(dat)));

% params
a = params(1);
b = params(2);
c = params(3);
d = params(4);

fast_start = 0;
slow_start = mean(data(1:10));

for i = 1:length(data) %rows/channels only? need another for loop for cols?
    % update the running averages
    if i > 1
        fast(i) = a*fast(i-1) + b*(data(i)-slow(i-1));
        slow(i) = c*slow(i-1) + d*data(i);
    else
        fast(i) = a*fast_start + b*(data(i)-slow_start);
        slow(i) = c*slow_start + d*data(i);
    end
    
    %calculate variance of moving average
    var(i) = ((data(i)-slow(i))*2)/(length(data)-1);
    
    % check for thresh
    %ind(i) = abs(fast(i))>=thresh;
    
end

% check for thresh, logical vector with 1s where shows blink (bad ch)
badChanFast = logical(abs(fast)>=fast_threshold);

% DIFFERENTIAL AVG is difference between a slower & faster moving avg
badChanDiff = logical((abs(fast)-abs(slow))>=diff_threshold);

% check for zero variance, set 0 var to be true so later we can reject all
% 1's in logical array
badChanVar = logical(var == 0);
    

% make 1 bad channel variable?
badChanLong = cat(2,badChanFast,badChanDiff,badChanVar);

% if a channel is bad for >20% of segments, reject ch from whole recording
y=zeros(size(badChanLong));

for i = 1:length(badChanLong)
   y(i) = sum(badChanLong==badChanLong(i));    
    
end






    % collect the bad channels
    if rejArt_artDetClassic
        if ~exist('badChanFast','var')
            badChanFast = [];
        end
        if ~exist('badChanDiff','var')
            badChanDiff = [];
        end
        if ~exist('badChanVar','var')
            badChanVar = [];
        end
        %want to add the 3 arrays instead of cat()?
        badChan = unique(cat(2,badChanFast,badChanDiff,badChanVar));
        badChan_str = cell(length(badChan,1));
        for i = 1:length(badChan)
            badChan_str{i} = sprintf('E%d',badChan(i));
        end
    end
    
    % reject bad channels
    if rejArt_artDetClassic && ~isempty(badChan)
        cfg_rv = [];
        cfg_rv.channel = eval(sprintf('{''all''%s};',sprintf(repmat(' ''-E%d''',1,length(badChan)),badChan)));
        cfg_rv.keepchannel = 'no';
        %cfg_rv.keepchannel = 'nan';
        
        data = ft_rejectvisual(cfg_rv,data); % change to something automatic
    elseif rejArt_artDetClassic && isempty(badChan)
        fprintf('No bad channels to reject!\n');
    end
    
    % ft_channelselection.m args to exclude chs & groups:
    % {'all', '-POz', '-Fp1', '-EOG'}
        
        %% Detect bad segments
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % look for threshold artifacts - only if using EGI's HydroCel GSN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(elecfile,'GSN-HydroCel-129.sfp') || strcmp(elecfile,'GSN-HydroCel-128.sfp')
            cfg = [];
            cfg.continuous = 'no';
            % get the trial definition for automated FT artifact rejection
            cfg.trl = ft_findcfg(data.cfg,'trl');
            
            % % exclude eye channels and neighbors - assumes we're using EGI's HCGSN
            % cfg.artfctdef.threshold.channel = {'all', ...
            %   '-E48', '-E128', '-E127', '-E126', '-E125', '-E119', ...
            %   '-E43', '-E32', '-E25', '-E21', '-E17', '-E14', '-E8', '-E1', '-E120', ...
            %   '-E26', '-E22', '-E15', '-E9', '-E2', ...
            %   '-E23', '-E18', '-E16', '-E10', '-E3', ...
            %   '-E19', '-E11', '-E4'};
            
            % exclude eye channels and neighbors and all periphery channels - assumes we're using EGI's HCGSN
            cfg.artfctdef.threshold.channel = {'all', '-E1', '-E8', '-E14', '-E17', '-E21', '-E25', '-E32', '-E38', '-E43', '-E44', '-E48', '-E49', '-E56', '-E57', '-E63', '-E64', '-E68', '-E69', '-E73', '-E74', '-E81', '-E82', '-E88', '-E89', '-E94', '-E95', '-E99', '-E100', '-E107', '-E113', '-E114', '-E119', '-E120', '-E121', '-E125', '-E126', '-E127', '-E128', ...
                '-E26', '-E22', '-E15', '-E9', '-E2', ...
                '-E23', '-E18', '-E16', '-E10', '-E3', ...
                '-E19', '-E11', '-E4'};
            
            cfg.artfctdef.threshold.bpfilter = 'yes';
            cfg.artfctdef.threshold.bpfreq = [0.3 30];
            cfg.artfctdef.threshold.bpfiltord = 4;
            
            cfg.artfctdef.threshold.min = ana.artifact.threshmin;
            cfg.artfctdef.threshold.max = ana.artifact.threshmax;
            cfg.artfctdef.threshold.range = ana.artifact.threshrange;
            
            fprintf('\nUsing EGI HydroCel GSN...\nChecking for voltages above %.1f uV and below %.1f uV, or peak-to-peak range of %.1f (excludes eye channels and neighbors and periphery channels)...\n',cfg.artfctdef.threshold.max,cfg.artfctdef.threshold.min,cfg.artfctdef.threshold.range);
            
            % auto mark zvalue artifacts
            [cfg, artifact_thresh] = ft_artifact_threshold(cfg, data);
        else
            warning('Not using EGI HydroCel GSN 128/129 electrode file! Threshold artifacts are not being assessed!!');
        end
        
        %% Reject automatically defined artifacts
        % do the actual rejection of artifact trials (complete or parial rejection)
if (rejArt_nsAuto || rejArt_zeroVar) && foundArt
  cfg.artfctdef.reject = ana.artifact.reject;
  data = ft_rejectartifact(cfg,data);
end
