%% load the analysis details

% keyboard; % I put this here for times I hit "Run" by mistake, so I can still stop the whole thing from running.

expName = 'EBIRD';

subDir = '';
% subDir = 'ica_test';
% subDir = 'new_art500';
% subDir = 'ft_data/data_art_ftManual_ftICA';

if ~isempty(subDir)
  warning('Loading data from subDir ''%s''...',subDir);
end
dataDir = fullfile(expName,'EEG','Sessions','ftpp',subDir);
% keyboard;
% Possible locations of the data files (dataroot)
serverDir = fullfile(filesep,'Volumes','curranlab','Data');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dreamDir = fullfile(filesep,'data','projects','curranlab');
localDir = fullfile(getenv('HOME'),'data');


% pick the right dataroot
if exist('serverDir','var') && exist(serverDir,'dir')
  dataroot = serverDir;
  %runLocally = 1;
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
  %runLocally = 1;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
  %runLocally = 0;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

procDir = fullfile(dataroot,dataDir,'ft_data/data_art_nsClassic/tla');
% procDir = fullfile(dataroot,dataDir,'ft_data/data_art_ftManual_ftICA/tla');

subjects = {
% %   'EBIRD049'; % Pilot. (due to short ses1 match, missing ses2 name)
% %   'EBIRD002'; % Pilot. (due to short ses1 match, missing ses2 name)
% %   'EBIRD003'; % Pilot. (due to missing ses7 name) - NB: LAST PILOT TO BE REPLACED
% %   'EBIRD004'; % DNF. Dropout. Last session: 8.
  'EBIRD005';
% %   'EBIRD006'; % DNF. Dropout. Last session: 2.
  'EBIRD007';
  'EBIRD008';
  'EBIRD009';
  'EBIRD010';
  'EBIRD011';
  'EBIRD012';
%   %'EBIRD013'; % DNF. Dropout. Last session: 5. Lost session 6 in HD crash.
%   %'EBIRD014'; % DNF. Rejected. Last session: 1.
%   %'EBIRD015'; % DNF. Lost in HD crash.
%   %'EBIRD016'; % DNF. Lost in HD crash.
%   %'EBIRD017'; % DNF. Lost in HD crash.
  'EBIRD018';
  'EBIRD019';
  'EBIRD020';
%   'EBIRD021'; %Bad behvior
%   %'EBIRD022'; % DNF. Dropout. Last session: 8.
%   %'EBIRD023'; % DNF. Dropout. Last session: 1.
  'EBIRD024';
  'EBIRD025';
  'EBIRD027';
  'EBIRD029';
  'EBIRD032';
  'EBIRD034';
  'EBIRD042';
  };

% only one cell, with all session names
% sesNames = {'session_1','session_8','session_9'};
sesNames = {'session_2','session_3','session_4', 'session_5','session_6','session_7'};
% sesNames = {'session_1','session_8'};
% sesNames = {'session_1'};

% replaceDataroot = {'/Users/matt/data','/Volumes/curranlab/Data'};
replaceDataroot = true;

[exper,ana,dirs,files] = mm_loadAD(procDir,subjects,sesNames,replaceDataroot,[],subDir);

%KLUDGE TO FIX FOLDER ORANIZATION PROBLEM. MOVE session_1 out of ica_test 
% dirs.dataDir = dataDir; 
% dirs.saveDirStem = dataDir;

% adFile = '/Users/matt/data/EBIRD/EEG/Sessions/ftpp/ft_data/match_stim_eq0_art_ftManual_ftICA/tla/analysisDetails.mat';
% % [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails(adFile,adFile5,true,true,true);
% 
% % % server_adFile = '';
% % if exist(server_adFile,'file')
% %   mm_mergeAnalysisDetails(adFile,server_adFile,true,false,false);
% % end
% 
% % [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);
% [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,false);

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% list the event values to analyze; specific to each experiment

% ana.trl_order.name_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'imgCond', 'isSubord', 'response', 'rt', 'acc'};

ana.eventValues = repmat({{'name_stim'}},1,length(exper.sessions));

ana.eventValuesSplit = repmat({...
  {{'basic', 'subord'} ...
  }},1,length(exper.sessions));

for ses = 1:length(exper.sessions)
  ana.trl_expr{ses} = {...
    {...
    sprintf('eventNumber == %d & isSubord == 0 ',find(ismember(exper.eventValues{ses},'name_stim'))) ...
    sprintf('eventNumber == %d & isSubord == 1 ',find(ismember(exper.eventValues{ses},'name_stim'))) ...
    }...
    };
end

% ana.trl_order.nametrain_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'block', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'imgCond', 'isSubord', 'response', 'rt', 'acc'};

% can only add nametrain stimuli if they were processed in initial script
addNametrain = false;
if addNametrain
  ana.eventValues{1} = cat(2,'nametrain_stim',ana.eventValues{1});
  
  ana.eventValuesSplit{1} = cat(2,{{'basic', 'subord'}},ana.eventValuesSplit{1});
  
  ana.trl_expr{1} = cat(2,...
    {{...
    sprintf('eventNumber == %d & isSubord == 0 ',find(ismember(exper.eventValues{1},'nametrain_stim'))) ...
    sprintf('eventNumber == %d & isSubord == 1 ',find(ismember(exper.eventValues{1},'nametrain_stim'))) ...
    }},...
    ana.trl_expr{1});
end

%% split by correct/incorrect

% ana.trl_order.name_stim = {'eventNumber', 'sesType', 'phaseType', 'phaseCount', 'trial', 'familyNum', 'speciesNum', 'exemplarNum', 'imgCond', 'isSubord', 'response', 'rt', 'acc'};

ana.eventValues = repmat({{'name_stim'}},1,length(exper.sessions));

% ana.eventValuesSplit = repmat({...
%   {{'basic_cor', 'basic_inc', 'subord_cor', 'subord_inc'} ...
%   }},1,length(exper.sessions));

ana.eventValuesSplit = repmat({...
  {{'basic_cor', 'subord_cor'} ...
  }},1,length(exper.sessions));

for ses = 1:length(exper.sessions)
  ana.trl_expr{ses} = {...
    {...
    sprintf('eventNumber == %d & isSubord == 0 & acc == 1',find(ismember(exper.eventValues{ses},'name_stim'))) ...
    %sprintf('eventNumber == %d & isSubord == 0 & acc == 0',find(ismember(exper.eventValues{ses},'name_stim'))) ...
    sprintf('eventNumber == %d & isSubord == 1 & acc == 1',find(ismember(exper.eventValues{ses},'name_stim'))) ...
    %sprintf('eventNumber == %d & isSubord == 1 & acc == 0',find(ismember(exper.eventValues{ses},'name_stim'))) ...
    }...
    };
end

%% load in the subject data

% keyboard;
keeptrials = false;
[data_tla,exper] = mm_loadSubjectData(exper,dirs,ana,'tla',keeptrials,'trialinfo');

% overwrite ana.eventValues with the new split events
ana.eventValues = ana.eventValuesSplit;

%% save data

saveDir = dirs.saveDirProc;

% save('/data/projects/curranlab/EBIRD/EEG/Sessions/ftpp/ft_data/data_art_nsClassic/erp/ebird_data_tla_train.mat','data_tla','exper','ana','dirs','files','subjects','-v7.3');
save(fullfile(saveDir,'ebird_data_tla_train.mat'),'data_tla','exper','ana','dirs','files','subjects','-v7.3');

%% load

subDir = '';
dataDir = fullfile('EBIRD','EEG','Sessions','ftpp',subDir);
% Possible locations of the data files (dataroot)
serverDir = fullfile(filesep,'Volumes','curranlab','Data');
serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dreamDir = fullfile(filesep,'data','projects','curranlab');
localDir = fullfile(getenv('HOME'),'data');

% pick the right dataroot
if exist('serverDir','var') && exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist('serverLocalDir','var') && exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
elseif exist('dreamDir','var') && exist(dreamDir,'dir')
  dataroot = dreamDir;
elseif exist('localDir','var') && exist(localDir,'dir')
  dataroot = localDir;
else
  error('Data directory not found.');
end

loadDir = fullfile(dataroot,dataDir,'ft_data/data_art_nsClassic/tla');

load(fullfile(loadDir,'ebird_data_tla_train.mat'));

replaceDataType = {};

[dirs] = mm_checkDirs(dirs,replaceDataType,subDir);

%% decide who to kick out based on trial counts

% Subjects with bad behavior
% exper.badBehSub = {{},{},{}};
% exper.badBehSub = {{'EBIRD021'},{'EBIRD021'},{'EBIRD021'}};

% exper.badBehSub = {repmat({'EBIRD021'},1,length(sesNames))};
exper.badBehSub = repmat({{}},1,length(sesNames));

evToCheck = ana.eventValues;
% evToCheck = repmat({...
%   {{'basic_cor', 'subord_cor'} ...
%   }},1,length(exper.sessions));

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,15,[],'vert', evToCheck);

%% let me know that it's done
emailme = 1;
if emailme
  subject = sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    };
  send_gmail(subject,mail_message);
end
keyboard;

%% Display blink info by condition
% 
% 
% tabLimit = 7;
% eventValues = ana.eventValues;
% printDirection = 'horiz';
% 
% tabchar_sub = '';
% if length(exper.subjects{1}) > tabLimit
%   for i = 1:floor(length(exper.subjects{1}) / tabLimit)
%     tabchar_sub = cat(2,tabchar_sub,'\t');
%   end
% end
% 
% fprintf('\n');
% 
% % print out the blink counts for each subject
% for ses = 1:length(exper.sessions)
%   fprintf('========================\n');
%   fprintf('Session: %s\n',exper.sesStr{ses});
%   fprintf('========================\n');
%   for evVal = 1:length(eventValues{ses})
%     if strcmp(printDirection,'horiz')
%     fprintf('Subject%s%s\n',sprintf(tabchar_sub),sprintf(repmat('\t%s',1,length(eventValues{ses}{evVal})),eventValues{ses}{evVal}{:}));
%     elseif strcmp(printDirection,'vert')
%       fprintf('Subject\n');
%     end
%     for sub = 1:length(exper.subjects)
%       subStr = exper.subjects{sub};
%       for es = 1:length(eventValues{ses}{evVal})
%         if strcmp(printDirection,'horiz')
%           tabchar_ev = '';
%           if length(eventValues{ses}{evVal}{es}) > tabLimit
%             for i = 1:floor(length(eventValues{ses}{evVal}{es}) / tabLimit)
%               tabchar_ev = cat(2,tabchar_ev,'\t');
%             end
%           end
%  
% %           keyboard;
%           
%           subStr = cat(2,subStr,sprintf('\t%d%s',sum(exper.artifacts.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es}).blink{sub}),sprintf(tabchar_ev)));
%         
% %         keyboard;
%         
%         elseif strcmp(printDirection,'vert')
%           if es == 1
%             subStr = sprintf('%s\n',subStr);
%           end
%           subStr = cat(2,subStr,sprintf('\t%s: %d\n',eventValues{ses}{evVal}{es},exper.nTrials.(exper.sesStr{ses}).(eventValues{ses}{evVal}{es}).blink.(sub)));
%         end
%       end
%       fprintf('%s\n',subStr);
%     end
%   end
% end

%% Calculate artifacts info across conditions

%                 blink: {18x1 cell}
%                eyemove: {18x1 cell}
%            manybadchan: {18x1 cell}
%        badneighborchan: {18x1 cell}
%           basic_ftAuto: {18x1 cell}
%            jump_ftAuto: {18x1 cell}
%       threshold_ftAuto: {18x1 cell}
%         basic_ftManual: {18x1 cell}
%          jump_ftManual: {18x1 cell}
%     threshold_ftManual: {18x1 cell}
%        visual_ftManual: {18x1 cell}
%            basic_ftICA: {18x1 cell}
%             jump_ftICA: {18x1 cell}
%        threshold_ftICA: {18x1 cell}
%           visual_ftICA: {18x1 cell}

artType = 'blink';

artCount = zeros(length(exper.subjects), length(exper.sessions));   
trialCount = zeros(length(exper.subjects), length(exper.sessions));  
for ses = 1:length(exper.sessions)
    for evVal = 1:length(ana.eventValues{ses})
         for sub = 1:length(exper.subjects)
%               for es = 1:length(eventValues{ses}{evVal})
                 artCount(sub, ses) =   artCount(sub, ses)+ sum(exper.artifacts.(exper.sesStr{ses}).match_stim.(artType){sub});  
                 trialCount(sub, ses) =   trialCount(sub, ses)+ length(exper.artifacts.(exper.sesStr{ses}).match_stim.(artType){sub});  
%               end
         end
    end
end

artPercents = artCount./trialCount;

%% Test plots to make sure data look ok

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.ylim = [-15 15];
cfg_ft.layout = ft_prepare_layout([],ana);
sub = 1;
ses = 1;
ev = 1;

for sub = 1:length(subjects)
% for sub = 1:4
% for ev = 1:length(ana.eventValues{ses}{1})
  figure
%   ft_multiplotER(cfg_ft,data_tla.(exper.sessions{ses}{1}).(ana.eventValues{ses}{1}{ev}).sub(sub).data); %1 condition
  ft_multiplotER(cfg_ft,data_tla.(exper.sessions{ses}{1}).(ana.eventValues{ses}{1}{ev}).sub(sub).data, data_tla.(exper.sessions{ses}{1}).(ana.eventValues{ses}{1}{ev+1}).sub(sub).data); %2 condition
%   title(strrep(ana.eventValues{1}{i},'_','-'));
end

% cfg_ft = [];
% cfg_ft.channel = {'E20'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% figure
% ft_singleplotER(cfg_ft,data_tla.(ana.eventValues{1}{1}).sub(1).ses(1).data,data_tla.(ana.eventValues{1}{2}).sub(1).ses(1).data);

% cfg_ft = [];
% cfg_ft.channel = {'E20'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% %figure
% %ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data,data_tla.(exper.eventValues{2}).sub(1).ses(1).data,data_tla.(exper.eventValues{3}).sub(1).ses(1).data);
% %legend(strrep(exper.eventValues,'_','-'),'Location','SouthEast');
% figure
% cfg_ft.graphcolor = 'b';
% ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% legend(strrep(exper.eventValues{1},'_','-'),'Location','SouthEast');
% hold on
% plot([cfg_ft.xlim(1) cfg_ft.xlim(2)],[0 0],'k--'); % horizontal
% plot([0 0],[-5 5],'k--'); % vertical
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cfg_ft = [];
% cfg_ft.channel = {'all'};
% cfg_ft.latency = [.5 .8];
% data_topo = ft_timelockanalysis(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% cfg_ft = [];
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% cfg_ft.interactive = 'yes';
% cfg_ft.colormap = 'jet';
% %cfg_ft.colormap = 'hot';
% cfg_ft.colorbar = 'yes';
% cfg_ft.xlim = 'maxmin';
% cfg_ft.layout = ft_prepare_layout([],ana);
% figure
% ft_topoplotER(cfg_ft,data_topo);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% lowpass filter and segment for ERPs

% data_tla_backup = data_tla;

lpfilt = true;

if lpfilt
  sampleRate = exper.sampleRate;
  lpfreq = 40;
  lofiltord = 3;
  lpfilttype = 'but';
end

% cfg_sel = [];
% cfg_sel.latency = [-0.2 1.0];

for ses = 1:length(exper.sesStr)
  for typ = 1:length(ana.eventValues{ses})
    for evVal = 1:length(ana.eventValues{ses}{typ})
      for sub = 1:length(exper.subjects)
        fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{typ}{evVal});
        
        if lpfilt
          % process individual trials
          if isfield(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,'trial')
            for i = 1:size(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial,1)
              data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial(i,:,:) = ft_preproc_lowpassfilter( ...
                squeeze(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.trial(i,:,:)), ...
                sampleRate,lpfreq,lofiltord,lpfilttype);
            end
          end
          
          % process subject average
          if isfield(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data,'avg')
            data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.avg = ft_preproc_lowpassfilter( ...
              data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data.avg, ...
              sampleRate,lpfreq,lofiltord,lpfilttype);
          end
        end
        
        % % select
        % data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data = ft_selectdata(cfg_sel,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).sub(sub).data);
        
      end
      
%       % grand average
%       if exist('ga_tla','var')
%         if lpfilt
%           if isfield(ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}),'avg')
%             ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).avg = ft_preproc_lowpassfilter( ...
%               ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}).avg, ...
%               sampleRate,lpfreq,lofiltord,lpfilttype);
%           end
%         end
%         
%         % select
%         ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = ft_selectdata_new(cfg_sel,ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}));
%       end
    end
  end
end

%% get the grand average

ga_tla = struct;

for ses = 1:length(exper.sesStr)
  % set up strings to put in grand average function
  cfg_ana = [];
  cfg_ana.is_ga = 0;
  %cfg_ana.conditions = ana.eventValues{ses};
  cfg_ana.data_str = 'data_tla';
  % cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);
  %cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
  
  cfg_ft = [];
  cfg_ft.keepindividual = 'no';
  %cfg_ft.parameter = 'trial';
  for typ = 1:length(ana.eventValues{ses})
    cfg_ana.conditions = ana.eventValues{ses}{typ};
    cfg_ana.sub_str = mm_catSubStr_multiSes(cfg_ana,exper,ses);
    
    for evVal = 1:length(ana.eventValues{ses}{typ})
      %tic
      fprintf('Running ft_timelockgrandaverage on %s...',ana.eventValues{ses}{typ}{evVal});
      %ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{evVal}){ses}));
      ga_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{typ}{evVal}) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{ses}{typ}{evVal})));
      fprintf('Done.\n');
      %toc
    end
  end
end

% % turn keeptrial data into average for statistical functions because proper
% % processing of dimord is currently broken
% 
% data_tla_avg = struct;
% 
% cfg = [];
% cfg.keeptrials = 'no';
% 
% for ses = 1:length(exper.sesStr)
%   for evVal = 1:length(ana.eventValues{ses})
%     for sub = 1:length(exper.subjects)
%       fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{ses}{evVal});
%       if isfield(data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}).sub(sub).data,'avg')
%         data_tla_avg.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}).sub(sub).data = ft_timelockanalysis(cfg,data_tla.(exper.sesStr{ses}).(ana.eventValues{ses}{evVal}).sub(sub).data);
%       end
%     end
%   end
% end

%% Multi plots of grand averages

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.ylim = [-15 15];
cfg_ft.layout = ft_prepare_layout([],ana);
ses = 1;
ev = 1;

% for ev = 1:length(ana.eventValues{ses}{1})
  figure
%   ft_multiplotER(cfg_ft,ga_tla.(exper.sessions{ses}{1}).(ana.eventValues{ses}{1}{ev})); %1 condition
  ft_multiplotER(cfg_ft,ga_tla.(exper.sessions{ses}{1}).(ana.eventValues{ses}{1}{ev}), ga_tla.(exper.sessions{ses}{1}).(ana.eventValues{ses}{1}{ev+1})); %2 conditions
%     title(strrep(ana.eventValues{1}{i},'_','-'));
% end

%% plot the conditions - simple

cfg_ft = [];
% cfg_ft.xlim = [0 1.0];
cfg_ft.xlim = [-0.05 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% %cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

% % cfg_plot.rois = {{'LAS'},{'LPS'}};
% cfg_plot.rois = {{'posterior'}};
% % cfg_plot.rois = {{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-8 8; -8 8];
% cfg_plot.legendlocs = {'SouthEast','NorthWest'};

% same as Scott et al. (2008)
cfg_plot.rois = {{'LPI2'},{'RPI2'}};
% cfg_plot.rois = {{'LPI2','RPI2'}};
cfg_plot.ylims = [-2 6.5; -2 6.5];
cfg_plot.legendlocs = {'NorthEast','NorthEast'};
% cfg_plot.legendlocs = {'SouthEast','SouthEast'};

% cfg_ft.xlim = [-0.2 1.0];
% cfg_plot.rois = {{'E70'},{'E83'}};
% cfg_plot.ylims = [-2 8; -2 8];
% cfg_plot.legendlocs = {'NorthEast','NorthEast'};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

%cfg_plot.condByTypeByROI = repmat({{{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}},size(cfg_plot.rois));

% cfg_plot.condByROI = repmat({{'cond1' 'cond2'}},size(cfg_plot.rois));

sesNum = [1, 2, 3, 4, 5, 6];

% cfg_plot.condByROI = repmat(ana.eventValues{sesNum},size(cfg_plot.rois));

cfg_plot.condByROI = repmat({{'basic', 'subord'}},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'basic_cor', 'subord_cor'}},size(cfg_plot.rois));

cfg_ft.graphcolor = 'mrcbgkmrcbgkmrcbgk';
cfg_ft.linestyle = {'-','-','-','-','-','-','--','--','--','--','--','--','-.','-.','-.','-.','-.','-.'};

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  

  %mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
  mm_ft_simpleplotER_multiSes(cfg_ft,cfg_plot,ana,exper,sesNum,ga_tla);
  %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
  
%   stimStr = '';
%   imgCond = '';
%   trainCond = '';
%   imgCond = 'Incongruent';
%   imgCond = 'Gray-AllSF';
%   imgCond = 'HiSF';
%   imgCond = 'LoSF';
  
  %legendLoc = 'NorthOutside';
  %legendLoc = 'NorthEast';
  legend({'T1 basic', 'T1 subord', 'T2 basic', 'T2 subord', 'T3 basic', 'T3 subord', 'T4 basic', 'T4 subord', 'T5 basic', 'T5 subord', 'T6 basic', 'T6 subord'},'Location',cfg_plot.legendlocs{r});
  %legend off
  
  saveFigs = false;
  
  publishfig(gcf,1,14);
  if saveFigs
    figFileName = sprintf('ERP_basic_subord_train');
    print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
  end
end

%% LINES: break pre/post tests into different graphs - training x img cond

% pretest only will collapse across training

% dataMeasure = 'dp';
% dataLabel = 'd''';
% ylimits = [0 4];

% dataMeasure = 'rt_hit';
% dataLabel = 'Response Time: Hits';
% ylimits = [0 5000];

% dataMeasure = 'hr';
% dataLabel = 'Hit Rate';
% ylimits = [0 1];

% dataMeasure = 'c';
% dataLabel = 'Response bias (criterion; c)';
% ylimits = [-0.6 0.6];
% % positive/conservative bias indicates a tendency to say 'new', whereas
% % negative/liberal bias indicates a tendency to say 'old'

% dataMeasure = 'Br';
% dataLabel = 'Response bias index (Br)';
% ylimits = [0 1];

% dataMeasure = 'Pr';
% dataLabel = 'Discrimination index (Pr)';
% ylimits = [0 1];

sessions = {'session_1', 'session_8', 'session_9'};
% sessions_str = {'Pretest', '1-Day Posttest','Delay'};
sessions_str = {'Pretest', '1-Day Posttest','1-Week Posttest'};
% sessions = {'session_8', 'session_9'};
% % sessions_str = {'Pretest', '1-Day Posttest','Delay'};
% sessions_str = {'1-Day Posttest','1-Week Posttest'};

% if collapsePhases
phases = {'match'};
% else
%   phases = {'match_1'};
% end
% training = {'trained','untrained'};

naming = {'basic','subord'};
naming_str = {'Basic','Subord'};
% naming = {'subord'};
% naming_str = {'Subord'};

% training = {'TT','UU','TU','UT'};
% trainStr = {'Trained','Untrained','TU','UT'};

% training = {'TT','UU'};
% trainingStr = {'Trained', 'Untrained'};

stimNum = {'1','2'};
% stimNum = {'1'};
% stimNum = {'2'};
% stimNum_str = {'Stim1','Stim2'};
stimNum_str = {'1','2'};

% % separate hemi
hemis = {'LPI2','RPI2'};
hemis_str = {'Left','Right'};

% % collapse hemi
% hemis = {{'LPI2','RPI2'}};
% hemis_str = {'LR'};

imgConds = {'norm','color','g'};
imgConds_str = {'Congruent','Incongruent','Gray'};
trainConds = {'un', 'tr'};
trainConds_str = {'Un', 'Tr'};
groupname = 'Color';

% imgConds = {'g','g_hi8','g_lo8'};
% % imgConds_str = {'Gray','Hi8','Lo8'};
% imgConds_str = {'AllSF','HiSF','LoSF'};
% groupname = 'SpatialFreq';

allBadSub = logical(sum(exper.badSub,2));

cfg = [];
% cfg.latency = [0.155 0.211]; % N170 (Scott)
% cfg.latency = [0.152 0.212]; % N170 (initial guess)
cfg.latency = [0.156 0.208]; % N170 (use this)
component_str = 'N170';

% cfg.latency = [0.234 0.334]; % N250
% component_str = 'N250';

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},naming{1},imgConds{1},trainConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},naming{1},imgConds{1},trainConds{1})).sub(1).data.time,cfg.latency(2));
% tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},naming{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
% tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},naming{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

%anovaData = [];
%theseData = [];

data_sub = nan(length(sessions),length(imgConds),length(naming),length(stimNum),length(trainConds), length(hemis),length(subjects));
% data_sub = nan(length(sessions),length(imgConds),length(naming),length(stimNum),length(hemis),length(subjects));

for ses = 1:length(sessions)
  for im = 1:length(imgConds)
    
    for n = 1:length(naming)
      for sn = 1:length(stimNum)
        for tn = 1:length(trainConds)

        condition = sprintf('stim%s_%s_%s_%s',stimNum{sn},naming{n},imgConds{im},trainConds{tn});
        
            for h = 1:length(hemis)
              for sub = 1:length(exper.subjects)
                if ~allBadSub(sub)
                  cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis{h})});
                  %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
                  %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));

                  % % collapse hemis
                  % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});

                  elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
                  data_sub(ses,im,n,sn,tn,h,sub) = mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2);

               end
            end
          end
        end
      end
    end
    %anovaData = cat(1,anovaData,theseData);
  end
end

% for s = 1:length(sessions)
%   for p = 1:length(phases)
%     for i = 1:length(imgConds)
%       
%       %collapseData = [];
%       
%       for t = 1:length(training)
%         for n = 1:length(naming)
%           data(t,i,s,p,:) = results.(sessions{s}).(phases{p}).(training{t}).(imgConds{i}).(naming{n}).(dataMeasure);
%         end
%       end
%       
%     end
%   end
% end

% % stats
% [h, p, ci, stats] = ttest(squeeze(data.(dataMeasure).(naming{1})(:,1,1,1,1)),squeeze(data.(dataMeasure).(naming{1})(:,2,1,1,1)));

%% Make line plots - separate session, hemi, stimNum, trainConds

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 7];
end

cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};

if strcmp(groupname,'Color')
  cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
elseif strcmp(groupname,'SpatialFreq')
  cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    for sn = 1:length(stimNum)
      for tn = 1:length(trainConds)

          data_mean = [];
          data_sem = [];
          %data_mean = nan(length(naming),length(imgConds));
          %data_sem = nan(length(naming),length(imgConds));

          thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s, %s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn},trainConds_str{tn});
          figFileName = sprintf('%s_%s_%s_%s_s%s_%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn},trainConds_str{tn});

          cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
          cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));

          poCount = 0;
          for im = 1:length(imgConds)
            for n = 1:length(naming)
                poCount = poCount + 1;
        
              condition = sprintf('%s_%s',imgConds{im},naming{n})
              cfg_plot.plot_order{poCount} = condition;
              %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
              cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));

              data_mean.(condition) = nanmean(data_sub(s,im,n,sn,tn,h,:),7);
%               keyboard
              data_sem.(condition) = nanstd(data_sub(s,im,n,sn,tn,h,:),0,7) ./ sqrt(sum(~allBadSub));
            end
          end
          
           hpm = nan(1,length(cfg_plot.plot_order));

          figure
          % plot the lines
          eval(sprintf('plot([%s],cfg_plot.linespec{1},''LineWidth'',cfg_plot.linewidth);',sprintf(repmat('data_mean.%s ',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:})));
  
          hold on
          
   
          for c = 1:length(cfg_plot.plot_order)
            % errorbars
            heb = errorbar(c,mean(data_mean.(cfg_plot.plot_order{c}),1),data_sem.(cfg_plot.plot_order{c}),cfg_plot.linespec{c},'LineWidth',cfg_plot.errwidth);
            % remove errorbar ends
            
            if cfg_plot.removeErrBarEnds
              chil = get(heb,'Children');
              xdata = get(chil(2),'XData');
              ydata = get(chil(2),'YData');
              xdata(cfg_plot.errBarEndMarkerInd) = NaN;
              ydata(cfg_plot.errBarEndMarkerInd) = NaN;
              set(chil(2),'XData',xdata);
              set(chil(2),'YData',ydata);
              set(heb,'Children',chil);
            end
            % plot the markers
            hpm(c) = plot(c,mean(data_mean.(cfg_plot.plot_order{c}),1),cfg_plot.linespec{c},'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor{c});

            text(c,  cfg_plot.ylim(1)-.2, cfg_plot.rename_conditions{c}, 'Rotation', 15, 'FontSize', 18, 'HorizontalAlignment', 'center');
          end
          if cfg_plot.plotLegend
            legend(hpm(1:length(cfg_plot.legendtext)),cfg_plot.legendtext);
          end

          cfg_plot.xlim = [0.5 (length(cfg_plot.rename_conditions) + .5)];

          plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal

          hold off

          set(gcf,'Name',sprintf('%s',thisTitle))

          title(thisTitle);

          % make it look good
          axis([cfg_plot.xlim(1) cfg_plot.xlim(2) cfg_plot.ylim(1) cfg_plot.ylim(2)])
          xlabel(cfg_plot.xlabel);
          ylabel(cfg_plot.ylabel);
          set(gca,'XTick',(1:length(cfg_plot.rename_conditions)))
          set(gca,'XTickLabel',[])
          %set(gca,'XTickLabel',strrep(cfg_plot.rename_conditions,'_',''))
          set(gca,'YTick',(cfg_plot.ylim(1):.5:cfg_plot.ylim(2)))
          axis square
          publishfig(gcf,0);

          if saveFigs
            print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
          end
      end
    end
  end
end

%% bar plots of subordinate - collapse stim order

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
  cfg_plot.ylim = [0 3];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 6];
end

% cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};
% 
% if strcmp(groupname,'Color')
%   cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
% elseif strcmp(groupname,'SpatialFreq')
%   cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
% end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

ylimits = cfg_plot.ylim;

cfg_plot.axisxy = false;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    %for sn = 1:length(stimNum)
      
      %thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      %figFileName = sprintf('%s_%s_%s_%s_s%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      thisTitle = sprintf('%s: %s: %s, %s hemi',groupname,component_str,sessions_str{s},hemis_str{h});
      figFileName = sprintf('bar_%s_%s_%s_%s',groupname,component_str,sessions_str{s},hemis_str{h});
      
      cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
      cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));
      
      data_mean = [];
      data_sem = [];
      %data_mean = nan(length(naming),length(imgConds));
      %data_sem = nan(length(naming),length(imgConds));
      
      poCount = 0;
      for im = 1:length(imgConds)
        for n = 1:length(naming)
          poCount = poCount + 1;
          
          condition = sprintf('%s_%s',imgConds{im},naming{n});
          cfg_plot.plot_order{poCount} = condition;
          %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
          cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));
          
          %data_mean.(condition) = nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6);
          %data_sem.(condition) = nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub));
          
          data_mean = cat(2,data_mean,nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6));
          data_sem = cat(2,data_sem,nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub)));
        end
      end
      
      bw_legend = imgConds_str;
      
      bw_title = thisTitle;
      %bw_title = sprintf('%s: %s, %s',groupname,sesStr{s},naming{n});
      %bw_title = sprintf('%s%s: %s: %s',upper(naming{n}(1)),naming{n}(2:end),strrep(phases{p},'_','\_'),strrep(imgConds{i},'_','\_'));
      %bw_groupnames = {'Pretest', 'Posttest', 'Delay'};
      if strcmp(sessions{s},'session_1')
        bw_groupnames = [];
        axis_x = 1.5;
      else
        %bw_groupnames = trainingStr;
        bw_groupnames = [];
        
        %axis_x = length(training) + 1.5;
        %axis_x = length(training) + 0.5;
        axis_x = 1.5;
      end
      %bw_xlabel = 'Test day';
      bw_xlabel = cfg_plot.xlabel;
      bw_ylabel = cfg_plot.ylabel;
      %if exist('linspecer','file')
      %  bw_colormap = 'linspecer';
      %else
      %  bw_colormap = 'gray';
      %end
      
      bw_colormap = flipud(jet);
      bw_data = data_mean;
      bw_errors = data_sem;
      bw_width = 0.75;
      hbp = barweb(bw_data,bw_errors,bw_width,bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend,[],'plot');
      set(hbp.legend,'Location','NorthEast');
      axis([0.5 axis_x ylimits(1) ylimits(2)]);
      
      %if length(sessions) == 1 && strcmp(sessions{s},'pretest')
      %  xlabel('Collapsed');
      %end
      publishfig(gcf,0);
      
      if cfg_plot.axisxy
        %axis xy;
        set(gca,'YDir','reverse');
        set(hbp.legend,'Location','SouthWest');
        figFileName = sprintf('%s_flip',figFileName);
      end
      
      figFileName = strrep(strrep(figFileName,' ','_'),'-','_');
      
      if saveFigs
        print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
        %print(gcf,figFormat,figRes,fullfile(figsDir,sprintf('prepost_trainUn_%s_%s_%s_%s_%s',phases{p},dataMeasure,sessions{s},naming{n},groupname)));
      end
      
  end
end

%% bar plots of subordinate - separate stim order

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
  cfg_plot.ylim = [0 3];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 6];
end

% cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};
% 
% if strcmp(groupname,'Color')
%   cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
% elseif strcmp(groupname,'SpatialFreq')
%   cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
% end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

ylimits = cfg_plot.ylim;

cfg_plot.axisxy = true;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    for sn = 1:length(stimNum)
      
      thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      figFileName = sprintf('bar_%s_%s_%s_%s_s%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      %thisTitle = sprintf('%s: %s: %s, %s hemi',groupname,component_str,sessions_str{s},hemis_str{h});
      %figFileName = sprintf('bar_%s_%s_%s_%s',groupname,component_str,sessions_str{s},hemis_str{h});
      
      cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
      cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));
      
      data_mean = [];
      data_sem = [];
      %data_mean = nan(length(naming),length(imgConds));
      %data_sem = nan(length(naming),length(imgConds));
      
      poCount = 0;
      for im = 1:length(imgConds)
        for n = 1:length(naming)
          poCount = poCount + 1;
          
          condition = sprintf('%s_%s',imgConds{im},naming{n});
          cfg_plot.plot_order{poCount} = condition;
          %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
          cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));
          
          %data_mean.(condition) = nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6);
          %data_sem.(condition) = nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub));
          
          data_mean = cat(2,data_mean,nanmean(data_sub(s,im,n,sn,h,:),6));
          data_sem = cat(2,data_sem,nanstd(data_sub(s,im,n,sn,h,:),0,6) ./ sqrt(sum(~allBadSub)));
        end
      end
      
      bw_legend = imgConds_str;
      
      bw_title = thisTitle;
      %bw_title = sprintf('%s: %s, %s',groupname,sesStr{s},naming{n});
      %bw_title = sprintf('%s%s: %s: %s',upper(naming{n}(1)),naming{n}(2:end),strrep(phases{p},'_','\_'),strrep(imgConds{i},'_','\_'));
      %bw_groupnames = {'Pretest', 'Posttest', 'Delay'};
      if strcmp(sessions{s},'session_1')
        bw_groupnames = [];
        axis_x = 1.5;
      else
        %bw_groupnames = trainingStr;
        bw_groupnames = [];
        
        %axis_x = length(training) + 1.5;
        %axis_x = length(training) + 0.5;
        axis_x = 1.5;
      end
      %bw_xlabel = 'Test day';
      bw_xlabel = cfg_plot.xlabel;
      bw_ylabel = cfg_plot.ylabel;
      %if exist('linspecer','file')
      %  bw_colormap = 'linspecer';
      %else
      %  bw_colormap = 'gray';
      %end
      
      bw_colormap = flipud(jet);
      bw_data = data_mean;
      bw_errors = data_sem;
      bw_width = 0.75;
      hbp = barweb(bw_data,bw_errors,bw_width,bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend,[],'plot');
      set(hbp.legend,'Location','NorthEast');
      axis([0.5 axis_x ylimits(1) ylimits(2)]);
      
      %if length(sessions) == 1 && strcmp(sessions{s},'pretest')
      %  xlabel('Collapsed');
      %end
      publishfig(gcf,0);
      
      if cfg_plot.axisxy
        %axis xy;
        set(gca,'YDir','reverse');
        set(hbp.legend,'Location','SouthWest');
        figFileName = sprintf('%s_flip',figFileName);
      end
      
      figFileName = strrep(strrep(figFileName,' ','_'),'-','_');
      
      if saveFigs
        print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
        %print(gcf,figFormat,figRes,fullfile(figsDir,sprintf('prepost_trainUn_%s_%s_%s_%s_%s',phases{p},dataMeasure,sessions{s},naming{n},groupname)));
      end
      
    end
  end
end

%% Make line plots - separate session (collapse hemi, stimNum)

% make some plots

saveFigs = true;

cfg_plot = [];
cfg_plot.linewidth = 2;
% cfg_plot.marksize = 10;
cfg_plot.marksize = 15;
cfg_plot.errwidth = 1;
cfg_plot.errBarEndMarkerInd = [4 5 7 8];
cfg_plot.removeErrBarEnds = 1;

cfg_plot.xlabel = 'Conditions';
cfg_plot.xlabel = '';
cfg_plot.ylabel = 'Voltage (\muV)';

if strcmp(component_str,'N170')
  cfg_plot.ylim = [-0.5 3.5];
elseif strcmp(component_str,'N250')
  cfg_plot.ylim = [3 7];
end

cfg_plot.linespec = {'ko','ks','ko','ks','ko','ks'};

if strcmp(groupname,'Color')
  cfg_plot.markcolor = {[1 1 1],[1 1 1],[0 0 0],[0 0 0],[0.6 0.6 0.6],[0.6 0.6 0.6]};
elseif strcmp(groupname,'SpatialFreq')
  cfg_plot.markcolor = {[0.6 0.6 0.6],[0.6 0.6 0.6],[1 1 1],[1 1 1],[0 0 0],[0 0 0]};
end

cfg_plot.plotLegend = true;

cfg_plot.legendtext = naming_str;

% for p = 1:length(phases)
for s = 1:length(sessions)
  for h = 1:length(hemis)
    %for sn = 1:length(stimNum)
      
      data_mean = [];
      data_sem = [];
      %data_mean = nan(length(naming),length(imgConds));
      %data_sem = nan(length(naming),length(imgConds));
      
      %thisTitle = sprintf('%s: %s: %s, %s hemi, Stim%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      %figFileName = sprintf('%s_%s_%s_%s_s%s',groupname,component_str,sessions_str{s},hemis_str{h},stimNum_str{sn});
      thisTitle = sprintf('%s: %s: %s, %s hemi',groupname,component_str,sessions_str{s},hemis_str{h});
      figFileName = sprintf('%s_%s_%s_%s',groupname,component_str,sessions_str{s},hemis_str{h});
      
      cfg_plot.plot_order = cell(1,prod([length(naming),length(imgConds)]));
      cfg_plot.rename_conditions = cell(1,prod([length(naming),length(imgConds)]));
      
      poCount = 0;
      for im = 1:length(imgConds)
        for n = 1:length(naming)
          poCount = poCount + 1;
          
          condition = sprintf('%s_%s',imgConds{im},naming{n});
          cfg_plot.plot_order{poCount} = condition;
          %cfg_plot.rename_conditions{poCount} = sprintf('%s %s',imgConds_str{im},naming_str{n});
          cfg_plot.rename_conditions{poCount} = sprintf('%s-%s',imgConds_str{im},naming_str{n}(1));
          
          data_mean.(condition) = nanmean(nanmean(data_sub(s,im,n,:,h,:),4),6);
          data_sem.(condition) = nanstd(nanmean(data_sub(s,im,n,:,h,:),4),0,6) ./ sqrt(sum(~allBadSub));
        end
      end
      
      hpm = nan(1,length(cfg_plot.plot_order));
      
      figure
      % plot the lines
      eval(sprintf('plot([%s],cfg_plot.linespec{1},''LineWidth'',cfg_plot.linewidth);',sprintf(repmat('data_mean.%s ',1,length(cfg_plot.plot_order)),cfg_plot.plot_order{:})));
      hold on
      for c = 1:length(cfg_plot.plot_order)
        % errorbars
        heb = errorbar(c,mean(data_mean.(cfg_plot.plot_order{c}),1),data_sem.(cfg_plot.plot_order{c}),cfg_plot.linespec{c},'LineWidth',cfg_plot.errwidth);
        % remove errorbar ends
        if cfg_plot.removeErrBarEnds
          chil = get(heb,'Children');
          xdata = get(chil(2),'XData');
          ydata = get(chil(2),'YData');
          xdata(cfg_plot.errBarEndMarkerInd) = NaN;
          ydata(cfg_plot.errBarEndMarkerInd) = NaN;
          set(chil(2),'XData',xdata);
          set(chil(2),'YData',ydata);
          set(heb,'Children',chil);
        end
        % plot the markers
        hpm(c) = plot(c,mean(data_mean.(cfg_plot.plot_order{c}),1),cfg_plot.linespec{c},'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor{c});
        
        text(c,  cfg_plot.ylim(1)-.2, cfg_plot.rename_conditions{c}, 'Rotation', 15, 'FontSize', 18, 'HorizontalAlignment', 'center');
      end
      if cfg_plot.plotLegend
        legend(hpm(1:length(cfg_plot.legendtext)),cfg_plot.legendtext);
      end
      
      cfg_plot.xlim = [0.5 (length(cfg_plot.rename_conditions) + .5)];
      
      plot([cfg_plot.xlim(1) cfg_plot.xlim(2)],[0 0],'k--'); % horizontal
      
      hold off
      
      set(gcf,'Name',sprintf('%s',thisTitle))
      
      title(thisTitle);
      
      % make it look good
      axis([cfg_plot.xlim(1) cfg_plot.xlim(2) cfg_plot.ylim(1) cfg_plot.ylim(2)])
      xlabel(cfg_plot.xlabel);
      ylabel(cfg_plot.ylabel);
      set(gca,'XTick',(1:length(cfg_plot.rename_conditions)))
      set(gca,'XTickLabel',[])
      %set(gca,'XTickLabel',strrep(cfg_plot.rename_conditions,'_',''))
      set(gca,'YTick',(cfg_plot.ylim(1):.5:cfg_plot.ylim(2)))
      axis square
      publishfig(gcf,0);
      
      if saveFigs
        print(gcf,'-dpng',fullfile(dirs.saveDirFigs,figFileName));
      end
      
    %end
  end
end

%% find the peak of N170
% 
% cfg = [];
% cfg.conditions = cellflat(ana.eventValues{1});
% cfg.data_str = 'ga_tla';
% cfg.is_ga = true;
% cfg.excludeBadSub = true;
% 
% sesNum = 1;
% 
% ana_str = mm_catSubStr_multiSes2(cfg,exper,sesNum);
% 
% cfg_ft = [];
% % %cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'posterior'})});
% % %cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS', 'RPS', 'LPI', 'RPI', 'LPS2', 'RPS2', 'LPI2', 'RPI2', 'PS', 'PI'})});
% cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPI', 'RPI', 'LPI2', 'RPI2', 'PI'})});
% cfg_ft.channel{end+1} = 'E81';
% 
% cfg_ft.latency = [0 0.5];
% 
% 
% % cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'RPI2'})});
% % cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'RPI3'})});
% 
% ga_allCond = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',ana_str));
% 
% figure
% ft_singleplotER(cfg_ft,ga_allCond);

%% find peak N170
% 
% cfg = [];
% % cfg.latency = [0.155 0.211];
% cfg.latency = [0.156 0.208];
% 
% ga_n170 = ft_selectdata_new(cfg,ga_allCond);
% 
% [y,i] = sort(mean(ga_n170.avg,2));
% 
% ga_n170.label(i)'
% y'

%% find peak N250
% 
% cfg = [];
% cfg.latency = [0.23 0.33];
% 
% ga_n250 = ft_selectdata_new(cfg,ga_allCond);
% 
% [y,i] = sort(mean(ga_n250.avg,2));
% 
% ga_n250.label(i)'
% y'

%% plot peaks
% 
% cfg_ft = [];
% cfg_ft.showlabels = 'yes';
% % cfg_ft.interactive = 'yes';
% cfg_ft.interactive = 'no';
% cfg_ft.showoutline = 'yes';
% cfg_ft.fontsize = 9;
% cfg_ft.ylim = [-2 5];
% cfg_ft.layout = ft_prepare_layout([],ana);
% 
% figure
% ft_multiplotER(cfg_ft,ga_allCond);

%% plot 90
% 
% cfg_ft = [];
% cfg_ft.channel = {'E90'};
% 
% figure
% ft_singleplotER(cfg_ft,ga_allCond);

%% peak to peak N170
% 
% pospeak = [0.1 0.15];
% 
% negpeak = [0.15 0.2];
% 
% cfg = [];
% cfg.latency = pospeak;
% ga_pos = ft_selectdata_new(cfg,ga_allCond);
% 
% 
% cfg = [];
% cfg.latency = negpeak;
% ga_neg = ft_selectdata_new(cfg,ga_allCond);
% 
% pos_max = max(ga_pos.avg,[],2);
% neg_min = min(ga_neg.avg,[],2);
% 
% posnegdiff = pos_max - abs(neg_min);
% 
% [y,i] = sort(posnegdiff);
% 
% ga_allCond.label(i)'
% y'
%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);


%% find peak = N710
cfg = [];

% images
% cfg.conditions = cellflat(ana.eventValues{1}{1});
% words
% cfg.conditions = cellflat(ana.eventValues{1}{2});

% % all together
% cfg.conditions = {'word_onePres','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2','word_RgH_rc_spac_p1','word_RgH_fo_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_mass_p1'};
% cfg.conditions = {'word_onePres','word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};
% cfg.conditions = {'word_RgH_rc_spac_p2','word_RgH_fo_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_mass_p2'};

cfg.conditions = {'stim1_basic_color_un', 'stim1_basic_g_un', 'stim1_basic_g_hi8_un', 'stim1_basic_g_lo8_un', 'stim1_basic_norm_un', ...
  'stim1_subord_color_un', 'stim1_subord_g_un', 'stim1_subord_g_hi8_un', 'stim1_subord_g_lo8_un', 'stim1_subord_norm_un', ...
  'stim2_basic_color_un', 'stim2_basic_g_un', 'stim2_basic_g_hi8_un', 'stim2_basic_g_lo8_un', 'stim2_basic_norm_un', ...
  'stim2_subord_color_un', 'stim2_subord_g_un', 'stim2_subord_g_hi8_un', 'stim2_subord_g_lo8_un', 'stim2_subord_norm_un', ...
  'stim1_basic_color_tr', 'stim1_basic_g_tr', 'stim1_basic_g_hi8_tr', 'stim1_basic_g_lo8_tr', 'stim1_basic_norm_tr', ...
  'stim1_subord_color_tr', 'stim1_subord_g_tr', 'stim1_subord_g_hi8_tr', 'stim1_subord_g_lo8_tr', 'stim1_subord_norm_tr', ...
  'stim2_basic_color_tr', 'stim2_basic_g_tr', 'stim2_basic_g_hi8_tr', 'stim2_basic_g_lo8_tr', 'stim2_basic_norm_tr', ...
  'stim2_subord_color_tr', 'stim2_subord_g_tr', 'stim2_subord_g_hi8_tr', 'stim2_subord_g_lo8_tr', 'stim2_subord_norm_tr' ...
    };

% % ses 1
% cfg.conditions = {'stim1_basic_color_un', 'stim1_basic_g_un', 'stim1_basic_g_hi8_un', 'stim1_basic_g_lo8_un', 'stim1_basic_norm_un', ...
%   'stim1_subord_color_un', 'stim1_subord_g_un', 'stim1_subord_g_hi8_un', 'stim1_subord_g_lo8_un', 'stim1_subord_norm_un', ...
%   'stim2_basic_color_un', 'stim2_basic_g_un', 'stim2_basic_g_hi8_un', 'stim2_basic_g_lo8_un', 'stim2_basic_norm_un', ...
%   'stim2_subord_color_un', 'stim2_subord_g_un', 'stim2_subord_g_hi8_un', 'stim2_subord_g_lo8_un', 'stim2_subord_norm_un', ...
%       };

% % N170
cfg.order = 'ascend'; % ascend = negative peaks first
cfg.roi = {'posterior'};
% cfg.roi = {'Cz'};
% cfg.latency = [0.2 0.6];
n170peak = 0.184;
% % % cfg.latency = [n400Peak-0.05 n400Peak+0.05]; % N400 - around GA peak (space+mass) +/- 50
cfg.latency = [n170peak-0.025 n170peak+0.025]; % N400 - around GA peak (space+mass) +/- 50
% cfg.latency = [n400Peak-0.1 n400Peak+0.1]; % N400 - around GA peak (space+mass) +/- 100

cfg.datadim = 'elec';
% cfg.datadim = 'time';
% % cfg.roi = {'Cz'};
% % cfg.roi = {'LPI3','RPI3'};
% cfg.roi = {'Pz'};
% % cfg.roi = {'PS2'};
% % cfg.roi = {'RPI3'};
% % cfg.roi = {'E84'}; % center of RPI3
% % cfg.roi = {'RPS2'};
% % cfg.roi = {'E85'}; % center of RPS2
% % cfg.roi = {'LPS2'};
% % cfg.latency = [0 1.0];
% cfg.latency = [0.2 0.9];

cfg.is_ga = false;
% cfg.is_ga = true;
% cfg.outputSubjects = true;
cfg.outputSubjects = false;
cfg.sesNum = 1; % How to do more than 1 session????

cfg.plotit = true;
cfg.voltlim = [-3 3]; % LPC
% cfg.voltlim = [-2 2]; % N400
% cfg.voltlim = [-1 5];

% % only for datadim='elec' and datadim='peak2peak'
% cfg.plottype = 'topo';
cfg.plottype = 'multi';
cfg_plot.interactive = 'yes';
cfg_plot.xlim=[-0.2 1.0];

% % % only for datadim='peak2peak'
cfg.datadim = 'peak2peak';
cfg.order = 'descend'; % descend = positive peaks first
% cfg.roi = {'posterior'};
cfg.pospeak = [0.09 0.16];
cfg.negpeak = [0.16 0.21];

% peakInfo = mm_findPeak(cfg,ana,exper,ga_tla);
% [peakInfo, ga_allCond] = mm_findPeak_tc(cfg,ana,exper,data_tla);
peakInfo = mm_findPeak(cfg,ana,exper,data_tla,cfg_plot);

%% save peak data

% % save('/Volumes/curranlab/Data/EBIRD/EEG/Sessions/ftpp/new_art500/ft_data/data_art_nsClassic_ftAuto/erp/ebird_data_tla.mat','data_tla','exper','ana','dirs','files','subjects','-v7.3');
% 
% save('peak_data.mat','peakInfo','ga_allCond','cfg','ana','exper','-v7.3');

%% BIG ANOVA

% test day (3) x image condition (5) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

cfg = [];

cfg.latency = [0.156 0.208]; % N170 (use this)
component_str = 'N170';

cfg.latency = [0.234 0.334]; % N250
component_str = 'N250';

% sessions = {'session_1', 'session_8', 'session_9'};
% lev_sessions = {'pre','post','delay'};

% sessions = {'session_1', 'session_8'};
% lev_sessions = {'pre','post'};

% sessions = {'session_1', 'session_9'};
% lev_sessions = {'pre','delay'};

sessions = {'session_1','session_8', 'session_9'};
lev_sessions = {'pre','post','delay'};

% groupname = 'Color';
% imgConds = {'norm','color','g'};
% lev_imgConds = {'cong','incon','gray'};

groupname = 'SpatialFreq';
imgConds = {'g','g_hi8','g_lo8'};
lev_imgConds = {'gray','g_hi8','g_lo8'};

naming = {'basic', 'subord'};

stimNum = {'1' '2'};
lev_stimNum = {'s1' 's2'};
% stimNum = {'2'};
% lev_stimNum = {'s2'};

trainConds = {'un', 'tr'};
lev_trainConds = {'un', 'tr'};

hemis = {'LPI2','RPI2'};
lev_hemis = {'left','right'};

allBadSub = logical(sum(exper.badSub,2));

fprintf('ANOVA: %s %s\n',groupname,component_str);

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},naming{1},imgConds{1},trainConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},naming{1},imgConds{1},trainConds{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];
lnDone = false;
vnDone = false;
for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    lnCount = 0;
    vnCount = 0;
    
    
    for ses = 1:length(sessions)
      for im = 1:length(imgConds)
        
        for sn = 1:length(stimNum)
          for n = 1:length(naming)
            for  tn = 1:length(trainConds) 
                
                condition = sprintf('stim%s_%s_%s_%s',stimNum{sn},naming{n},imgConds{im},trainConds{tn});

                for h = 1:length(hemis)
                    
                  vnCount = vnCount + 1; 
                  if ~vnDone
                        variableNames{vnCount} = sprintf('Y%d',vnCount);
                  end
          
                          if ~lnDone
                            lnCount = lnCount + 1;
                            levelNames{lnCount,1} = sessions{ses};
                            levelNames{lnCount,2} = imgConds{im};
                            levelNames{lnCount,3} = stimNum{sn};
                            levelNames{lnCount,4} = naming{n};
                            levelNames{lnCount,5} = trainConds{tn};
                            levelNames{lnCount,6} = hemis{h};
                          end
 
                  cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis(h))});
 
                  % % collapse hemis
                  % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});

                  elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
%                   theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
               
                  anovaData(sub,vnCount) = mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2);
% keyboard

                end
            end
          end
        end
      end
    end
%     anovaData = cat(1,anovaData,theseData);
  end
  
  nVariables = vnCount;
  vnDone = true;
  lnDone = true;
  
end
% factorNames = {'sessions', 'imgConds', 'stimNum', 'naming', 'trainConds', 'hemis'};
factorNames = {'ses', 'img', 'stm', 'nam', 'trn', 'hem'};
anovaTable = array2table(anovaData,'VariableNames',variableNames);
within = cell2table(levelNames,'VariableNames',factorNames);
rm = fitrm(anovaTable,'Y1-Y144~1','WithinDesign',within);
[ranovatbl,A,C,D] = ranova(rm, 'WithinModel','ses*img*stm*nam*trn*hem')


% varnames = {'TestDay', 'ImgCond','StimNum','Basic/Subord','Hemisphere'};
% 
% levelnames = {lev_sessions lev_imgConds lev_stimNum training lev_hemis};
% 
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% % O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);
% 
% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord'};
% levelnames = {sessions imgConds stimNum training};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training)], varnames);

%% BIG ANOVA - stim1 only, collapse hemi

% test day (3) x image condition (5) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

cfg = [];

cfg.latency = [0.156 0.208]; % N170 (use this)
component_str = 'N170';

cfg.latency = [0.234 0.334]; % N250
component_str = 'N250';

sessions = {'session_1', 'session_8', 'session_9'};
lev_sessions = {'pre','post','delay'};

% sessions = {'session_1', 'session_8'};
% lev_sessions = {'pre','post'};

% sessions = {'session_1', 'session_9'};
% lev_sessions = {'pre','delay'};

% sessions = {'session_8', 'session_9'};
% lev_sessions = {'post','delay'};

groupname = 'Color';
imgConds = {'norm','color','g'};
lev_imgConds = {'cong','incon','gray'};

groupname = 'SpatialFreq';
imgConds = {'g','g_hi8','g_lo8'};
lev_imgConds = {'gray','g_hi8','g_lo8'};

training = {'basic', 'subord'};

% stimNum = {'1' '2'};
% lev_stimNum = {'s1' 's2'};
stimNum = {'1'};
lev_stimNum = {'s1'};

% hemis = {'LPI2','RPI2'};
% lev_hemis = {'left','right'};

hemis = {{'LPI2','RPI2'}};
lev_hemis = {'LR'};

allBadSub = logical(sum(exper.badSub,2));

fprintf('ANOVA: %s %s\n',groupname,component_str);

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},naming{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},naming{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
      for im = 1:length(imgConds)
        
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%s_%s_%s',stimNum{sn},training{tr},imgConds{im});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis{h})});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

% varnames = {'TestDay', 'ImgCond','StimNum','Basic/Subord','Hemisphere'};
% 
% levelnames = {lev_sessions lev_imgConds lev_stimNum training lev_hemis};
% 
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% % O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);

varnames = {'testDay', 'ImgCond','Basic/Subord'};
levelnames = {lev_sessions lev_imgConds training};
O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training)], varnames,[],[],[],[],[],[],levelnames);

%% normal condition only ANOVA

% test day (3) x image condition (5) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

sessions = {'session_1', 'session_8', 'session_9'};
% sessions = {'session_1', 'session_8'};
% sessions = {'session_1', 'session_9'};
% sessions = {'session_8', 'session_9'};

% imgConds = {'norm'};
% imgConds = {'g'};
imgConds = {'color'};
% imgConds = {'norm','g','color'};
% imgConds = {'g','g_hi8','g_lo8'};

training = {'basic', 'subord'};

stimNum = {'1' '2'};
% stimNum = {'2'};

hemis = {'LPI2','RPI2'};

allBadSub = logical(sum(exper.badSub,2));

cfg = [];
% cfg.latency = [0.155 0.211]; % N170 (Scott)
% cfg.latency = [0.152 0.212]; % N170 (initial guess)
cfg.latency = [0.156 0.208]; % N170 (use this)

% cfg.latency = [0.234 0.334]; % N250

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
      for im = 1:length(imgConds)
        
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%s_%s_%s',stimNum{sn},training{tr},imgConds{im});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis(h))});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'testDay','StimNum','Basic/Subord','Hemisphere'};

lev_sessions = {'pre','post','delay'};
% lev_sessions = {'pre','post'};
% lev_sessions = {'pre','delay'};
% lev_sessions = {'post','delay'};

% lev_imgConds = {'cong','gray','incon'};
% lev_imgConds = {'gray','g_hi8','g_lo8'};

lev_stimNum = {'s1' 's2'};
% lev_stimNum = {'s2'};

lev_hemis = {'left','right'};

levelnames = {lev_sessions lev_stimNum training lev_hemis};

O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord'};
% levelnames = {sessions imgConds stimNum training};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training)], varnames);

%% Pretest ANOVA

% test day (3) x image condition (5) x stim1/stim2 (2) x training basic/subord (2) x hemi (2)

sessions = {'session_1'};

imgConds = {'norm','g','color'};
% imgConds = {'g','g_hi8','g_lo8'};

training = {'basic', 'subord'};

stimNum = {'1' '2'};
% stimNum = {'2'};

hemis = {'LPI2','RPI2'};

allBadSub = logical(sum(exper.badSub,2));

cfg = [];
% cfg.latency = [0.155 0.211]; % N170 (Scott)
% cfg.latency = [0.152 0.212]; % N170 (initial guess)
cfg.latency = [0.156 0.208]; % N170 (use this)

cfg.latency = [0.234 0.334]; % N250

tbeg = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(1));
tend = nearest(data_tla.(sessions{1}).(sprintf('stim%s_%s_%s',stimNum{1},training{1},imgConds{1})).sub(1).data.time,cfg.latency(2));

anovaData = [];

for sub = 1:length(exper.subjects)
  if ~allBadSub(sub)
    theseData = [];
    for ses = 1:length(sessions)
      for im = 1:length(imgConds)
        
        for sn = 1:length(stimNum)
          for tr = 1:length(training)
            
            condition = sprintf('stim%s_%s_%s',stimNum{sn},training{tr},imgConds{im});
            
            for h = 1:length(hemis)
              cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis(h))});
              %dat = ft_selectdata_new(cfg,data_tla.(sessions{ses}).(condition).sub(sub).data);
              %theseData = cat(2,theseData,mean(mean(dat.avg,1),2));
              
              % % collapse hemis
              % cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,hemis)});
              
              elecInd = ismember(data_tla.(sessions{ses}).(condition).sub(sub).data.label,cfg.channel);
              theseData = cat(2,theseData,mean(mean(data_tla.(sessions{ses}).(condition).sub(sub).data.avg(elecInd,tbeg:tend),1),2));
              
            end
          end
        end
      end
    end
    anovaData = cat(1,anovaData,theseData);
  end
end

varnames = {'ImgCond','StimNum','Basic/Subord','Hemisphere'};

%lev_sessions = {'pre'};

lev_imgConds = {'cong','gray','incon'};
% lev_imgConds = {'gray','g_hi8','g_lo8'};

lev_stimNum = {'s1' 's2'};
% lev_stimNum = {'s2'};

lev_hemis = {'left','right'};

levelnames = {lev_imgConds lev_stimNum training lev_hemis};

O = teg_repeated_measures_ANOVA(anovaData, [length(imgConds) length(stimNum) length(training) length(hemis)], varnames,[],[],[],[],[],[],levelnames);
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','Basic/Subord','Hemisphere'};
% levelnames = {sessions imgConds training hemis};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(training) length(hemis)], varnames);

% varnames = {'testDay', 'ImgCond','StimNum','Basic/Subord'};
% levelnames = {sessions imgConds stimNum training};
% O = teg_repeated_measures_ANOVA(anovaData, [length(sessions) length(imgConds) length(stimNum) length(training)], varnames);

%% subplots of each subject's ERPs

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.rois = {{'LAS'},{'LPS'}};
% cfg_plot.rois = {{'E70'},{'E83'}};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-0.2 1.0];
cfg_plot.ylim = [-10 10];
cfg_plot.parameter = 'avg';

% cfg_plot.rois = {{'E83'}};
% cfg_plot.xlim = [-0.2 1.0];
% cfg_plot.ylim = [-10 10];

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% % abbreviations for the condition types
% cfg_plot.typesByROI = {...
%   {'C2','C6'},...
%   {'C2','C6'}};

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotER(cfg_plot,ana,exper,data_tla);
end

% %% plot the conditions
% 
% cfg_ft = [];
% cfg_ft.xlim = [-0.2 1.5];
% cfg_ft.parameter = 'avg';
% 
% cfg_plot = [];
% 
% %cfg_plot.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -1 6; -1 6; -1 6];
% % vertical solid lines to plot
% cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
% cfg_plot.plotLegend = 0;
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
% cfg_plot.plotTitle = 0;
% 
% cfg_plot.is_ga = 1;
% cfg_plot.excludeBadSub = 1;
% 
% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% 
% % cfg_plot.condByTypeByROI = {...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   };
% % cfg_plot.typesByROI = repmat({{'C2','C6'}},size(cfg_plot.condByTypeByROI));
% 
% %cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));
% 
% for r = 1:length(cfg_plot.rois)
%   cfg_plot.roi = cfg_plot.rois{r};
%   %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
%   cfg_plot.conditions = cfg_plot.condByROI{r};
%   cfg_ft.ylim = cfg_plot.ylims(r,:);
%   cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
%   if cfg_plot.plotLegend
%     cfg_plot.legendloc = cfg_plot.legendlocs{r};
%   end
%   
%   mm_ft_plotERP(cfg_ft,cfg_plot,ana,exper,files,dirs,ga_tla);
% end

%% make some GA plots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [-0.2 1.0]; % time
% cfg_ft.xlim = [-0.2 1.5]; % time
%cfg_ft.xlim = [-0.2 2.0]; % time
cfg_ft.xlim = [-0.5 0.6]; % time

cfg_ft.parameter = 'avg';

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_plot.rois = {{'LPI2'},{'RPI2'}};
cfg_plot.rois = {{'LPI2','RPI2'}};
cfg_plot.ylims = [-0.5 6; -0.5 6];
% cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.x_bounds = [];
cfg_plot.plotLegend = 0;
cfg_plot.legendlocs = {'SouthEast','SouthEast'};

cfg_plot.xlabel = 'Time (s)';
cfg_plot.ylabel = 'Voltage (\muV)';
% cfg_plot.xlabel = '';
% cfg_plot.ylabel = '';

% cfg_plot.ftFxn = 'ft_topoplotER';
% cfg_plot.ylims = [-6 6];
% %cfg_plot.ylims = 'maxmin';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% %cfg_ft.comment = 'no';
% % cfg_plot.rois = {'all'};
% % cfg_ft.xlim = [0 1.2]; % time
% % cfg_plot.subplot = 1;
% cfg_plot.rois = {{'LAS'}};
% cfg_ft.xlim = [0.3 0.5]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [0.5 0.8]; % time
% %cfg_plot.rois = {{'LPS'}};
% %cfg_ft.xlim = [1.0 1.5]; % time

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% %cfg_plot.rois = {{'FS'},{'LAS','RAS'},{'LPS','RPS'}};
% %cfg_plot.rois = {{'FS'},{'PS'}};
% %cfg_plot.rois = {'E71'};
% cfg_plot.rois = {'all'};
% cfg_plot.ylims = [-5 5]; % voltage in multiplot
% %cfg_plot.ylims = repmat('maxmin',size(cfg_plot.rois,2),1); % voltage in multiplot

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
% cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
cfg_plot.condByROI = repmat({{{'FSC','FSI','N'}}},size(cfg_plot.rois));
cfg_plot.rename_condByROI = repmat({{{'FSC','FSI','CR'}}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));
%cfg_plot.rename_condByROI = repmat({{{'SC','SI','CR'}}},size(cfg_plot.rois));

% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% cfg_plot.condByTypeByROI = repmat({{{'HSC2','HSI2','CR2'},{'HSC6','HSI6','CR6'}}},size(cfg_plot.rois));
% % cfg_plot.condByTypeByROI = {...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}...
% %   };
% cfg_plot.typesByROI = repmat({{'C2','C6'}},size(cfg_plot.condByTypeByROI));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.types = cfg_plot.typesByROI{r};
  cfg_plot.rename_conditions = cfg_plot.rename_condByROI{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  
  if strcmp(cfg_plot.ftFxn,'ft_singleplotER')
    cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
    if cfg_plot.plotLegend
      cfg_plot.legendloc = cfg_plot.legendlocs{r};
    end
  end
  
  mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_ft = [];
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
cfg_ft.colormap = 'jet';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'no';

%cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';
%cfg_plot.conditions = {{'HSC2','CR2'},{'HSI2','CR2'},{'HSC2','HSI2'},{'HSC6','CR6'},{'HSI6','CR6'},{'HSC6','HSI6'}}; % {'H2','CR2'}, {'H6','CR6'},
%cfg_plot.conditions = {{'RH','RCR'},{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}};
% cfg_plot.conditions = {{'SC','CR'},{'SC','SI'},{'SI','CR'}};
%cfg_plot.conditions = {{'RHSC','RHSI'}};
%cfg_plot.conditions = {{'FSC','RSSI'}};
cfg_plot.conditions = {{'Face','House'}};


cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-5.5 5.5]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
% cfg_ft.comment = 'no';
cfg_ft.comment = 'xlim';
cfg_ft.commentpos = 'middletop';

cfg_plot.roi = {'E73'};
%cfg_plot.roi = {'LAS'};
% cfg_ft.xlim = [0.01 0.8]; % time

% cfg_plot.roi = {'LPS','RPS'};
% cfg_ft.xlim = [0.5 0.8]; % time

%cfg_plot.roi = {'RAS'};
%cfg_ft.xlim = [1.1 1.9]; % time
% cfg_ft.xlim = [0 1.5]; % time
% cfg_plot.roi = {'all'};

% cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time
cfg_ft.xlim = (0:0.05:1.0); % time

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.xlim = [-0.2 1.5]; % time
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.ylim = [-1 1]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_ft.xlim = [-0.2 1.5]; % time
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.ylim = [-1 1]; % volt

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
% and the times that correspond to each set of ROIs

%cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'},{'LPS','RPS'},{'LPS','RPS'},{'LAS','RAS'},{'LPS','RPS'},{'LPS','RPS'},{'LPS','RPS'}};
%cfg_ana.latencies = [0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];

% cfg_ana.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_ana.latencies = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];

% cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_ana.rois = {{'FC'}};

cfg_ana.rois = {{'LPI2','RPI2'},{'LPI2','RPI2'},{'LPI2','RPI2'}};

cfg_ana.latencies = [0.156 0.208; 0.156 0.208; 0.156 0.208];

% cfg.latency = [0.156 0.208]; % N170 (use this)
% 
% cfg.latency = [0.234 0.334]; % N250


% % LF O/N
% cfg_ana.rois = {{'RAS'},{'RAS'},{'RAI'},{'RAI'}};
% cfg_ana.latencies = [1.1 1.4; 1.4 1.9; 1.1 1.4; 1.4 1.9];

% % LPN
% cfg_ana.rois = {{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_ana.latencies = [1.2 1.8; 1.2 1.8; 1.2 1.8];


% cfg_ana.conditions = {'all'};
%cfg_ana.conditions = {{'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},{'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'}};
%cfg_ana.conditions = {{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}}; % {'RH','RCR'},
%cfg_ana.conditions = {{'SC','CR'},{'SI','CR'},{'SC','SI'}}; % {'RH','CR'},
% cfg_ana.conditions = {{'Face','House'}};
% cfg_ana.conditions = {{'RgH','CR'}};

cfg_ana.conditions = {{'stim2_basic_norm' 'stim2_subord_norm'} {'stim2_basic_g' 'stim2_subord_g'} {'stim2_basic_color' 'stim2_subord_color'}};


%cfg_ana.conditions = {{'all'}};
%cfg_ana.conditions = {{'all_within_types'}};
%cfg_ana.conditions = {{'all_across_types'}};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.parameter = 'avg';
cfg_ft.correctm = 'fdr';

% line plot parameters
cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 1;
%cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 1 4; 1 4; 1 4; -4 -1; 1 4; 1 4; 1 4];
%cfg_plot.ylims = [-4 -1; 1 4; 1 4];
% cfg_plot.ylims = [-5 -2; -5 -2; -5 -2; 1.5 4.5; 1.5 4.5;];
cfg_plot.ylims = [-2 6; -2 6; -2 6];
%cfg_plot.plot_order = {'CR2','H2','HSC2','HSI2','CR6','H6','HSC6','HSI6'};
%cfg_plot.plot_order = {'RHSC','RHSI','RCR'};
% cfg_plot.plot_order = {'SC','SI','CR'};
%cfg_plot.rename_conditions = {'SC','SI','CR'};
% cfg_plot.xlabel = 'Condition';
cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
% cfg_plot.ylabel = '';

sesNum = 1;

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  % use data_tla_avg because FieldTrip doesn't deal with dimord properly
  % when single-trials exist
  %mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla_avg,sesNum);
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla,sesNum);
end

%% output some values

cfg = [];

cfg.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg.latencies = [0.3 0.5; 0.5 0.8];
% cfg.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg.latencies = [0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];
cfg.condByROI = repmat({{'SC','SI','CR'}},size(cfg.rois));

% cfg.rois = {{'FC'}};
% cfg.latencies = [0.3 0.5];
% cfg.condByROI = repmat({{'FSC','FSI','N'}},size(cfg.rois));

cfg.parameter = 'avg';

cfg.excludeBadSub = false;

%cfg.direction = 'columns';
cfg.direction = 'rows';

mm_printDataToText(cfg,exper,ana,dirs,data_tla);

%% 3-way ANOVA: Hemisphere x Block Type x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% IV2: abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% IV3: outermost cell holds one cell for each ROI; each ROI cell holds one
% cell for each event type; each event type cell holds strings for its
% conditions
cfg_ana.condByTypeByROI = {...
  %{{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
  {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}},...
  {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

% For each ROI, what's common among the conditions in each type
cfg_ana.condCommonByROI = {...
  %{'CR','H','HSC','HSI'},...
  {'CR','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Block Type','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByTypeByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov33ER(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Hemisphere x Condition - separate colors

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% % IV2: define the conditions tested for each set of ROIs
% cfg_ana.condByROI = {...
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};
% 
% % For each ROI, what's common among the conditions in each type
% cfg_ana.condCommonByROI = {...
%   {'CR','H','HSC','HSI'},...
%   {'CR','HSC','HSI'}};
% 
% % abbreviations for the condition types
% cfg_ana.typesByROI = {...
%   {'C2','C6'},...
%   {'C2','C6'}};

% IV2: define the conditions tested for each set of ROIs
%cfg_ana.condByROI = {{'RCR','RH'},{'RCR','RHSC','RHSI'}};
cfg_ana.condByROI = {{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
% for the one-way ANOVAs
cfg_ana.condCommonByROI = {...
  %{'CR','H'},...
  {'CR','HSC','HSI'},...
  {'CR','HSC','HSI'}};

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  %cfg_ana.types = cfg_ana.typesByROI{r};
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% IV2: define the conditions tested for each set of ROIs
%cfg_ana.condByROI = {{'RH','RCR'},{'RCR','RHSC','RHSI'}};
%cfg_ana.condByROI = {{'RCR','RHSC','RHSI'},{'RCR','RHSC','RHSI'}};
cfg_ana.condByROI = {exper.eventValues,exper.eventValues};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% cfg_ana.condCommonByROI = {...
%   {'CR','HSC','HSI'},...
%   {'CR','HSC','HSI'}};
cfg_ana.condCommonByROI = {...
  exper.eventValues,...
  exper.eventValues};

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  
  mm_ft_rmaov2ER(cfg_ana,exper,ana,data_tla);
end

%% 1-way ANOVA: Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 0;

% define which regions to average across for the test
cfg_ana.rois = {{'FC'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5];

cfg_ana.condByROI = repmat({{'FSC', 'FSI', 'N'}},size(cfg_ana.rois));
% Define the IVs (type: event, roi, latency)
cfg_ana.IV1.name = 'FSC/FSI/CR';
cfg_ana.IV1.cond = {'FSC', 'FSI', 'N'};
cfg_ana.IV1.type = 'event';

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  
  mm_ft_rmaov1ER_spec(cfg_ana,exper,ana,data_tla);
end

%% cluster statistics

cfg_ft = [];
%cfg_ft.avgovertime = 'no';
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'no';

cfg_ft.parameter = 'avg';

cfg_ana = [];
cfg_ana.roi = 'all';
%cfg_ana.latencies = [0 1.0; 1.0 2.0];
%cfg_ana.latencies = [0.2 0.6];
cfg_ana.latencies = [0.3 0.5];

cfg_ft.numrandomization = 500;
% cfg_ft.clusteralpha = 0.05;
cfg_ft.clusteralpha = 0.1;
cfg_ft.alpha = 0.05;

% extra directory info
cfg_ana.dirStr = '';
if strcmp(cfg_ft.avgovertime,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgT'];
end
if strcmp(cfg_ft.avgoverchan,'yes')
  cfg_ana.dirStr = [cfg_ana.dirStr,'_avgC'];
end

% cfg_ana.conditions = {...
%   {'CR2','H2'},{'CR2','HSC2'},{'CR2','HSI2'},{'HSC2','HSI2'},...
%   {'CR6','H6'},{'CR6','HSC6'},{'CR6','HSI6'},{'HSC6','HSI6'},...
%   {'CR2','CR6'},{'H2','H6'},{'HSC2','HSC6'},{'HSI2','HSI6'}};
cfg_ana.conditions = {'all_within_types'};
%cfg_ana.conditions = {{'RCR','RH'},{'RCR','RHSC'},{'RCR','RHSI'},{'RHSC','RHSI'}};

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  
  stat_clus = mm_ft_clusterstatER(cfg_ft,cfg_ana,exper,ana,dirs,data_tla);
end

%% plot the cluster statistics

%files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = 0.1;

cfg_plot = [];
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.dirStr = cfg_ana.dirStr;

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  
  mm_ft_clusterplotER(cfg_ft,cfg_plot,ana,files,dirs);
end

%% let me know that it's done
emailme = 1;
if emailme
  subject = sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with %s tla:%s',exper.name,sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    };
  send_gmail(subject,mail_message);
end

%% correlations

cfg_ana = [];

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.dpTypesByROI = {...
  {'Item','Source'},...
  {'Item','Source'}};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds two cells, one for each
% d' type; each d' cell contains strings for its conditions
cfg_ana.condByROI = {...
  {{{'CR2','H2'},{'HSC2','HSI2'}}, {{'CR6','H6'},{'HSC6','HSI6'}}}...
  {{{'CR2','H2'},{'HSC2','HSI2'}}, {{'CR6','H6'},{'HSC6','HSI6'}}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'C2','C6'},...
  {'C2','C6'}};

% C2 d' values
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_item =  abs([]);
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_source =  abs([]);
% C6 d' values
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_item =  abs([]);
cfg_ana.(cfg_ana.typesByROI{1}{2}).d_source =  abs([]);

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data_tla);
end
