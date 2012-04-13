% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'COSI';

exper.sampleRate = 250;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 2.0];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
%exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
exper.eventValues = sort({'CCR','CHSC','CHSI','SCR','SHSC','SHSI'});
%exper.eventValues = sort({'CF','SF','CN','SN','CRO','SRO','CRS','SRS'});

% combine some events into higher-level categories
exper.eventValuesExtra.toCombine = {{'CHSC','CHSI'},{'SHSC','SHSI'}};
exper.eventValuesExtra.newValue = {{'CH'},{'SH'}};
%exper.eventValuesExtra.toCombine = {{'CCR','SCR'},{'CHSC','CHSI','SHSC','SHSI'},{'CHSC','SHSC'},{'CHSI','SHSI'}};
%exper.eventValuesExtra.newValue = {{'RCR'},{'RH'},{'RHSC'},{'RHSI'}};
%exper.eventValuesExtra.toCombine = {{'CF','SF'},{'CN','SN'},{'CRO','SRO'},{'CRS','SRS'}};
%exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 0;
exper.eventValuesExtra.equateExtrasSeparately = 0;

exper.subjects = {
  'COSI001';
  'COSI002';
  'COSI003';
  'COSI004';
  'COSI005';
  'COSI006';
  'COSI007';
%   'COSI008';
%   'COSI009';
%   'COSI010';
  'COSI011';
  'COSI012';
  'COSI013';
  'COSI014';
  'COSI015';
  'COSI016';
  'COSI017';
  'COSI018';
  'COSI019';
  'COSI020';
%   'COSI021';
  'COSI022';
  'COSI023';
  'COSI024';
  'COSI025';
  'COSI026';
  'COSI027';
  'COSI028';
  'COSI029';
  'COSI030';
  'COSI031';
  'COSI032';
  'COSI033';
  'COSI034';
  'COSI035';
  };
% COSI009 blinked a lot (> 90%); did not do second session

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.

%exper.sessions = {'session_0'};
%exper.sessions = {'session_1'};
exper.sessions = {{'session_0','session_1'}};

%% set up file and directory handling parameters

% directory where the data to read is located
%dirs.subDir = 'RK';
dirs.subDir = '';
dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
%dirs.dataDir = fullfile(exper.name,'eeg','nspp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile('/Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile('/Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile('/data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'data');

% pick the right dirs.dataroot
if exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
  %runLocally = 1;
elseif exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
  %runLocally = 1;
elseif exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
  %runLocally = 0;
elseif exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

% Use the FT chan locs file
files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

% figure printing options - see mm_ft_setSaveDirs for other options
files.saveFigs = 1;
files.figFontName = 'Helvetica';
files.figPrintFormat = 'epsc2';
files.figPrintRes = 150;

%files.figPrintFormat = 'tiff';
%files.figPrintRes = 1000;

% %% add NS's artifact information to the event structure
% nsEvFilters.eventValues = exper.eventValues;
% % CCR
% nsEvFilters.CCR.type = 'LURE_PRES';
% nsEvFilters.CCR.testType = 'side';
% nsEvFilters.CCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % CHSC
% nsEvFilters.CHSC.type = 'TARG_PRES';
% nsEvFilters.CHSC.testType = 'side';
% nsEvFilters.CHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % CHSI
% nsEvFilters.CHSI.type = 'TARG_PRES';
% nsEvFilters.CHSI.testType = 'side';
% nsEvFilters.CHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% % SCR
% nsEvFilters.SCR.type = 'LURE_PRES';
% nsEvFilters.SCR.testType = 'task';
% nsEvFilters.SCR.filters = {'rec_isTarg == 0', 'rec_correct == 1'};
% % SHSC
% nsEvFilters.SHSC.type = 'TARG_PRES';
% nsEvFilters.SHSC.testType = 'task';
% nsEvFilters.SHSC.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 1'};
% % SHSI
% nsEvFilters.SHSI.type = 'TARG_PRES';
% nsEvFilters.SHSI.testType = 'task';
% nsEvFilters.SHSI.filters = {'rec_isTarg == 1', 'rec_correct == 1', 'src_correct == 0'};
% 
% for sub = 1:length(exper.subjects)
%   for ses = 1:length(exper.sessions)
%     overwriteArtFields = 1;
%     ns_addArtifactInfo(dirs.dataroot,exper.subjects{sub},exper.sessions{ses},nsEvFilters,overwriteArtFields);
%   end
% end

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';
ana.artifact.type = {'zeroVar'};
%ana.artifact.type = {'nsAuto'};
ana.overwrite.raw = 1;

% process the data
ana.ftFxn = 'ft_timelockanalysis';
% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'tla';
ana.overwrite.proc = 1;

% any preprocessing?
cfg_pp = [];

% single precision to save space
cfg_pp.precision = 'single';

% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];

cfg_proc = [];
cfg_proc.keeptrials = 'no';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);
process_ft_data(ana,cfg_proc,exper,dirs);

%% save the analysis details

% overwrite if it already exists
saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
%if ~exist(saveFile,'file')
fprintf('Saving %s...',saveFile);
save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
fprintf('Done.\n');
%else
%  error('Not saving! %s already exists.\n',saveFile);
%end

%% let me know that it's done
emailme = 1;
if emailme
  subject = sprintf('Done with%s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:}));
  mail_message = {...
    sprintf('Done with%s %s',sprintf(repmat(' %s',1,length(exper.eventValues)),exper.eventValues{:})),...
    sprintf('%s',saveFile),...
    };
  send_gmail(subject,mail_message);
end

%% load the analysis details

%adFile = '/Volumes/curranlab/Data/COSI/eeg/nspp/-1000_2000/ft_data/CCR_CH_CHSC_CHSI_SCR_SH_SHSC_SHSI_eq0_art_nsAuto/tla_-1000_2000_avg/analysisDetails.mat';
adFile = '/Volumes/curranlab/Data/COSI/eeg/eppp/-1000_2000/ft_data/CCR_CH_CHSC_CHSI_SCR_SH_SHSC_SHSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';
[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

files.figFontName = 'Helvetica';
files.figPrintFormat = 'epsc2';
files.figPrintRes = 150;

%files.figPrintFormat = 'tiff';
%files.figPrintRes = 1000;

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% list the event values to analyze; specific to each experiment

% this is useful for when there are multiple types of event values, for
% example, hits and CRs in two conditions. You don't have to enter anything
% if you just want all events from exper.eventValues together in a single
% cell because it will get set to {exper.eventValues}, but it needs to be a
% cell containing a cell of eventValue strings

% this is only used by mm_ft_checkCondComps to create pairwise combinations
% either within event types {'all_within_types'} or across all event types
% {'all_across_types'}; mm_ft_checkCondComps is called within subsequent
% analysis functions

% list the values separated by types: Color, Side
ana.eventValues = {{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}};
%ana.eventValues = {{'RCR','RH','RHSC','RHSI'}};
%ana.eventValues = {exper.eventValues};
%ana.eventValues = {{'F','N','RO','RS'}};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end

%% load in the subject data

[data_tla] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'tla');

%% Test plots to make sure data look ok

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.layout = ft_prepare_layout([],ana);
figure
sub=1;
ses=1;
for i = 1:2
  figure
  ft_multiplotER(cfg_ft,data_tla.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

% cfg_ft = [];
% cfg_ft.channel = {'E20'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.xlim = [-0.2 1.0];
% %figure
% %ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data,data_tla.(exper.eventValues{2}).sub(1).ses(1).data,data_tla.(exper.eventValues{3}).sub(1).ses(1).data);
% %legend(strrep(exper.eventValues,'_',''),'Location','SouthEast');
% figure
% cfg_ft.graphcolor = 'b';
% ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% legend(strrep(exper.eventValues{1},'_',''),'Location','SouthEast');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Remove the dof field. I'm not sure what degrees of freedom refers to, but the tutorial says to.
% for evVal = 1:length(exper.eventValues)
%   for sub = 1:length(exper.subjects)
%     for ses = 1:length(exper.sessions)
%       if isfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof')
%         data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data = rmfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof');
%       end
%     end
%   end
% end

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {};
exper.badBehSub = {'COSI031','COSI002','COSI035','COSI007','COSI030'};

% COSI031 had too many false alarms (too many "old" responses)

% noisy subjects
%
% COSI002, COSI035

% weird/nonexistent P1/N1 complex
%
% COSI002,'COSI013','COSI017','COSI027'

% COSI017 had lots of "shorted" channels (ID'd by EP), though I'm not sure
% I should reject them based on this

% weird FN400 voltage patterns: 007, 030

% exclude subjects with low event counts
[exper] = mm_threshSubs(exper,ana,20);

%% get the grand average

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = ana.eventValues;
cfg_ana.data_str = 'data_tla';
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

ga_tla = struct;

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sessions)
  for typ = 1:length(ana.eventValues)
    for evVal = 1:length(ana.eventValues{typ})
      %tic
      fprintf('Running ft_timelockgrandaverage on %s...',ana.eventValues{typ}{evVal});
      ga_tla.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
      fprintf('Done.\n');
      %toc
    end
  end
end

%% plot the conditions - simple

cfg_ft = [];
cfg_ft.xlim = [-.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];
cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
%cfg_plot.rois = {{'RAS'},{'LPS'}};
cfg_plot.ylims = [-5 2; -5 2; -5 2; -2 5; -2 5];
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% cfg_ft.xlim = [-.2 1.0];
% cfg_plot.rois = {{'E83'}};
% cfg_plot.ylims = [-10 10];
% cfg_plot.legendlocs = {'NorthEast'};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
%   {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}}};
cfg_plot.condByTypeByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,ga_tla);
  %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
end

%% subplots of each subject's ERPs

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.rois = {{'LAS'},{'LPS'}};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-0.2 1.0];
cfg_plot.ylim = [-10 10];
cfg_plot.parameter = 'avg';

% %cfg_plot.rois = {{'E70'}}; % left
% %cfg_plot.rois = {{'E75'}}; % center
% cfg_plot.rois = {{'E83'}}; % right
% %cfg_plot.xlim = [-0.2 1.0];
% cfg_plot.xlim = [0 0.3];
% cfg_plot.ylim = [-10 10];

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
%   {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}}};
% % abbreviations for the condition types
% cfg_plot.typesByROI = {...
%   {'Color','Side'},...
%   {'Color','Side'}};

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotER(cfg_plot,ana,exper,data_tla);
end

%% individual headplots

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.xlim = [0 0.2];
cfg_ft.fontsize = 9;
cfg_ft.layout = ft_prepare_layout([],ana);
sub=3;
ses=1;
for i = 3
  figure
  ft_multiplotER(cfg_ft,data_tla.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

%% plot the conditions

cfg_ft = [];
cfg_ft.xlim = [-0.2 1.5];
cfg_ft.parameter = 'avg';

cfg_plot = [];

%cfg_plot.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
cfg_plot.rois = {{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -1 6; -1 6; -1 6];
% vertical solid lines to plot
cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.plotLegend = 1;
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
cfg_plot.plotTitle = 1;

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByTypeByROI = {...
%   {{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
%   {{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
%   {{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
%   {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}}...
%   {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}}...
%   {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}}...
%   };
% cfg_plot.typesByROI = repmat({{'Color','Side'}},size(cfg_plot.condByTypeByROI));

%cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'RHSC','RHSI','RCR'}},size(cfg_plot.rois));

cfg_plot.condByTypeByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.conditions = cfg_plot.condByROI{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
  if cfg_plot.plotLegend
    cfg_plot.legendloc = cfg_plot.legendlocs{r};
  end
  
  mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);
end

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_ft = [];
cfg_ft.xlim = [-0.2 1]; % time
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'no';

%cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';
cfg_plot.conditions = {{'CHSC','CCR'},{'CHSI','CCR'},{'CHSC','CHSI'},{'SHSC','SCR'},{'SHSI','SCR'},{'SHSC','SHSI'}}; % {'CH','CCR'},{'SH','SCR'},

cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-1 1]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;

% % cfg_plot.roi = {'LAS','RAS'};
% % cfg_ft.xlim = [0.3 0.5]; % time
% cfg_plot.roi = {'LPS','RPS'};
% cfg_ft.xlim = [0.5 0.8]; % time
% cfg_ft.comment = 'no';

cfg_ft.xlim = [0 1.0]; % time
%cfg_ft.xlim = (0:0.05:1.0); % time
cfg_plot.subplot = 1;

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.ylim = [-1 1]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.ylim = [-1 1]; % volt

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
%cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_ana.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% define the times that correspond to each set of ROIs
%cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
cfg_ana.latencies = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];

%cfg_ana.conditions = {{'CCR','CH'},{'CCR','CHSC'},{'CCR','CHSI'},{'CHSC','CHSI'},{'SCR','SH'},{'SCR','SHSC'},{'SCR','SHSI'},{'SHSC','SHSI'}};
%cfg_ana.conditions = {{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}}; % {'RH','RCR'},
cfg_ana.conditions = {{'all_within_types'}};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.parameter = 'avg';
cfg_ft.correctm = 'fdr';

% line plot parameters
cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
%cfg_plot.ylims = [-4 -1; 1 4];
cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 1 4; 1 4];
%cfg_plot.plot_order = {'CCR','CH','CHSC','CHSI','SCR','SH','SHSC','SHSI'};

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla);
end

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
  {'Color','Side'},...
  {'Color','Side'}};

% IV3: outermost cell holds one cell for each ROI; each ROI cell holds one
% cell for each event type; each event type cell holds strings for its
% conditions
cfg_ana.condByTypeByROI = {...
  %{{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
  {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}},...
  {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}}};

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

% IV2: define the conditions tested for each set of ROIs
cfg_ana.condByROI = {...
  {{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
  {{'CCR','CHSC','CHSI'},{'SCR','SHSC','SHSI'}}};

% For each ROI, what's common among the conditions in each type
cfg_ana.condCommonByROI = {...
  {'CR','H','HSC','HSI'},...
  {'CR','HSC','HSI'}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'Color','Side'},...
  {'Color','Side'}};

cfg_ana.IV_names = {'ROI','Condition'};

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.condCommon = cfg_ana.condCommonByROI{r};
  cfg_ana.types = cfg_ana.typesByROI{r};
  
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

%% cluster statistics

cfg_ft = [];
cfg_ft.avgovertime = 'no';
cfg_ft.avgoverchan = 'no';

cfg_ft.parameter = 'avg';

cfg_ana = [];
cfg_ana.roi = 'all';
cfg_ana.latencies = [0 1.0; 1.0 2.0];

% cfg_ana.conditions = {...
%   {'CCR','CH'},{'CCR','CHSC'},{'CCR','CHSI'},{'CHSC','CHSI'},...
%   {'SCR','SH'},{'SCR','SHSC'},{'SCR','SHSI'},{'SHSC','SHSI'},...
%   {'CCR','SCR'},{'CH','SH'},{'CHSC','SHSC'},{'CHSI','SHSI'}};
%cfg_ana.conditions = {'all_within_types'};
cfg_ana.conditions = {{'CHSC','CCR'},{'CHSI','CCR'},{'CHSC','CHSI'},...
  {'SHSC','SCR'},{'SHSI','SCR'},{'SHSC','SHSI'}}; % {'CH','CCR'},{'SH','SCR'},

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  
  stat_clus = mm_ft_clusterstatER(cfg_ft,cfg_ana,exper,ana,dirs,data_tla);
end

%% plot the cluster statistics

files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = .1;

cfg_plot = [];
cfg_plot.latencies = cfg_ana.latencies;
cfg_plot.conditions = cfg_ana.conditions;

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
  {{{'CCR','CH'},{'CHSC','CHSI'}}, {{'SCR','SH'},{'SHSC','SHSI'}}}...
  {{{'CCR','CH'},{'CHSC','CHSI'}}, {{'SCR','SH'},{'SHSC','SHSI'}}}};

% abbreviations for the condition types
cfg_ana.typesByROI = {...
  {'Color','Side'},...
  {'Color','Side'}};

% Color d' values
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_item =  abs([]);
cfg_ana.(cfg_ana.typesByROI{1}{1}).d_source =  abs([]);
% Side d' values
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
