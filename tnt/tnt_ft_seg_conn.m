% Connectivity analyses on time-frequency data

% See Maris & Oostenveld (2007) for a info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'TNT';

exper.sampleRate = 250;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 1.7];

% equate the number of trials across event values?
exper.equateTrials = 1;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
exper.eegFileExt = 'egis';
%exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files
exper.eventValues = sort({'B1of2','B2of2','NT1of2For','NT1of2Rec','NT2of2For','NT2of2Rec','T1of2For','T1of2Rec','T2of2For','T2of2Rec'});
%exper.eventValues = sort({'NT1of2For','NT1of2Rec','NT2of2For','NT2of2Rec','T1of2For','T1of2Rec','T2of2For','T2of2Rec'});
%exper.eventValues = sort({'NT2of2Rec','T2of2Rec'});

% combine some events into higher-level categories
exper.eventValuesExtra.toCombine = {...
  {'B1of2','B2of2'},...
  {'NT1of2For','NT2of2For','NT1of2Rec','NT2of2Rec'},...
  {'T1of2For','T2of2For','T1of2Rec','T2of2Rec'}...
%   {'NT1of2For','NT1of2Rec'},{'NT2of2For','NT2of2Rec'}...
%   {'T1of2For','T1of2Rec'},{'T2of2For','T2of2Rec'}...
%   {'NT1of2For','NT2of2For'},{'NT1of2Rec','NT2of2Rec'}...
%   {'T1of2For','T2of2For'},{'T1of2Rec','T2of2Rec'}...
  };
exper.eventValuesExtra.newValue = {...
  {'B'},...
  {'NT'},...
  {'TH'}...
%   {'NT1'},{'NT2'}...
%   {'TH1'},{'TH2'}...
%   {'NTF'},{'NTR'}...
%   {'THF'},{'THR'}...
  };

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 1;

exper.subjects = {
  'TNT 06';
  'TNT 07';
  'TNT 08';
  'TNT 09';
  'TNT 11';
  'TNT 13';
  'TNT 14';
  'TNT 15';
  'TNT 17';
  'TNT 19';
  'TNT 20';
  'TNT 21';
  'TNT 22';
  'TNT 23';
  'TNT 25';
  'TNT 26';
  'TNT 27';
  'TNT 28';
  'TNT 30';
  'TNT 32';
  'TNT 33';
  'TNT 35';
  'TNT 39';
  'TNT 41';
  'TNT 42';
  'TNT 44';
  'TNT 45';
  'TNT 46';
  'TNT 47';
  'TNT 48';
  'TNT 49';
  'TNT 50';
  'TNT 51';
  'TNT 52';
  'TNT 53';
  'TNT 54';
  };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'ses1'};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.dataDir = fullfile(exper.name,'TNT_matt','eeg',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));

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
files.figPrintFormat = 'png';
%files.figPrintFormat = 'epsc2';

%% Convert the data to FieldTrip structs - excludes NS artifact trials

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_freqanalysis';
ana.artifact.type = 'nsAuto';

% % otherFxn - other processing; must be numbered sequentially in {index}
% ana.otherFxn = {};

% ana.otherFxn{1} = 'ft_scalpcurrentdensity';
% ana.cfg_other = [];
% ana.cfg_other{1}.elecfile = files.elecfile;
% ana.cfg_other{1}.method = 'spline';
% % % set the type to go in the file name
% ana.cfg_other{1}.ftype = 'scd';

% ana.otherFxn{2} = 'ft_resampledata';
% ana.cfg_other = [];
% % set the type to go in the file name
% ana.cfg_other{2}.ftype = 'resamp';
% ana.cfg_other{2}.resamplefs = 100;
% ana.cfg_other{2}.detrend = 'no';

% any preprocessing?
cfg_pp = [];
% single precision to save space
cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.pad = 'maxperlen';

%cfg_proc.precison = 'single';

cfg_proc.output = 'fourier';
cfg_proc.channelcmb = {'all','all'};
% need to keep trials for fourier; not for powandcsd
cfg_proc.keeptrials = 'yes';
cfg_proc.keeptapers = 'yes';

%cfg_proc.output = 'powandcsd';
% % channelcmb is set up as {'channel','refchan'} pairs
%cfg_proc.channelcmb = {'E11','E62';'E11','E52';'E11','E92';'E11','E45';'E11','E108'}; % Fz, Pz; Fz, P3/P4; Fz, T7/T8
% % cfg_proc.channelcmb = {'E3','E60';'E3','E62'};
% % cfg_proc.channelcmb = {'E3','E60';'E3','E62';'E11','E60';'E11','E62'};
%cfg_proc.channel = unique(cfg_proc.channelcmb);
% % do not need to keep trials for powandcsd
%cfg_proc.keeptrials = 'no';
%cfg_proc.keeptapers = 'no';

% % MTM FFT
% cfg_proc.method = 'mtmfft';
% cfg_proc.taper = 'dpss';
% %cfg_proc.foilim = [3 50];
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% %cfg_proc.foi = 3:freqstep:50;
% cfg_proc.foi = 3:freqstep:9;
% cfg_proc.tapsmofrq = 5;
% cfg_proc.toi = -0:0.04:1.0;

% multi-taper method
cfg_proc.method = 'mtmconvol';
cfg_proc.taper = 'hanning';
%cfg_proc.taper = 'dpss';
%cfg_proc.toi = -0.8:0.04:3.0;
cfg_proc.toi = -0.5:0.04:1.0;
freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
%cfg_proc.foi = 3:freqstep:50;
cfg_proc.foi = 3:freqstep:9;
%cfg_proc.foi = 3:1:9;
%cfg_proc.foi = 2:2:30;
cfg_proc.t_ftimwin = 4./cfg_proc.foi;
% tapsmofrq is not used for hanning taper; it is used for dpss
%cfg_proc.tapsmofrq = 0.4*cfg_proc.foi;

% % wavelet
% cfg_proc.method = 'wavelet';
% cfg_proc.width = 5;
% %cfg_proc.toi = -0.8:0.04:3.0;
% cfg_proc.toi = -0.3:0.04:1.0;
% % evenly spaced frequencies, but not as many as foilim makes
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% %cfg_proc.foi = 3:freqstep:50;
% cfg_proc.foi = 3:freqstep:9;
% %cfg_proc.foilim = [3 9];

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'conn');

% set ftype to name the output file
ana.ftype = cfg_proc.output;

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

%% connectivity analysis: should I calculate coherence or phase-locking values?

% input is individual subject data with powandcsd

cfg_ft = [];
cfg_ft.method = 'wpli';
%cfg_ft.method = 'coh';
%cfg_ft.method = 'plv';
%cfg_ft.channelcmb = {'E11','E62';'E11','E52';'E11','E92';'E11','E45';'E11','E108'}; % Fz, Pz; Fz, P3/P4; Fz, T7/T8
%cfg_ft.channelcmb = {'E3','E60';'E3','E62'};
%cfg_ft.channelcmb = {'E3','E60';'E3','E62';'E11','E60';'E11','E62'};

% do some pairwise combinations
chans_to_pair = {'E1'; 'E6'; 'E9'; 'E11'; 'E15'; 'E22'; 'E24'; 'E29'; 'E32'; 'E33'; 'E34'; 'E37'; 'E41'; 'E43'; 'E45'; 'E51'; 'E52'; 'E53'; 'E55'; 'E62'; 'E64'; 'E66'; 'E69'; 'E72'; 'E75'; 'E81'; 'E84'; 'E86'; 'E87'; 'E89'; 'E92'; 'E95'; 'E97'; 'E103'; 'E108'; 'E111'; 'E116'; 'E120'; 'E122'; 'E124'};
cfg_ft.channelcmb = nchoosek(chans_to_pair,2);

% % do all pairwise combinations
% chans_to_pair = cell(129,1);
% for i = 1:length(chans_to_pair)
%   chans_to_pair{i} = sprintf('E%d',i);
% end
% chans_to_pair{end} = 'Cz';
% cfg_ft.channelcmb = nchoosek(chans_to_pair,2);

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    if iscell(exper.sessions{ses}) && length(exper.sessions{ses}) > 1
      sesStr = sprintf(repmat('%s ',1,length(exper.sessions{ses})),exper.sessions{ses}{:});
    elseif ~iscell(exper.sessions{ses}) || (iscell(exper.sessions{ses}) && length(exper.sessions{ses}) == 1)
      sesStr = exper.sessions{ses};
    end
    for evVal = 1:length(exper.eventValues)
      
      cfg_ft.inputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',ana.ftype,exper.eventValues{evVal}));
      cfg_ft.outputfile = fullfile(dirs.saveDirProc,exper.subjects{sub},sprintf('ses%d',ses),sprintf('data_%s_%s.mat',cfg_ft.method,exper.eventValues{evVal}));
      
      if ~exist(cfg_ft.outputfile,'file')
        fprintf('Processing %s, %s, %s...\n',exper.subjects{sub},sesStr,exper.eventValues{evVal});
        ft_connectivityanalysis(cfg_ft);
        fprintf('Done.\n');
      else
        fprintf('ALREADY EXISTS: %s\n',cfg_ft.outputfile);
      end
    end
  end
end

%% load the analysis details

adFile = '/Volumes/curranlab/Data/TNT/TNT_matt/eeg/-1000_1700/ft_data/B_NT_TH_eq1/conn_scd_mtmconvol_hanning_fourier_-500_980_3_9/analysisDetails.mat';

[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,false);
%[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,{'/data/projects/curranlab','/Volumes/curranlab/Data/TNT'});

%% set up channel groups

ana = mm_ft_elecGroups(ana);

%% set up the event values to analyze; specific to each experiment

% this is useful for when there are multiple types of event values, for
% example, hits and CRs in two conditions. You don't have to enter anything
% if you just want all events from exper.eventValues together in a single
% cell because it will get set to {exper.eventValues}, but it needs to be a
% cell containing a cell of eventValue strings

% this is only used by mm_ft_checkCondComps to create pairwise combinations
% either within event types {'all_within_types'} or across all event types
% {'all_across_types'}; mm_ft_checkCondComps is called within subsequent
% analysis functions

ana.eventValues = {exper.eventValues};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end

%% load some data

[data_conn] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'wpli');

%[data_conn] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'coh');

%[data_conn] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'plv');

% rename cohspctrm and plvspctrm into powspctrm because ft_freqstatistics
% has a special case for creating the labels when labelcmb exists in
% cohspctrm data, and I'd rather just stick with the non-special case

%orig_field = 'plvspctrm';
orig_field = 'wplispctrm';
%new_field = 'cohspctrm';
new_field = 'powspctrm';
if isfield(data_conn.(exper.eventValues{1}).sub(1).ses(1).data,orig_field)
  fprintf('Renaming %s into %s...',orig_field,new_field);
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      for typ = 1:length(ana.eventValues)
        for evVal = 1:length(ana.eventValues{typ})
          if isfield(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,orig_field)
            data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(new_field) = data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(orig_field);
            data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,orig_field);
          end
        end
      end
    end
  end
  fprintf('Done.\n');
end

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {};

% exclude subjects with low event counts
[exper] = mm_threshSubs(exper,ana,15);

%% plot coherence - test a single subject

% Newlabel = cell(size(data_conn.(exper.eventValues{1}).sub(1).ses(1).data.labelcmb,1),1);
% for j = 1:size(data_conn.(exper.eventValues{1}).sub(1).ses(1).data.labelcmb,1)
%   Newlabel{j} = [data_conn.(exper.eventValues{1}).sub(1).ses(1).data.labelcmb{j,1},'_',data_conn.(exper.eventValues{1}).sub(1).ses(1).data.labelcmb{j,2}];
% end
% %data_conn.(exper.eventValues{1}).sub(1).ses(1).data.labelold = data_conn.(exper.eventValues{1}).sub(1).ses(1).data.label;
% data_conn.(exper.eventValues{1}).sub(1).ses(1).data.label = Newlabel;

% one channel, time x freq x coherence
cfg_ft = [];
cfg_ft.xparam = 'time';
cfg_ft.yparam = 'freq';
%cfg_ft.parameter = 'cohspctrm';
%cfg_ft.parameter = 'plvspctrm';
cfg_ft.parameter = 'powspctrm';

% I don't know if you're supposed to baseline correct coherence
%cfg_ft.baseline = [-0.3 -0.1];
%cfg_ft.baselinetype = 'absolute';

cfg_ft.refchannel = {'E11'};
%cfg_ft.refchannel = {'E3'};

% cfg_ft.xlim = 'maxmin'; % time
% cfg_ft.ylim = 'maxmin'; % freq
cfg_ft.zlim = [0 1]; % spctrm
cfg_ft.channel = {'E62'};
%cfg_ft.channel = {'E60'};
sub=1;
ses=1;
for i = 1:3
  figure
  ft_singleplotTFR(cfg_ft,data_conn.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

% cfg_cc = [];
% cfg_cc.layout = ft_prepare_layout([],ana);
% cfg_cc.foi = [8];
% ft_topoplotCC(cfg_cc,data_conn.(exper.eventValues{1}).sub(1).ses(1).data);

% % one channel, freq x coherence (MTMFFT, no time dimension)
% cfg_ft = [];
% cfg_ft.xparam = 'freq';
% cfg_ft.parameter = 'cohspctrm';
% %cfg_ft.xlim = [3 50];
% cfg_ft.refchannel = 'E11';
% cfg_ft.channel = 'E62';
% figure
% ft_singleplotER(cfg_ft,data_conn.(exper.eventValues{1}).sub(1).ses(1).data);

cfg_ft = [];
cfg_ft.xparam = 'time';
cfg_ft.yparam = 'freq';
%cfg_ft.parameter = 'cohspctrm';
%cfg_ft.parameter = 'plvspctrm';
cfg_ft.parameter = 'powspctrm';
cfg_ft.refchannel = {'E6'};
cfg_ft.channel = 'all';
cfg_ft.xlim = [0.5 0.8]; % time
cfg_ft.ylim = [3 8]; % freq
cfg_ft.zlim = [0 1]; % spctrm
cfg_ft.layout = ft_prepare_layout([],ana);
cfg_ft.marker = 'labels';
cfg_ft.colorbar = 'yes';
sub=1;
ses=1;
for i = 1:2
  figure
  ft_singleplotTFR(cfg_ft,data_conn.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

cfg_ft = [];
cfg_ft.xparam = 'time';
cfg_ft.yparam = 'freq';
%cfg_ft.parameter = 'cohspctrm';
cfg_ft.parameter = 'plvspctrm';
cfg_ft.refchannel = {'E3'};
cfg_ft.channel = 'all';
cfg_ft.xlim = [0.5 0.8]; % time
cfg_ft.ylim = [3 8]; % freq
cfg_ft.zlim = [0 1]; % spctrm
cfg_ft.elecfile='GSN-HydroCel-129.sfp';
cfg_ft.layout = ft_prepare_layout(cfg_ft);
cfg_ft.showlabels = 'yes';
cfg_ft.colorbar = 'yes';
cfg_ft.showoutline = 'yes';
sub=1;
ses=1;
for i = 1:2
  figure
  ft_singleplotTFR(cfg_ft,data_conn.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

% cfg_ft.xlim = 'maxmin'; % time
% cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = [0 1]; % spctrm
% cfg_ft.layout = ft_prepare_layout([],data_freq);
% cfg_ft.showlabels = 'yes';
% figure
% ft_multiplotTFR(cfg_ft,data_conn.(exper.eventValues{1}).sub(1).ses(1).data);

cfg_ft = [];
cfg_ft.xparam = 'time';
cfg_ft.yparam = 'freq';
cfg_ft.parameter = 'cohspctrm';
%cfg_ft.parameter = 'plvspctrm';
%cfg_ft.parameter = 'powspctrm';
cfg_ft.refchannel = 'E6';

cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
cfg_ft.ylim = [3 8]; % freq
%cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.zlim = [0 1]; % pow

cfg_plot = [];
cfg_plot.conditions = {{'NT','TH'}};
%cfg_plot.conditions = {{'NT2of2Rec'}};
cfg_plot.plotTitle = 1;
cfg_plot.ftFxn = 'ft_topoplotTFR';
cfg_plot.subplot = 1;
cfg_ft.xlim = (0:0.1:1.0); % time

mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,data_conn.(ana.eventValues{1}{1}).sub(1).ses(1).data);

% % won't work with data_conn full structre because mm_ft_plotTFR wants
% % ga_conn
% mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,data_conn);

% statfun_indepsamplesZcoh - single subject only

%% label-field cration hack

% 1) make label, 2) grand average, 3) put labelcmb back for GA, 4) if using
% fieldname cohspctrm, put labelcmb field back for data, otherwise keep
% label

%labelcmb_backup = data_conn.(ana.eventValues{1}{1}).sub(1).ses(1).data.labelcmb;

fprintf('Data: Creating the ''label'' field {''E1 - E1''} from ''labelcmb'' {''E1'', ''E2''}, removing ''labelcmb''...');
for sub = 1:length(exper.subjects)
  ses=1;
  for typ = 1:length(ana.eventValues)
    for evVal = 1:length(ana.eventValues{typ})
      for lab =1:size(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.labelcmb,1)
        data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label{lab} = sprintf('%s - %s',...
          data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.labelcmb{lab,1},...
          data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.labelcmb{lab,2});
      end
      data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'labelcmb');
      data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label = data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label';
    end
  end
end
fprintf('Done.\n');

%% get the grand average

ga_conn = struct;

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = ana.eventValues;
cfg_ana.data_str = 'data_conn';
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sessions)
  for typ = 1:length(ana.eventValues)
    for evVal = 1:length(ana.eventValues{typ})
      %tic
      fprintf('Running ft_freqgrandaverage on %s...\n',ana.eventValues{typ}{evVal});
      ga_conn.(ana.eventValues{typ}{evVal})(ses) = eval(sprintf('ft_freqgrandaverage(cfg_ft,%s);',cfg_ana.sub_str.(ana.eventValues{typ}{evVal}){ses}));
      fprintf('Done.\n');
      %toc
    end
  end
end

%% put the labelcmb field back

% fprintf('Data: Creating the ''labelcmb'' {''E1'', ''E2''} field from ''label'' {''E1 - E1''}, removing ''label''\n');
% for sub = 1:length(exper.subjects)
%   ses=1;
%   for typ = 1:length(ana.eventValues)
%     for evVal = 1:length(ana.eventValues{typ})
%       for lab =1:size(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label,1)
%         chan1 = data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label{lab}(1:strfind(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label{lab},' - ')-1);
%         chan2 = data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label{lab}(strfind(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label{lab},' - ')+3:end);
%         data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.labelcmb{lab,1} = chan1;
%         data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.labelcmb{lab,2} = chan2;
%       end
%       data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = rmfield(data_conn.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'label');
%     end
%   end
% end
% fprintf('Done.\n');

fprintf('GA: Creating the ''labelcmb'' {''E1'', ''E2''} field from ''label'' {''E1 - E1''}, removing ''label''...');
for typ = 1:length(ana.eventValues)
  for evVal = 1:length(ana.eventValues{typ})
    for lab =1:size(ga_conn.(ana.eventValues{typ}{evVal}).label,1)
      chan1 = ga_conn.(ana.eventValues{typ}{evVal}).label{lab}(1:strfind(ga_conn.(ana.eventValues{typ}{evVal}).label{lab},' - ')-1);
      chan2 = ga_conn.(ana.eventValues{typ}{evVal}).label{lab}(strfind(ga_conn.(ana.eventValues{typ}{evVal}).label{lab},' - ')+3:end);
      ga_conn.(ana.eventValues{typ}{evVal}).labelcmb{lab,1} = chan1;
      ga_conn.(ana.eventValues{typ}{evVal}).labelcmb{lab,2} = chan2;
    end
    ga_conn.(ana.eventValues{typ}{evVal}) = rmfield(ga_conn.(ana.eventValues{typ}{evVal}),'label');
  end
end
fprintf('Done.\n');

%% make some GA plots

files.saveFigs = 0;

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [.5 1.0]; % time
cfg_ft.ylim = [3 8]; % freq
%cfg_ft.ylim = [8 12]; % freq
%cfg_ft.ylim = [12 28]; % freq
%cfg_ft.ylim = [28 50]; % freq
cfg_ft.zlim = [-0 1]; % pow

cfg_ft.parameter = 'powspctrm';

cfg_plot = [];
cfg_plot.plotTitle = 1;

% % ft_singleplotTFR
% cfg_plot.rois = {{'E62'},{'E52'}};
% cfg_plot.refchans = {'E11'};

% % ft_topoplotTFR
cfg_plot.refchans = {'E1'; 'E6'; 'E9'; 'E11'; 'E15'; 'E22'; 'E24'; 'E29'; 'E32'; 'E33'; 'E34'; 'E37'; 'E41'; 'E43'; 'E45'; 'E51'; 'E52'; 'E53'; 'E55'; 'E62'; 'E64'; 'E66'; 'E69'; 'E72'; 'E75'; 'E81'; 'E84'; 'E86'; 'E87'; 'E89'; 'E92'; 'E95'; 'E97'; 'E103'; 'E108'; 'E111'; 'E116'; 'E120'; 'E122'; 'E124'};
cfg_plot.rois = {{'all'}};

cfg_plot.is_ga = 1;
% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions
cfg_plot.condByROI = repmat({{exper.eventValues}},size(cfg_plot.refchans));
%cfg_plot.condByROI = repmat({{exper.eventValues}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'B','NT','TH'}},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'NTF','NTR','THF','THR'}},size(cfg_plot.rois));

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

% cfg_plot.ftFxn = 'ft_singleplotTFR';

cfg_plot.ftFxn = 'ft_topoplotTFR';
cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.comment = 'auto';
cfg_ft.xlim = [0.5 0.8]; % time
%cfg_plot.subplot = 1;
%cfg_ft.xlim = [0 1.0]; % time

% cfg_plot.ftFxn = 'ft_multiplotTFR';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';

for c = 1:length(cfg_plot.refchans)
  cfg_ft.refchannel = cfg_plot.refchans{c};
  cfg_plot.conditions = cfg_plot.condByROI{c};
  for r = 1:length(cfg_plot.rois)
    cfg_plot.roi = cfg_plot.rois{r};
    
    mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_conn);
    keyboard
    close all
  end
end

% %% plot the conditions
% 
% cfg_ft = [];
% cfg_ft.xlim = [-0.2 1.0];
% cfg_ft.parameter = 'powspctrm';
% 
% cfg_plot = [];
% 
% %cfg_plot.rois = {{'E11 - E62'},{'E11 - E52'}};
% cfg_plot.rois = {{'E62'},{'E52'}};
% cfg_ft.refchannel = {'E11'};
% cfg_plot.ylims = repmat([0 1],size(cfg_plot.rois));
% % vertical solid lines to plot
% cfg_plot.plotLegend = 1;
% cfg_plot.legendlocs = repmat({'NorthWest'},size(cfg_plot.rois));
% cfg_plot.plotTitle = 0;
% 
% cfg_plot.is_ga = 1;
% cfg_plot.excludeBadSub = 1;
% 
% % outermost cell holds one cell for each ROI; each ROI cell holds one cell
% % for each event type; each event type cell holds strings for its
% % conditions
% % cfg_plot.condByROI = {...
% %   {{'RCR','RH','RHSC','RHSI'}},...
% %   {{'RCR','RH','RHSC','RHSI'}},...
% %   {{'RCR','RH','RHSC','RHSI'}},...
% %   {{'RCR','RHSC','RHSI'}},...
% %   {{'RCR','RHSC','RHSI'}},...
% %   {{'RCR','RHSC','RHSI'}}};
% 
% cfg_plot.condByROI = repmat({{'NT2of2Rec','T2of2Rec'}},size(cfg_plot.rois));
% 
% for r = 1:length(cfg_plot.rois)
%   cfg_plot.roi = cfg_plot.rois{r};
%   cfg_plot.legendloc = cfg_plot.legendlocs{r};
%   cfg_ft.ylim = cfg_plot.ylims(r,:);
%   %cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
%   cfg_plot.conditions = cfg_plot.condByROI{r};
%   
%   mm_ft_plotTFR(cfg_ft,cfg_plot,ana,files,dirs,ga_conn);
% end

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
%cfg_ana.rois = {{'PS'},{'MF'},{'LPS','RPS'},{'PS'},{'PS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.5 0.8];
%cfg_ana.latencies = [0.8 1.0];
% define the frequencies that correspond to each set of ROIs
cfg_ana.frequencies = [3 8];

%cfg_ana.conditions = {{'TH','NT'}};
cfg_ana.conditions = {{'all'}};

% coherence analysis need this format: {'channel - refchannel'}
cfg_ana.rois = {{'E11 - E62'}};
%cfg_ana.rois = {{'E11 - E52'}};
%cfg_ana.rois = {{'E11 - E92'}};
%cfg_ana.rois = {{'E11 - E45'}};
%cfg_ana.rois = {{'E11 - E108'}};

% only data with plvspctrm can use this
%cfg_ana.refchannel = {'E62'};
%cfg_ana.refchannel = {'E52'};
%cfg_ana.refchannel = {'E92'};
%cfg_ana.refchannel = {'E45'};
%cfg_ana.refchannel = {'E108'};

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.avgoverfreq = 'yes';
cfg_ft.correctm = 'fdr';

% % cohspctrm needs the label field because there is a special case that
% % requires it that takes care of creating labelcmb
% cfg_ft.parameter = 'cohspctrm';

% powspctrm needs the labelcmb field
cfg_ft.parameter = 'powspctrm';

% % plvspctrm needs the labelcmb field
% cfg_ft.parameter = 'plvspctrm';

%cfg_ft.labelcmb = {'E3','E60'};

% line plot parameters
cfg_plot = [];
cfg_plot.individ_plots = 0;
cfg_plot.line_plots = 0;
%cfg_plot.ylims = [-1 1; -1 1; -1 1];
cfg_plot.ylims = [0 1];

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_ft.frequency = cfg_ana.frequencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestTFR(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_conn);
end

%% cluster statistics

cfg_ft = [];
cfg_ft.avgoverchan = 'no';
cfg_ft.avgovertime = 'no';
cfg_ft.avgoverfreq = 'yes';
%cfg_ft.avgoverfreq = 'no';

cfg_ft.parameter = 'powspctrm';

% debugging
%cfg_ft.numrandomization = 100;

cfg_ft.numrandomization = 500;

cfg_ana = [];
cfg_ana.roi = 'all';
%cfg_ana.conditions = {{'TH','NT'}};
cfg_ana.conditions = {'all'};
%cfg_ana.conditions = {{'TH','NT'},{'TRec','NTRec'},{'TFor','NTFor'},{'TRec','TFor'},{'NTRec','NTFor'}};
%cfg_ana.conditions = {{'THR','THF'},{'NTR','NTF'},{'THR','NTR'},{'THF','NTF'}};
%cfg_ana.conditions = {{'TH1','TH2'},{'NT1','NT2'},{'TH1','NT1'},{'TH2','NT2'}};

cfg_ana.frequencies = [3 8];
%cfg_ana.frequencies = [3 8; 8 12; 12 28; 28 50; 50 100];
cfg_ana.latencies = [0 1.0];
%cfg_ana.latencies = [0 0.5; 0.5 1.0];

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  for fr = 1:size(cfg_ana.frequencies,1)
    cfg_ft.frequency = cfg_ana.frequencies(fr,:);
    
    [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_conn);
  end
end

%% plot the cluster statistics

files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = .1;

cfg_plot = [];
cfg_plot.conditions = cfg_ana.conditions;
cfg_plot.frequencies = cfg_ana.frequencies;
cfg_plot.latencies = cfg_ana.latencies;

for lat = 1:size(cfg_plot.latencies,1)
  cfg_ft.latency = cfg_plot.latencies(lat,:);
  for fr = 1:size(cfg_plot.frequencies,1)
    cfg_ft.frequency = cfg_plot.frequencies(fr,:);
    
    mm_ft_clusterplotTFR(cfg_ft,cfg_plot,ana,files,dirs);
  end
end

%% old stuff

% %% try coh stats
% 
% cfg_ana = [];
% 
% % exclude the bad subjects from the subject count
% cfg_ana.numSub = length(exper.subjects) - sum(exper.badSub);
% cfg_ana.numConds = 2;
% 
% cfg_ft = [];
% cfg_ft.avgovertime = 'yes';
% cfg_ft.avgoverchan = 'yes';
% cfg_ft.avgoverfreq = 'yes';
% cfg_ft.parameter = 'cohspctrm';
% %cfg_ft.parameter = 'plvspctrm';
% 
% cfg_ft.latency = [0.5 1.0];
% % define the frequencies that correspond to each set of ROIs
% cfg_ft.frequency = [3 8];
% 
% %cfg_ft.correctm = 'fdr';
% cfg_ft.method    = 'analytic';
% cfg_ft.statistic = 'depsamplesT';
% cfg_ft.alpha     = 0.05;
% cfg_ft.correctm  = 'no';
% 
% % make the design matrix
% cfg_ft.design = zeros(2,cfg_ana.numSub*cfg_ana.numConds);
% % set the unit and independent variables
% cfg_ft.uvar = 1; % the 1st row in cfg_ft.design contains the units of observation (subject number)
% cfg_ft.ivar = 2; % the 2nd row in cfg_ft.design contains the independent variable (data/condition)
% for i = 1:cfg_ana.numSub
%   for j = 1:cfg_ana.numConds
%     subIndex = i + ((j - 1)*cfg_ana.numSub);
%     cfg_ft.design(1,subIndex) = i; % UO/DV (subject #s)
%     cfg_ft.design(2,subIndex) = j; % IV (condition #s)
%   end
% end
% 
% %cfg_ft.channel = {'E11 - E60'};
% cfg_ft.channel = {'E11 - E62'};
% %cfg_ft.channel = {'E11'};
% 
% cfg_ana.is_ga = 0;
% cfg_ana.conditions = ana.eventValues{1};
% cfg_ana.data_str = 'data_conn';
% cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);
% 
% stat1 = eval(sprintf('ft_freqstatistics(cfg_ft,%s,%s);',cell2mat(cfg_ana.sub_str.(ana.eventValues{1}{1})),cell2mat(cfg_ana.sub_str.(ana.eventValues{1}{1}))));
% 
% %% try single subject coh stats
% 
% % cfg_ana = [];
% % 
% % % exclude the bad subjects from the subject count
% % cfg_ana.numSub = length(exper.subjects) - sum(exper.badSub);
% 
% cfg_ft = [];
% 
% cfg_ft.parameter = 'fourierspctrm';
% %cfg_ft.parameter = 'cohspctrm';
% %cfg_ft.parameter = 'plvspctrm';
% cfg_ft.method    = 'montecarlo';
% cfg_ft.statistic = 'indepsamplesZcoh';
% cfg_ft.alpha     = 0.05;
% cfg_ft.correctm  = 'cluster';
% cfg_ft.numrandomization = 100;
% 
% cfg_ft.latency = [0 1.0];
% % define the frequencies that correspond to each set of ROIs
% cfg_ft.frequency = [3 8];
% 
% cfg_ft.avgovertime = 'no';
% cfg_ft.avgoverchan = 'no';
% cfg_ft.avgoverfreq = 'yes';
% 
% cfg_ft.elec = ana.elec;
% 
% % set up the design for the statistical test
% cfg_ft.ivar = 1; % the 1st row in cfg_ft.design contains the independent variable (data/condition)
% %cfg_ft.uvar = 2; % the 2nd row in cfg_ft.design contains the units of observation (subject number)
% %cfg_ft.design = [ones(1,size(data_conn.(exper.eventValues{1}).sub(1).ses(1).data.(cfg_ft.parameter),1)) 2*ones(1,size(data_conn.(exper.eventValues{2}).sub(1).ses(1).data.(cfg_ft.parameter),1))]; % IV
% cfg_ft.design = [ones(1,size(freqNT.(cfg_ft.parameter),1)) 2*ones(1,size(freqT.(cfg_ft.parameter),1))]; % IV
% %cfg_ft.design(2,1:2*cfg_ana.numSub) = [1:cfg_ana.numSub 1:cfg_ana.numSub]; % DV
% 
% %cfg_ft.label = {'E11 - E62'};
% cfg_ft.label = {'E11','E62'};
% 
% %eval(sprintf('ft_freqstatistics(cfg_ft,%s);',ana.sub_str.NT{1},ana.sub_str.TH{1}));
% 
% %stat2 = ft_freqstatistics(cfg_ft,data_conn.(exper.eventValues{1}).sub(1).ses(1).data,data_conn.(exper.eventValues{2}).sub(1).ses(1).data);
% %stat2 = ft_freqstatistics(cfg_ft,data_conn1,data_conn2);
% stat2 = ft_freqstatistics(cfg_ft,freqNT,freqT);
