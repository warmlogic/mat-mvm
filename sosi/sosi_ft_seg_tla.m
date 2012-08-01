% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SOSI';

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
%exper.eventValues = sort({'CR','SC','SI'});
%exper.eventValues = sort({'F','N','RO','RS'});

exper.eventValues = sort({'FSC','FSI','NM','NS','ROSC','ROSI','RSSC','RSSI'});

% combine some events into higher-level categories

% exper.eventValuesExtra.toCombine = {{'SC','SI'}};
% exper.eventValuesExtra.newValue = {{'RH'}};

% exper.eventValuesExtra.toCombine = {{'FSC','FSI'},{'NS','NM'},{'ROSC','ROSI'},{'RSSC','RSSI'},...
%   {'ROSC','RSSC'},{'ROSI','RSSI'}};
% exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'},...
%   {'RSC'},{'RSI'}};

%exper.eventValuesExtra.toCombine = {{'FSC','FSI'},{'NS','NM'},{'ROSC','ROSI'},{'RSSC','RSSI'}};
%exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'}};

exper.eventValuesExtra.toCombine = {{'NS','NM'},{'FSC','ROSC','RSSC'},{'FSI','ROSI','RSSI'},{'FSC','ROSC','RSSC','FSI','ROSI','RSSI'}};
exper.eventValuesExtra.newValue = {{'CR'},{'SC'},{'SI'},{'H'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 1;
exper.eventValuesExtra.equateExtrasSeparately = 0;

exper.subjects = {
  'SOSI001';
  'SOSI002';
  'SOSI003';
  'SOSI004';
  'SOSI005';
  'SOSI006';
  'SOSI007';
  'SOSI008';
  'SOSI009';
  'SOSI010';
  'SOSI011';
  'SOSI012';
  'SOSI013';
  'SOSI014';
  'SOSI015';
  'SOSI016';
  'SOSI017';
  'SOSI018';
  'SOSI019';
  'SOSI020';
  'SOSI021';
  'SOSI022';
  'SOSI023';
  'SOSI024';
  'SOSI025';
  'SOSI026';
  'SOSI027';
  'SOSI028';
  'SOSI029';
  'SOSI030';
  };
% original SOSI019 was replaced because the first didn't finish

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_0'};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.subDir = '';
%dirs.subDir = 'RK';
%dirs.subDir = 'RKSCSI';
dirs.dataDir = fullfile(exper.name,'eeg','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
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
exper.refChan = {'Cz'};

% figure printing options - see mm_ft_setSaveDirs for other options
files.saveFigs = 1;
files.figFontName = 'Helvetica';

files.figPrintFormat = 'epsc2';
files.figPrintRes = 150;

%files.figPrintFormat = 'png';
%files.figPrintFormat = 'tiff';
%files.figPrintRes = 1000;

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';
%ana.artifact.type = {'zeroVar','badChanManual','badChanEP'};
ana.artifact.type = {'zeroVar'};
ana.overwrite.raw = 1;

% process the data
ana.ftFxn = 'ft_timelockanalysis';
% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'tla';
ana.overwrite.proc = 1;

% any preprocessing?
cfg_pp = [];
% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];
% single precision to save space
%cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.keeptrials = 'no';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);
process_ft_data(ana,cfg_proc,exper,dirs);

% %% get the bad channel information
% cfg = [];
% cfg.badChanManual = false;
% cfg.badChanEP = true;
% [exper] = mm_getBadChan(cfg,exper,dirs);

% save the analysis details
saveFile = fullfile(dirs.saveDirProc,sprintf('analysisDetails.mat'));
% if ~exist(saveFile,'file')
fprintf('Saving analysis details: %s...',saveFile);
save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
fprintf('Done.\n');
% else
%   error('Not saving! %s already exists.\n',saveFile);
% end

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

% % F/N/RO/RS x SC/SI all combos
% adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/F_FSC_FSI_N_NM_NS_RO_ROSC_ROSI_RS_RSC_RSI_RSSC_RSSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';

% CR SC SI (with bad chan info), use this
adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_SC_SI_eq0_art_zeroVar_badChanManual_badChanEP/tla_-1000_2000_avg/analysisDetails.mat';

% % CR SC SI H
% adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/CR_H_SC_SI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';

%adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/RCR_RH_RHSC_RHSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';
%adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/F_N_RO_RS_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';
%adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/FSC_FSI_N_ROSC_ROSI_RSSC_RSSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';

%adFile = '/Volumes/curranlab/Data/SOSI/eeg/eppp/-1000_2000/ft_data/FSC_FSI_NM_NS_ROSC_ROSI_RSSC_RSSI_eq0_art_zeroVar_badChanManual_badChanEP/tla_-1000_2000_avg/analysisDetails.mat';

[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

% files.figFontName = 'Helvetica';
% %files.figPrintFormat = 'epsc2';
% files.figPrintFormat = 'png';
% files.figPrintRes = 150;

%files.figPrintFormat = 'tiff';
%files.figPrintRes = 1000;

%% set up channel groups

% pre-defined in this function
ana = mm_ft_elecGroups(ana);

%% create fields for analysis and GA plotting; specific to each experiment

% this is useful for when there are multiple types of event values, for
% example, hits and CRs in two conditions. You don't have to enter anything
% if you just want all events from exper.eventValues together in a single
% cell because it will get set to {exper.eventValues}, but it needs to be a
% cell containing a cell of eventValue strings

% this is only used by mm_ft_checkCondComps to create pairwise combinations
% either within event types {'all_within_types'} or across all event types
% {'all_across_types'}; mm_ft_checkCondComps is called within subsequent
% analysis functions

% ana.eventValues = {exper.eventValues};

ana.eventValues = {{'CR','SC','SI'}};
%ana.eventValues = {{'FSC','FSI','N','RSSC','RSSI','ROSC','ROSI'}};
%ana.eventValues = {{'FSC','FSI','N','RSSC','RSSI'}};
%ana.eventValues = {{'FSC','FSI','N','RSSC','ROSC'}};
%ana.eventValues = {{'FSC','FSI','N','RSSC'}};
%ana.eventValues = {{'F','N','RO','RS'}};
%ana.eventValues = {{'FSC','RSSI'}};

%ana.eventValues = {{'FSC','FSI','N','RSC','RSI'}};
% ana.eventValues = {{'FSC','FSI','N'}};
% ana.eventValues = {{'FSC','FSI'}};

% make sure ana.eventValues is set properly
if ~iscell(ana.eventValues{1})
  ana.eventValues = {ana.eventValues};
end
if ~isfield(ana,'eventValues') || isempty(ana.eventValues{1})
  ana.eventValues = {exper.eventValues};
end

%% load in the subject data

[data_tla] = mm_ft_loadSubjectData(exper,dirs,ana.eventValues,'tla');

% %% get rid of the bad channels
% 
% cfg = [];
% cfg.printRoi = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% [data_tla] = mm_rmBadChan(cfg,exper,ana,data_tla);

%% Test plots to make sure data look ok

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.layout = ft_prepare_layout([],ana);
sub=2;
ses=1;
for i = 1:2
  figure
  ft_multiplotER(cfg_ft,data_tla.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

% cfg_ft = [];
% cfg_ft.showlabels = 'yes';
% cfg_ft.interactive = 'yes';
% cfg_ft.showoutline = 'yes';
% cfg_ft.fontsize = 9;
% cfg_ft.layout = ft_prepare_layout([],ana);
% figure
% ft_multiplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data);
% 
% cfg_ft = [];
% cfg_ft.channel = {'E26'};
% %cfg_ft.linewidth = 2;
% cfg_ft.graphcolor = 'rbk';
% %cfg_ft.linestyle = {'-','--','-.'};
% figure
% ft_singleplotER(cfg_ft,data_tla.(exper.eventValues{1}).sub(1).ses(1).data,data_tla.(exper.eventValues{2}).sub(1).ses(1).data,data_tla.(exper.eventValues{3}).sub(1).ses(1).data);
% legend((exper.eventValues{1}),(exper.eventValues{2}),(exper.eventValues{3}),'Location','SouthEast');
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
% cfg_ft.layout = ft_prepare_layout([],data_tla);
% figure
% ft_topoplotER(cfg_ft,data_topo);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Remove the dof field. I'm not sure what degrees of freedom refers to, but the tutorial says to.
% for evVal = 1:length(exper.eventValues)
%   for sub = 1:length(exper.subjects)
%     for ses = 1:length(exper.sesStr)
%       if isfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof')
%         data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data = rmfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof');
%       end
%     end
%   end
% end

%% decide who to kick out based on trial counts

% Subjects with bad behavior
%exper.badBehSub = {};
%exper.badBehSub = {'SOSI011','SOSI030'};
exper.badBehSub = {'SOSI011','SOSI030','SOSI007'}; % for publication; also 001 has low trial counts

% % equating source d' 15 subjects
% exper.badBehSub = {'SOSI011','SOSI030','SOSI007','SOSI004','SOSI015','SOSI013','SOSI025','SOSI029','SOSI008','SOSI023','SOSI028','SOSI027','SOSI018','SOSI009'};
% % equating source d' 13 subjects
% exper.badBehSub = {'SOSI011','SOSI030','SOSI007','SOSI004','SOSI015','SOSI013','SOSI025','SOSI029','SOSI008','SOSI023','SOSI028','SOSI027','SOSI018','SOSI009','SOSI026','SOSI017'};



% 007 had poor accuracy (and 003, 005, 024?)

% 011 and 030 had no F responses; 001 has weird voltages

% exper.p1n1_good = {'SOSI001','SOSI006','SOSI007','SOSI008','SOSI010','SOSI012','SOSI013','SOSI017','SOSI021','SOSI023','SOSI025','SOSI030'};
% exper.p1n1_ok = {'SOSI002','SOSI004','SOSI016','SOSI019','SOSI024','SOSI026','SOSI028','SOSI029'};
% exper.p1n1_bad = {'SOSI003','SOSI005','SOSI009','SOSI011','SOSI014','SOSI015','SOSI018','SOSI020','SOSI022','SOSI027'};

% plot good
%exper.badBehSub = unique(cat(2,exper.p1n1_ok,exper.p1n1_bad));

% plot ok
%exper.badBehSub = unique(cat(2,exper.p1n1_good,exper.p1n1_bad));

% plot bad
%exper.badBehSub = unique(cat(2,exper.p1n1_good,exper.p1n1_ok));

% % subsample of subjects (thresh=14) to remove for Familiar contrasts;
% % should have 17 subjects left over
% exper.badBehSub = {
%   'SOSI001'
%   'SOSI003'
%   'SOSI005'
%   'SOSI007'
%   'SOSI009'
%   'SOSI011'
%   'SOSI012'
%   'SOSI017'
%   'SOSI018'
%   'SOSI022'
%   'SOSI023'
%   'SOSI024'
%   'SOSI030'
%   };

% % subsample of subjects (thresh=15) to remove for Familiar contrasts;
% % should have 16 subjects left over
% exper.badBehSub = {
%   'SOSI001'
%   'SOSI003'
%   'SOSI005'
%   'SOSI007'
%   'SOSI009'
%   'SOSI011'
%   'SOSI012'
%   'SOSI017'
%   'SOSI018'
%   'SOSI022'
%   'SOSI023'
%   'SOSI024'
%   'SOSI029'
%   'SOSI030'
%   };

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs(exper,ana,15);

% thresholds = [5:18];
% sumOfAvgs = zeros(1,length(thresholds));
% avgOfSum = zeros(1,length(thresholds));
% for i = 1:length(thresholds)
%   % exclude subjects with low event counts
%   [exper,ana] = mm_threshSubs(exper,ana,thresholds(i));
%   
%   sumOfAvgs(i) = sum([mean(exper.nTrials.FSC(~exper.badSub)), mean(exper.nTrials.FSI(~exper.badSub))]);
%   fprintf('Side: sum([mean(FSC),mean(FSI)]) = %.1f\n',sum([mean(exper.nTrials.FSC(~exper.badSub)), mean(exper.nTrials.FSI(~exper.badSub))]));
%   
%   avgOfSum(i) = sum([exper.nTrials.FSC(~exper.badSub); exper.nTrials.FSI(~exper.badSub)]) / length(~exper.badSub);
%   fprintf('Side: sum([(FSC),(FSI)]) / #sub = %.1f\n',sum([exper.nTrials.FSC(~exper.badSub); exper.nTrials.FSI(~exper.badSub)]) / length(~exper.badSub));
% end
% 
% figure;
% plot(thresholds,sumOfAvgs);
% title('sum of averages');
% figure;
% plot(thresholds,avgOfSum);
% title('average of sum');

%% get the grand average ERPs

% set up strings to put in grand average function
cfg_ana = [];
cfg_ana.is_ga = 0;
cfg_ana.conditions = ana.eventValues;
cfg_ana.data_str = 'data_tla';
cfg_ana.sub_str = mm_ft_catSubStr(cfg_ana,exper);

ga_tla = struct;

cfg_ft = [];
cfg_ft.keepindividual = 'no';
for ses = 1:length(exper.sesStr)
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

%% Time-Frequency: evoked (run on average ERP)

cfg_ft = [];
cfg_ft.pad = 'maxperlen';
% cfg_ft.output = 'pow';
% % cfg_ft.output = 'powandcsd';
% cfg_ft.keeptrials = 'no';
% cfg_ft.keeptapers = 'no';
cfg_ft.output = 'fourier';
cfg_ft.keeptrials = 'yes';
cfg_ft.keeptapers = 'yes';

% wavelet
cfg_ft.method = 'wavelet';
cfg_ft.width = 6;
%cfg_ft.toi = -0.8:0.04:3.0;
cfg_ft.toi = -0.5:0.04:1.0;
% % evenly spaced frequencies, but not as many as foilim makes
% freqstep = (exper.sampleRate/(diff(exper.prepost)*exper.sampleRate)) * 2;
% % cfg_ft.foi = 3:freqstep:9;
% cfg_ft.foi = 3:freqstep:60;
cfg_ft.foi = 4:1:30;
%cfg_ft.foi = 4:1:100;

baseline = [-0.4 -0.2];

data_evoked = struct;
for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
        if isfield(data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'avg')
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqanalysis(cfg_ft,data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
          
%           time = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.time;
%           
%           % power
%           param = 'powspctrm';
%           data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = abs(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm).^2;
%           blt = time >= baseline(1) & time <= baseline(2);
%           nSmp = size(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param),4);
%           blm = repmat(nanmean(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,:,:,blt),4),[1,1,1,nSmp]);
%           
%           % % divide by the baseline (relative), then log transform
%           % data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) ./ blm);
%           % % subtract and divide by the baseline (relative change), then log transform
%           % data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10((data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) - blm) ./ blm);
%           % do a log transform, then an absolute change subtraction
%           data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)) - log10(blm);
%           
%           % get rid of any zeros
%           data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) == 0) = eps(0);
%           
%           % phase
%           param = 'phasespctrm';
%           data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = angle(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm);
        else
          warning([mfilename,':erp2freq'],'NO DATA for %s! Monving on!\n',ana.eventValues{typ}{evVal});
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data;
        end
      end
    end
  end
end

%% Time-Frequency: save it

% data_evoked.SC = data_evoked.RHSC;
% data_evoked = rmfield(data_evoked,'RHSC');
% data_evoked.SI = data_evoked.RHSI;
% data_evoked = rmfield(data_evoked,'RHSI');
% data_evoked.CR = data_evoked.RCR;
% data_evoked = rmfield(data_evoked,'RCR');

save(fullfile(dirs.saveDirProc,'data_evoked.mat'),'data_evoked');

%% Time-Frequency: plots

%chan=11; % Fz
%chan=62; % Pz
chan=20; % LAS middle
%chan=118; % RAS middle
%chan=53; % LPS middle
%chan=86; % RPS middle

sub=1;
ses=1;
typ=1;
evVal = 3;

freq = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.freq;
time = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.time;

figure

subplot(3,1,1)
plot(data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.time,data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.avg(chan,:));
title(sprintf('ERP, %s, chan %d',ana.eventValues{typ}{evVal},chan));
xlabel('Time (s)');
ylabel('Voltage (\muV)');

param = 'powspctrm';
% %pow = squeeze(abs(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm).^2);
% blt = time >= baseline(1) & time <= baseline(2);
% %nTrl = size(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param),1);
% nSmp = size(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param),4);
% blm = repmat(nanmean(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,:,:,blt),4),[1,1,1,nSmp]);
% data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) - blm);
% 
% %blm = repmat(nanmean(log10(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,:,:,blt)),4),[1,1,1,nSmp]);
% 
% % % blstd = repmat(nanmean(nanstd(log10(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,:,:,blt)),0,1),4),[nTrl,1,1,nSmp]);
% 
% %data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = (log10(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)) - blm);
% %data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) - blm);
% 
% %data_pow_blc = (log10(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)) - blm);
% data_pow_blc = log10(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) ./ blm);
% 
% %data_pow_blc = (10.^(data_pow_blc/10)-1)*100;

%figure
subplot(3,1,2)
surf(time,freq,squeeze(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,chan,:,:)));
shading interp;view([0,90]);axis tight;
%imagesc(time,freq,squeeze(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,chan,:,:)),[-3 3]);
%axis xy
colorbar
title(sprintf('Power, %s, chan %d',ana.eventValues{typ}{evVal},chan));
xlabel('Time (s)');
ylabel('Frequency');

%phase = squeeze(angle(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm));
%figure
param = 'phasespctrm';
subplot(3,1,3)
surf(time,freq,squeeze(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,chan,:,:)));
shading interp;view([0,90]);axis tight;
colorbar
title(sprintf('Phase, %s, chan %d',ana.eventValues{typ}{evVal},chan));
xlabel('Time (s)');
ylabel('Frequency');

% coh = squeeze(data_pow.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm ./ abs(data_pow.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm));
% figure;
% surf(time,freq,abs(squeeze(coh(chan,:,:))));
% shading interp;view([0,90]);axis tight;

%% plot the conditions - simple

cfg_ft = [];
cfg_ft.xlim = [-0.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

% cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-5 2; -5 2; -5 2; -1 6; -1 6];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'LAS'},{'LPS'}};
% cfg_plot.ylims = [-5 2; -1 6];
cfg_plot.rois = {{'FC'}};
cfg_plot.ylims = [-5 2; -5 2];
cfg_plot.legendlocs = {'SouthEast','NorthWest'};

%cfg_plot.rois = {{'FS'},{'FI'}};
%cfg_plot.ylims = [-5 2; -5 2];
%cfg_plot.legendlocs = {'SouthEast','SouthEast'};

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

% cfg_ft.xlim = [-.2 1.0];
% cfg_plot.rois = {{'E83'}};
% cfg_plot.ylims = [-10 10];
% cfg_plot.legendlocs = {'NorthEast'};

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByROI = {...
%   {{'CR','RH','SC','SI'}},...
%   {{'CR','SC','SI'}}};

%cfg_plot.condByROI = {'all','all'};
cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
%cfg_plot.condByROI = repmat({{'SC','SI','CR'}},size(cfg_plot.rois));

% cfg_ft.linestyle = {'-','--','-.'};
% %cfg_ft.linestyle = {'-','--','-.'};
% cfg_ft.linewidth = 3;

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_simpleplotER(cfg_ft,cfg_plot,ana,exper,ga_tla);
  %print(gcf,'-dpng',sprintf('~/Desktop/%s_good_%d',exper.name,length(exper.subjects) - length(exper.badBehSub)));
end

%% subplots of each subject's ERPs

cfg_plot = [];
%cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
%cfg_plot.rois = {{'Cz'}};
cfg_plot.excludeBadSub = 0;
cfg_plot.numCols = 5;
cfg_plot.xlim = [-0.2 1.0];
cfg_plot.ylim = [-10 10];
cfg_plot.parameter = 'avg';

% cfg_plot.rois = {{'E83'}};
% cfg_plot.xlim = [-0.2 1.0];
% cfg_plot.xlim = [0 0.5];
% cfg_plot.ylim = [-10 10];

% outermost cell holds one cell for each ROI; each ROI cell holds one cell
% for each event type; each event type cell holds strings for its
% conditions

% cfg_plot.condByROI = {...
%   {'CR','RH','SC','SI'},...
%   {'CR','SC','SI'}};
%cfg_plot.condByROI = {'all','all'};
cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
  mm_ft_subjplotER(cfg_plot,ana,exper,data_tla);
end

%% make some GA plots

cfg_ft = [];
cfg_ft.colorbar = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showlabels = 'yes';
%cfg_ft.xlim = 'maxmin'; % time
%cfg_ft.ylim = 'maxmin'; % freq
% cfg_ft.zlim = 'maxmin'; % pow
%cfg_ft.xlim = [-0.2 1.0]; % time
%cfg_ft.xlim = [-0.2 2.0]; % time
cfg_ft.xlim = [-0.2 1.5]; % time

cfg_ft.parameter = 'avg';

%cfg_ft.linewidth = 2;

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

cfg_plot.ftFxn = 'ft_singleplotER';
cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -4.5 2.5; -1 6; -1 6; -1 6];
% cfg_plot.rois = {{'FC'}};
% cfg_plot.ylims = [-5 2];
cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.plotLegend = 0;
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};

% cfg_plot.xlabel = 'Time (s)';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

% cfg_plot.ftFxn = 'ft_topoplotER';
% cfg_plot.ylims = [-6 6];
% %cfg_plot.ylims = 'maxmin';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% %cfg_ft.comment = 'no';
% %cfg_plot.rois = {'all'};
% cfg_plot.subplot = 0;
% cfg_plot.rois = {{'FS'}};
% cfg_ft.xlim = [0.3 0.5]; % time
% % cfg_plot.rois = {{'LPS'}};
% % cfg_ft.xlim = [0.5 0.8]; % time
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

% cfg_plot.condByROI = repmat({{{'FSC','FSI','N'}}},size(cfg_plot.rois));
% cfg_plot.rename_condByROI = repmat({{{'FSC','FSI','CR'}}},size(cfg_plot.rois));

cfg_plot.condByROI = repmat({{{'SC','SI','CR'}}},size(cfg_plot.rois));
cfg_plot.rename_condByROI = repmat({{{'SC','SI','CR'}}},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
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

%% find the max difference between conditions

cfg = [];

cfg.contrast = {'FSC','FSI'};

% cfg.contrast = {'SC','SI'};
% cfg.contrast = {'SC','CR'};
% cfg.contrast = {'SI','CR'};
cfg.latency = [0.3 0.5];

%cfg.roi = {'LAS'};
%cfg.roi = {'FS'};
cfg.roi = {'all129'};

cfg.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg.roi)});
cfg.chansel = ismember(ft_channelselection({'all','-Fid*'},ana.elec.label),cfg.channel);
cfg.timesel = find(ga_tla.(cfg.contrast{1}).time >= cfg.latency(1) & ga_tla.(cfg.contrast{1}).time <= cfg.latency(2));

chanDiff = mean(ga_tla.(cfg.contrast{1}).avg(cfg.chansel,cfg.timesel),2) - mean(ga_tla.(cfg.contrast{2}).avg(cfg.chansel,cfg.timesel),2);

[vSort,cSort] = sort(chanDiff);
cfg.channel(cSort)

%% plot the contrasts

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_ft = [];
%cfg_ft.xlim = [-0.2 1.5]; % time
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'no';

% comparisons to make
%cfg_plot.conditions = {'all'};
%cfg_plot.conditions = {{'RH','CR'},{'SC','CR'},{'SC','SI'},{'SI','CR'}};
cfg_plot.conditions = {{'SC','CR'},{'SC','SI'},{'SI','CR'}};
%cfg_plot.conditions = {{'SC','CR'}}; % {'RH','CR'},
%cfg_plot.conditions = {{'FSC','RSSI'}};
% cfg_plot.conditions = {{'H','CR'}};

cfg_plot.ftFxn = 'ft_topoplotER';
%cfg_ft.zlim = 'maxmin'; % volt
cfg_ft.zlim = [-1.5 1.5]; % volt
%cfg_ft.zlim = [-1.0 1.0]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.comment = 'no';

 cfg_plot.roi = {'LAS','RAS'};
%cfg_plot.roi = {'LAS'};
%cfg_plot.roi = {'FS'};
cfg_ft.xlim = [0.3 0.5]; % time

% cfg_plot.roi = {'LPS','RPS'};
% cfg_ft.xlim = [0.5 0.8]; % time

%cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time
%cfg_ft.xlim = (0:0.05:1.0); % time
%cfg_plot.roi = {'PS'};

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.ylim = [-1 1]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.ylim = [-2 2]; % volt

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test and
% the times that correspond to each set of ROIs

% cfg_ana.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg_ana.latencies = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];

%cfg_ana.rois = {{'LAS'},{'RAS'}};
%cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% cfg_ana.rois = {{'LAS'}};
% cfg_ana.rois = {{'FC'}};
% cfg_ana.latencies = [0.3 0.5]; % t=13: marg
%cfg_ana.latencies = [0.375 0.45]; % t=13: sig cluster
%cfg_ana.latencies = [0.312 0.436]; % t=14: sig cluster
%cfg_ana.latencies = [0.325 0.475]; % t=14: sig

% % LF O/N
% cfg_ana.rois = {{'RAS'},{'RAS'},{'RAI'},{'RAI'}};
% cfg_ana.latencies = [1.1 1.4; 1.4 1.9; 1.1 1.4; 1.4 1.9];

% % LPN
% cfg_ana.rois = {{'LPS'},{'RPS'},{'LPS','RPS'}};
% cfg_ana.latencies = [1.2 1.8; 1.2 1.8; 1.2 1.8];

%cfg_ana.conditions = {'all'};
%cfg_ana.conditions = {{'RH','CR'},{'SC','CR'},{'SI','CR'},{'SC','SI'}};
%cfg_ana.conditions = {{'F','N'},{'RS','N'},{'RS','RO'},{'RS','F'},{'RO','N'},{'RO','F'}};
%cfg_ana.conditions = {{'SC','CR'},{'SI','CR'},{'SC','SI'}}; % {'RH','CR'},
%cfg_ana.conditions = {{'RH','CR'}}; % {'RH','CR'},
cfg_ana.conditions = {{'FSC','N'},{'FSI','N'},{'FSC','FSI'}};

% % late right frontal old/new
% cfg_ana.conditions = {{'SC','CR'}};
% cfg_ana.rois = {{'RAS'}};
% cfg_ana.latencies = [1.2 1.4];

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
% cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 2.5 5.5; 2.5 5.5];
cfg_plot.ylims = [-4 -1; 2.5 5.5];
%cfg_plot.plot_order = {'CR','RH','SC','SI'};
cfg_plot.plot_order = {'SC','SI','CR'};
%cfg_plot.rename_conditions = {'SC','SI','CR'};
% cfg_plot.xlabel = 'Condition';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

cfg_plot.linespec = 'k--o';
cfg_plot.markcolor = 'k';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  if cfg_plot.individ_plots || cfg_plot.line_plots
    cfg_plot.ylim = cfg_plot.ylims(r,:);
  end
  
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla);
end

%% output some values

cfg = [];

% cfg.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg.latencies = [0.3 0.5; 0.5 0.8];
cfg.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
cfg.latencies = [0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];
cfg.condByROI = repmat({{'SC','SI','CR'}},size(cfg.rois));

% cfg.rois = {{'FC'}};
% cfg.latencies = [0.3 0.5];
% cfg.condByROI = repmat({{'FSC','FSI','N'}},size(cfg.rois));

cfg.parameter = 'avg';

cfg.excludeBadSub = false;

cfg.direction = 'columns';
%cfg.direction = 'rows';

mm_printDataToText(cfg,exper,ana,dirs,data_tla);

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% IV1: define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

% IV2: define the conditions tested for each set of ROIs
%cfg_ana.condByROI = {{'CR','RH'},{'CR','SC','SI'}};
cfg_ana.condByROI = {{'CR','SC','SI'},{'CR','SC','SI'}};
%cfg_ana.condByROI = {exper.eventValues,exper.eventValues};

% for the one-way ANOVAs
cfg_ana.condCommonByROI = {...
  %{'CR','H'},...
  {'CR','SC','SI'},...
  {'CR','SC','SI'}};
% cfg_ana.condCommonByROI = {...
%   exper.eventValues,...
%   exper.eventValues};

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

%% run the cluster statistics

cfg_ft = [];
%cfg_ft.avgovertime = 'no';
cfg_ft.avgovertime = 'no';
cfg_ft.avgoverchan = 'yes';

cfg_ft.parameter = 'avg';

cfg_ana = [];
%cfg_ana.roi = 'all';
cfg_ana.roi = 'FC';
%cfg_ana.latencies = [0 1.0; 1.0 2.0];
%cfg_ana.latencies = [0.2 0.6];
cfg_ana.latencies = [0.3 0.5];
%cfg_ana.latencies = [0.3 0.4];

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

%cfg_ana.conditions = {'all'};
%cfg_ana.conditions = {{'CR','RH'},{'CR','SC'},{'CR','SI'},{'SC','SI'}};

cfg_ana.conditions = {{'FSC','FSI'}};
%cfg_ana.conditions = {{'FSC','N'}};

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  
  stat_clus = mm_ft_clusterstatER(cfg_ft,cfg_ana,exper,ana,dirs,data_tla);
end

%% plot the cluster statistics

files.saveFigs = 1;

cfg_ft = [];
cfg_ft.alpha = 0.2;

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
  {{{'CR','RH'},{'SC','SI'}}}...
  {{{'CR','RH'},{'SC','SI'}}}...
  };

% C2 d' values
cfg_ana.d_item =  abs([2.4097 1.0515 1.96 1.3415 0.8847 1.2262 0.185 1.9627 1.9172 1.2042 2.4869 0.9108 1.0307 1.7126 1.5148 0.9685 0.8161 2.01 1.4023 1.1812 1.1075 1.158 2.297 0.4566 1.4412 1.5941 1.4036 1.7561 1.5984 1.2559]);
cfg_ana.d_source =  abs([2.6106 1.2027 1.4935 1.6098 1.3775 1.0969 -0.0091 1.8163 2.6303 1.1994 2.3238 1.3773 1.6098 1.558 1.6074 0.8615 1.5947 2.5703 1.1201 1.1413 1.4435 0.8991 2.1249 1.0527 1.6722 1.5655 2.2779 2.1516 1.6759 1.6462]);

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.dpTypes = cfg_ana.dpTypesByROI{r};
  
  mm_ft_corr_dprimeER(cfg_ana,ana,exper,files,dirs,data_tla);
end

%% Make contrast plots (with culster stat info) - old function

% set up contrast
cfg_ana = [];
cfg_ana.include_clus_stat = 1;
cfg_ana.timeS = (0:0.05:1.0);
cfg_ana.timeSamp = round(linspace(1,exper.sampleRate,length(cfg_ana.timeS)));

cfg_plot = [];
cfg_plot.minMaxVolt = [-1 1];
cfg_plot.numRows = 4;

cfg_ft = [];
cfg_ft.interactive = 'no';
cfg_ft.elec = ana.elec;
cfg_ft.highlight = 'on';
if cfg_ana.include_clus_stat == 0
  cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS','RPS','LPS'})});
end
cfg_ft.comment = 'xlim';
cfg_ft.commentpos = 'title';

% create contrast
cont_topo = [];
cont_topo.RHvsCR = ga_tla.RH;
cont_topo.RHvsCR.avg = ga_tla.RH.avg - ga_tla.CR.avg;
cont_topo.RHvsCR.individual = ga_tla.RH.individual - ga_tla.CR.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.RHvsCR.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.RHvsCR);
end
set(gcf,'Name','H - CR')

% create contrast
cont_topo.SCvsSI = ga_tla.SC;
cont_topo.SCvsSI.avg = ga_tla.SC.avg - ga_tla.SI.avg;
cont_topo.SCvsSI.individual = ga_tla.SC.individual - ga_tla.SI.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.SCvsSI.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.SCvsSI);
end
set(gcf,'Name','SC - SI')

% create contrast
cont_topo.SCvsCR = ga_tla.SC;
cont_topo.SCvsCR.avg = ga_tla.SC.avg - ga_tla.CR.avg;
cont_topo.SCvsCR.individual = ga_tla.SC.individual - ga_tla.CR.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.SCvsCR.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.SCvsCR);
end
set(gcf,'Name','SC - CR')

% % mecklinger plot
% cfg_ft.xlim = [1200 1800];
% cfg_ft.zlim = [-2 2];
% cfg_ft.colorbar = 'yes';
% figure
% ft_topoplotER(cfg_ft,cont_topo.SCvsCR);

% create contrast
cont_topo.SIvsCR = ga_tla.SI;
cont_topo.SIvsCR.avg = ga_tla.SI.avg - ga_tla.CR.avg;
cont_topo.SIvsCR.individual = ga_tla.SI.individual - ga_tla.CR.individual;
if cfg_ana.include_clus_stat == 1
  pos = stat_clus.SIvsCR.posclusterslabelmat==1;
end
% make a plot
figure
for k = 1:length(cfg_ana.timeS)-1
  subplot(cfg_plot.numRows,(length(cfg_ana.timeS)-1)/cfg_plot.numRows,k);
  cfg_ft.xlim = [cfg_ana.timeS(k) cfg_ana.timeS(k+1)];
  cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];
  if cfg_ana.include_clus_stat == 1
    pos_int = mean(pos(:,cfg_ana.timeSamp(k):cfg_ana.timeSamp(k+1)),2);
    cfg_ft.highlightchannel = find(pos_int==1);
  end
  ft_topoplotER(cfg_ft,cont_topo.SIvsCR);
end
set(gcf,'Name','SI - CR')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % because of a bug (might be fixed now)
% if ~isfield(stat_clus.SCvsSIvsCR,'negclusters') && isfield(stat_clus.SCvsSIvsCR,'posclusters')
%   fprintf('No neg clusters found\n');
%   stat_clus.SCvsSIvsCR.negclusters.prob = .5;
%   stat_clus.SCvsSIvsCR.negclusters.clusterstat = 0;
%   stat_clus.SCvsSIvsCR.negclusterslabelmat = zeros(size(stat_clus.SCvsSIvsCR.posclusterslabelmat));
%   stat_clus.SCvsSIvsCR.negdistribution = zeros(size(stat_clus.SCvsSIvsCR.posdistribution));
% end
% if ~isfield(stat_clus.SCvsSIvsCR,'posclusters') && isfield(stat_clus.SCvsSIvsCR,'negclusters')
%   fprintf('No pos clusters found\n');
%   stat_clus.SCvsSIvsCR.posclusters.prob = 1;
%   stat_clus.SCvsSIvsCR.posclusters.clusterstat = 0;
%   stat_clus.SCvsSIvsCR.posclusterslabelmat = zeros(size(stat_clus.SCvsSIvsCR.negclusterslabelmat));
%   stat_clus.SCvsSIvsCR.posdistribution = zeros(size(stat_clus.SCvsSIvsCR.negdistribution));
% end
%
% cfg_ft = [];
% % p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
% cfg_ft.highlightsymbolseries = ['*','*','.','.','.'];
% cfg_ft.layout = ft_prepare_layout(cfg_ft,ga_tla);
% cfg_ft.contournum = 0;
% cfg_ft.emarker = '.';
% cfg_ft.alpha  = 0.05;
% cfg_ft.parameter = 'stat';
% cfg_ft.zlim = [-5 5];
% ft_clusterplot(cfg_ft,stat_clus.SCvsSIvsCR);
