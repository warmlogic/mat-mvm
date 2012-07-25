% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'COSI2';

exper.sampleRate = 500;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 2.0];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
%exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% NB: exporting to raw because the EGIS tool won't export reference chan

% types of events to find in the NS file; these must be the same as the
% events in the NS files

% exper.eventValues = sort({'CCR','CHSC','CHSI','SCR','SHSC','SHSI'});
% exper.eventValues = sort({'CF','SF','CN','SN','CRO','SRO','CRS','SRS'});

exper.eventValues = sort({...
  'CFSC','CFSI','CNM','CNS','CROSC','CROSI','CRSSC','CRSSI',...
  'SFSC','SFSI','SNM','SNS','SROSC','SROSI','SRSSC','SRSSI'});


% combine some events into higher-level categories
% exper.eventValuesExtra.toCombine = {{'CHSC','CHSI'},{'SHSC','SHSI'}};
% exper.eventValuesExtra.newValue = {{'CH'},{'SH'}};
%exper.eventValuesExtra.toCombine = {{'CCR','SCR'},{'CHSC','CHSI','SHSC','SHSI'},{'CHSC','SHSC'},{'CHSI','SHSI'}};
%exper.eventValuesExtra.newValue = {{'RCR'},{'RH'},{'RHSC'},{'RHSI'}};
%exper.eventValuesExtra.toCombine = {{'CF','SF'},{'CN','SN'},{'CRO','SRO'},{'CRS','SRS'}};
%exper.eventValuesExtra.newValue = {{'F'},{'N'},{'RO'},{'RS'}};
%exper.eventValuesExtra.toCombine = {{'CNS','CNM'},{'SNS','SNM'}};
%exper.eventValuesExtra.newValue = {{'CN'},{'SN'}};

% exper.eventValuesExtra.toCombine = {...
%   {'CFSC','CFSI'},{'CNS','CNM'},{'CROSC','CROSI'},{'CRSSC','CRSSI'},{'CROSC','CRSSC'},{'CROSI','CRSSI'},...
%   {'SFSC','SFSI'},{'SNS','SNM'},{'SROSC','SROSI'},{'SRSSC','SRSSI'},{'SROSC','SRSSC'},{'SROSI','SRSSI'}};
% exper.eventValuesExtra.newValue = {...
%   {'CF'},{'CN'},{'CRO'},{'CRS'},{'CRSC'},{'CRSI'},...
%   {'SF'},{'SN'},{'SRO'},{'SRS'},{'SRSC'},{'SRSI'}};

exper.eventValuesExtra.toCombine = {...
  {'CNS','CNM'},{'CFSC','CROSC','CRSSC'},{'CFSI','CROSI','CRSSI'},{'CFSC','CROSC','CRSSC','CFSI','CROSI','CRSSI'},...
  {'SNS','SNM'},{'SFSC','SROSC','SRSSC'},{'SFSI','SROSI','SRSSI'},{'SFSC','SROSC','SRSSC','SFSI','SROSI','SRSSI'}};
exper.eventValuesExtra.newValue = {...
  {'CCR'},{'CSC'},{'CSI'},{'CH'},...
  {'SCR'},{'SSC'},{'SSI'},{'SH'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 1;
exper.eventValuesExtra.equateExtrasSeparately = 0;

exper.subjects = {
%   'COSI2001';
%   'COSI2002';
%   'COSI2003';
%   'COSI2004';
%   'COSI2005';
%   'COSI2006';
%   'COSI2007';
  'COSI2008';
  'COSI2009';
  'COSI2010';
%   'COSI2011'; % will not have a session_1, didn't like EEG
  'COSI2012';
  'COSI2013';
%   'COSI2014'; % no session_1, didn't perform well in session_0
  'COSI2015';
  'COSI2016';
  'COSI2017';
  'COSI2018';
  'COSI2019';
  'COSI2020';
  'COSI2021';
  'COSI2022';
  'COSI2023';
  'COSI2024';
  'COSI2025';
  'COSI2026';
  'COSI2027';
  'COSI2028';
  'COSI2029';
  'COSI2030';
%   'COSI2031'; % Thought reference electrode messed up. No session_1.
  'COSI2032';
  'COSI2033';
  'COSI2034';
  'COSI2035';
  'COSI2036';
  'COSI2037';
  'COSI2038'; % COSI2038: potentially bad session_1 (bathroom, sick)
  'COSI2039';
  'COSI2040'; % EEG is pretty noisy overall
%   'COSI2041'; % COSI2041: no-show, no session_1
  'COSI2042';
  'COSI2043';
  'COSI2044';
  'COSI2045';
  };

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

% figure printing options - see mm_ft_setSaveDirs for other options
files.saveFigs = 1;
files.figFontName = 'Helvetica';
%files.figPrintFormat = 'epsc2';
files.figPrintFormat = 'png';
files.figPrintRes = 150;

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';

% eppp
%ana.artifact.type = {'zeroVar','badChanManual','badChanEP'};
ana.artifact.type = {'zeroVar'};
% % nspp
% ana.artifact.type = {'nsAuto'};

ana.overwrite.raw = 1;

% process the data
ana.ftFxn = 'ft_timelockanalysis';
% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'tla';
ana.overwrite.proc = 1;

% any preprocessing?
cfg_pp = [];

% % single precision to save space
% cfg_pp.precision = 'single';

% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];

% do a lowpass filter at 40 Hz because this is currently 100 Hz
cfg_pp.lpfilter = 'yes';
cfg_pp.lpfreq = 40;

cfg_proc = [];
cfg_proc.keeptrials = 'no';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'tla');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);
process_ft_data(ana,cfg_proc,exper,dirs);

% %% get the bad channel information
% 
% cfg = [];
% cfg.badChanManual = false;
% cfg.badChanEP = true;
% [exper] = mm_getBadChan(cfg,exper,dirs);

%% save the analysis details

% overwrite if it already exists
saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
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

% % nspp
% adFile = '/Volumes/curranlab/Data/COSI2/eeg/nspp/-1000_2000/ft_data/CCR_CSC_CSI_SCR_SSC_SSI_eq0_art_nsAuto/tla_-1000_2000_avg/analysisDetails.mat';

% % old - CR SC SI
% adFile = '/Volumes/curranlab/Data/COSI2/eeg/eppp/-1000_2000/ft_data/CCR_CSC_CSI_SCR_SSC_SSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';

% C/S x CR SC SI - with badChan info (use this)
adFile = '/Volumes/curranlab/Data/COSI2/eeg/eppp/-1000_2000/ft_data/CCR_CSC_CSI_SCR_SSC_SSI_eq0_art_zeroVar_badChanManual_badChanEP/tla_-1000_2000_avg/analysisDetails.mat';

% % C/S x CR SC SI H
% adFile = '/Volumes/curranlab/Data/COSI2/eeg/eppp/-1000_2000/ft_data/CCR_CH_CSC_CSI_SCR_SH_SSC_SSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';

% % F/N/RO/RS x SC/SI all combos
% adFile = '/Volumes/curranlab/Data/COSI2/eeg/eppp/-1000_2000/ft_data/CF_CFSC_CFSI_CN_CNM_CNS_CRO_CROSC_CROSI_CRS_CRSC_CRSI_CRSSC_CRSSI_SF_SFSC_SFSI_SN_SNM_SNS_SRO_SROSC_SROSI_SRS_SRSC_SRSI_SRSSC_SRSSI_eq0_art_zeroVar/tla_-1000_2000_avg/analysisDetails.mat';

[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

files.figFontName = 'Helvetica';
files.figPrintFormat = 'epsc2';
%files.figPrintFormat = 'png';
files.figPrintRes = 150;

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

% % list the values separated by types: Color, Side

ana.eventValues = {{'CCR','CSC','CSI'},{'SCR','SSC','SSI'}};
% ana.eventValues = {{'CCR','CSC','CSI','CH'},{'SCR','SSC','SSI','SH'}};

%ana.eventValues = {{'RCR','RH','RHSC','RHSI'}};
%ana.eventValues = {exper.eventValues};
%ana.eventValues = {{'F','N','RO','RS'}};

% % all possible categories (F N RO RS x S/I)
% ana.eventValues = {...
%   {'CFSC','CFSI','CNS','CNM','CROSC','CROSI','CRSSC','CRSSI'},...
%   {'SFSC','SFSI','SNS','SNM','SROSC','SROSI','SRSSC','SRSSI'}};
% 
% ana.eventValues = {...
%   {'CFSC','CFSI','CNS','CNM','CRSSC','CRSSI'},...
%   {'SFSC','SFSI','SNS','SNM','SRSSC','SRSSI'}};
% 
% ana.eventValues = {...
%   {'CFSC','CFSI','CN','CRSSC'},...
%   {'SFSC','SFSI','SN','SRSSC'}};
% 
% ana.eventValues = {...
%   {'CFSC','CFSI','CN','CRSC','CRSI'},...
%   {'SFSC','SFSI','SN','SRSC','SRSI'}};

% ana.eventValues = {...
%   {'CFSC','CFSI','CN'},...
%   {'SFSC','SFSI','SN'}};

% ana.eventValues = {...
%   {'CFSC','CFSI'},...
%   {'SFSC','SFSI'}};

% % F N RO RS
% ana.eventValues = {...
%   {'CF','CN','CRO','CRS'},...
%   {'SF','SN','SRO','SRS'}};

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
% [data_tla,ana] = mm_rmBadChan(cfg,exper,ana,data_tla);

%% Test plots to make sure data look ok

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.interactive = 'yes';
cfg_ft.showoutline = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.layout = ft_prepare_layout([],ana);
%figure
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
exper.badBehSub = {'COSI2008','COSI2009','COSI2020','COSI2025','COSI2038','COSI2011','COSI2014','COSI2031','COSI2041'}; % ,'COSI2035','COSI2037'

% 8, 9, 20, 25: no F responses in one color/side SC/SI bin

% 11, 14, 31, 41: no session_1

% 16, 29 have fewer than 15 trials for Side-SI;
% 39 has fewer than 15 trials for Color-SI

% 38: bad session_1 (puker)

%%%%%
% not using as a basis for exclusion
%%%%%

% 37 has too many bad channels (25) - a relevant ROI has >=20% bad channels

% 35 and 38 have noisy ERPs

% 40 session_0 has noisy EEG (NS file)


% % subsample of subjects (thresh=14) to remove for Familiar contrasts;
% % should have 13 subjects left over
% exper.badBehSub = {
%   'COSI2008'
%   'COSI2009'
%   'COSI2012'
%   'COSI2013'
%   'COSI2015'
%   'COSI2016'
%   'COSI2018'
%   'COSI2020'
%   'COSI2021'
%   'COSI2022'
%   'COSI2023'
%   'COSI2025'
%   'COSI2027'
%   'COSI2029'
%   'COSI2032'
%   'COSI2035'
%   'COSI2038'
%   'COSI2039'
%   'COSI2040'
%   'COSI2044'
%   'COSI2045'
%   };

% % subsample of subjects (thresh=15) to remove for Familiar contrasts;
% % should have 12 subjects left over
% exper.badBehSub = {
%   'COSI2008'
%   'COSI2009'
%   'COSI2012'
%   'COSI2013'
%   'COSI2015'
%   'COSI2016'
%   'COSI2018'
%   'COSI2020'
%   'COSI2021'
%   'COSI2022'
%   'COSI2023'
%   'COSI2024'
%   'COSI2025'
%   'COSI2027'
%   'COSI2029'
%   'COSI2032'
%   'COSI2035'
%   'COSI2038'
%   'COSI2039'
%   'COSI2040'
%   'COSI2044'
%   'COSI2045'
%   };

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs(exper,ana,15);

% fprintf('Color: sum([mean(FSC),mean(FSI)]) = %.1f\n',sum([mean(exper.nTrials.CFSC(~exper.badSub)), mean(exper.nTrials.CFSI(~exper.badSub))]));
% fprintf('Side: sum([mean(FSC),mean(FSI)]) = %.1f\n',sum([mean(exper.nTrials.SFSC(~exper.badSub)), mean(exper.nTrials.SFSI(~exper.badSub))]));

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
%cfg_ft.foi = 4:1:100;
cfg_ft.foi = 4:1:100;

baseline = [-0.4 -0.2];

data_evoked = struct;
for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    for typ = 1:length(ana.eventValues)
      for evVal = 1:length(ana.eventValues{typ})
        fprintf('%s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
        if isfield(data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'avg')
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_freqanalysis(cfg_ft,data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
          
          time = data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.time;
          
          % power
          param = 'powspctrm';
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = abs(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm).^2;
          blt = time >= baseline(1) & time <= baseline(2);
          nSmp = size(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param),4);
          blm = repmat(nanmean(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,:,:,blt),4),[1,1,1,nSmp]);
          
          % % divide by the baseline (relative), then log transform
          % data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) ./ blm);
          % % subtract and divide by the baseline (relative change), then log transform
          % data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10((data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) - blm) ./ blm);
          % do a log transform, then an absolute change subtraction
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = log10(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)) - log10(blm);
          
          % get rid of any zeros
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) == 0) = eps(0);
          
          % phase
          param = 'phasespctrm';
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param) = angle(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm);
        else
          warning([mfilename,':erp2freq'],'NO DATA for %s! Monving on!\n',ana.eventValues{typ}{evVal});
          data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = data_tla.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data;
        end
      end
    end
  end
end

%% Time-Frequency: save it

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
%cfg_ft.xlim = [-0.5 1.5];
cfg_ft.xlim = [-0.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

% cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.4; -4.5 2.4; -2 5; -2 5];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

%cfg_plot.rois = {{'LAS'},{'LPS'}};
%cfg_plot.ylims = [-5.5 2.5; -2 5];
cfg_plot.rois = {{'LAS'},{'FC'},{'C'}};
%cfg_plot.rois = {{'C'}};
%cfg_plot.rois = {{'LAS'}};
cfg_plot.ylims = [-5.5 2.5; -5.5 2.5; -5.5 2.5];
cfg_plot.legendlocs = {'SouthEast','NorthWest','NorthWest'};

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
%cfg_plot.condByTypeByROI = repmat({{{'CSC','CSI','CCR'},{'SSC','SSI','SCR'}}},size(cfg_plot.rois));

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
cfg_plot.excludeBadSub = 0;
cfg_plot.parameter = 'avg';
cfg_plot.numCols = 5;

cfg_plot.rois = {{'LAS'},{'LPS'}};
cfg_plot.xlim = [-0.2 1.0];
cfg_plot.ylim = [-10 10];

% %cfg_plot.rois = {{'E70'}}; % left
% cfg_plot.rois = {{'E75'}}; % center
% %cfg_plot.rois = {{'E83'}}; % right
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
for i = 1:length(ana.eventValues{1})
  figure
  ft_multiplotER(cfg_ft,data_tla.(ana.eventValues{1}{i}).sub(sub).ses(ses).data);
  title(ana.eventValues{1}{i});
end

%% make some GA plots

cfg_ft = [];
cfg_ft.xlim = [-0.2 1.5];
% cfg_ft.xlim = [-0.2 1.0];
cfg_ft.parameter = 'avg';

cfg_plot = [];

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

cfg_plot.ftFxn = 'ft_singleplotER';
cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5; -2 5];
cfg_plot.rois = {{'FC'}};
cfg_plot.ylims = [-4.5 2.5];
cfg_plot.x_bounds = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8; 0.5 0.8];
cfg_plot.plotLegend = 0;
cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','SouthEast','NorthWest','NorthWest','NorthWest'};
cfg_plot.plotTitle = 0;

% cfg_plot.xlabel = 'Time (s)';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

% cfg_plot.ftFxn = 'ft_topoplotER';
% cfg_plot.plotLegend = 0;
% cfg_plot.ylims = [-6 6];
% %cfg_plot.ylims = 'maxmin';
% %cfg_ft.marker = 'on';
% cfg_ft.marker = 'labels';
% cfg_ft.markerfontsize = 9;
% %cfg_ft.comment = 'no';
% % cfg_plot.rois = {'all'};
% % cfg_ft.xlim = [0 1.2]; % time
% % cfg_plot.subplot = 1;
% %cfg_plot.rois = {{'LAS'}};
% cfg_plot.rois = {{'FS'}};
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

% INSTRUCTIONS: outermost cell holds one cell for each ROI; each ROI cell
% holds one cell for each event type; each event type cell holds strings
% for its conditions

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
% cfg_plot.rename_condByROI = repmat({ana.eventValues},size(cfg_plot.rois));
cfg_plot.typesByROI = repmat({{'Color','Side'}},size(cfg_plot.condByTypeByROI));
cfg_plot.rename_condByROI = repmat({{{'FSC','FSI','CR'},{'FSC','FSI','CR'}}},size(cfg_plot.rois));

% cfg_plot.condByTypeByROI = repmat({{{'CSC','CSI','CCR'},{'SSC','SSI','SCR'}}},size(cfg_plot.rois));
% cfg_plot.typesByROI = repmat({{'Color','Side'}},size(cfg_plot.condByTypeByROI));
% cfg_plot.rename_condByROI = repmat({{{'SC','SI','CR'}}},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  %cfg_plot.conditions = cfg_plot.condByROI{r};
  cfg_plot.types = cfg_plot.typesByROI{r};
  cfg_plot.rename_conditions = cfg_plot.rename_condByROI{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  if strcmp(cfg_plot.ftFxn,'ft_singleplotER')
    cfg_plot.x_bound = cfg_plot.x_bounds(r,:);
  end
  if cfg_plot.plotLegend
    cfg_plot.legendloc = cfg_plot.legendlocs{r};
  end
  
  mm_ft_plotER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);
end

%% find the max difference between conditions

cfg = [];

% cfg.contrast = {'CFSC','CFSI'};

% cfg.contrast = {'CSC','CSI'};
% cfg.contrast = {'CSC','CCR'};
% cfg.contrast = {'CSI','CCR'};

cfg.contrast = {'SFSC','SFSI'};

% cfg.contrast = {'SSC','SSI'};
% cfg.contrast = {'SSC','SCR'};
% cfg.contrast = {'SSI','SCR'};

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
cfg_ft.parameter = 'avg';
cfg_ft.interactive = 'no';
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'no';

%cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';
cfg_plot.conditions = {{'CSC','CCR'},{'CSI','CCR'},{'CSC','CSI'},{'SSC','SCR'},{'SSI','SCR'},{'SSC','SSI'}};
%cfg_plot.conditions = {{'CH','CCR'},{'SH','SCR'}};

cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-1.5 1.5]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;

cfg_ft.comment = 'no';

% cfg_plot.roi = {'LAS','RAS'};
% cfg_ft.xlim = [0.3 0.5]; % time

cfg_plot.roi = {'LPS','RPS'};
cfg_ft.xlim = [0.5 0.8]; % time

% cfg_ft.xlim = [0 1.0]; % time
% %cfg_ft.xlim = (0:0.05:1.0); % time
% cfg_plot.subplot = 1;

% cfg_plot.ftFxn = 'ft_multiplotER';
% cfg_ft.showlabels = 'yes';
% cfg_ft.comment = '';
% cfg_ft.xlim = [-0.2 1]; % time
% cfg_ft.ylim = [-1 1]; % volt

% cfg_plot.ftFxn = 'ft_singleplotER';
% cfg_plot.roi = {'LPS'};
% cfg_ft.showlabels = 'yes';
% cfg_ft.xlim = [-0.2 1]; % time
% cfg_ft.ylim = [-1 1]; % volt

mm_ft_contrastER(cfg_ft,cfg_plot,ana,files,dirs,ga_tla);

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
% and define the times that correspond to each set of ROIs

%cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_ana.rois = {{'FC'}};
cfg_ana.latencies = [0.3 0.5];
%cfg_ana.rois = {{'E12'}};
%cfg_ana.rois = {{'LAS'},{'LPS'}};
%cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
%cfg_ana.rois = {{'FS'},{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
%cfg_ana.latencies = [0.3 0.5; 0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];

%cfg_ana.conditions = {{'CSC','CCR'},{'CSI','CCR'},{'CSC','CSI'},{'SSC','SCR'},{'SSI','SCR'},{'SSC','SSI'}};
% cfg_ana.conditions = {{'CSC','CCR'},{'CSI','CCR'},{'CSC','CSI'}};
% cfg_ana.conditions = {{'SSC','SCR'},{'SSI','SCR'},{'SSC','SSI'}};
cfg_ana.conditions = {{'all_within_types'}};
%cfg_ana.conditions = {{'CFSC','CFSI'},{'SFSC','SFSI'},{'CFSC','CNS'},{'CFSI','CNS'},{'SFSC','SNS'},{'SFSI','SNS'}};

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
%cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 1 4; 1 4];
cfg_plot.ylims = [-4 -1; 1 4];
%cfg_plot.plot_order = {'CSC','CSI','CCR','SSC','SSI','SCR'};
cfg_plot.plot_order = {'CSC','CSI','CCR'};
% cfg_plot.plot_order = {'SSC','SSI','SCR'};
cfg_plot.rename_conditions = {'SC','SI','CR'};
% cfg_plot.xlabel = 'Condition';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

cfg_plot.linespec = 'k--o';
cfg_plot.plotLegend = 1;
cfg_plot.legendtext = {'Color'};
cfg_plot.markcolor = 'w';
% cfg_plot.legendtext = {'Location'};
% cfg_plot.markcolor = 'k';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla);
end

%% output some values

cfg = [];

% % cfg.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% % cfg.latencies = [0.3 0.5; 0.5 0.8];
% cfg.rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};
% cfg.latencies = [0.3 0.5; 0.3 0.5; 0.5 0.8; 0.5 0.8];
% cfg.condByROI = repmat({{'CSC','CSI','CCR','SSC','SSI','SCR'}},size(cfg.rois));

cfg.rois = {{'FC'}};
cfg.latencies = [0.3 0.5];
cfg.condByROI = repmat({{'CFSC','CFSI','CN','SFSC','SFSI','SN'}},size(cfg.rois));

cfg.excludeBadSub = false;

cfg.parameter = 'avg';

cfg.direction = 'columns';
% cfg.direction = 'rows';

mm_printDataToText(cfg,exper,ana,dirs,data_tla);

%% 3-way ANOVA: Hemisphere x Block Type x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 0;

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

% cfg_ana.condByTypeByROI = {...
%   %{{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
%   {{'CCR','CSC','CSI'},{'SCR','SSC','SSI'}},...
%   {{'CCR','CSC','CSI'},{'SCR','SSC','SSI'}}};

cfg_ana.condByTypeByROI = {...
  %{{'CCR','CH','CHSC','CHSI'},{'SCR','SH','SHSC','SHSI'}},...
  {{'CSC','CSI'},{'SSC','SSI'}},...
  {{'CSC','CSI'},{'SSC','SSI'}}};

% cfg_ana.condByTypeByROI = {...
%   {{'CF','CN','CRO','CRS'},{'SF','SN','SRO','SRS'}},...
%   {{'CF','CN','CRO','CRS'},{'SF','SN','SRO','SRS'}}};

% For each ROI, what's common among the conditions in each type

% cfg_ana.condCommonByROI = {...
%   %{'CR','H','HSC','HSI'},...
%   {'CR','SC','SI'},...
%   {'CR','SC','SI'}};

cfg_ana.condCommonByROI = {...
  {'SC','SI'},...
  {'SC','SI'}};

% cfg_ana.condCommonByROI = {...
%   {'F','N','RO','RS'},...
%   {'F','N','RO','RS'}};

cfg_ana.IV_names = {'ROI','Source Condition','Trial Condition'};

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
  {{'CCR','CSC','CSI'},{'SCR','SSC','SSI'}},...
  {{'CCR','CSC','CSI'},{'SCR','SSC','SSI'}}};

% For each ROI, what's common among the conditions in each type
cfg_ana.condCommonByROI = {...
  {'CR','SC','SI'},...
  {'CR','SC','SI'}};

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

%% 2-way ANOVA: Source condition x Trial condition - critical

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 0;

% define which regions to average across for the test
cfg_ana.rois = {{'FC'}};

% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5];

%cfg_ana.condByROI = repmat({{'NT1of2', 'NT2of2', 'T1of2', 'T2of2', 'B1of2', 'B2of2'}},size(cfg_ana.rois));
% cfg_ana.condByROI = repmat({{'CFSC', 'CFSI', 'SFSC', 'SFSI'}},size(cfg_ana.rois));
cfg_ana.condByROI = repmat({{'CFSC', 'CFSI', 'CN', 'SFSC', 'SFSI', 'SN'}},size(cfg_ana.rois));
% Define the IVs (type: event, roi, latency)
cfg_ana.IV1.name = 'source condition';
cfg_ana.IV1.cond = {'C','S'};
cfg_ana.IV1.type = 'event';
cfg_ana.IV2.name = 'trial condition';
% cfg_ana.IV2.cond = {'FSC','FSI'};
cfg_ana.IV2.cond = {'FSC','FSI','N'};
cfg_ana.IV2.type = 'event';

cfg_ana.parameter = 'avg';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ana.conditions = cfg_ana.condByROI{r};
  cfg_ana.latency = cfg_ana.latencies(r,:);
  
  mm_ft_rmaov2ER_spec(cfg_ana,exper,ana,data_tla);
end

%% cluster statistics

cfg_ft = [];
%cfg_ft.avgovertime = 'no';
cfg_ft.avgovertime = 'no';
cfg_ft.avgoverchan = 'yes';

cfg_ft.parameter = 'avg';

cfg_ana = [];
cfg_ana.roi = 'all';
%cfg_ana.latencies = [0 1.0; 1.0 2.0];
%cfg_ana.latencies = [0 1.0];
cfg_ana.latencies = [0.2 0.5];
% cfg_ana.latencies = [0.3 0.5];

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
%   {'CCR','CH'},{'CCR','CHSC'},{'CCR','CHSI'},{'CHSC','CHSI'},...
%   {'SCR','SH'},{'SCR','SHSC'},{'SCR','SHSI'},{'SHSC','SHSI'},...
%   {'CCR','SCR'},{'CH','SH'},{'CHSC','SHSC'},{'CHSI','SHSI'}};

% cfg_ana.conditions = {...
%   {'CSC','CCR'},{'CSI','CCR'},{'CSC','CSI'},...
%   {'SSC','SCR'},{'SSI','SCR'},{'SSC','SSI'}};

% cfg_ana.conditions = {'all_within_types'};

cfg_ana.conditions = {...
  {'CFSC','CFSI'},...
  {'SFSC','SFSI'}};

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
