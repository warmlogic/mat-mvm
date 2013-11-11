% Make plots and do analyses for timelocked EEG (ERPs)

% See Maris & Oostenveld (2007) for info on nonparametric statistics

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'SPACE';

exper.sampleRate = 250;

% pre- and post-stimulus times to read, in seconds (pre is negative)
exper.prepost = [-1.0 1.0];

% equate the number of trials across event values?
exper.equateTrials = 0;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
exper.eegFileExt = 'raw';

% types of events to find in the NS file; these must be the same as the
% events in the NS files; or space_trialfun.m must be set up to find the
% corrct events
exper.eventValues = sort({'Face', 'House'});

% exper.eventValues = sort({'Face VA','Face SA','Face SU','Face VU','House VA','House SA','House SU','House VU'});
% 
% % combine some events into higher-level categories
% exper.eventValuesExtra.toCombine = {{'Face VA','Face SA','Face SU','Face VU'},{'House VA','House SA','House SU','House VU'}};
% exper.eventValuesExtra.newValue = {{'Face'},{'House'}};

% keep only the combined (extra) events and throw out the original events?
exper.eventValuesExtra.onlyKeepExtras = 0;
exper.eventValuesExtra.equateExtrasSeparately = 0;

exper.subjects = {
%   'SPACE001';
%   'SPACE002';
%   'SPACE003';
%   'SPACE004';
  'SPACE005';
%   'SPACE006';
%   'SPACE007';
  };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {'session_1'};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.subDir = '';
%dirs.subDir = 'RKSCSI';
%dirs.subDir = 'RK';
% dirs.dataDir = fullfile(exper.name,'EEG','Sessions','face_house_ratings','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
dirs.behDir = fullfile(exper.name,'Behavioral','Sessions',dirs.subDir);
dirs.dataDir = fullfile(exper.name,'EEG','Sessions','face_house_ratings','ftpp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
% Possible locations of the data files (dataroot)
dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'data');

% pick the right dirs.dataroot
if isfield(dirs,'serverDir') && exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
  %runLocally = 1;
elseif isfield(dirs,'serverLocalDir') && exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
  %runLocally = 1;
elseif isfield(dirs,'dreamDir') && exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
  %runLocally = 0;
elseif isfield(dirs,'localDir') && exist(dirs.localDir,'dir')
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

%files.figPrintFormat = 'png';
%files.figPrintFormat = 'tiff';
%files.figPrintRes = 1000;

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';

% ana.continuous = 'no';
% ana.trialFxn = 'seg_trialfun';

ana.continuous = 'yes';
ana.trialFxn = 'space_trialfun';

ana.artifact.type = {'ftManual', 'ftICA'};
%ana.artifact.type = {'zeroVar', 'badChanManual', 'badChanEP'};
% ana.artifact.type = {'zeroVar'};
ana.overwrite.raw = 1;

% process the data
ana.ftFxn = 'ft_timelockanalysis';
% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'tla';
ana.overwrite.proc = 1;

% any preprocessing? (run after processing artifacts)
cfg_pp = [];
% average rereference
cfg_pp.reref = 'yes';
cfg_pp.refchannel = 'all'; 
cfg_pp.implicitref = 'Cz'; 
% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];
% single precision to save space
%cfg_pp.precision = 'single';

cfg_proc = [];
cfg_proc.keeptrials = 'yes';

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

% save the analysis details

% overwrite if it already exists
saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
if ~exist(saveFile,'file')
  fprintf('Saving analysis details: %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  fprintf('Done.\n');
else
  % combine them
  fprintf('Combining these analysis details with previous runs: %s...',saveFile);
  
  old_anaDetails = load(saveFile);
  
  if all(ismember(old_anaDetails.exper.eventValues,exper.eventValues))
    subjects_new = ~ismember(exper.subjects,old_anaDetails.exper.subjects);
    
    if any(subjects_new)
      subjects_all = cat(1,old_anaDetails.exper.subjects,exper.subjects(subjects_new));
      
      newSubInd = find(subjects_new);
      
      nTrials_all = old_anaDetails.exper.nTrials;
      nTr_fn = fieldnames(nTrials_all);
      for i = 1:length(nTr_fn)
        for j = 1:length(newSubInd)
          nTr_thisSub = exper.nTrials.(nTr_fn{i});
          nTrials_all.(nTr_fn{i}) = cat(1,nTrials_all.(nTr_fn{i}),nTr_thisSub(newSubInd(j),:));
        end
      end
      
      badChan_all = cat(1,old_anaDetails.exper.badChan,exper.badChan(newSubInd,:));
      badEv_all = cat(1,old_anaDetails.exper.badEv,exper.badEv(newSubInd,:));
      
      % add the combined info into the struct we want to save
      exper.subjects = subjects_all;
      exper.nTrials = nTrials_all;
      exper.badChan = badChan_all;
      exper.badEv = badEv_all;
      
      save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
      %error('Not saving! %s already exists.\n',saveFile);
      fprintf('Done.\n');
    else
      fprintf('There are no new subjects in exper.subjects (everyone is already in the old analysisDetails.mat).  Not saving new analysisDetails.mat.\n');
    end
  else
    fprintf('exper.eventValues does not match up between old and new subjects! Not saving new analysisDetails.mat.\n');
  end
end

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

adFile = '/Users/matt/data/SPACE/EEG/Sessions/face_house_ratings/ftpp/-1000_1000/ft_data/Face_Face_SA_Face_SU_Face_VA_Face_VU_House_House_SA_House_SU_House_VA_House_VU_eq0_art_ftManual_ftICA/tla_-1000_1000/analysisDetails.mat';

% adFile = '/Volumes/curranlab/Data/SPACE/EEG/Sessions/face_house_ratings/eppp/-1000_1000/ft_data/Face_Face_SA_Face_SU_Face_VA_Face_VU_House_House_SA_House_SU_House_VA_House_VU_eq0_art_zeroVar/tla_-1000_1000/analysisDetails.mat';

[exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_ft_loadAD(adFile,true);

files.figFontName = 'Helvetica';
files.figPrintFormat = 'epsc2';
%files.figPrintFormat = 'png';
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

% list the values separated by types: 2Colors, 6Colors
%ana.eventValues = {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}};
%ana.eventValues = {{'HSC2','HSI2','CR2'},{'HSC6','HSI6','CR6'}};
%ana.eventValues = {{'RCR','RH','RHSC','RHSI'}};
%ana.eventValues = {{'RCR','RHSC','RHSI'}};

% ana.eventValues = {exper.eventValues};

ana.eventValues = {{'Face','House'}};
%ana.eventValues = {{'FSC','FSI','N','RSSC','RSSI','ROSC','ROSI'}};
%ana.eventValues = {{'FSC','FSI','N','RSSC','RSSI'}};
%ana.eventValues = {{'FSC','FSI','N','RSSC','ROSC'}};
%ana.eventValues = {{'F','N','RO','RS'}};
%ana.eventValues = {{'FSC','RSSI'}};

% ana.eventValues = {{'FSC','FSI','N','RSC','RSI'}};
% ana.eventValues = {{'FSC','FSI','N'}};

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
sub=1;
ses=1;
for i = 1:length(ana.eventValues{1})
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
exper.badBehSub = {};

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs(exper,ana,15);

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

data_evoked.SC = data_evoked.RHSC;
data_evoked = rmfield(data_evoked,'RHSC');
data_evoked.SI = data_evoked.RHSI;
data_evoked = rmfield(data_evoked,'RHSI');
data_evoked.CR = data_evoked.RCR;
data_evoked = rmfield(data_evoked,'RCR');

save(fullfile(dirs.saveDirProc,'data_evoked.mat'),'data_evoked');

%% Time-Frequency: plots

%chan=11; % Fz
%chan=62; % Pz
%chan=20; % LAS middle
%chan=118; % RAS middle
chan=53; % LPS middle
%chan=86; % RPS middle

sub=1;
ses=1;
typ=1;
evVal = 2;

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
%surf(time,freq,squeeze(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,chan,:,:)));
%shading interp;view([0,90]);axis tight;
imagesc(time,freq,squeeze(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,chan,:,:)),[-2 2]);
axis xy
%colorbar
title(sprintf('Power, %s, chan %d',ana.eventValues{typ}{evVal},chan));
xlabel('Time (s)');
ylabel('Frequency');

%phase = squeeze(angle(data_freq.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.fourierspctrm));
%figure
param = 'phasespctrm';
subplot(3,1,3)
surf(time,freq,squeeze(data_evoked.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(:,chan,:,:)));
shading interp;view([0,90]);axis tight;
%colorbar
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

% %cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% cfg_plot.rois = {{'LAS'},{'RAS'},{'FS'},{'LPS'},{'RPS'}};
% cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5];
% cfg_plot.legendlocs = {'SouthEast','SouthEast','SouthEast','NorthWest','NorthWest'};

cfg_plot.rois = {{'FS'},{'LPS'}};
% cfg_plot.rois = {{'FC'}};
cfg_plot.ylims = [-7 2; -1 8];
cfg_plot.legendlocs = {'SouthEast','NorthWest'};

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
%   {{'CR2','H2','HSC2','HSI2'},{'CR6','H6','HSC6','HSI6'}},...
%   {{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}};

%cfg_plot.condByTypeByROI = repmat({{{'CR2','HSC2','HSI2'},{'CR6','HSC6','HSI6'}}},size(cfg_plot.rois));

cfg_plot.condByROI = repmat({ana.eventValues},size(cfg_plot.rois));

for r = 1:length(cfg_plot.rois)
  cfg_plot.roi = cfg_plot.rois{r};
  cfg_plot.legendloc = cfg_plot.legendlocs{r};
  cfg_ft.ylim = cfg_plot.ylims(r,:);
  %cfg_plot.conditions = cfg_plot.condByTypeByROI{r};
  cfg_plot.conditions = cfg_plot.condByROI{r};
  
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
cfg_ft.xlim = [-0.2 1.5]; % time
%cfg_ft.xlim = [-0.2 2.0]; % time

cfg_ft.parameter = 'avg';

cfg_plot = [];
cfg_plot.plotTitle = 0;

cfg_plot.is_ga = 1;
cfg_plot.excludeBadSub = 1;

%%%%%%%%%%%%%%%
% Type of plot
%%%%%%%%%%%%%%%

cfg_plot.ftFxn = 'ft_singleplotER';
cfg_plot.rois = {{'FS'},{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
cfg_plot.ylims = [-4.5 2.5; -4.5 2.5; -4.5 2.5; -4.5 2.5; -2 5; -2 5; -2 5];
cfg_plot.rois = {{'FC'}};
cfg_plot.ylims = [-7.5 2];
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
%cfg_ft.colormap = 'hot';
cfg_ft.colorbar = 'no';

%cfg_plot.conditions = {{'all_within_types'}};
%cfg_plot.conditions = {{'all_across_types'}};
%cfg_plot.condMethod = 'pairwise';
%cfg_plot.conditions = {{'HSC2','CR2'},{'HSI2','CR2'},{'HSC2','HSI2'},{'HSC6','CR6'},{'HSI6','CR6'},{'HSC6','HSI6'}}; % {'H2','CR2'}, {'H6','CR6'},
%cfg_plot.conditions = {{'RH','RCR'},{'RHSC','RCR'},{'RHSI','RCR'},{'RHSC','RHSI'}};
cfg_plot.conditions = {{'SC','CR'},{'SC','SI'},{'SI','CR'}};
%cfg_plot.conditions = {{'RHSC','RHSI'}};
%cfg_plot.conditions = {{'FSC','RSSI'}};


cfg_plot.ftFxn = 'ft_topoplotER';
cfg_ft.zlim = [-1.5 1.5]; % volt
cfg_ft.marker = 'on';
%cfg_ft.marker = 'labels';
cfg_ft.markerfontsize = 9;
cfg_ft.comment = 'no';

cfg_plot.roi = {'LAS','RAS'};
%cfg_plot.roi = {'LAS'};
cfg_ft.xlim = [0.3 0.5]; % time

% cfg_plot.roi = {'LPS','RPS'};
% cfg_ft.xlim = [0.5 0.8]; % time

%cfg_plot.roi = {'RAS'};
%cfg_ft.xlim = [1.1 1.9]; % time
% cfg_ft.xlim = [0 1.5]; % time
% cfg_plot.roi = {'all'};

% cfg_plot.subplot = 1;
% cfg_ft.xlim = [0 1.0]; % time
% %cfg_ft.xlim = (0:0.05:1.0); % time

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

cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
cfg_ana.rois = {{'FC'}};
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

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
cfg_ana.conditions = {{'Face','House'}};

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
cfg_plot.line_plots = 0;
%cfg_plot.ylims = [-4 -1; -4 -1; -4 -1; 1 4; 1 4; 1 4; -4 -1; 1 4; 1 4; 1 4];
%cfg_plot.ylims = [-4 -1; 1 4; 1 4];
% cfg_plot.ylims = [-5 -2; -5 -2; -5 -2; 1.5 4.5; 1.5 4.5;];
cfg_plot.ylims = [-5 -2; 1.5 4.5];
%cfg_plot.plot_order = {'CR2','H2','HSC2','HSI2','CR6','H6','HSC6','HSI6'};
%cfg_plot.plot_order = {'RHSC','RHSI','RCR'};
% cfg_plot.plot_order = {'SC','SI','CR'};
%cfg_plot.rename_conditions = {'SC','SI','CR'};
% cfg_plot.xlabel = 'Condition';
% cfg_plot.ylabel = 'Voltage (\muV)';
cfg_plot.xlabel = '';
cfg_plot.ylabel = '';

for r = 1:length(cfg_ana.rois)
  cfg_ana.roi = cfg_ana.rois{r};
  cfg_ft.latency = cfg_ana.latencies(r,:);
  cfg_plot.ylim = cfg_plot.ylims(r,:);
  
  mm_ft_ttestER(cfg_ft,cfg_ana,cfg_plot,exper,ana,files,dirs,data_tla);
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
