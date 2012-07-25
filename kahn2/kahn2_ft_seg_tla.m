exper.name = 'KAHN2matt';

% types of events to save; must be the same as the events in the NS files
exper.eventValues = sort({'CR__','CoTa','InTa'});
exper.eventNames = {'Correct Rejections','Correct Source','Incorrect Source'};
%exper.eventValues = sort({'CR__','CoTa','InTa','PP__','PR__','PTa_','RP__','RR__','RTa_'});
%exper.eventValues = sort({'CR__','PP__','PR__','RP__','RR__'});

if ~isfield(exper,'eventValuesExtra')
  exper.eventValuesExtra = {};
end

% pre- and post-stimulus times to save, in seconds
%exper.prepost = [-1.0 2.0];
exper.prepost = [-0.8 1.5];

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw/sbin must be put in
% dirs.dataroot/ns_raw or egis must be put in dirs.dataroot/ns_egis
exper.eegFileExt = 'set';

% name of the folder to save the FT data in
dirs.saveDirName = 'ft_data';

exper.subjects = {
  'KAHN2 01';
  'KAHN2 03';
  'KAHN2 04';
  'KAHN2 05';
  'KAHN2 07';
  'KAHN2 08';
  'KAHN2 09';
  'KAHN2 10';
  'KAHN2 11';
  'KAHN2 12';
  'KAHN2 13';
  'KAHN2 14';
  'KAHN2 16';
  'KAHN2 17';
  'KAHN2 18';
  'KAHN2 19';
  'KAHN2 20';
  'KAHN2 21';
  'KAHN2 22';
  'KAHN2 23';
  'KAHN2 24';
  'KAHN2 25';
  'KAHN2 26';
  'KAHN2 27';
  'KAHN2 28';
  'KAHN2 29';
  'KAHN2 31';
  'KAHN2 32';
  'KAHN2 34';
  'KAHN2 38';
  'KAHN2 47';
  'KAHN2 94';
  };
% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, the ns_bci directory. They are not
% necessarily the session directory names where the FieldTrip data is saved
% for each subject because of the option to combine sessions. See 'help
% create_ft_struct' for more information.
exper.sessions = {{'_1', '_2'}};

%% set up parameters

% directory where the data to read is located
dirs.dataDir = fullfile(exp.name,exp.name);

% directory to save the FT data; if undefined, set to dirs.dataDir
dirs.saveDirStem = fullfile('KAHN2matt','eeg','eegpp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000));

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

%% Convert the data to FieldTrip structs

ana.segFxn = 'seg2ft';
ana.ftFxn = 'ft_timelockanalysis';

cfg_pp = [];

cfg_proc = [];
cfg_proc.keeptrials = 'no';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'pow');

ana.ftype = 'tla';

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,cfg_proc,exper,dirs,files);

%% if already saved and not yet loaded, load the ft_timelockanalysis files

if ~exist('data_tla','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('data_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000)));
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDirProc,savedFiles(sf).name));
    fprintf('Done.\n');
  end
end

%% fix the labels in the data becaues the EEGLAB data labels are different

for evVal = 1:length(exper.eventValues)
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data.label = data_tla.elec.label;
    end
  end
end

%% Test plots to make sure data look ok

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
% cfg_ft.channel = {'e20'};
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
% cfg_ft.layout = ft_prepare_layout([],data_tla);
% figure
% ft_topoplotER(cfg_ft,data_topo);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Remove the dof field. I'm not sure what degrees of freedom refers to, but the tutorial says to.

for evVal = 1:length(exper.eventValues)
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      if isfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof')
        data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data = rmfield(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data,'dof');
      end
    end
  end
end

%% initialize the analysis struct

ana = struct;

%% set up channel groups

ana.elecGroups = {...
  {'e32','e33','e38','e39','e43','e44','e128'},... % LAI
  {'e1','e114','e115','e120','e121','e122','e125'},... % RAI
  {'e12','e13','e19','e20','e24','e28','e29'},... % LAS
  {'e4','e5','e111','e112','e117','e118','e124'},... % RAS
  {'e37','e42','e52','e53','e54','e60','e61'},... % LPS
  {'e78','e79','e85','e86','e87','e92','e93'},... % RPS
  {'e57','e58','e63','e64','e65','e68','e69'},... %LPI
  {'e89','e90','e94','e95','e96','e99','e100'},... % RPI
  {'e4','e5','e6','e11','e12','e13','e19','e112'},... % Mid-frontal
  {'e52','e53','e60','e61','e59','e66','e67'},... % LPS2
  {'e77','e78','e84','e85','e86','e91','e92'},... % RPS2
  {'e75'},... % P1 (central)
  {'e64'},... % N1 (lateral)
  };
ana.elecGroupsStr = {'LAI','RAI','LAS','RAS','LPS','RPS','LPI','RPI','MF','LPS2','RPS2','P1','N1'};

%% create fields for plotting and analysis (pa); modify for each experiment

ana.pa_ga_plot{1} = sprintf('ga_tla.%s',exper.eventValues{1});
ana.pa_evVal{1} = exper.eventValues(1);
ana.pa_evNames{1} = exper.eventNames(1);
ana.pa_type{1} = '';
for evVal = 2:length(exper.eventValues)
  ana.pa_ga_plot{1} = cat(2,ana.pa_ga_plot{1},sprintf(',ga_tla.%s',exper.eventValues{evVal}));
  ana.pa_evVal{1} = cat(2,ana.pa_evVal{1},exper.eventValues{evVal});
  ana.pa_evNames{1} = cat(2,ana.pa_evNames{1},exper.eventNames{evVal});
end
if length(ana.pa_type) == 1
  ana.pa_evUnique = exper.eventValues;
else
  ana.pa_evUnique = {'1','2','3'};
end

%% decide who to kick out based on trial counts

for evVal = 1:length(exper.eventValues)
  numEv.(exper.eventValues{evVal}) = zeros(length(exper.subjects),length(exper.sessions));
end

for evVal = 1:length(exper.eventValues)
  for sub = 1:length(exper.subjects)
    for ses = 1:length(exper.sessions)
      numEv.(exper.eventValues{evVal})(sub,ses) = size(ft_findcfg(data_tla.(exper.eventValues{evVal}).sub(sub).ses(ses).data.cfg,'trl'),1);
    end
  end
end

fprintf('\t%s\n',sprintf(repmat('\t%s',1,length(exper.eventValues)),exper.eventValues{:}));
for ses = 1:length(exper.sessions)
  for sub = 1:length(exper.subjects)
    subStr = exper.subjects{sub};
    for evVal = 1:length(exper.eventValues)
      subStr = cat(2,subStr,sprintf('\t%d',numEv.(exper.eventValues{evVal})(sub,ses)));
    end
    fprintf('%s\n',subStr);
  end
end

numEv.thresh = 15;
numEv.lowNum = zeros(length(exper.subjects),length(exper.sessions));
for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    for evVal = 1:length(exper.eventValues)
      if numEv.(exper.eventValues{evVal})(sub,ses) < numEv.thresh
        numEv.lowNum(sub,ses) = 1;
      end
    end
  end
end
numEv.lowNum = logical(numEv.lowNum);

% who we're rejecting
fprintf('Subjects with low trial counts:\n');
exper.subjects(numEv.lowNum)

exper.badBehSub = {};

fprintf('Subjects with bad behavior:\n');
exper.badBehSub

exper.badSub = zeros(length(exper.subjects),length(exper.sessions));
for ses = 1:length(exper.sessions)
  exper.badSub(:,ses) = sign(sum([ismember(exper.subjects,exper.badBehSub) numEv.lowNum(:,ses)],2));
end
exper.badSub = logical(exper.badSub);

fprintf('Number of events included in EEG analyses (%d subjects; threshold: %d events):\n',sum(~exper.badSub),numEv.thresh);
for evVal = 1:length(exper.eventValues)
  numEv.meanTrials = mean(numEv.(exper.eventValues{evVal})(~exper.badSub));
  numEv.sem = std(numEv.(exper.eventValues{evVal})(~exper.badSub),0,1)/sqrt(sum(~exper.badSub));
  numEv.sd = std(numEv.(exper.eventValues{evVal})(~exper.badSub),0,1);
  fprintf('%s:\tM=%.3f,\tSEM=%.3f,\tSD=%.3f\n',exper.eventValues{evVal},numEv.meanTrials,numEv.sem,numEv.sd);
end

%% set up strings to put in ft_timelockgrandaverage

count = 1;
for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sessions)
    if exper.badSub(sub,ses)
      fprintf('BadSubject: %s\n',exper.subjects{sub});
      continue
    else
      if count == 1
        for evVal = 1:length(exper.eventValues)
          ana.str.(exper.eventValues{evVal}){ses} = sprintf('data_tla.%s.sub(%d).ses(%d).data',exper.eventValues{evVal},sub,ses);
        end
        count = count + 1;
      elseif count > 1
        for evVal = 1:length(exper.eventValues)
          ana.str.(exper.eventValues{evVal}){ses} = sprintf('%s,data_tla.%s.sub(%d).ses(%d).data',ana.str.(exper.eventValues{evVal}){ses},exper.eventValues{evVal},sub,ses);
        end
        count = count + 1;
      end
    end
  end
end

%% get the grand average and add the electrode locations again

cfg_ft = [];
cfg_ft.keepindividual = 'yes';
for ses = 1:length(exper.sessions)
  for evVal = 1:length(exper.eventValues)
    tic
    fprintf('Running ft_timelockgrandaverage on %s...',exper.eventValues{evVal});
    ga_tla.(exper.eventValues{evVal})(ses) = eval(sprintf('ft_timelockgrandaverage(cfg_ft,%s);',ana.str.(exper.eventValues{evVal}){ses}));
    fprintf('Done.\n');
    toc
  end
end
ga_tla.elec = data_tla.elec;

%% save

saveFile = fullfile(dirs.saveDirProc,sprintf('ga_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000));
if ~exist(saveFile,'file')
  save(saveFile,'ga_tla');
end

%% if already saved and not yet loaded, load the ft_timelockgrandaverage files

if ~exist('ga_tla','var')
  savedFiles = dir(fullfile(dirs.saveDirProc,sprintf('ga_tla_%d_%d.mat',exper.prepost(1)*1000,exper.prepost(2)*1000)));
  for sf = 1:length(savedFiles)
    fprintf('Loading %s...',savedFiles(sf).name);
    load(fullfile(dirs.saveDirProc,savedFiles(sf).name));
    fprintf('Done.\n');
  end
end

%% plot the conditions - simple

% ft_singleplotER seems to want individual subject data and not grand
% average data. we're about to give it GA data. when it gets this, it will
% give 2 warnings, but it should be ok to proceed

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.xlim = [-.2 1.0];
cfg_ft.linewidth = 2;
cfg_ft.graphcolor = 'rbk';
cfg_ft.linestyle = {'-','--','-.'};

for pa = 1:length(ana.pa_evVal)
  % FN400
  cfg_ft.ylim = [-4 3];
  cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS'})});
  cfg_ft.channel = 'e20';
  figure
  eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana.pa_ga_plot{pa}));
  legend(strrep(ana.pa_evVal{pa},'_',''),'Location','SouthEast');
  % P O/N
  cfg_ft.ylim = [-1 6];
  cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LPS','RPS'})});
  figure
  eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana.pa_ga_plot{pa}));
  legend(strrep(ana.pa_evVal{pa},'_',''),'Location','NorthEast');
end

%% subplot with all subjects' ERPs

ses = 1;
cfg_plot = [];
%cfg_plot.roi = {'LAS','RAS'};
cfg_plot.roi = {'RAS'};
%cfg_plot.roi = {'LPS','RPS'};
%cfg_plot.roi = {'LPS'};
cfg_plot.chan_str = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
cfg_plot.chan = str2double(strrep(cfg_plot.chan_str,'e',''));
cfg_plot.numCols = 5;
cfg_plot.excludeBadSub = 0;
if cfg_plot.excludeBadSub
  cfg_plot.numRows = ceil((length(exper.subjects) - sum(exper.badSub)) / cfg_plot.numCols);
else
  cfg_plot.numRows = ceil((length(exper.subjects)) / cfg_plot.numCols);
end
cfg_plot.graphcolor = 'rbkgcmy';

for pa = 1:length(ana.pa_evVal)
  figure
  for sub = 1:length(exper.subjects)
    subplot(cfg_plot.numRows,cfg_plot.numCols,sub);
    evStr = [];
    for evVal = 1:length(ana.pa_evVal{pa})
      evStr = cat(2,evStr,sprintf(' %s:%d',strrep(ana.pa_evVal{pa}{evVal},'_',''),numEv.(ana.pa_evVal{pa}{evVal})(sub)));
      plot(data_tla.(ana.pa_evVal{pa}{evVal}).sub(sub).ses(ses).data.time,mean(data_tla.(ana.pa_evVal{pa}{evVal}).sub(sub).ses(ses).data.avg(cfg_plot.chan,:),1),cfg_plot.graphcolor(evVal));
      hold on
    end
    if exper.badSub(sub,ses)
      fprintf('BadSubject: %s\n',exper.subjects{sub});
      title(sprintf('*BAD*: %s;%s',exper.subjects{sub},evStr));
    else
      title(sprintf('%s;%s',exper.subjects{sub},evStr));
    end
    %ylim([0 1.9e-13]);
    axis([-0.2 1.0 -10 10]);
    hold off
  end
  subplot(cfg_plot.numRows,cfg_plot.numCols,length(exper.subjects) + 1);
  text(0.5,0.9,sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}),'color','k');
  ycoord = 0.7;
  for evVal = 1:length(ana.pa_evVal{pa})
    text(0.5,ycoord,strrep(ana.pa_evVal{pa}{evVal},'_',''),'color',cfg_plot.graphcolor(evVal));
    ycoord = ycoord - 0.2;
  end
  axis off
end

%% plot the conditions

cfg_ft = [];
cfg_ft.showlabels = 'yes';
cfg_ft.fontsize = 9;
cfg_ft.elec = ga_tla.elec;
cfg_ft.linewidth = 2;
cfg_ft.graphcolor = 'rbk';
cfg_ft.linestyle = {'-','--','-.'};

cfg_plot = [];
%cfg_plot.fillcolor = [.8,.8,.8];
%cfg_plot.fillalpha = 0.3;
%cfg_plot.filledge = 'none';
cfg_plot.minMaxVolt = [-5 2; 0 7];
cfg_plot.latency = [-0.2 1.0];
cfg_plot.plotLegend = 1;
cfg_plot.plot_rois = {{'MF'},{'LAS'},{'RAS'},{'LAS','RAS'},{'LPS'},{'RPS'},{'LPS','RPS'}};
%cfg_plot.plot_rois = {{'LAS'},{'RAS'},{'LPS'},{'RPS'}};

for pa = 1:length(ana.pa_ga_plot)
  for i = 1:length(cfg_plot.plot_rois)
  cfg_plot.roi = cfg_plot.plot_rois{i};
  cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
  figure
  %hold on
  eval(sprintf('ft_singleplotER(cfg_ft,%s);',ana.pa_ga_plot{pa}));
  hold on
  plot([cfg_plot.latency(1) cfg_plot.latency(2)],[0 0],'k--'); % horizontal
  if ismember('LAS',cfg_plot.roi) || ismember('RAS',cfg_plot.roi)
    plot([0 0],[cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)],'k--'); % vertical
    if cfg_plot.plotLegend
      legend(strrep(ana.pa_evVal{pa},'_',''),'Location','SouthEast');
    end
    
    cfg_plot.vert_latency = [0.3 0.5];
    %h = fill([cfg_plot.vert_latency(1),cfg_plot.vert_latency(2),cfg_plot.vert_latency(2),cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(1,1),cfg_plot.minMaxVolt(1,1),cfg_plot.minMaxVolt(1,2),cfg_plot.minMaxVolt(1,2)],fillcolor);
    %set(h,'FaceAlpha',fillalpha);
    %set(h,'EdgeColor',filledge)
    plot([cfg_plot.vert_latency(1) cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)],'k'); % vertical
    plot([cfg_plot.vert_latency(2) cfg_plot.vert_latency(2)],[cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)],'k'); % vertical
    
    axis([cfg_plot.latency(1) cfg_plot.latency(2) cfg_plot.minMaxVolt(1,1) cfg_plot.minMaxVolt(1,2)])
  elseif ismember('LPS',cfg_plot.roi) || ismember('RPS',cfg_plot.roi)
    plot([0 0],[cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)],'k--'); % vertical
    if cfg_plot.plotLegend
      legend(strrep(ana.pa_evVal{pa},'_',''),'Location','NorthEast');
    end
    
    cfg_plot.vert_latency = [0.5 0.8];
    %h = fill([cfg_plot.vert_latency(1),cfg_plot.vert_latency(2),cfg_plot.vert_latency(2),cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(2,1),cfg_plot.minMaxVolt(2,1),cfg_plot.minMaxVolt(2,2),cfg_plot.minMaxVolt(2,2)],fillcolor);
    %set(h,'FaceAlpha',fillalpha);
    %set(h,'EdgeColor',filledge);
    plot([cfg_plot.vert_latency(1) cfg_plot.vert_latency(1)],[cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)],'k'); % vertical
    plot([cfg_plot.vert_latency(2) cfg_plot.vert_latency(2)],[cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)],'k'); % vertical
    
    axis([cfg_plot.latency(1) cfg_plot.latency(2) cfg_plot.minMaxVolt(2,1) cfg_plot.minMaxVolt(2,2)])
  end
  xlabel('Time (s)');
  ylabel('Voltage (\muV)');
  set(gcf,'Name',sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}))
  %title(sprintf(repmat('%s ',1,length(cfg_plot.roi)),cfg_plot.roi{:}))
  %axis ij % negative up
  publishfig(gcf,1);
  if cfg_plot.plotLegend
    cfg_plot.legend_str = '_legend';
  else
    cfg_plot.legend_str = '';
  end
  if files.saveFigs
    if ~isempty(ana.pa_type{pa})
      cfg_plot.figfilename = sprintf('tla_erp_ga_%s_%s%d_%d%s.%s',ana.pa_type{pa},sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_plot.latency(1)*1000,cfg_plot.latency(2)*1000,cfg_plot.legend_str,files.figFileExt);
    else
      cfg_plot.figfilename = sprintf('tla_erp_ga_%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_plot.latency(1)*1000,cfg_plot.latency(2)*1000,cfg_plot.legend_str,files.figFileExt);
    end
    print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
  end
  end
end

%% Plot average time contrast topoplots

cfg_plot = [];
cfg_plot.minMaxVolt = [-1 1];
cfg_plot.plotTitle = 1;
cfg_plot.plotColorbar = 1;
if cfg_plot.plotColorbar
  cfg_plot.colorbar_str = '_cb';
else
  cfg_plot.colorbar_str = '';
end

% define which regions to highlight in the plot
cfg_plot.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of rois
cfg_plot.latencies = [0.3 0.5; 0.5 0.8];

cfg_ft = [];
cfg_ft.interactive = 'no';
%cfg_ft.colorbar = 'yes';
cfg_ft.layout = ft_prepare_layout([],ga_tla);
cfg_ft.highlight = 'on';
cfg_ft.highlightsize = 10;
%cfg_ft.comment = 'xlim';
%cfg_ft.commentpos = 'title';
cfg_ft.comment = 'no';
%cfg_ft.colormap = colormap('jet'); % default; blue to red
%cfg_ft.colormap = colormap('hot'); % dark to light; better for b&w printers
cfg_ft.zlim = [cfg_plot.minMaxVolt(1) cfg_plot.minMaxVolt(2)];

% initialize for storing the contrast topoplots
cont_topo = [];

for pa = 1:length(ana.pa_evVal)
  % find all the pairwise event combinations for these contrasts
  cfg_plot.comboEv = nchoosek(ana.pa_evVal{pa},2);
  % get the names for plotting
  cfg_plot.comboName = nchoosek(ana.pa_evNames{pa},2);
  
  for r = 1:length(cfg_plot.rois)
    cfg_ft.xlim = cfg_plot.latencies(r,:);
    cfg_plot.roi = cfg_plot.rois{r};
    cfg_ft.highlightchannel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_plot.roi)});
    for c = 1:length(cfg_plot.comboEv)
      cfg_plot.cond = cfg_plot.comboEv(c,:);
      cfg_plot.condNames = cfg_plot.comboName(c,:);
      
      % create contrast
      cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})) = ga_tla.(cfg_plot.cond{1});
      if isfield(ga_tla.(cfg_plot.cond{1}),'avg')
        cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})).avg = ga_tla.(cfg_plot.cond{1}).avg - ga_tla.(cfg_plot.cond{2}).avg;
      end
      if isfield(ga_tla.(cfg_plot.cond{1}),'individual')
        cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})).individual = ga_tla.(cfg_plot.cond{1}).individual - ga_tla.(cfg_plot.cond{2}).individual;
      end
      
      % make a plot
      figure
      ft_topoplotER(cfg_ft,cont_topo.(sprintf('%svs%s',cfg_plot.cond{1},cfg_plot.cond{2})));
      set(gcf,'Name',sprintf('%s - %s',cfg_plot.cond{1},cfg_plot.cond{2}))
      if cfg_plot.plotColorbar
        h = colorbar;
        set(get(h,'YLabel'),'string','Voltage (\muV)');
      end
      if cfg_plot.plotTitle
        title(sprintf('%s - %s',cfg_plot.condNames{1},cfg_plot.condNames{2}));
      end
      publishfig(gcf,0);
      if files.saveFigs
        if ~isempty(ana.pa_type{pa})
          cfg_plot.figfilename = sprintf('tla_topo_ga_%s_%s%s%d_%d%s.%s',ana.pa_type{pa},sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
        else
          cfg_plot.figfilename = sprintf('tla_topo_ga_%s%s%d_%d%s.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_plot.roi)),cfg_plot.roi{:}),cfg_ft.xlim(1)*1000,cfg_ft.xlim(2)*1000,cfg_plot.colorbar_str,files.figFileExt);
        end
        print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
      end
    end % c
  end % r
end % pa

%% descriptive statistics: ttest

cfg_ana = [];
% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];
% exclude the bad subjects from the subject count
cfg_ana.numSub = length(exper.subjects) - sum(exper.badSub);

cfg_ana.make_plots = 1;

% set parameters for the statistical test
cfg_ft = [];
cfg_ft.avgovertime = 'yes';
cfg_ft.avgoverchan = 'yes';
cfg_ft.parameter   = 'individual';
cfg_ft.method      = 'analytic';
cfg_ft.statistic   = 'depsamplesT';
cfg_ft.alpha       = 0.05;
cfg_ft.correctm    = 'no';

% set up the design for the statistical test
cfg_ft.design(1,1:2*cfg_ana.numSub) = [ones(1,cfg_ana.numSub) 2*ones(1,cfg_ana.numSub)];
cfg_ft.design(2,1:2*cfg_ana.numSub) = [1:cfg_ana.numSub 1:cfg_ana.numSub];
cfg_ft.ivar = 1; % the 1st row in cfg_ft.design contains the independent variable
cfg_ft.uvar = 2; % the 2nd row in cfg_ft.design contains the units of observation (subject number)

for pa = 1:length(ana.pa_evVal)
  % find all the pairwise event combinations for t-tests
  cfg_ana.comboEv = nchoosek(ana.pa_evVal{pa},2);
  
  % do t-tests at each ROI
  for r = 1:length(cfg_ana.rois)
    cfg_ana.roi = cfg_ana.rois{r};
    cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi)});
    cfg_ft.latency = cfg_ana.latencies(r,:);
    
    % get times, voltage, SEM
    cfg_ana.chan = unique(str2double(strrep(cfg_ft.channel,'e','')));
    for evVal = 1:length(ana.pa_evVal{pa})
      cfg_ana.timesel.(ana.pa_evVal{pa}{evVal}) = find(ga_tla.(ana.pa_evVal{pa}{evVal}).time >= cfg_ft.latency(1) & ga_tla.(ana.pa_evVal{pa}{evVal}).time <= cfg_ft.latency(2));
      cfg_ana.values.(ana.pa_evVal{pa}{evVal}) = mean(mean(ga_tla.(ana.pa_evVal{pa}{evVal}).individual(:,cfg_ana.chan,cfg_ana.timesel.(ana.pa_evVal{pa}{evVal})),2),3);
      cfg_ana.sem.(ana.pa_evVal{pa}{evVal}) = std(cfg_ana.values.(ana.pa_evVal{pa}{evVal}))/sqrt(length(cfg_ana.values.(ana.pa_evVal{pa}{evVal})));
    end
    
    fprintf('\n-------------------------------------\n');
    fprintf('%s: ROI:%s; Times: %.3f--%.3f\n',ana.pa_type{pa},sprintf(repmat(' %s',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ana.latencies(r,1),cfg_ana.latencies(r,2));
    fprintf('-------------------------------------\n\n');

    % print GA and sub avg voltages
    cfg_ana.goodSub = exper.subjects(~exper.badSub);
    % ga
    fprintf('\t%s\n',sprintf(repmat('\t%s',1,length(ana.pa_evVal{pa})),ana.pa_evVal{pa}{:}));
    gaStr = sprintf('GA\t');
    for evVal = 1:length(ana.pa_evVal{pa})
      gaStr = cat(2,gaStr,sprintf('\t%.3f',mean(cfg_ana.values.(ana.pa_evVal{pa}{evVal}))));
    end
    fprintf('%s\n',gaStr);
    % sub avg
    fprintf('Subject\t%s\n',sprintf(repmat('\t%s',1,length(ana.pa_evVal{pa})),ana.pa_evVal{pa}{:}));
    for sub = 1:length(cfg_ana.goodSub)
      subStr = exper.subjects{sub};
      for evVal = 1:length(ana.pa_evVal{pa})
        subStr = cat(2,subStr,sprintf('\t%.3f',cfg_ana.values.(ana.pa_evVal{pa}{evVal})(sub)));
      end
      fprintf('%s\n',subStr);
    end
    
    % do the t-test
    for c = 1:length(cfg_ana.comboEv)
      cfg_ana.cond = cfg_ana.comboEv(c,:);
      cfg_ana.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})) = ft_timelockstatistics(cfg_ft,ga_tla.(cfg_ana.cond{1}),ga_tla.(cfg_ana.cond{2}));
      % % matlab dependent samples ttest
      % sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2}) = cfg_ana.values.(cfg_ana.cond{1}) - cfg_ana.values.(cfg_ana.cond{2});
      % [h,p,ci,stats] = ttest(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2}), 0, 0.05); % H0: mean = 0, alpha 0.05
    end
    % print the results
    for c = 1:length(cfg_ana.comboEv)
      cfg_ana.cond = cfg_ana.comboEv(c,:);
      if ~isempty(ana.pa_type{pa})
        fprintf('%s: %s (M=%.3f; SEM=%.3f) vs %s (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',ana.pa_type{pa},cfg_ana.cond{1},mean(cfg_ana.values.(cfg_ana.cond{1}),1),cfg_ana.sem.(cfg_ana.cond{1}),cfg_ana.cond{2},mean(cfg_ana.values.(cfg_ana.cond{2}),1),cfg_ana.sem.(cfg_ana.cond{2}),cfg_ana.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})).df,cfg_ana.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})).stat,cfg_ana.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})).prob);
      else
        fprintf('%s (M=%.3f; SEM=%.3f) vs %s (M=%.3f; SEM=%.3f):\tt(%d)=%.4f, p=%.10f\n',cfg_ana.cond{1},mean(cfg_ana.values.(cfg_ana.cond{1}),1),cfg_ana.sem.(cfg_ana.cond{1}),cfg_ana.cond{2},mean(cfg_ana.values.(cfg_ana.cond{2}),1),cfg_ana.sem.(cfg_ana.cond{2}),cfg_ana.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})).df,cfg_ana.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})).stat,cfg_ana.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})).prob);
      end
    end
    
    if cfg_ana.make_plots == 1
      % plot with all subjects on it to see the effect in each subject
      for c = 1:length(cfg_ana.comboEv)
        cfg_ana.cond = cfg_ana.comboEv(c,:);
        figure
        plot([cfg_ana.values.(cfg_ana.cond{1}) cfg_ana.values.(cfg_ana.cond{2})]','o-');
        xlim([0.5 2.5])
        title(sprintf('%.1fs--%.1fs, %s',cfg_ft.latency(1),cfg_ft.latency(2),sprintf(repmat('%s ',1,length(cfg_ana.roi)),cfg_ana.roi{:})));
        ylabel('Microvolts (\muV)');
        set(gca,'XTickLabel',{'',cfg_ana.cond{1},'',cfg_ana.cond{2},''})
      end
      
      % do the mean amplitude line plots
      cfg_plot = [];
      
      if ismember('LAS',cfg_ana.roi) || ismember('RAS',cfg_ana.roi)
        cfg_plot.conditions = ana.pa_evVal(pa);
        %cfg_plot.yminmax = [-5 -1];
      elseif ismember('LPS',cfg_ana.roi) || ismember('RPS',cfg_ana.roi)
        cfg_plot.conditions = ana.pa_evVal(pa);
        %cfg_plot.yminmax = [2 6];
      end
      cfg_plot.yminmax = eval(sprintf('[floor(min([%s])) ceil(max([%s]))]',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(ana.pa_evVal{pa})),ana.pa_evVal{pa}{:}),sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(ana.pa_evVal{pa})),ana.pa_evVal{pa}{:})));
      
      % set up how the lines will look
      cfg_plot.linewidth = 2;
      cfg_plot.marksize = 10;
      cfg_plot.linespec = 'k--o';
      cfg_plot.markcolor = 'w';
      cfg_plot.errwidth = 1;
      cfg_plot.errBarEndMarkerInd = [4 5 7 8];
      cfg_plot.removeErrBarEnds = 1;
      
      for cnd = 1:length(cfg_plot.conditions)
        cfg_plot.cond = cfg_plot.conditions{cnd};
        
        figure
        % plot the lines
        eval(sprintf('plot([%s],cfg_plot.linespec,''LineWidth'',cfg_plot.linewidth);',sprintf(repmat('mean(cfg_ana.values.%s,1) ',1,length(cfg_plot.cond)),cfg_plot.cond{:})));
        hold on
        for c = 1:length(cfg_plot.cond)
          % errorbars
          h = errorbar(c,mean(cfg_ana.values.(cfg_plot.cond{c}),1),cfg_ana.sem.(cfg_plot.cond{c}),cfg_plot.linespec,'LineWidth',cfg_plot.errwidth);
          % remove errorbar ends
          if cfg_plot.removeErrBarEnds
            chil = get(h,'Children');
            xdata = get(chil(2),'XData');
            ydata = get(chil(2),'YData');
            xdata(cfg_plot.errBarEndMarkerInd) = NaN;
            ydata(cfg_plot.errBarEndMarkerInd) = NaN;
            set(chil(2),'XData',xdata);
            set(chil(2),'YData',ydata);
            set(h,'Children',chil);
          end
          % plot the markers
          plot(c,mean(cfg_ana.values.(cfg_plot.cond{c}),1),cfg_plot.linespec,'LineWidth',cfg_plot.linewidth,'MarkerSize',cfg_plot.marksize,'MarkerFaceColor',cfg_plot.markcolor);
        end
        hold off
        
        % make it look good
        axis([.5 (length(cfg_plot.cond) + .5) cfg_plot.yminmax(1) cfg_plot.yminmax(2)])
        xlabel('Condition');
        ylabel('Voltage (\muV)');
        set(gca,'XTick',(1:length(cfg_plot.cond)))
        set(gca,'XTickLabel',strrep(cfg_plot.cond,'_',''))
        set(gca,'YTick',(cfg_plot.yminmax(1):.5:cfg_plot.yminmax(2)))
        axis square
        publishfig(gcf,0);
        if files.saveFigs
          if ~isempty(ana.pa_type{pa})
            cfg_plot.figfilename = sprintf('tla_line_ga_%s_%s%s%d_%d.%s',ana.pa_type{pa},sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000,files.figFileExt);
          else
            cfg_plot.figfilename = sprintf('tla_line_ga_%s%s%d_%d.%s',sprintf(repmat('%s_',1,length(cfg_plot.cond)),cfg_plot.cond{:}),sprintf(repmat('%s_',1,length(cfg_ana.roi)),cfg_ana.roi{:}),cfg_ft.latency(1)*1000,cfg_ft.latency(2)*1000,files.figFileExt);
          end
          print(gcf,files.figPrintFormat,fullfile(dirs.saveDirFigs,cfg_plot.figfilename));
        end
      end % cnd; line plots
    end % cfg_ana.make_plots
  end % r
end % pa

%% 2-way ANOVA: Hemisphere x Condition

cfg_ana = [];
cfg_ana.alpha = 0.05;
cfg_ana.showtable = 1;
cfg_ana.printTable_tex = 1;

% define which regions to average across for the test
cfg_ana.rois = {{'LAS','RAS'},{'LPS','RPS'}};
% define the times that correspond to each set of ROIs
cfg_ana.latencies = [0.3 0.5; 0.5 0.8];

for pa = 1:length(ana.pa_evVal)
  % define the conditions that correspond to each set of ROIs
  cfg_ana.conditions = {ana.pa_evVal{pa}, ana.pa_evVal{pa}};
  
  for r = 1:length(cfg_ana.rois)
    cfg_ana.roi = cfg_ana.rois{r};
    cfg_ana.latency = cfg_ana.latencies(r,:);
    cfg_ana.cond = cfg_ana.conditions{r};
    
    % need to do Greenhouse-Geisser or Huynh-Feldt correction?
    if length(cfg_ana.roi) > 2 || length(cfg_ana.cond) > 2
      cfg_ana.calcGGHF = 1;
    else
      cfg_ana.calcGGHF = 0;
    end
    
    fprintf('====================== RMAOV2_mod =========================\n');
    if ~isempty(ana.pa_type{pa})
      fprintf('%s: IV1: ROI (%d;%s), IV2: Condition (%d;%s)\n',ana.pa_type{pa},length(cfg_ana.roi),sprintf(repmat(' %s',1,length(cfg_ana.roi)),cfg_ana.roi{:}),length(cfg_ana.cond),sprintf(repmat(' %s',1,length(cfg_ana.cond)),cfg_ana.cond{:}));
    else
      fprintf('IV1: ROI (%d;%s), IV2: Condition (%d;%s)\n',length(cfg_ana.roi),sprintf(repmat(' %s',1,length(cfg_ana.roi)),cfg_ana.roi{:}),length(cfg_ana.cond),sprintf(repmat(' %s',1,length(cfg_ana.cond)),cfg_ana.cond{:}));
    end
    
    cfg_ana.anovamat = [];
    for reg = 1:length(cfg_ana.roi)
      cfg_ana.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,cfg_ana.roi(reg))});
      cfg_ana.chan = unique(str2double(strrep(cfg_ana.channel,'e','')));
      
      for c = 1:length(cfg_ana.cond)
        cfg_ana.timesel.(cfg_ana.cond{c}) = find(ga_tla.(cfg_ana.cond{c}).time >= cfg_ana.latency(1) & ga_tla.(cfg_ana.cond{c}).time <= cfg_ana.latency(2));
        cfg_ana.values.(cfg_ana.cond{c}) = mean(mean(ga_tla.(cfg_ana.cond{c}).individual(:,cfg_ana.chan,cfg_ana.timesel.(cfg_ana.cond{c})),2),3);
        
        %index = reg + (c - 1) + (reg - 1);
        for sub = 1:size(cfg_ana.values.(cfg_ana.cond{c}),1)
          
          % format for RMAOV2_mod: [data roi cond subNum]
          cfg_ana.anovamat = cat(1,cfg_ana.anovamat,[cfg_ana.values.(cfg_ana.cond{c})(sub) reg c sub]);
        end
      end
    end
    [cfg_ana.p] = RMAOV2_mod(cfg_ana.anovamat,cfg_ana.alpha,cfg_ana.showtable,cfg_ana.calcGGHF,cfg_ana.printTable_tex);
  end % r
end % pa

%% cluster statistics

cfg_ft = [];
%start and stop time of analysis, in sec
cfg_ft.latency = [0.0 1.0];
cfg_ft.avgovertime = 'no';

cfg_ft.elec = ga_tla.elec;
cfg_ft.channel = 'all';
%cfg_ft.channel = cat(2,ana.elecGroups{ismember(ana.elecGroupsStr,{'LAS','RAS'})});
%cfg_ft.avgoverchan = 'yes';

cfg_ft.parameter = 'individual';

% use the Monte Carlo Method to calculate the significance probability
cfg_ft.method = 'montecarlo';
cfg_ft.correctm = 'cluster';
% alpha level of the sample-specific test statistic that will be used for
% thresholding
cfg_ft.clusteralpha = 0.05;
% test statistic that will be evaluated under the permutation distribution
cfg_ft.clusterstatistic = 'maxsum';
% minimum number of neighborhood channels that is required for a selected
% sample to be included in the clustering algorithm (default = 0)
cfg_ft.minnbchan = 2;
% alpha level of the permutation test
cfg_ft.alpha = 0.025;
% number of draws from the permutation distribution, should be 1000
cfg_ft.numrandomization = 1000;

% set the unit and independent variables
cfg_ft.uvar = 1; % row of design matrix containing subject numbers (DV)
cfg_ft.ivar = 2; % row of design matrix containing conditions (IV)

cfg_ana = [];
cfg_ana.numConds = 2;
% use the dependent samples T-statistic as a measure to evaluate the
% effect at the sample level
if cfg_ana.numConds == 2
  cfg_ft.statistic = 'depsamplesT';
  % test tails: -1 = left tail, 0 = two tail, 1 = right tail?
  cfg_ft.tail = 0;
  cfg_ft.clustertail = 0;
elseif cfg_ana.numConds > 2
  cfg_ft.statistic = 'depsamplesF';
  % test tails: -1 = left tail, 0 = two tail, 1 = right tail?
  cfg_ft.tail = 1;
  cfg_ft.clustertail = 1;
end

% make the design matrix
cfg_ana.numSub = length(exper.subjects) - sum(exper.badSub);
cfg_ft.design = zeros(2,cfg_ana.numSub*cfg_ana.numConds);
for i = 1:cfg_ana.numConds
  for j = 1:cfg_ana.numSub
    cfg_ft.design(1,1+((i - 1)*cfg_ana.numSub) + (j - 1)) = j; % subject #s
    cfg_ft.design(2,1+((i - 1)*cfg_ana.numSub) + (j - 1)) = i; % condition #s
  end
end

for pa = 1:length(ana.pa_evVal)
  % find all the pairwise event combinations for t-tests
  cfg_ana.comboEv = nchoosek(ana.pa_evVal{pa},2);
  
  % run the nonparametric cluster statistics
  stat_clus = [];
  for c = 1:length(cfg_ana.comboEv)
    cfg_ana.cond = cfg_ana.comboEv(c,:);
    stat_clus.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})) = ft_timelockstatistics(cfg_ft,ga_tla.(cfg_ana.cond{1}),ga_tla.(cfg_ana.cond{2}));
  end
end

% save(fullfile(dirs.saveDirProc,'stat_clus.mat'),'stat_clus');
% load(fullfile(dirs.saveDirProc,'stat_clus.mat'));

% % run the nonparametric cluster statistics
% stat_clus.(sprintf('%svs%svs%s',cfg_ana.cond{1},cfg_ana.cond{2},cfg_ana.cond{3})) = ft_timelockstatistics(cfg_ft,ga_tla.(cfg_ana.cond{1}),ga_tla.(cfg_ana.cond{2}),ga_tla.(cfg_ana.cond{3}));

%% plot the cluster statistics

cfg_ft = [];
% p-val markers; default ['*','x','+','o','.'], p < [0.01 0.05 0.1 0.2 0.3]
cfg_ft.highlightsymbolseries = ['*','+','.','.','.'];
cfg_ft.highlightcolorpos = [0.5 0 1];
cfg_ft.highlightcolorneg = [0 0.5 0];
cfg_ft.elec = ga_tla.elec;
cfg_ft.contournum = 0;
cfg_ft.emarker = '.';
cfg_ft.alpha  = 0.05;
cfg_ft.parameter = 'stat';
cfg_ft.zlim = [-5 5];
for c = 1:length(cfg_ana.comboEv)
  cfg_ana.cond = cfg_ana.comboEv(c,:);
  ft_clusterplot(cfg_ft,stat_clus.(sprintf('%svs%s',cfg_ana.cond{1},cfg_ana.cond{2})));
end
