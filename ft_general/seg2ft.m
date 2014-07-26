function [ft_raw,badChanAllSes,badEvEvVals,artifacts,trialinfo_allEv] = seg2ft(dataroot,subject,session,sesNum_orig,eventValue,eventValue_orig,prepost,elecfile,ana,exper,dirs)
%SEG2FT: take segmented EEG data and put it in FieldTrip format
%
% [ft_raw,badChan,badEv,artifacts,trialinfo_allEv] = seg2ft(dataroot,subject,session,sesNum_orig,eventValue,eventValue_orig,prepost,elecfile,ana,exper,dirs)
%
% Output:
%   ft_raw  = struct with one field for each event value
%   badChan = bad channel information
%   badEv   = bad event information
%   artifacts = artifact information
%   trialinfo_allEv = full set of trialinfo, includes trials rejected for
%                     artifacts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP:
%
% Export Net Station data as either EGIS (129 channels, where the final
% channel [Cz] is the reference; extension: 'egis') or NS Simple Binary
% (extension: 'raw' or 'sbin') calibrated, including the reference channel.
% These options are in the File Export tool. create_ft_struct, the function
% that calls seg2ft, expects EGIS files to be stored in a 'ns_egis'
% directory at the level of dirs.dataDir. If using raw files, they should
% be in 'ns_raw' instead.
%
% This function can deal with event values that have zero events. It is
% probably better to check on your event count and exclude those subjects
% with zero events for any of the event values you're trying to keep before
% trying to run this script. Nonetheless, it will insert an empty event
% entry for an empty eventValue, and the subjects will be excluded when
% using mm_threshSubs.m.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFACT INFORMATION:
%
% ana.artifact.type can be 'none', 'nsAuto', 'zeroVar', 'preRejManual',
% 'ftManual', 'ftICA', 'badChanManual', 'badChanEP', 'rmBadChan', 'nsClassic'.
% Default is {'none'}.
%
% It can be one of those strings or a cell array of multiple strings (e.g.,
% {'nsAuto','preRejManual'} to do both Net Station artifact rejection and
% FieldTrip manual ("visual") rejection). 'ftICA' also includes manual
% rejection after manually assessing components. Though it is not prohibited,
% 'nsAuto' and 'zeroVar' should not be used together because there is
% probably no good reason to find bad trials with both EP Toolkit and NS.
%
% 'nsAuto', 'zeroVar', and 'preRejManual' are processed first, then 'ftManual',
% then 'ftICA'. Subquent processing will not include earlier rejected artifacts.
% Note: any FT artifact processing requires manual intervention (as does 'preRejManual'),
% while 'nsAuto' and 'zeroVar' artifact processing does not. 'preRejManual' is for
% inspecting the artifacts that have been previously identified by other
% software (NS, EP Toolkit, etc.). After manual inspection, both 'preRejManual'
% and 'ftManual' give the option to repair individual channels (for all
% trials) using FT_CHANNELREPAIR, so be sure to keep track of any channels
% that you want to repair (i.e., instead of rejecting them as artifacts).
%
% If using NS artifacts ('nsAuto'), this function expects to find a Net Station
% segment info file with a .bci extension; this contains artifact information.
% It is exported from Net Station using the File Export tool. To set up the tool,
% export format is metadata and check the segment information option. Run
% the tool on the file that was exported to egis/raw (i.e., the baseline
% correction file or the average rereference file). The bci files should be
% stored in a 'ns_bci' directory in fullfile(dirs.dataroot,dirs.dataDir).
%
% Rejecting trials with zero variance ('zeroVar') should be used when using
% bad trial detection done by EP Toolkit because it flattens all channels
% of any bad trials.
%
% If 'ftManual', a visualization of all channels for each event will appear,
% where each trial is shown one-by-one.
%
% If 'ftICA', ICA will run on all trials across all event values.
% Individual components can be rejected after this.  Finally, a
% visualization of all channels for each event will appear, where each
% trial is shown one-by-one.
%
% 'badChanManual' requires a tab-delimited file titled
% [exper.name,'_badChan.txt'] to reside in
% fullfile(dirs.dataroot,dirs.dataDir). The three tab columns are subject
% name (e.g., EXPER001), session name (e.g., session_0), and bad channel
% numbers listed as integers in brackets (e.g., [56 93]). Using this option
% does not modify the data.
%
% 'badChanEP' requires the Artifact_Correction_Log output from EP Toolkit
% artifact processing, and must reside in a directory labeled with the
% session name (from exper.sessions) which is in a directory called
% 'ep_art' in fullfile(dirs.dataroot,dirs.dataDir). This will only look for
% channels listed as being globally bad. Using this option does not modify
% the data.
%
% For the badChan methods, 'rmBadChan' gives the option to delete those
% channels from the data using ft_rejectvisual. Using this option will
% return the data without the bad channels.
% NB: If you edit mm_ft_artifact you can turn those channels into NaNs, or
% do that later on your own.
%
% !!!EXTREMELY IMPORTANT!!! (Disclaimer: I'm not 100% about this)
% Do not reject ICA components from data that has already had
% ICA components rejected. Also, be very wary about rejecting ICA
% components if you want to do phase analyses; I think ICA screws up phase
% information, but I need to gather more information on this. See this PDF
% for more details:
% http://www.appliedneuroscience.com/Tutorial%20on%20ICA%20Phase%20Adulteration.pdf
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB DATA
%
% If using eeglab data, no artifact detection is done and no bci file is
% expected to exist. Also, the directory structure is different and can be
% gleaned by examining the code here, but right now it is only set up to
% process Erika Nyhus's KAHN2 data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also: CREATE_FT_STRUCT, MM_FT_ARTIFACT, PROCESS_FT_DATA,
% FT_CHANNELREPAIR, FT_REJECTVISUAL
%

%% set the artifact processing parameters

if ischar(ana.artifact.type)
  ana.artifact.type = {ana.artifact.type};
end

artifactOpts = {'none','nsClassic','nsAuto','zeroVar','badChanManual','badChanEP','rmBadChan','repairBadChan','preRejManual','ftAuto','ftManual','ftICA'};

if any(~ismember(ana.artifact.type,artifactOpts))
  wrongType = ana.artifact.type(~ismember(ana.artifact.type,artifactOpts));
  error('An artifact option was not set correctly. Incorrect type(s):%s',sprintf(repmat(' %s',1,length(wrongType)),wrongType{:}));
end

% set artifact defaults
if any(ismember(ana.artifact.type,artifactOpts)) && ~ismember('none',ana.artifact.type)
  rejArt = true;
else
  rejArt = false;
end

if strcmp(ana.continuous,'yes')
  if ~isfield(ana.artifact,'continuousReject') || isempty(ana.artifact.continuousReject)
    ana.artifact.continuousReject = false;
  end
  if ~isfield(ana.artifact,'continuousRepair') || isempty(ana.artifact.continuousRepair)
    ana.artifact.continuousRepair = false;
  end
  if ~isfield(ana.artifact,'continuousICA') || isempty(ana.artifact.continuousICA)
    ana.artifact.continuousICA = false;
  end
end

% initialize
badChan = {};
badEv = [];
artfctdefEv = struct;
% maintain a list of all artifact types
artfctdefEv.types = {};
artfctdefSamp = [];

%% set up some processing parameters

% make sure eventValue is set up correctly
if ~iscell(eventValue)
  eventValue = {eventValue};
end
% if length(eventValue) > 1
%   error('Expecting only one eventValue.');
% end

if ~iscell(session)
  session = {session};
end
if length(session) > 1
  append_data = struct;
end

if strcmpi(exper.eegFileExt,'raw') || strcmpi(exper.eegFileExt,'sbin')
  ftype = 'egi_sbin';
  nsDir = 'ns_raw';
elseif strcmpi(exper.eegFileExt,'egis')
  ftype = 'egi_egis';
  nsDir = 'ns_egis';
elseif strcmpi(exper.eegFileExt,'mff')
  % see MFF info here: http://fieldtrip.fcdonders.nl/getting_started/egi
  ftype = 'egi_mff_v2';
  nsDir = 'ns_mff';
elseif strcmpi(exper.eegFileExt,'set')
  ftype = 'eeglab_set';
  nsDir = subject;
else
  error('ftype not set because extension was not properly set.');
end

% make sure the chan locs file exists
%
% if this is an EGI electrode location file included with FieldTrip, the 3
% Fid (fiduciary) points are included, meaning there are 3 non-electrodes
% in elec.label
if ~exist(elecfile,'file')
  error('Cannot find channel locs file at %s',elecfile);
else
  [cpath,cname,cext] = fileparts(elecfile);
  if strcmpi(cext,'.sfp')
    locsFormat = 'besa_sfp';
  else
    locsFormat = [];
  end
  elec = ft_read_sens(elecfile,'fileformat',locsFormat);
end
% get rid of the fiduciary channels
elec.label = ft_channelselection({'all','-Fid*'},elec.label);
nChan_elecfile = size(elec.label,1);

% find reference channel index
if isnumeric(exper.refChan)
  refChanInd = false(size(elec.label));
  refChanInd(exper.refChan) = true;
else
  refChanInd = ismember(elec.label,exper.refChan);
end
if isempty(find(refChanInd,1))
  error('Could not find reference channel.');
end

badChanAllSes = {};
badEvAllSes = [];

%% for each session, read in the EEG file

% initialize to save all trial metadata
allTrialinfo = [];

% event number column comes from the trialfun function when the trl gets
% created
trialinfo_eventNumCol = 1;

for ses = 1:length(session)
  sesName = session{ses};
  
  % set sesStr to make sure it starts with a character, not a #, etc.
  sesStr = sprintf('ses_%s',sesName);
  
  if strcmpi(exper.eegFileExt,'sbin') || strcmpi(exper.eegFileExt,'raw') || strcmpi(exper.eegFileExt,'egis') || strcmpi(exper.eegFileExt,'mff')
    % make sure the EEG file exists
    nsfile = dir(fullfile(dataroot,sesName,nsDir,[subject,'*.',exper.eegFileExt]));
    if isempty(nsfile)
      error('Cannot find %s*.%s file in %s. Make sure the subject number matches!',subject,exper.eegFileExt,fullfile(dataroot,sesName,nsDir));
    elseif length(nsfile) > 1
      error('More than one %s*.%s file found in %s',subject,exper.eegFileExt,fullfile(dataroot,sesName,nsDir));
    elseif length(nsfile) == 1
      infile_ns = fullfile(dataroot,sesName,nsDir,nsfile.name);
    end
    
  elseif strcmpi(exper.eegFileExt,'set')
    % this is really just set up to analyze Erika Nyhus's KAHN2 data
    
    %isclean = 1;
    clean_str = 'clean';
    
    %if isclean
    %  clean_str = 'clean';
    %else
    %  clean_str = '';
    %end
    
    nsfile = dir(fullfile(dataroot,nsDir,[subject,sprintf('%s%s%s.',sesName,cell2mat(eventValue),clean_str),exper.eegFileExt]));
    if isempty(nsfile)
      error('Cannot find %s file in %s',[subject,sprintf('%s%s%s.',sesName,cell2mat(eventValue),clean_str),exper.eegFileExt],fullfile(dataroot,nsDir));
    elseif length(nsfile) > 1
      error('More than one %s file found in %s',[subject,sprintf('%s%s%s.',sesName,cell2mat(eventValue),clean_str),exper.eegFileExt],fullfile(dataroot,nsDir));
    elseif length(nsfile) == 1
      infile_ns = fullfile(dataroot,nsDir,nsfile.name);
    end
  end
  
  % % debug
  % hdr = ft_read_header(infile_ns,'dataformat',ftype,'headerformat',ftype);
  % data = ft_read_data(infile_ns,'dataformat',ftype,'headerformat',ftype);
  % event = ft_read_event(infile_ns,'eventformat',ftype,'dataformat',ftype,'headerformat',ftype);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initial parameters for reading the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cfg = [];
  cfg.dataset = infile_ns;
  cfg.headerfile = infile_ns;
  if ~isempty(ftype)
    cfg.dataformat = ftype;
    cfg.headerformat = ftype;
  end
  cfg.continuous = ana.continuous;
  cfg.checksize = ana.checksize;
  
  ana.flatChan = {};
  ana.flatRef = false;
  if strcmp(ana.continuous,'yes')
    % do some initial processing of raw data
    
    % combine the initial cfg and the continuous data preprocessing
    % settings from ana.cfg_cont
    if isfield(ana,'cfg_cont')
      fprintf('%s.m: The field ana.cfg_cont exists. Using these user-defined settings for preprocessing continuous data instead of defaults.\n',mfilename);
      
      fn_cfg = fieldnames(cfg);
      if ~isempty(ana.cfg_cont)
        fn_cfg_pre = fieldnames(ana.cfg_cont);
      else
        fn_cfg_pre = [];
      end
      fn = [fn_cfg; fn_cfg_pre];
      
      if length(fn) == unique(length(fn))
        c1 = struct2cell(cfg);
        if ~isempty(ana.cfg_cont)
          c2 = struct2cell(ana.cfg_cont);
        else
          c2 = [];
        end
        c = [c1; c2];
        cfg_cont = cell2struct(c,fn,1);
      else
        error('Non-unique field names in ana.cfg_cont!\n');
        %fprintf('Non-unique field names in ana.cfg_cont!\n');
        %keyboard
      end
    else
      warning('%s: The field ana.cfg_cont does not exist. Using defaults for preprocessing continuous data!',mfilename);
      
      cfg_cont = cfg;
      % set some reasoable defaults
      cfg_cont.lpfilter = 'yes';
      cfg_cont.lpfreq = 100;
      cfg_cont.hpfilter = 'yes';
      cfg_cont.hpfreq = 0.1;
      cfg_cont.hpfilttype = 'but';
      cfg_cont.hpfiltord = 4;
      cfg_cont.bsfilter = 'yes';
      cfg_cont.bsfreq = [59 61];
      % % cfg_cont.dftfilter = 'yes';
      % % cfg_cont.dftfreq = [60 120 180];
    end
    fprintf('Preprocessing continuous data: %s...\n',infile_ns);
    data = ft_preprocessing(cfg_cont);
    fprintf('Done.\n');
    
    % find out how many channels are in the data
    nChan_data = length(data.label);
    
    %% Check on channel information
    
    % check on whether we have the reference channel
    %
    % The channel files included with FieldTrip have 3 "extra" (fiduciary)
    % channels defined, so we need to also check using an extra 3 chans
    % subtracted off
    if strcmpi(exper.eegFileExt,'mff') && (nChan_data == nChan_elecfile || (nChan_data == nChan_elecfile - 3)) && ...
        strcmp(data.label{refChanInd},'REF')
      warning('Renaming reference channel ''REF'' to ''%s'' so it matches the elecfile',exper.refChan);
      % rename the ref channel so average reference puts data in chan Cz
      data.label{refChanInd} = exper.refChan;
    end
    
    fprintf('Looking for flat channels...')
    flatChannel = false(size(data.trial{1},1),1);
    for i = 1:length(flatChannel)
      if nanvar(data.trial{1}(i,:),0,2) == 0
        flatChannel(i) = true;
      end
    end
    if any(flatChannel)
      fprintf('Done. Found %d flat channels:%s.\n.',sum(flatChannel),sprintf(repmat(' %s',1,length(data.label(flatChannel))),data.label{flatChannel}));
      flatChans = data.label(flatChannel)';
      ana.flatChan = cell(1,sum(flatChannel));
      for i = 1:length(flatChans)
        ana.flatChan{i} = sprintf('-%s',flatChans{i});
      end
    else
      fprintf('Done. Found none.\n');
    end
    if ismember(exper.refChan,flatChans)
      ana.flatRef = true;
      % MFF files has reference channel but it has zero variance
      warning('This dataset is not rereferenced (flat reference channel was included). Let''s hope you are rereferencing in FieldTrip!');
    end
    
    %% channel repair, data rejection, run ICA on continuous data
    
    if ana.artifact.continuousRepair
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Channel repair
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      rejArt_repair = [];
      while isempty(rejArt_repair) || (rejArt_repair ~= 0 && rejArt_repair ~= 1)
        rejArt_repair = input('\nDo you want to see whether there are channels to repair? (1 or 0, then press ''return''):\n\n');
      end
      if rejArt_repair
        [data,badChan] = mm_ft_artifact_repairChan(data,badChan,elecfile,'yes',[-1000 1000],30);
      end
      
      cfg_manArt = [];
      cfg_manArt.continuous = 'yes';
      cfg_manArt.blocksize = 120;
      %cfg_manArt.viewmode = 'butterfly';
      cfg_manArt.viewmode = 'vertical';
      cfg_manArt.elecfile = elecfile;
      cfg_manArt.plotlabels = 'some';
      cfg_manArt.ylim = [-100 100];
    end
      
    if ana.artifact.continuousReject
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Data rejection
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % use cursor drag and click to mark artifacts;
      % use arrows to advance to next trial;
      % use the q key to quit the data browser
      
      % manual rejection
      fprintf('Processing continuous data...\n');
      fprintf('\n\nManual artifact rejection:\n');
      fprintf('\tDrag mouse to select artifact area; click area to mark an artifact.\n');
      fprintf('\tUse arrows to move to next trial.\n');
      if strcmp(cfg_manArt.viewmode,'butterfly')
        fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
      end
      fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
      fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n');
      
      cfg_manArt = ft_databrowser(cfg_manArt,data);
      % bug when calling rejectartifact right after databrowser, pause first
      pause(1);
      % rename the visual field
      if isfield(cfg_manArt.artfctdef,'visual')
        cfg_manArt.artfctdef.visual_continuous = cfg_manArt.artfctdef.visual;
        cfg_manArt.artfctdef = rmfield(cfg_manArt.artfctdef,'visual');
        cfg_manArt.artfctdef.reject = 'partial';
        data = ft_rejectartifact(cfg_manArt, data);
      end
    end
      
    if ana.artifact.continuousICA
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % ICA on continuous data
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cfg_ica = [];
      cfg_ica.method = 'runica';
      %cfg_ica.demean = 'no';
      
      cfg_ica.channel = [{'all'}, ana.flatChan];
      
      if ~isempty(badChan)
        %       % Method 1: exclude repaired channels
        %       %
        %       % http://sccn.ucsd.edu/pipermail/eeglablist/2010/003490.html
        %
        %       fprintf('\n\tIMPORTANT! You have repaired channels:%s\n\n',sprintf(repmat(' %s',1,length(badChan_str)),badChan_str{:}));
        %       ica_chanNum = 0;
        %       fprintf('\tTherefore, you must run ICA on a subset of channels.\n');
        %     else
        %       ica_chanNum = [];
        %       fprintf('\nWe believe that you have NOT repaired any channels. Thus, you can run ICA on all channels (option ''1'').\n');
        %       fprintf('\tBut if that somehow is not the case, you must run ICA on a subset of channels (option ''0'').\n');
        
        % Method 2: run only on a given number of principle components
        %
        % http://mailman.science.ru.nl/pipermail/fieldtrip/2013-June/006656.html
        %
        % but be careful!: http://sccn.ucsd.edu/pipermail/eeglablist/2010/003339.html
        
        fprintf('Determining rank of data...');
        cfg_ica.(cfg_ica.method).pca = rank(data.trial{1});
        fprintf('Done.\n');
        
        fprintf('Running ICA after doing PCA on %d components.\n',cfg_ica.(cfg_ica.method).pca);
        fprintf('Inspired by this post: http://mailman.science.ru.nl/pipermail/fieldtrip/2013-June/006656.html\n');
        fprintf('\tBut be careful!: http://sccn.ucsd.edu/pipermail/eeglablist/2010/003339.html\n');
        fprintf('The alternative is to exclude bad channels from ICA: http://sccn.ucsd.edu/pipermail/eeglablist/2010/003490.html\n');
        fprintf('\tHowever, you are NOT doing that! You would need to uncomment some parts of %s (and comment others) to do so.\n',mfilename);
      end
      
      comp = ft_componentanalysis(cfg_ica,data);
      
      keepChoosingICAcomps = true;
      while keepChoosingICAcomps
        cfg_ica = [];
        cfg_ica.viewmode = 'component';
        cfg_ica.continuous = 'yes';
        % number of seconds to display
        cfg_ica.blocksize = 30;
        %cfg_ica.blocksize = 10;
        %cfg_ica.channels = 1:nComponents;
        cfg_ica.plotlabels = 'yes';
        cfg_ica.layout = elecfile;
        cfg_ica.elecfile = elecfile;
        cfg_ica.ylim = [-2000 2000];
        ft_databrowser(cfg_ica,comp);
        % bug when calling rejectartifact right after databrowser, pause first
        pause(1);
        
        fprintf('Processing%s...\n',sprintf(repmat(' ''%s''',1,length(eventValue_orig)),eventValue_orig{:}));
        %fprintf('\n\nViewing the first %d components.\n',nComponents);
        fprintf('ICA component browsing:\n');
        fprintf('\t1. Look for patterns that are indicative of artifacts.\n');
        fprintf('\t\tPress the ''channel >'' button to see the next set of components.\n');
        fprintf('\t\tComponents may not be numbered, so KEEP TRACK of where you are (top component has the lowest number). Write down component numbers for rejection.\n');
        fprintf('\t2. Manually close the components window when finished browsing.\n');
        
        rej_comp = [];
        while isempty(rej_comp) || (rej_comp ~= 0 && rej_comp ~= 1)
          rej_comp = input('Were there components to reject? (1 or 0, then press ''return''):\n\n');
        end
        if rej_comp
          % prompt the user for the component numbers to reject
          componentsToReject = input(sprintf('\t3. Type component numbers to reject (on a single line) and press ''return'',\n\teven if these instructions move up due to output while browsing components (e.g., ''1, 4, 11'' without quotes):\n\n'),'s');
          
          % reject the bad components
          if ~isempty(componentsToReject)
            cfg_ica = [];
            cfg_ica.component = str2double(regexp(componentsToReject,'\d*','match')');
            data_ica_rej = ft_rejectcomponent(cfg_ica, comp, data);
          end
        end
        
        cfg_browse = [];
        %cfg_browse.viewmode = 'butterfly';
        cfg_browse.viewmode = 'vertical';
        cfg_browse.continuous = 'yes';
        cfg_browse.blocksize = 120;
        cfg_browse.elecfile = elecfile;
        cfg_browse.plotlabels = 'some';
        cfg_browse.ylim = [-100 100];
        
        fprintf('\nViewing continuous data...\n');
        fprintf('\tUse arrows to move to next trial.\n');
        if strcmp(cfg_browse.viewmode,'butterfly')
          fprintf('\tUse the ''i'' key and mouse to identify channels in the data browser.\n');
        end
        fprintf('\tUse the ''q'' key to quit the data browser when finished.\n');
        fprintf('\tPress / (or any key besides q, t, i, h, c, v, or a number) to view the help screen.\n\n');
        
        ft_databrowser(cfg_browse, data_ica_rej);
        % bug when calling rejectartifact right after databrowser, pause first
        pause(1);
        
        done_with_ica = [];
        while isempty(done_with_ica) || (done_with_ica ~= 0 && done_with_ica ~= 1)
          done_with_ica = input('\nAre you happy with your post-ICA results? 1 to move on to next step, 0 to redo ICA component rejection. (1 or 0, then press ''return''):\n\n');
        end
        
        if done_with_ica
          keepChoosingICAcomps = false;
          data = data_ica_rej;
          clear data_ica_rej comp
        end
      end
    end
  end
  
  % % debug
  % data = ft_preprocessing(cfg);
  
  %% Select events
  
  % % debug
  % cfg = [];
  % cfg.dataset = infile_ns;
  % cfg.trialdef.eventtype = '?';
  % allEv = ft_definetrial(cfg);
  
  % find out which events are in infile_ns and throw an error if eventValue
  % is not one of these
  if strcmp(cfg.continuous,'no')
    cfg_noEv = [];
    cfg_noEv.dataset = infile_ns;
    cfg_noEv.trialdef.eventtype = '?';
    allEv = ft_definetrial(cfg_noEv);
    evVals = cell(size(allEv.event));
    for i = 1:length(allEv.event)
      evVals{i} = allEv.event(i).value;
    end
    evVals = unique(evVals);
    if ~ismember(eventValue_orig,evVals)
      fprintf('The available event values in %s are: %s\n',infile_ns,sprintf(repmat('''%s'' ',1,length(evVals)),evVals{:}));
      error('%s is not in the EEG file. You should redefine exper.eventValues.',sprintf(repmat('''%s'' ',1,length(eventValue_orig)),eventValue_orig{:}));
    elseif ismember(eventValue_orig,evVals)
      fprintf('You can safely ignore the warning about ''no trialfun was specified''.\n')
    end
  end
  
  % set up for defining the trials based on file type
  %if strcmp(cfg.continuous,'no')
  cfg.trialdef.eventvalue = eventValue_orig;
  %end
  
  if strcmpi(exper.eegFileExt,'sbin') || strcmpi(exper.eegFileExt,'raw') || strcmpi(exper.eegFileExt,'egis') || strcmpi(exper.eegFileExt,'mff')
    cfg.trialfun = ana.trialFxn;
    if strcmp(cfg.continuous,'no')
      cfg.trialdef.eventtype = 'trial';
      
      cfg.trialdef.prestim = abs(prepost(1)); % in seconds; must be positive
      cfg.trialdef.poststim = prepost(2); % in seconds; must be positive
      
      %       if ana.useEvents
      %         cfg.eventinfo.events = subEvents.events;
      %       end
      %       if ana.useExpParam
      %         cfg.eventinfo.expParam = expParam;
      %       end
      %       if ana.useNsEvt
      %         cfg.eventinfo.ns_evt = ns_evt;
      %       end
      %       if ana.useExpInfo
      %         cfg.eventinfo.trl_order = ana.trl_order;
      %         cfg.eventinfo.sessionNames = ana.sessionNames;
      %         cfg.eventinfo.sessionNum = ses;
      %         cfg.eventinfo.phaseNames = ana.phaseNames;
      %         cfg.eventinfo.eventValues = eventValue_orig;
      %         cfg.eventinfo.prepost = prepost;
      %       end
      %
      %       cfg.trialdef.eventtype = 'trigger';
    elseif strcmp(cfg.continuous,'yes')
      % put in extra info used for creating trialinfo into cfg.eventinfo
      % field, but be sure to remove it or the cfg will take up way too
      % much disk space (this happens automatically after running
      % ft_definetrial)
      cfg.eventinfo.useMetadata = ana.useMetadata;
      if cfg.eventinfo.useMetadata
        % need to know what kind of metadata to access
        cfg.eventinfo.metadata.types = ana.metadata.types;
        % add metadata specific to this run
        cfg.eventinfo.metadata.dataroot = dataroot;
        cfg.eventinfo.metadata.subject = subject;
        cfg.eventinfo.metadata.sesName = sesName;
        % add the basic directories struct
        cfg.eventinfo.metadata.dirs = dirs;
      end
      
      if ana.useExpInfo
        cfg.eventinfo.trl_order = ana.trl_order;
        cfg.eventinfo.sessionNames = ana.sessionNames;
        cfg.eventinfo.sessionNum = sesNum_orig;
        cfg.eventinfo.phaseNames = ana.phaseNames;
        cfg.eventinfo.eventValues = eventValue_orig;
        cfg.eventinfo.prepost = prepost;
        
        if isfield(ana,'evtToleranceMS')
          cfg.eventinfo.evtToleranceMS = ana.evtToleranceMS;
        end
        
        cfg.eventinfo.usePhotodiodeDIN = ana.usePhotodiodeDIN;
        if ana.usePhotodiodeDIN
          cfg.eventinfo.photodiodeDIN_toleranceMS = ana.photodiodeDIN_toleranceMS;
          cfg.eventinfo.photodiodeDIN_str = ana.photodiodeDIN_str;
        end
      end
      
      if strcmpi(exper.eegFileExt,'mff')
        cfg.trialdef.eventtype = 'ECI_TCPIP_55513';
      else
        cfg.trialdef.eventtype = 'trigger';
      end
    end
  elseif strcmpi(exper.eegFileExt,'set')
    cfg.trialfun = 'trialfun_general';
    cfg.trialdef.eventtype = 'trigger';
  end
  % define the trials
  try
    fprintf('Searching for %strials...\n',sprintf(repmat('''%s'' ',1,length(eventValue_orig)),eventValue_orig{:}));
    cfg = ft_definetrial(cfg);
  catch ME
    fprintf('\n');
    warning('An error occurred!!!');
    fprintf('\tFunction: %s\n',ME.stack(1).name);
    fprintf('\tLine number: %d\n',ME.stack(1).line);
    fprintf('\tReason: %s\n',ME.message);
    
    % if there were zero trials for this event type
    if strfind(ME.message,'no trials were defined')
      fprintf('No %s events found!\n',sprintf(repmat('''%s'' ',1,length(eventValue_orig)),eventValue_orig{:}));
    end
    warning('Returning an empty dataset for %s. This will save an error file when running the ft_*analysis function.',...
      sprintf(repmat('''%s'' ',1,length(eventValue_orig)),eventValue_orig{:}));
    
    % something is wrong in your trialfun, figure it out; or dbcont
    %keyboard
    %error('something is wrong in your trialfun, figure it out');
    
    % set an empty cell and return to the calling function
    ft_raw.data.trial = {};
    return
  end
  
  % remove these fields or the cfg will take up way too much disk space
  if isfield(cfg,'eventinfo')
    cfg = rmfield(cfg,'eventinfo');
  end
  if isfield(cfg.callinfo.usercfg,'eventinfo')
    cfg.callinfo.usercfg = rmfield(cfg.callinfo.usercfg,'eventinfo');
  end
  
  %% Get the data and process it if necessary
  
  % collect the event type numbers
  allTrialinfo = cat(1,allTrialinfo,cfg.trl(:,4:end));
  
  % get the actual data
  if strcmp(cfg.continuous,'no')
    data_seg = ft_preprocessing(cfg);
  elseif strcmp(cfg.continuous,'yes')
    
    if exist('cfg_manArt','var') && isfield(cfg_manArt.artfctdef,'visual_continuous') && ...
        isfield(cfg_manArt.artfctdef.visual_continuous,'artifact') && ~isempty(cfg_manArt.artfctdef.visual_continuous.artifact)
      % initialize to store whether there was an artifact for each trial
      if isempty(badEv)
        combineArtLists = false;
        %badEv = [(1:size(data.sampleinfo,1))', zeros(size(data.sampleinfo,1), 1)];
        badEv = false(size(cfg.trl,1), 1);
      else
        combineArtLists = true;
      end
      
      % find out what kind of artifacts we're dealing with
      fn = fieldnames(cfg_manArt.artfctdef);
      theseArt = {};
      for i = 1:length(fn)
        if isstruct(cfg_manArt.artfctdef.(fn{i})) && isfield(cfg_manArt.artfctdef.(fn{i}),'artifact')
          theseArt = cat(2,theseArt,fn{i});
          if ~ismember(fn{i},artfctdefEv.types)
            artfctdefEv.types = cat(2,artfctdefEv.types,fn{i});
          else
            warning('''%s'' is already in the artifact types, if you continue you will overwrite the previous artifact information for this type.',fn{i});
            keyboard
          end
          artfctdefEv.(fn{i}) = false(size(badEv));
        end
      end
      
      % store artifacts for the only events we are checking
      foundArtEv = false(size(cfg.trl,1),length(theseArt));
      
      % find out which samples were marked as artifacts
      if ~isempty(theseArt)
        artSamp = false(cfg.trl(end,2),length(theseArt));
        for i = 1:length(theseArt)
          for j = 1:size(cfg_manArt.artfctdef.(theseArt{i}).artifact,1)
            % mark that it was a particular type of artifact
            artSamp(cfg_manArt.artfctdef.(theseArt{i}).artifact(j,1):cfg_manArt.artfctdef.(theseArt{i}).artifact(j,2),ismember(theseArt,theseArt{i})) = true;
          end
        end
        % save a list of trials with artifact status
        for k = 1:size(cfg.trl,1)
          if any(any(artSamp(cfg.trl(k,1):cfg.trl(k,2),:),1),2)
            foundArtEv(k,any(artSamp(cfg.trl(k,1):cfg.trl(k,2),:),1)) = true;
          end
        end
        clear artSamp
      end
      
      if combineArtLists
        % put the new artifacts into the old list
        rCount = 0;
        for i = 1:size(badEv,1)
          %if badEv(i,2) == 0
          if badEv(i,1) == 0
            rCount = rCount + 1;
            if any(foundArtEv(rCount,:))
              %badEv(i,2) = 1;
              badEv(i,1) = 1;
              setTheseArt = theseArt(foundArtEv(rCount,:));
              for a = 1:length(setTheseArt)
                artfctdefEv.(setTheseArt{a})(i) = true;
              end
            end
          end
        end
        if ~isempty(theseArt)
          if isempty(artfctdefSamp)
            artfctdefSamp = cfg_manArt.artfctdef;
          else
            for i = 1:length(theseArt)
              if isfield(artfctdefSamp,theseArt{i})
                artfctdefSamp.(theseArt{i}).artifact = cat(1,artfctdefSamp.(theseArt{i}).artifact,cfg_manArt.artfctdef.(theseArt{i}).artifact);
              else
                artfctdefSamp.(theseArt{i}).artifact = cfg_manArt.artfctdef.(theseArt{i}).artifact;
              end
            end
          end
        end
      else
        badEv = logical(sum(foundArtEv,2));
        artfctdefSamp = cfg_manArt.artfctdef;
        for a = 1:length(theseArt)
          artfctdefEv.(theseArt{a}) = foundArtEv(:,a);
        end
      end
      
      if any(badEv)
        % remove trials from rejected regions
        fprintf('Removing %d trials from rejected data regions.\n',sum(badEv));
        cfg.trl = cfg.trl(~badEv,:);
      end
    end
    
    fprintf('Segmenting continuous data...');
    data_seg = ft_redefinetrial(cfg, data);
    fprintf('Done.\n');
    
    % save memory
    clear data
  end
  
  % collect the event type numbers - do this with cfg  up above instead
  %allTrialinfo = cat(1,allTrialinfo,data_seg.trialinfo);
  
  % hack: renumber samples so they don't overlap, or throw an error
  overlap = false;
  for i = 2:size(data_seg.sampleinfo,1)
    if data_seg.sampleinfo(i,1) < data_seg.sampleinfo(i-1,2)
      overlap = true;
      fprintf('Found trials with overlapping sample indices! This could be problematic during artifact rejection.\n');
      break
    end
  end
  if overlap && ana.allowTrialOverlap
    if ana.renumberSamplesContiguous
      % overwrite the sample numbers so trials do not overlap
      fprintf('Renumbering sample indices so trials are contiguous....');
      
      nSamp = length(data_seg.sampleinfo(1,1):data_seg.sampleinfo(1,2));
      % set the beginning sample
      data_seg.sampleinfo(1,1) = 1;
      % set the ending sample
      data_seg.sampleinfo(1,2) = nSamp;
      for i = 2:size(data_seg.sampleinfo,1)
        nSamp = length(data_seg.sampleinfo(i,1):data_seg.sampleinfo(i,2));
        % set the beginning sample
        data_seg.sampleinfo(i,1) = data_seg.sampleinfo(i-1,2) + 1;
        % set the ending sample
        data_seg.sampleinfo(i,2) = data_seg.sampleinfo(i,1) + nSamp - 1;
      end
      % and overwrite the trl stored in cfg
      data_seg.cfg.trl(:,1:2) = data_seg.sampleinfo(:,1:2);
      fprintf('Done.\n');
    end
  elseif overlap && ~ana.allowTrialOverlap
    error('Found overlapping trials (e.g., trials %d and %d, but there are probably more). You have set ana.allowTrialOverlap=false, so this is not allowed!',i-1,i);
  elseif ~overlap && ana.renumberSamplesContiguous
    fprintf('No trials with overlapping samples found. Renumbering sample indices anyway so trials are contiguous....');
    nSamp = length(data_seg.sampleinfo(1,1):data_seg.sampleinfo(1,2));
    % set the beginning sample
    data_seg.sampleinfo(1,1) = 1;
    % set the ending sample
    data_seg.sampleinfo(1,2) = nSamp;
    for i = 2:size(data_seg.sampleinfo,1)
      nSamp = length(data_seg.sampleinfo(i,1):data_seg.sampleinfo(i,2));
      % set the beginning sample
      data_seg.sampleinfo(i,1) = data_seg.sampleinfo(i-1,2) + 1;
      % set the ending sample
      data_seg.sampleinfo(i,2) = data_seg.sampleinfo(i,1) + nSamp - 1;
    end
    % and overwrite the trl stored in cfg
    data_seg.cfg.trl(:,1:2) = data_seg.sampleinfo(:,1:2);
    fprintf('Done.\n');
  end
  
  % find out how many channels are in the segmented data
  nChan_data_seg = length(data_seg.label);
  
  %% Check on channel information
  
  % check on whether we have the reference channel
  %
  % The channel files included with FieldTrip have 3 "extra" (fiduciary)
  % channels defined, so we need to also check using an extra 3 chans
  % subtracted off
  if (nChan_data_seg == nChan_elecfile - 1) || (nChan_data_seg == nChan_elecfile - 4)
    % one less channel because we're checking to see if the reference
    % channel is missing
    
    %error('This dataset is either not rereferenced or the reference channel was not exported. Go back and rereference or export the reference channel in Net Station before running this script!');
    warning('This dataset is either not rereferenced or the reference channel was not exported. Let''s hope you are rereferencing in FieldTrip!');
  elseif strcmp(cfg.continuous,'no') && ~strcmpi(exper.eegFileExt,'mff') && (nChan_data_seg == nChan_elecfile || nChan_data_seg == nChan_elecfile - 3)
    % grab data from all of the trials
    trialData = cat(3,data_seg.trial{:});
    
    % check the variance across time for the reference channel
    if sum(nanvar(trialData(refChanInd,:,:),0,2) ~= 0) == 0
      % if none of trials have a non-zero variance reference channel, then
      % it has not been rereferenced. Some trials may have zero variance
      % because of how bad trial rejection works in EP Toolkit (it zeros
      % out all channels for bad trials).
      %
      % var=0 means that the final (reference) electrode is flat and this
      % data set has not been (average) rereferenced
      error('This dataset is not rereferenced. Go back and rereference in Net Station before running this script!');
    else
      
      % has full number of channels and is already rereferenced (ref channel
      % is not flat); check multiple trials
      fprintf('Channels are already (average) rereferenced, as they should be.\n');
      
      % depending on whether the channel string was capitalized or lowercase
      % in the electrode template, make the data elec label match. This is
      % actually important for how FieldTrip deals with electrode numbers.
      %
      % TODO: We now always want to use capital letters, so this should
      % probably be changed.
      if strcmp(elec.label{ceil(nChan_data_seg/2)}(1),'E')
        isCapital = 1;
      elseif strcmp(elec.label{ceil(nChan_data_seg/2)}(1),'e')
        isCapital = 0;
      else
        warning([mfilename,':electrodeCapitalization'],'There is no ''E'' or ''e'' at the start of the electrode number! Going with uppercase.')
        isCapital = 1;
      end
      
      if isCapital
        % capitalize the E for each electrode, or add it in if it's not there
        for c = 1:nChan_data_seg
          if strcmp(data_seg.label{c}(1),'e')
            data_seg.label{c} = upper(data_seg.label{c});
          elseif ~strcmp(data_seg.label{c}(1),'e') && ~strcmp(data_seg.label{c}(1),'E')
            data_seg.label{c} = ['E' data_seg.label{c}];
          end
        end
      elseif ~isCapital
        % make sure the e for each electrode is lowercase, or add it in if
        % it's not there
        for c = 1:nChan_data_seg
          if strcmp(data_seg.label{c}(1),'E')
            data_seg.label{c} = lower(data_seg.label{c});
          elseif ~strcmp(data_seg.label{c}(1),'e') && ~strcmp(data_seg.label{c}(1),'E')
            data_seg.label{c} = ['e' data_seg.label{c}];
          end
        end
      end
      
      % set the last channel name to 'Cz' if that's what was set in
      % elec.label (e.g., instead of 'E129')
      if strcmp(elec.label{end},'Cz')
        if isCapital
          lastChanStr = sprintf('E%d',nChan_data_seg);
        elseif ~isCapital
          lastChanStr = sprintf('e%d',nChan_data_seg);
        end
        %lastChanStr = 'Cz';
        chanindx = find(strcmpi(data_seg.label,lastChanStr));
        if ~isempty(chanindx)
          % set the label for the reference channel
          %data_seg.label{chanindx} = elec.label{chanindx};
          data_seg.label{chanindx} = elec.label{end};
        end
      end
    end
  elseif ~strcmpi(exper.eegFileExt,'mff')
    error('Not sure what to do about rereferencing!');
  end
  
  %% artifact rejection
  
  if ~rejArt
    fprintf('Not performing any artifact rejection.\n');
  else
%     % save samples of visual_continuous artifacts (might have been
%     % renumbered)
%     if exist('cfg_manArt','var') && isfield(cfg_manArt.artfctdef,'visual_continuous') && ...
%         isfield(cfg_manArt.artfctdef.visual_continuous,'artifact') && ~isempty(cfg_manArt.artfctdef.visual_continuous.artifact)
%       cfg_manArt.artfctdef.visual_continuous.artifact = data_seg.sampleinfo(badEv,:);
%     end
    
    % detect/reject other artifacts
    [data_seg,badChan,badEv,artfctdefEv] = mm_ft_artifact(dataroot,subject,sesName,eventValue_orig,ana,exper,elecfile,data_seg,dirs,badChan,badEv,artfctdefEv,artfctdefSamp);
    badChanAllSes = unique(cat(1,badChanAllSes,badChan));
    % Concatenate sessions together if they're getting combined (appended).
    % Otherwise cat() won't make any difference.
    badEvAllSes = cat(1,badEvAllSes,badEv);
    
%     foundArt = false(size(badEv));
%     for i = 1:length(artfctdefEv.types)
%       if ~strcmp(artfctdefEv.types{i},'visual_continuous')
%         foundArt = logical(foundArt + artfctdefEv.(artfctdefEv.types{i}));
%       end
%     end
    
%     % make sure we get rid of those visual_continuous artifacts
%     if exist('cfg_manArt','var') && isfield(cfg_manArt.artfctdef,'visual_continuous') && ...
%         isfield(cfg_manArt.artfctdef.visual_continuous,'artifact') && ~isempty(cfg_manArt.artfctdef.visual_continuous.artifact)
%         
%       cfg_manArt.artfctdef.reject = 'complete';
%       data_seg = ft_rejectartifact(cfg_manArt, data_seg);
%     end
    
%     if ~exist('artfctdefSampAllSes','var')
%       artfctdefSampAllSes = artfctdefSamp;
%     else
%       fn = fieldnames(artfctdefSamp);
%       for f = 1:length(fn)
%         if isfield(artfctdefSampAllSes,fn{f})
%           artfctdefSampAllSes.(fn{f}) = cat(1,artfctdefSampAllSes.(fn{f}),artfctdefSamp.(fn{f}));
%         else
%           artfctdefSampAllSes.(fn{f}) = artfctdefSamp.(fn{f});
%         end
%       end
%     end
    if ~exist('artfctdefEvAllSes','var')
      artfctdefEvAllSes = artfctdefEv;
    else
      fn = fieldnames(artfctdefEv);
      for f = 1:length(fn)
        if strcmp(fn{f},'types')
          if isfield(artfctdefEvAllSes,(fn{f}))
            artfctdefEvAllSes.(fn{f}) = cat(2,artfctdefEvAllSes.(fn{f}),artfctdefEv.(fn{f}));
          else
            artfctdefEvAllSes.(fn{f}) = artfctdefEv.(fn{f});
          end
        else
          if isfield(artfctdefEvAllSes,fn{f})
            artfctdefEvAllSes.(fn{f}) = cat(1,artfctdefEvAllSes.(fn{f}),artfctdefEv.(fn{f}));
          else
            artfctdefEvAllSes.(fn{f}) = artfctdefEv.(fn{f});
          end
        end
      end
    end
  end
  
  %% if we're combining multiple sessions, add the data to the append struct
  if length(session) > 1
    append_data.(sesStr) = data_seg;
  end
end % ses

%% Append sessions, if necessary
  
% run ft_appenddata if we're combining multiple sessions
if length(session) > 1
  sesStr = sprintf('ses_%s',session{1});
  append_str = sprintf('append_data.%s',sesStr);
  
  for ses = 2:length(session)
    sesStr = sprintf('ses_%s',session{ses});
    append_str = cat(2,append_str,sprintf(',append_data.%s',sesStr));
  end
  
  data_seg = eval(sprintf('ft_appenddata([],%s);',append_str));
end

%% Separate the event values

% initialize the struct to return
ft_raw = struct;

badEvEvVals = struct;
artifacts = [];
if rejArt
  artTypes = fieldnames(artfctdefEvAllSes);
  artTypes = artTypes(~ismember(artTypes,'types'));
end
trialinfo_allEv = [];

for evVal = 1:length(eventValue)
  badEvEvVals.(eventValue{evVal}) = [];
  if ana.useExpInfo
    trialinfo_allEv.(eventValue{evVal}) = [];
  end
  if rejArt
    for at = 1:length(artTypes)
      artifacts.(eventValue{evVal}).(artTypes{at}) = [];
    end
  end
end

% if length(eventValue) > 1
for evVal = 1:length(eventValue)
  
  cfg_split = [];
  % select the correct trials for this event value
  cfg_split.trials = data_seg.trialinfo(:,trialinfo_eventNumCol) == evVal;
  
  if sum(cfg_split.trials) > 0
    fprintf('Selecting %d trials for %s...\n',sum(cfg_split.trials),eventValue{evVal});
    % get the data for only this event value
    ft_raw.(eventValue{evVal}) = ft_redefinetrial(cfg_split,data_seg);
    
    if rejArt
      %badEvEvVals.(eventValue{evVal}) = badEvAllSes(:,trialinfo_eventNumCol) == evVal;
      %badEvEvVals.(eventValue{evVal}) = badEvAllSes(cfg_split.trials);
      badEvEvVals.(eventValue{evVal}) = badEvAllSes(allTrialinfo(:,trialinfo_eventNumCol) == evVal);
      
      for at = 1:length(artTypes)
        %artifacts.(eventValue{evVal}).(artTypes{at}) = artfctdefEvAllSes.(artTypes{at})(:,trialinfo_eventNumCol) == evVal;
        %artifacts.(eventValue{evVal}).(artTypes{at}) = artfctdefEvAllSes.(artTypes{at})(cfg_split.trials);
        artifacts.(eventValue{evVal}).(artTypes{at}) = artfctdefEvAllSes.(artTypes{at})(allTrialinfo(:,trialinfo_eventNumCol) == evVal);
      end
    end
    
    if ana.useExpInfo
      % remove the buffer trialinfo -1s; those were set because the cfg.trl
      % matrix needed to hold all eventValues
      ft_raw.(eventValue{evVal}).trialinfo = ft_raw.(eventValue{evVal}).trialinfo(:,1:length(ana.trl_order.(eventValue{evVal})));
      trialinfo_allEv.(eventValue{evVal}) = allTrialinfo(allTrialinfo(:,trialinfo_eventNumCol) == evVal,1:length(ana.trl_order.(eventValue{evVal})));
    end
    
    fprintf('Done.\n');
  else
    warning('No trials found for %s!',eventValue{evVal});
    ft_raw.(eventValue{evVal}).trial = {};
    
    % something is wrong, figure it out; or dbcont
    %keyboard
  end
end
% elseif length(eventValue) == 1
%   ft_raw.(eventValue{1}) = data_seg;
%   % remove the buffer trialinfo -1s; those were set because the cfg.trl
%   % matrix needed to hold all eventValues
%   ft_raw.(eventValue{1}).trialinfo = ft_raw.(eventValue{1}).trialinfo(:,1:length(ana.trl_order.(eventValue{1})));
% end

end
