function cosi_prepData_summary(rejectArt,saveFiles,averageSes,dataroot)
% function cosi_prepData_summary(rejectArt,saveFiles,averageSes,dataroot)
%
% Save column-formatted data as a csv file
%
% Looks for each subject's events.mat and saves summary file(s) in dataroot
%

expName = 'COSI';

subsets = {'_color','_side'};

if nargin < 4
  serverDir = fullfile('/Volumes/curranlab/Data',expName,'eeg/behavioral');
  serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data',expName,'eeg/behavioral');
  if exist(serverDir,'dir')
    dataroot = serverDir;
  elseif exist(serverLocalDir,'dir')
    dataroot = serverLocalDir;
  else
    uroot = getenv('HOME');
    dataroot = fullfile(uroot,'data',expName,'eeg/behavioral');
  end
  if nargin < 3
    averageSes = 0;
    if nargin < 2
      saveFiles = 1;
      if nargin < 1
        rejectArt = 0;
      end
    end
  end
end

subjects = {
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
  'COSI020';
  'COSI019';
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
% original COSI019 was replaced because the first didn't finish

% sessions = {'session_0'};
sessions = {'session_0','session_1'};

%% Set up the headers
if rejectArt == 0
  tableHeader = {'sub','ses',... % subject info
    'R targ','R lure','S targ','S lure',... % raw numbers of targs, lures
    'Item d''','Source d''','Item RespBias (c)','Source RespBias (c)',... % accuracy
    'RH','RM','RCR','RFA','SH','SM','SCR','SFA',... % raw numbers for accuracy
    'RHR','RMR','RCRR','RFAR','SHR','SMR','SCRR','SFAR',... % accuracy rates
    'RH-SCor','RH-SInc',... % raw numbers for recognition hits
    'RH-SCorR','RH-SIncR',... % accuracy rates for recognition hits
    'RH-RS','RH-RO','RH-Fam',... % raw numbers for recognition hit response collapsed across source accuracy
    'RH-SCor RS','RH-SCor RO','RH-SCor Fam','RH-SInc RS','RH-SInc RO','RH-SInc Fam','RM Sure','RM Maybe','RCR Sure','RCR Maybe','RFA RS','RFA RO','RFA Fam',... % raw numbers by confidence
    'RH-SCor RS R','RH-SCor RO R','RH-SCor Fam R','RH-SInc RS R','RH-SInc RO R','RH-SInc Fam R','RM Sure R','RM Maybe R','RCR Sure R','RCR Maybe R','RFA RS R','RFA RO R','RFA Fam R',... % rates by confidence
    'IRK Fam','IRK Corr Fam','IRK Inc Fam',... % Independent Remember--Know Familiarity rates
    'RS-SCor WIR','RS-SInc WIR','RO-SCor WIR','RO-SInc WIR','F-SCor WIR','F-SInc WIR','RM Sure WIR','RM Maybe WIR','RCR Sure WIR','RCR Maybe WIR','RFA RS WIR','RFA RO WIR','RFA Fam WIR',... % within-response rates
    'RH rt','RM rt','RCR rt','RFA rt',... % rec reaction times
    'SH rt','SM rt','SCR rt','SFA rt',... % src reaction times
    'RH-SCor rt','RH-SInc rt',... % recognition hits reaction times
    'RH-SCor RS rt','RH-SCor RO rt','RH-SCor Fam rt','RH-SInc RS rt','RH-SInc RO rt','RH-SInc Fam rt','RM Sure rt','RM Maybe rt','RCR Sure rt','RCR Maybe rt','RFA RS rt','RFA RO rt','RFA Fam rt',... % rt by confidence
    'Rec Discrimination (Pr)','Rec Response Bias (Br)',... % more accuracy
    };
elseif rejectArt == 1
  tableHeader = {'sub','ses',... % subject info
    'R targ','R lure','S targ','S lure',... % raw numbers of targs, lures
    'Source d''','Source RespBias (c)',... % accuracy
    'RH','RCR','SH','SM','SCR','SFA',... % raw numbers for accuracy
    'SHR','SMR','SCRR','SFAR',... % accuracy rates
    'RH-SCor','RH-SInc',... % raw numbers for recognition hits
    'RH-SCorR','RH-SIncR',... % accuracy rates for recognition hits
    'RH-RS','RH-RO','RH-Fam',... % raw numbers for recognition hit response collapsed across source accuracy
    'RH-SCor RS','RH-SCor RO','RH-SCor Fam','RH-SInc RS','RH-SInc RO','RH-SInc Fam','RCR Sure','RCR Maybe',... % raw numbers by confidence
    'RH-SCor RS R','RH-SCor RO R','RH-SCor Fam R','RH-SInc RS R','RH-SInc RO R','RH-SInc Fam R','RCR Sure R','RCR Maybe R',... % rates by confidence
    'IRK Fam','IRK Corr Fam','IRK Inc Fam',... % Independent Remember--Know Familiarity rates
    'RS-SCor WIR','RS-SInc WIR','RO-SCor WIR','RO-SInc WIR','F-SCor WIR','F-SInc WIR','RCR Sure WIR','RCR Maybe WIR',... % within-response rates
    'RH rt','RCR rt',... % rec reaction times
    'SH rt','SM rt','SCR rt','SFA rt',... % src reaction times
    'RH-SCor rt','RH-SInc rt',... % recognition hits reaction times
    'RH-SCor RS rt','RH-SCor RO rt','RH-SCor Fam rt','RH-SInc RS rt','RH-SInc RO rt','RH-SInc Fam rt','RCR Sure rt','RCR Maybe rt',... % rt by confidence
    };
end

if rejectArt == 1
  artSegStr = sprintf('Segments marked bad for eye artifacts');
  tableHeader{end+1} = artSegStr;
end

%matlabpool open 3

%% for each subset
for s = 1:length(subsets)
  if saveFiles == 1
    % include artifacts in event count
    eventsTable = fullfile(dataroot,[expName,'_summary',subsets{s}]);
    if averageSes == 1
      eventsTable = sprintf('%s_avgSes',eventsTable);
    end
    if rejectArt == 1
      eventsTable = sprintf('%s_rejArt',eventsTable);
    end
    eventsTable = sprintf('%s.csv',eventsTable);
    
    if ~exist(eventsTable,'file')
      outfile = fopen(eventsTable,'wt');
      if isempty(subsets{s})
        fprintf(outfile,'%s\n','All together');
      else
        fprintf(outfile,'%s\n',strrep(subsets{s},'_',''));
      end
      fprintf(outfile,'%s',tableHeader{1});
      for headCol = 2:length(tableHeader)
        fprintf(outfile,',%s',tableHeader{headCol});
      end
      fprintf(outfile,'\n');
    else
      fprintf('%s already exists, not overwriting! Moving on to next file.\n',eventsTable);
      continue
      %saveFiles = 0;
    end
  else
    fprintf('saveFiles = 0, not writing out any files\n')
  end
  
  %% initialize
  numEv = struct;
  rates = struct;
  rates_wir = struct;
  rt = struct;
  acc = struct;
  
  % number of events
  numEv.rec_targEv = nan(length(subjects),length(sessions));
  numEv.rec_lureEv = nan(length(subjects),length(sessions));
  numEv.src_targEv = nan(length(subjects),length(sessions));
  numEv.src_lureEv = nan(length(subjects),length(sessions));
  % raw numbers for accuracy
  numEv.rec_h = nan(length(subjects),length(sessions));
  numEv.rec_cr = nan(length(subjects),length(sessions));
  if rejectArt == 0
    numEv.rec_m = nan(length(subjects),length(sessions));
    numEv.rec_fa = nan(length(subjects),length(sessions));
  end
  numEv.src_h = nan(length(subjects),length(sessions));
  numEv.src_cr = nan(length(subjects),length(sessions));
  numEv.src_m = nan(length(subjects),length(sessions));
  numEv.src_fa = nan(length(subjects),length(sessions));
  % raw numbers for recognition hits
  numEv.rec_h_srcCor = nan(length(subjects),length(sessions));
  numEv.rec_h_srcInc = nan(length(subjects),length(sessions));
  % raw numbers by confidence
  numEv.rec_h_srcCor_rs = nan(length(subjects),length(sessions));
  numEv.rec_h_srcCor_ro = nan(length(subjects),length(sessions));
  numEv.rec_h_srcCor_k = nan(length(subjects),length(sessions));
  numEv.rec_h_srcInc_rs = nan(length(subjects),length(sessions));
  numEv.rec_h_srcInc_ro = nan(length(subjects),length(sessions));
  numEv.rec_h_srcInc_k = nan(length(subjects),length(sessions));
  numEv.rec_m_sure = nan(length(subjects),length(sessions));
  numEv.rec_m_maybe = nan(length(subjects),length(sessions));
  numEv.rec_cr_sure = nan(length(subjects),length(sessions));
  numEv.rec_cr_maybe = nan(length(subjects),length(sessions));
  numEv.rec_fa_rs = nan(length(subjects),length(sessions));
  numEv.rec_fa_ro = nan(length(subjects),length(sessions));
  numEv.rec_fa_k = nan(length(subjects),length(sessions));
  
  % rates
  rates.rec_h = nan(length(subjects),length(sessions));
  rates.rec_cr = nan(length(subjects),length(sessions));
  if rejectArt == 0
    rates.rec_m = nan(length(subjects),length(sessions));
    rates.rec_fa = nan(length(subjects),length(sessions));
  end
  rates.src_h = nan(length(subjects),length(sessions));
  rates.src_m = nan(length(subjects),length(sessions));
  rates.src_cr = nan(length(subjects),length(sessions));
  rates.src_fa = nan(length(subjects),length(sessions));
  rates.rec_h_srcCor = nan(length(subjects),length(sessions));
  rates.rec_h_srcInc = nan(length(subjects),length(sessions));
  % rates by confidence
  rates.rec_h_srcCor_rs = nan(length(subjects),length(sessions));
  rates.rec_h_srcCor_ro = nan(length(subjects),length(sessions));
  rates.rec_h_srcCor_k = nan(length(subjects),length(sessions));
  rates.rec_h_srcInc_rs = nan(length(subjects),length(sessions));
  rates.rec_h_srcInc_ro = nan(length(subjects),length(sessions));
  rates.rec_h_srcInc_k = nan(length(subjects),length(sessions));
  rates.rec_h_rs = nan(length(subjects),length(sessions));
  rates.rec_h_ro = nan(length(subjects),length(sessions));
  rates.rec_h_k = nan(length(subjects),length(sessions));
  rates.rec_cr_sure = nan(length(subjects),length(sessions));
  rates.rec_cr_maybe = nan(length(subjects),length(sessions));
  if rejectArt == 0
    rates.rec_m_sure = nan(length(subjects),length(sessions));
    rates.rec_m_maybe = nan(length(subjects),length(sessions));
    rates.rec_fa_rs = nan(length(subjects),length(sessions));
    rates.rec_fa_ro = nan(length(subjects),length(sessions));
    rates.rec_fa_k = nan(length(subjects),length(sessions));
  end
  % Independent Remember--Know Familiarity rates
  rates.irk_Fam = nan(length(subjects),length(sessions));
  %I think this is wrong
  rates.irk_correct_Fam = nan(length(subjects),length(sessions));
  rates.irk_incorrect_Fam = nan(length(subjects),length(sessions));
  
  % accuracy
  if rejectArt == 0
    acc.item_dp = nan(length(subjects),length(sessions));
  end
  acc.source_dp = nan(length(subjects),length(sessions));
  % response bias (criterion) from Macmillan & Creelman p. 29
  if rejectArt == 0
    acc.item_c = nan(length(subjects),length(sessions));
  end
  acc.source_c = nan(length(subjects),length(sessions));
  if rejectArt == 0
    % From Mecklinger et al. (2007) and Corwin (1994)
    %
    % discrimination index
    acc.Pr = nan(length(subjects),length(sessions));
    % response bias index
    acc.Br = nan(length(subjects),length(sessions));
  end
  
  % within-response rates
  rates_wir.rec_h_srcCor_rs = nan(length(subjects),length(sessions));
  rates_wir.rec_h_srcCor_ro = nan(length(subjects),length(sessions));
  rates_wir.rec_h_srcCor_k = nan(length(subjects),length(sessions));
  rates_wir.rec_h_srcInc_rs = nan(length(subjects),length(sessions));
  rates_wir.rec_h_srcInc_ro = nan(length(subjects),length(sessions));
  rates_wir.rec_h_srcInc_k = nan(length(subjects),length(sessions));
  rates_wir.rec_cr_sure = nan(length(subjects),length(sessions));
  rates_wir.rec_cr_maybe = nan(length(subjects),length(sessions));
  if rejectArt == 0
    rates_wir.rec_m_sure = nan(length(subjects),length(sessions));
    rates_wir.rec_m_maybe = nan(length(subjects),length(sessions));
    rates_wir.rec_fa_rs = nan(length(subjects),length(sessions));
    rates_wir.rec_fa_ro = nan(length(subjects),length(sessions));
    rates_wir.rec_fa_k = nan(length(subjects),length(sessions));
  end
  
  % reaction times
  rt.rec_h = nan(length(subjects),length(sessions));
  rt.rec_cr = nan(length(subjects),length(sessions));
  if rejectArt == 0
    rt.rec_m = nan(length(subjects),length(sessions));
    rt.rec_fa = nan(length(subjects),length(sessions));
  end
  rt.src_h = nan(length(subjects),length(sessions));
  rt.src_m = nan(length(subjects),length(sessions));
  rt.src_cr = nan(length(subjects),length(sessions));
  rt.src_fa = nan(length(subjects),length(sessions));
  rt.rec_h_srcCor = nan(length(subjects),length(sessions));
  rt.rec_h_srcInc = nan(length(subjects),length(sessions));
  rt.rec_h_srcCor_rs = nan(length(subjects),length(sessions));
  rt.rec_h_srcCor_ro = nan(length(subjects),length(sessions));
  rt.rec_h_srcCor_k = nan(length(subjects),length(sessions));
  rt.rec_h_srcInc_rs = nan(length(subjects),length(sessions));
  rt.rec_h_srcInc_ro = nan(length(subjects),length(sessions));
  rt.rec_h_srcInc_k = nan(length(subjects),length(sessions));
  rt.rec_cr_sure = nan(length(subjects),length(sessions));
  rt.rec_cr_maybe = nan(length(subjects),length(sessions));
  if rejectArt == 0
    rt.rec_m_sure = nan(length(subjects),length(sessions));
    rt.rec_m_maybe = nan(length(subjects),length(sessions));
    rt.rec_fa_rs = nan(length(subjects),length(sessions));
    rt.rec_fa_ro = nan(length(subjects),length(sessions));
    rt.rec_fa_k = nan(length(subjects),length(sessions));
  end
  
  %% for each subject
  for sub = 1:length(subjects)
    fprintf('Getting data for %s...',subjects{sub});
    
    %% for each session
    for ses = 1:length(sessions)
      fprintf('%s...\n',sessions{ses});
      % set the subject events directory
      eventsDir_sub = fullfile(dataroot,subjects{sub},sessions{ses},'events');
      if exist(fullfile(eventsDir_sub,'events.mat'),'file')
        events = loadEvents(fullfile(eventsDir_sub,'events.mat'));
      else
        fprintf('events.mat for %s %s does not exist. Moving on.\n',subjects{sub},sessions{ses});
        continue
      end
      
      % include events with defined study color name and that was tested
      events = filterStruct(events,'src_correct ~= -1');
      
      if rejectArt == 1
        % remove artifacts
        subSesEv_all = events;
        events = filterStruct(subSesEv_all,'nsArt == 0');
      end
      
      % get only the events for this session
      %subSesEv = filterStruct(events,'ismember(subject,varargin{1}) & ismember(session,varargin{2})',subs_str{sub},sess(ses));
      
      % get only the events for this subset
      if strcmp(subsets{s},'_color')
        % get only the size study events
        fprintf('Color\n');
        subSesEv = filterStruct(events,'ismember(cond,varargin{1})',{'color'});
      elseif strcmp(subsets{s},'_side')
        % get only the life study events
        fprintf('Side\n');
        subSesEv = filterStruct(events,'ismember(cond,varargin{1})',{'side'});
      else
        fprintf('All\n');
        subSesEv = events;
      end
      
      %% RAW NUMBERS
      
      % When analyzing counts/rates when rejecting artifacts:
      %
      % Don't currently have all event information if using rejectArt; only have
      % info for recognition hits and correct rejections. Source hits, CR,
      % missed, and FAs are ok because they come from recognition hits. This is
      % because we're not segmenting misses and false alarms in NS.
      
      %         % get all the study events (excluding buffers)
      %         studyEv = filterStruct(subSesEv,'ismember(type,varargin{1})',{'STUDY_TARGET'});
      %
      %         study_rh = filterStruct(studyEv,'rec_correct == 1');
      %         study_rm = filterStruct(studyEv,'rec_correct == 0');
      %
      %         study_rh_srcCor = filterStruct(study_rh,'src_correct == 1');
      %         study_rh_srcInc = filterStruct(study_rh,'src_correct == 0');
      
      % get all the targets (old study items presented at test)
      rec_targEv = filterStruct(subSesEv,'rec_isTarg == 1 & ismember(type,varargin{1})',{'TEST_TARGET'});
      % get all the lures (new test items)
      rec_lureEv = filterStruct(subSesEv,'rec_isTarg == 0 & ismember(type,varargin{1})',{'TEST_LURE'});
      
      % all source reponses (targets and lures) are conditionalized on
      % getting a recognition hit
      %
      % source targets (for Hs and Ms) and lures (for CRs and FAs) come
      % from the source coding of target items, as explained in
      % cosi_createEvents.m
      src_targEv = filterStruct(rec_targEv,'src_isTarg == 1 & rec_correct == 1');
      src_lureEv = filterStruct(rec_targEv,'src_isTarg == 0 & rec_correct == 1');
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % TEST-PRESENTATION EVENTS %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % recognition events only
      rec_h = filterStruct(rec_targEv,'rec_correct == 1');
      rec_cr = filterStruct(rec_lureEv,'rec_correct == 1');
      if rejectArt == 0
        rec_m = filterStruct(rec_targEv,'rec_correct == 0');
        rec_fa = filterStruct(rec_lureEv,'rec_correct == 0');
      end
      
      % source events only
      src_h = filterStruct(src_targEv,'src_correct == 1');
      src_m = filterStruct(src_targEv,'src_correct == 0');
      % source CRs and FAs come from the source coding of target items, as
      % explained in sosi_createEvents.m
      src_cr = filterStruct(src_lureEv,'src_correct == 1');
      src_fa = filterStruct(src_lureEv,'src_correct == 0');
      
      % recognition hit - old items called "old"
      %
      % collapsed across source hits and CRs
      rec_h_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1 & ismember(rkn_resp,varargin{1})',{'REMEMBER_SOURCE','REMEMBER_OTHER','KNOW'});
      rec_h_srcCor_rs = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1 & ismember(rkn_resp,varargin{1})',{'REMEMBER_SOURCE'});
      rec_h_srcCor_ro = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1 & ismember(rkn_resp,varargin{1})',{'REMEMBER_OTHER'});
      rec_h_srcCor_k = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1 & ismember(rkn_resp,varargin{1})',{'KNOW'});
      % collapsed across source misses and FAs
      rec_h_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0 & ismember(rkn_resp,varargin{1})',{'REMEMBER_SOURCE','REMEMBER_OTHER','KNOW'});
      rec_h_srcInc_rs = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0 & ismember(rkn_resp,varargin{1})',{'REMEMBER_SOURCE'});
      rec_h_srcInc_ro = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0 & ismember(rkn_resp,varargin{1})',{'REMEMBER_OTHER'});
      rec_h_srcInc_k = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0 & ismember(rkn_resp,varargin{1})',{'KNOW'});
      
      % recognition numbers collapsed across source accuracy
      numEv.rec_h_rs(sub,ses) = length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs);
      numEv.rec_h_ro(sub,ses) = length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro);
      numEv.rec_h_k(sub,ses) = length(rec_h_srcCor_k) + length(rec_h_srcInc_k);
      
      % recognition correct rejection - new items called "new"
      rec_cr_sure = filterStruct(rec_cr,'src_correct == 1 & ismember(rkn_resp,varargin{1})',{'SURE'});
      rec_cr_maybe = filterStruct(rec_cr,'src_correct == 1 & ismember(rkn_resp,varargin{1})',{'MAYBE'});
      
      if rejectArt == 0
        % recognition miss - old items called "new"
        rec_m_sure = filterStruct(rec_m,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'SURE'});
        rec_m_maybe = filterStruct(rec_m,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'MAYBE'});
        
        % recognition false alarm - new items called "old"
        rec_fa_rs = filterStruct(rec_fa,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'REMEMBER_SOURCE'});
        rec_fa_ro = filterStruct(rec_fa,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'REMEMBER_OTHER'});
        rec_fa_k = filterStruct(rec_fa,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'KNOW'});
      end
      
      %% RATES
      
      % accuracy rates
      %
      % TODO: for some reason I wasn't calculating rates.rec_h and
      % rates.rec_cr when rejecting artifacts. Not sure why.
      rates.rec_h(sub,ses) = length(rec_h) / length(rec_targEv);
      rates.rec_cr(sub,ses) = length(rec_cr) / length(rec_lureEv);
      if rejectArt == 0
        rates.rec_m(sub,ses) = length(rec_m) / length(rec_targEv);
        rates.rec_fa(sub,ses) = length(rec_fa) / length(rec_lureEv);
      end
      rates.src_h(sub,ses) = length(src_h) / length(src_targEv);
      rates.src_m(sub,ses) = length(src_m) / length(src_targEv);
      rates.src_cr(sub,ses) = length(src_cr) / length(src_lureEv);
      rates.src_fa(sub,ses) = length(src_fa) / length(src_lureEv);
      
      % correct for rates = 1 or 0
      %
      % rec
      if rates.rec_h(sub,ses) == 1
        rates.rec_h(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      elseif rates.rec_h(sub,ses) == 0
        rates.rec_h(sub,ses) = 1/(2*length(rec_targEv));
      end
      if rates.rec_cr(sub,ses) == 1
        rates.rec_cr(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
      elseif rates.rec_cr(sub,ses) == 0
        rates.rec_cr(sub,ses) = 1/(2*length(rec_lureEv));
      end
      if rejectArt == 0
        if rates.rec_m(sub,ses) == 1
          rates.rec_m(sub,ses) = 1 - (1/(2*length(rec_targEv)));
        elseif rates.rec_m(sub,ses) == 0
          rates.rec_m(sub,ses) = 1/(2*length(rec_targEv));
        end
        if rates.rec_fa(sub,ses) == 1
          rates.rec_fa(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
        elseif rates.rec_fa(sub,ses) == 0
          rates.rec_fa(sub,ses) = 1/(2*length(rec_lureEv));
        end
      end
      % src
      if rates.src_h(sub,ses) == 1
        rates.src_h(sub,ses) = 1 - (1/(2*length(src_targEv)));
      elseif rates.src_h(sub,ses) == 0
        rates.src_h(sub,ses) = 1/(2*length(src_targEv));
      end
      if rates.src_m(sub,ses) == 1
        rates.src_m(sub,ses) = 1 - (1/(2*length(src_targEv)));
      elseif rates.src_m(sub,ses) == 0
        rates.src_m(sub,ses) = 1/(2*length(src_targEv));
      end
      if rates.src_cr(sub,ses) == 1
        rates.src_cr(sub,ses) = 1 - (1/(2*length(src_lureEv)));
      elseif rates.src_cr(sub,ses) == 0
        rates.src_cr(sub,ses) = 1/(2*length(src_lureEv));
      end
      if rates.src_fa(sub,ses) == 1
        rates.src_fa(sub,ses) = 1 - (1/(2*length(src_lureEv)));
      elseif rates.src_fa(sub,ses) == 0
        rates.src_fa(sub,ses) = 1/(2*length(src_lureEv));
      end
      
      % d primes
      if rejectArt == 0
        acc.item_dp(sub,ses) = norminv(rates.rec_h(sub,ses),0,1) - norminv(rates.rec_fa(sub,ses),0,1);
      end
      acc.source_dp(sub,ses) = norminv(rates.src_h(sub,ses),0,1) - norminv(rates.src_fa(sub,ses),0,1);
      
      % response bias (criterion) from Macmillan & Creelman p. 29
      if rejectArt == 0
        acc.item_c(sub,ses) = -0.5 * (norminv(rates.rec_h(sub,ses),0,1) + norminv(rates.rec_fa(sub,ses),0,1));
      end
      acc.source_c(sub,ses) = -0.5 * (norminv(rates.src_h(sub,ses),0,1) + norminv(rates.src_fa(sub,ses),0,1));
      
      if rejectArt == 0
        % From Mecklinger et al. (2007) and Corwin (1994)
        %
        % discrimination index
        acc.Pr(sub,ses) = rates.rec_h(sub,ses) - rates.rec_fa(sub,ses);
        % response bias index
        acc.Br(sub,ses) = rates.rec_fa(sub,ses) / (1 - acc.Pr(sub,ses));
      end
      
      % accuracy rates for recognition hits
      rates.rec_h_srcCor(sub,ses) = length(rec_h_srcCor) / length(rec_h);
      rates.rec_h_srcInc(sub,ses) = length(rec_h_srcInc) / length(rec_h);
      
      % rates by confidence
      rates.rec_h_srcCor_rs(sub,ses) = length(rec_h_srcCor_rs) / length(rec_targEv);
      rates.rec_h_srcCor_ro(sub,ses) = length(rec_h_srcCor_ro) / length(rec_targEv);
      rates.rec_h_srcCor_k(sub,ses) = length(rec_h_srcCor_k) / length(rec_targEv);
      rates.rec_h_srcInc_rs(sub,ses) = length(rec_h_srcInc_rs) / length(rec_targEv);
      rates.rec_h_srcInc_ro(sub,ses) = length(rec_h_srcInc_ro) / length(rec_targEv);
      rates.rec_h_srcInc_k(sub,ses) = length(rec_h_srcInc_k) / length(rec_targEv);
      rates.rec_h_rs(sub,ses) = numEv.rec_h_rs(sub,ses) / length(rec_targEv);
      rates.rec_h_ro(sub,ses) = numEv.rec_h_ro(sub,ses) / length(rec_targEv);
      rates.rec_h_k(sub,ses) = numEv.rec_h_k(sub,ses) / length(rec_targEv);
      rates.rec_cr_sure(sub,ses) = length(rec_cr_sure) / length(rec_lureEv);
      rates.rec_cr_maybe(sub,ses) = length(rec_cr_maybe) / length(rec_lureEv);
      if rejectArt == 0
        rates.rec_m_sure(sub,ses) = length(rec_m_sure) / length(rec_targEv);
        rates.rec_m_maybe(sub,ses) = length(rec_m_maybe) / length(rec_targEv);
        rates.rec_fa_rs(sub,ses) = length(rec_fa_rs) / length(rec_lureEv);
        rates.rec_fa_ro(sub,ses) = length(rec_fa_ro) / length(rec_lureEv);
        rates.rec_fa_k(sub,ses) = length(rec_fa_k) / length(rec_lureEv);
      end
      
      %% fix for when rates are 1 or 0 (Macmillan and Creelman, 1991)
      %
      % srcCor
      if rates.rec_h_srcCor_rs(sub,ses) == 1
        rates.rec_h_srcCor_rs(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      elseif rates.rec_h_srcCor_rs(sub,ses) == 0
        rates.rec_h_srcCor_rs(sub,ses) = 1/(2*length(rec_targEv));
      end
      if rates.rec_h_srcCor_ro(sub,ses) == 1
        rates.rec_h_srcCor_ro(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      elseif rates.rec_h_srcCor_ro(sub,ses) == 0
        rates.rec_h_srcCor_ro(sub,ses) = 1/(2*length(rec_targEv));
      end
      if rates.rec_h_srcCor_k(sub,ses) == 1
        rates.rec_h_srcCor_k(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      elseif rates.rec_h_srcCor_k(sub,ses) == 0
        rates.rec_h_srcCor_k(sub,ses) = 1/(2*length(rec_targEv));
      end
      % srcInc
      if rates.rec_h_srcInc_rs(sub,ses) == 1
        rates.rec_h_srcInc_rs(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      elseif rates.rec_h_srcInc_rs(sub,ses) == 0
        rates.rec_h_srcInc_rs(sub,ses) = 1/(2*length(rec_targEv));
      end
      if rates.rec_h_srcInc_ro(sub,ses) == 1
        rates.rec_h_srcInc_ro(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      elseif rates.rec_h_srcInc_ro(sub,ses) == 0
        rates.rec_h_srcInc_ro(sub,ses) = 1/(2*length(rec_targEv));
      end
      if rates.rec_h_srcInc_k(sub,ses) == 1
        rates.rec_h_srcInc_k(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      elseif rates.rec_h_srcInc_k(sub,ses) == 0
        rates.rec_h_srcInc_k(sub,ses) = 1/(2*length(rec_targEv));
      end
      % cr
      if rates.rec_cr_sure(sub,ses) == 1
        rates.rec_cr_sure(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
      elseif rates.rec_cr_sure(sub,ses) == 0
        rates.rec_cr_sure(sub,ses) = 1/(2*length(rec_lureEv));
      end
      if rates.rec_cr_maybe(sub,ses) == 1
        rates.rec_cr_maybe(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
      elseif rates.rec_cr_maybe(sub,ses) == 0
        rates.rec_cr_maybe(sub,ses) = 1/(2*length(rec_lureEv));
      end
      if rejectArt == 0
        % miss
        if rates.rec_m_sure(sub,ses) == 1
          rates.rec_m_sure(sub,ses) = 1 - (1/(2*length(rec_targEv)));
        elseif rates.rec_m_sure(sub,ses) == 0
          rates.rec_m_sure(sub,ses) = 1/(2*length(rec_targEv));
        end
        if rates.rec_m_maybe(sub,ses) == 1
          rates.rec_m_maybe(sub,ses) = 1 - (1/(2*length(rec_targEv)));
        elseif rates.rec_m_maybe(sub,ses) == 0
          rates.rec_m_maybe(sub,ses) = 1/(2*length(rec_targEv));
        end
        % fa
        if rates.rec_fa_rs(sub,ses) == 1
          rates.rec_fa_rs(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
        elseif rates.rec_fa_rs(sub,ses) == 0
          rates.rec_fa_rs(sub,ses) = 1/(2*length(rec_lureEv));
        end
        if rates.rec_fa_ro(sub,ses) == 1
          rates.rec_fa_ro(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
        elseif rates.rec_fa_ro(sub,ses) == 0
          rates.rec_fa_ro(sub,ses) = 1/(2*length(rec_lureEv));
        end
        if rates.rec_fa_k(sub,ses) == 1
          rates.rec_fa_k(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
        elseif rates.rec_fa_k(sub,ses) == 0
          rates.rec_fa_k(sub,ses) = 1/(2*length(rec_lureEv));
        end
      end
      
      %% Independent Remember--Know Familiarity rate
      rates.irk_Fam(sub,ses) = rates.rec_h_k(sub,ses) / (1 - (rates.rec_h_rs(sub,ses) + rates.rec_h_ro(sub,ses)));
      
      %% Independent Remember--Know Familiarity rates - I think this is wrong
      rates.irk_correct_Fam(sub,ses) = rates.rec_h_srcCor_k(sub,ses) / (1 - (rates.rec_h_srcCor_rs(sub,ses) + rates.rec_h_srcCor_ro(sub,ses)));
      rates.irk_incorrect_Fam(sub,ses) = rates.rec_h_srcInc_k(sub,ses) / (1 - (rates.rec_h_srcInc_rs(sub,ses) + rates.rec_h_srcInc_ro(sub,ses)));
      
      % Independent Remember--Know Recollect Other rates
      %irk_correct_RcO(sub,ses) = rates.rec_h_srcCor_ro(sub,ses) / (1 - (rates.rec_h_srcCor_rs(sub,ses)));
      %irk_incorrect_RcO(sub,ses) = rates.rec_h_srcInc_ro(sub,ses) / (1 - (rates.rec_h_srcInc_rs(sub,ses)));
      
      %% within-response rates
      %
      % rs
      if isempty(rec_h_srcCor_rs)
        rates_wir.rec_h_srcCor_rs(sub,ses) = 1/(2*length(rec_targEv));
      else
        rates_wir.rec_h_srcCor_rs(sub,ses) = length(rec_h_srcCor_rs) / (length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs));
      end
      if rates_wir.rec_h_srcCor_rs(sub,ses) == 1
        rates_wir.rec_h_srcCor_rs(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      end
      if isempty(rec_h_srcInc_rs)
        rates_wir.rec_h_srcInc_rs(sub,ses) = 1/(2*length(rec_targEv));
      else
        rates_wir.rec_h_srcInc_rs(sub,ses) = length(rec_h_srcInc_rs) / (length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs));
      end
      if rates_wir.rec_h_srcInc_rs(sub,ses) == 1
        rates_wir.rec_h_srcInc_rs(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      end
      % ro
      if isempty(rec_h_srcCor_ro)
        rates_wir.rec_h_srcCor_ro(sub,ses) = 1/(2*length(rec_targEv));
      else
        rates_wir.rec_h_srcCor_ro(sub,ses) = length(rec_h_srcCor_ro) / (length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro));
      end
      if rates_wir.rec_h_srcCor_ro(sub,ses) == 1
        rates_wir.rec_h_srcCor_ro(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      end
      if isempty(rec_h_srcInc_ro)
        rates_wir.rec_h_srcInc_ro(sub,ses) = 1/(2*length(rec_targEv));
      else
        rates_wir.rec_h_srcInc_ro(sub,ses) = length(rec_h_srcInc_ro) / (length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro));
      end
      if rates_wir.rec_h_srcInc_ro(sub,ses) == 1
        rates_wir.rec_h_srcInc_ro(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      end
      % k
      if isempty(rec_h_srcCor_k)
        rates_wir.rec_h_srcCor_k(sub,ses) = 1/(2*length(rec_targEv));
      else
        rates_wir.rec_h_srcCor_k(sub,ses) = length(rec_h_srcCor_k) / (length(rec_h_srcCor_k) + length(rec_h_srcInc_k));
      end
      if rates_wir.rec_h_srcCor_k(sub,ses) == 1
        rates_wir.rec_h_srcCor_k(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      end
      if isempty(rec_h_srcInc_k)
        rates_wir.rec_h_srcInc_k(sub,ses) = 1/(2*length(rec_targEv));
      else
        rates_wir.rec_h_srcInc_k(sub,ses) = length(rec_h_srcInc_k) / (length(rec_h_srcCor_k) + length(rec_h_srcInc_k));
      end
      if rates_wir.rec_h_srcInc_k(sub,ses) == 1
        rates_wir.rec_h_srcInc_k(sub,ses) = 1 - (1/(2*length(rec_targEv)));
      end
      % cr
      if isempty(rec_cr_sure)
        rates_wir.rec_cr_sure(sub,ses) = 1/(2*length(rec_lureEv));
      else
        rates_wir.rec_cr_sure(sub,ses) = length(rec_cr_sure) / (length(rec_cr_sure) + length(rec_cr_maybe));
      end
      if rates_wir.rec_cr_sure(sub,ses) == 1
        rates_wir.rec_cr_sure(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
      end
      if isempty(rec_cr_maybe)
        rates_wir.rec_cr_maybe(sub,ses) = 1/(2*length(rec_lureEv));
      else
        rates_wir.rec_cr_maybe(sub,ses) = length(rec_cr_maybe) / (length(rec_cr_sure) + length(rec_cr_maybe));
      end
      if rates_wir.rec_cr_maybe(sub,ses) == 1
        rates_wir.rec_cr_maybe(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
      end
      if rejectArt == 0
        % miss
        if isempty(rec_m_sure)
          rates_wir.rec_m_sure(sub,ses) = 1/(2*length(rec_targEv));
        else
          rates_wir.rec_m_sure(sub,ses) = length(rec_m_sure) / (length(rec_m_sure) + length(rec_m_maybe));
        end
        if rates_wir.rec_m_sure(sub,ses) == 1
          rates_wir.rec_m_sure(sub,ses) = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_m_maybe)
          rates_wir.rec_m_maybe(sub,ses) = 1/(2*length(rec_targEv));
        else
          rates_wir.rec_m_maybe(sub,ses) = length(rec_m_maybe) / (length(rec_m_sure) + length(rec_m_maybe));
        end
        if rates_wir.rec_m_maybe(sub,ses) == 1
          rates_wir.rec_m_maybe(sub,ses) = 1 - (1/(2*length(rec_targEv)));
        end
        % fa
        if isempty(rec_fa_rs)
          rates_wir.rec_fa_rs(sub,ses) = 1/(2*length(rec_lureEv));
        else
          rates_wir.rec_fa_rs(sub,ses) = length(rec_fa_rs) / (length(rec_fa_rs) + length(rec_fa_ro) + length(rec_fa_k));
        end
        if rates_wir.rec_fa_rs(sub,ses) == 1
          rates_wir.rec_fa_rs(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
        end
        if isempty(rec_fa_ro)
          rates_wir.rec_fa_ro(sub,ses) = 1/(2*length(rec_lureEv));
        else
          rates_wir.rec_fa_ro(sub,ses) = length(rec_fa_ro) / (length(rec_fa_rs) + length(rec_fa_ro) + length(rec_fa_k));
        end
        if rates_wir.rec_fa_ro(sub,ses) == 1
          rates_wir.rec_fa_ro(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
        end
        if isempty(rec_fa_k)
          rates_wir.rec_fa_k(sub,ses) = 1/(2*length(rec_lureEv));
        else
          rates_wir.rec_fa_k(sub,ses) = length(rec_fa_k) / (length(rec_fa_rs) + length(rec_fa_ro) + length(rec_fa_k));
        end
        if rates_wir.rec_fa_k(sub,ses) == 1
          rates_wir.rec_fa_k(sub,ses) = 1 - (1/(2*length(rec_lureEv)));
        end
      end
      
      %% REACTION TIMES
      
      % test items
      rt.rec_h(sub,ses) = mean(getStructField(rec_h,'src_rt'));
      rt.rec_cr(sub,ses) = mean(getStructField(rec_cr,'src_rt'));
      if rejectArt == 0
        rt.rec_m(sub,ses) = mean(getStructField(rec_m,'src_rt'));
        rt.rec_fa(sub,ses) = mean(getStructField(rec_fa,'src_rt'));
      end
      
      rt.src_h(sub,ses) = mean(getStructField(src_h,'src_rt'));
      rt.src_m(sub,ses) = mean(getStructField(src_m,'src_rt'));
      rt.src_cr(sub,ses) = mean(getStructField(src_cr,'src_rt'));
      rt.src_fa(sub,ses) = mean(getStructField(src_fa,'src_rt'));
      
      rt.rec_h_srcCor(sub,ses) = mean(getStructField(rec_h_srcCor,'src_rt'));
      rt.rec_h_srcInc(sub,ses) = mean(getStructField(rec_h_srcInc,'src_rt'));
      
      rt.rec_h_srcCor_rs(sub,ses) = mean(getStructField(rec_h_srcCor_rs,'src_rt'));
      rt.rec_h_srcCor_ro(sub,ses) = mean(getStructField(rec_h_srcCor_ro,'src_rt'));
      rt.rec_h_srcCor_k(sub,ses) = mean(getStructField(rec_h_srcCor_k,'src_rt'));
      rt.rec_h_srcInc_rs(sub,ses) = mean(getStructField(rec_h_srcInc_rs,'src_rt'));
      rt.rec_h_srcInc_ro(sub,ses) = mean(getStructField(rec_h_srcInc_ro,'src_rt'));
      rt.rec_h_srcInc_k(sub,ses) = mean(getStructField(rec_h_srcInc_k,'src_rt'));
      
      rt.rec_cr_sure(sub,ses) = mean(getStructField(rec_cr_sure,'src_rt'));
      rt.rec_cr_maybe(sub,ses) = mean(getStructField(rec_cr_maybe,'src_rt'));
      if rejectArt == 0
        rt.rec_m_sure(sub,ses) = mean(getStructField(rec_m_sure,'src_rt'));
        rt.rec_m_maybe(sub,ses) = mean(getStructField(rec_m_maybe,'src_rt'));
        rt.rec_fa_rs(sub,ses) = mean(getStructField(rec_fa_rs,'src_rt'));
        rt.rec_fa_ro(sub,ses) = mean(getStructField(rec_fa_ro,'src_rt'));
        rt.rec_fa_k(sub,ses) = mean(getStructField(rec_fa_k,'src_rt'));
      end
      
      %% print some info
      fprintf('\n\tRecognition accuracy info:\n');
      fprintf('\tHits: %.4f\n',rates.rec_h(sub,ses));
      if rejectArt == 0
        fprintf('\tMisses: %.4f\n',rates.rec_m(sub,ses));
      end
      fprintf('\tCorrect rejections: %.4f\n',rates.rec_cr(sub,ses));
      if rejectArt == 0
        fprintf('\tFalse alarms: %.4f\n',rates.rec_fa(sub,ses));
      end
      
      fprintf('\n\tSource accuracy info:\n');
      fprintf('\tHits: %.4f\n',rates.src_h(sub,ses));
      fprintf('\tMisses: %.4f\n',rates.src_m(sub,ses));
      fprintf('\tCorrect rejections: %.4f\n',rates.src_cr(sub,ses));
      fprintf('\tFalse alarms: %.4f\n\n',rates.src_fa(sub,ses));
      
      
      % ===================================================================
      
      if rejectArt == 1
        % count the artifacts
        
        % get all the targets (old study items presented at test)
        rec_targEv_art = filterStruct(subSesEv_all,'nsArt == 1 & rec_isTarg == 1 & rec_correct == 1 & ismember(type,varargin{1})',{'TEST_TARGET'});
        % get all the lures (new test items)
        rec_lureEv_art = filterStruct(subSesEv_all,'nsArt == 1 & rec_isTarg == 0 & rec_correct == 1 & ismember(type,varargin{1})',{'TEST_LURE'});
        
        % get all the targets (old study items presented at test)
        rec_targEv_all = filterStruct(subSesEv_all,'rec_isTarg == 1 & rec_correct == 1 & ismember(type,varargin{1})',{'TEST_TARGET'});
        % get all the lures (new test items)
        rec_lureEv_all = filterStruct(subSesEv_all,'rec_isTarg == 0 & rec_correct == 1 & ismember(type,varargin{1})',{'TEST_LURE'});
        
        numEv.art(sub,ses) = length(rec_targEv_art) + length(rec_lureEv_art);
        
        numEv.total(sub,ses) = length(rec_targEv_all) + length(rec_lureEv_all);
      end
      
      % raw numbers of targs, lures
      numEv.rec_targEv(sub,ses) = length(rec_targEv);
      numEv.rec_lureEv(sub,ses) = length(rec_lureEv);
      numEv.src_targEv(sub,ses) = length(src_targEv);
      numEv.src_lureEv(sub,ses) = length(src_lureEv);
      % raw numbers for accuracy
      numEv.rec_h(sub,ses) = length(rec_h);
      numEv.rec_cr(sub,ses) = length(rec_cr);
      if rejectArt == 0
        numEv.rec_m(sub,ses) = length(rec_m);
        numEv.rec_fa(sub,ses) = length(rec_fa);
      end
      numEv.src_h(sub,ses) = length(src_h);
      numEv.src_cr(sub,ses) = length(src_cr);
      numEv.src_m(sub,ses) = length(src_m);
      numEv.src_fa(sub,ses) = length(src_fa);
      % raw numbers for recognition hits
      numEv.rec_h_srcCor(sub,ses) = length(rec_h_srcCor);
      numEv.rec_h_srcInc(sub,ses) = length(rec_h_srcInc);
      % raw numbers by confidence
      numEv.rec_h_srcCor_rs(sub,ses) = length(rec_h_srcCor_rs);
      numEv.rec_h_srcCor_ro(sub,ses) = length(rec_h_srcCor_ro);
      numEv.rec_h_srcCor_k(sub,ses) = length(rec_h_srcCor_k);
      numEv.rec_h_srcInc_rs(sub,ses) = length(rec_h_srcInc_rs);
      numEv.rec_h_srcInc_ro(sub,ses) = length(rec_h_srcInc_ro);
      numEv.rec_h_srcInc_k(sub,ses) = length(rec_h_srcInc_k);
      numEv.rec_m_sure(sub,ses) = length(rec_m_sure);
      numEv.rec_m_maybe(sub,ses) = length(rec_m_maybe);
      numEv.rec_cr_sure(sub,ses) = length(rec_cr_sure);
      numEv.rec_cr_maybe(sub,ses) = length(rec_cr_maybe);
      numEv.rec_fa_rs(sub,ses) = length(rec_fa_rs);
      numEv.rec_fa_ro(sub,ses) = length(rec_fa_ro);
      numEv.rec_fa_k(sub,ses) = length(rec_fa_k);
      
    end % ses
  end % sub
  
  if averageSes == 1
    
    numEvFN = fieldnames(numEv);
    ratesFN = fieldnames(rates);
    rates_wirFN = fieldnames(rates_wir);
    accFN = fieldnames(acc);
    rtFN = fieldnames(rt);
    
    % sum or average across sessions
    for f = 1:length(numEvFN)
      numEv.(numEvFN{f}) = sum(numEv.(numEvFN{f}),2);
    end
    for f = 1:length(ratesFN)
      rates.(ratesFN{f}) = mean(rates.(ratesFN{f}),2);
    end
    for f = 1:length(rates_wirFN)
      rates_wir.(rates_wirFN{f}) = mean(rates_wir.(rates_wirFN{f}),2);
    end
    for f = 1:length(accFN)
      acc.(accFN{f}) = mean(acc.(accFN{f}),2);
    end
    for f = 1:length(rtFN)
      rt.(rtFN{f}) = mean(rt.(rtFN{f}),2);
    end
  end
  
  for sub = 1:length(subjects)
    for ses = 1:size(numEv.rec_targEv,2)
      if rejectArt == 0
        tableData = [...
          numEv.rec_targEv(sub,ses),numEv.rec_lureEv(sub,ses),numEv.src_targEv(sub,ses),numEv.src_lureEv(sub,ses),... % raw numbers of targs, lures
          acc.item_dp(sub,ses),acc.source_dp(sub,ses),acc.item_c(sub,ses),acc.source_c(sub,ses),... % accuracy
          numEv.rec_h(sub,ses),numEv.rec_m(sub,ses),numEv.rec_cr(sub,ses),numEv.rec_fa(sub,ses),numEv.src_h(sub,ses),numEv.src_m(sub,ses),numEv.src_cr(sub,ses),numEv.src_fa(sub,ses),... % raw numbers for accuracy
          rates.rec_h(sub,ses),rates.rec_m(sub,ses),rates.rec_cr(sub,ses),rates.rec_fa(sub,ses),rates.src_h(sub,ses),rates.src_m(sub,ses),rates.src_cr(sub,ses),rates.src_fa(sub,ses),... % accuracy rates
          numEv.rec_h_srcCor(sub,ses),numEv.rec_h_srcInc(sub,ses),... % raw numbers for recognition hits
          rates.rec_h_srcCor(sub,ses),rates.rec_h_srcInc(sub,ses),... % accuracy rates for recognition hits
          numEv.rec_h_rs(sub,ses),numEv.rec_h_ro(sub,ses),numEv.rec_h_k(sub,ses),... % raw numbers for recognition hit response collapsed across source accuracy
          numEv.rec_h_srcCor_rs(sub,ses),numEv.rec_h_srcCor_ro(sub,ses),numEv.rec_h_srcCor_k(sub,ses),numEv.rec_h_srcInc_rs(sub,ses),numEv.rec_h_srcInc_ro(sub,ses),numEv.rec_h_srcInc_k(sub,ses),numEv.rec_m_sure(sub,ses),numEv.rec_m_maybe(sub,ses),numEv.rec_cr_sure(sub,ses),numEv.rec_cr_maybe(sub,ses),numEv.rec_fa_rs(sub,ses),numEv.rec_fa_ro(sub,ses),numEv.rec_fa_k(sub,ses),... % raw numbers by confidence
          rates.rec_h_srcCor_rs(sub,ses),rates.rec_h_srcCor_ro(sub,ses),rates.rec_h_srcCor_k(sub,ses),rates.rec_h_srcInc_rs(sub,ses),rates.rec_h_srcInc_ro(sub,ses),rates.rec_h_srcInc_k(sub,ses),rates.rec_m_sure(sub,ses),rates.rec_m_maybe(sub,ses),rates.rec_cr_sure(sub,ses),rates.rec_cr_maybe(sub,ses),rates.rec_fa_rs(sub,ses),rates.rec_fa_ro(sub,ses),rates.rec_fa_k(sub,ses),... % rates by confidence
          rates.irk_Fam(sub,ses),rates.irk_correct_Fam(sub,ses),rates.irk_incorrect_Fam(sub,ses),... % Independent Remember--Know Familiarity rates
          rates_wir.rec_h_srcCor_rs(sub,ses),rates_wir.rec_h_srcInc_rs(sub,ses),rates_wir.rec_h_srcCor_ro(sub,ses),rates_wir.rec_h_srcInc_ro(sub,ses),rates_wir.rec_h_srcCor_k(sub,ses),rates_wir.rec_h_srcInc_k(sub,ses),rates_wir.rec_m_sure(sub,ses),rates_wir.rec_m_maybe(sub,ses),rates_wir.rec_cr_sure(sub,ses),rates_wir.rec_cr_maybe(sub,ses),rates_wir.rec_fa_rs(sub,ses),rates_wir.rec_fa_ro(sub,ses),rates_wir.rec_fa_k(sub,ses),... % within-response rates
          rt.rec_h(sub,ses),rt.rec_m(sub,ses),rt.rec_cr(sub,ses),rt.rec_fa(sub,ses),... % rec reaction times
          rt.src_h(sub,ses),rt.src_m(sub,ses),rt.src_cr(sub,ses),rt.src_fa(sub,ses),... % src reaction times
          rt.rec_h_srcCor(sub,ses),rt.rec_h_srcInc(sub,ses),... % recognition hits reaction times
          rt.rec_h_srcCor_rs(sub,ses),rt.rec_h_srcCor_ro(sub,ses),rt.rec_h_srcCor_k(sub,ses),rt.rec_h_srcInc_rs(sub,ses),rt.rec_h_srcInc_ro(sub,ses),rt.rec_h_srcInc_k(sub,ses),rt.rec_m_sure(sub,ses),rt.rec_m_maybe(sub,ses),rt.rec_cr_sure(sub,ses),rt.rec_cr_maybe(sub,ses),rt.rec_fa_rs(sub,ses),rt.rec_fa_ro(sub,ses),rt.rec_fa_k(sub,ses),... % rt by confidence
          acc.Pr(sub,ses),acc.Br(sub,ses),... % more accuracy
          ];
      elseif rejectArt == 1
        tableData = [...
          numEv.rec_targEv(sub,ses),numEv.rec_lureEv(sub,ses),numEv.src_targEv(sub,ses),numEv.src_lureEv(sub,ses),... % raw numbers of targs, lures
          acc.source_dp(sub,ses),acc.source_c(sub,ses),... % accuracy
          numEv.rec_h(sub,ses),numEv.rec_cr(sub,ses),numEv.src_h(sub,ses),numEv.src_m(sub,ses),numEv.src_cr(sub,ses),numEv.src_fa(sub,ses),... % raw numbers for accuracy
          rates.src_h(sub,ses),rates.src_m(sub,ses),rates.src_cr(sub,ses),rates.src_fa(sub,ses),... % accuracy rates
          numEv.rec_h_srcCor(sub,ses),numEv.rec_h_srcInc(sub,ses),... % raw numbers for recognition hits
          rates.rec_h_srcCor(sub,ses),rates.rec_h_srcInc(sub,ses),... % accuracy rates for recognition hits
          numEv.rec_h_rs(sub,ses),numEv.rec_h_ro(sub,ses),numEv.rec_h_k(sub,ses),... % raw numbers for recognition hit response collapsed across source accuracy
          numEv.rec_h_srcCor_rs(sub,ses),numEv.rec_h_srcCor_ro(sub,ses),numEv.rec_h_srcCor_k(sub,ses),numEv.rec_h_srcInc_rs(sub,ses),numEv.rec_h_srcInc_ro(sub,ses),numEv.rec_h_srcInc_k(sub,ses),numEv.rec_cr_sure(sub,ses),numEv.rec_cr_maybe(sub,ses),... % raw numbers by confidence
          rates.rec_h_srcCor_rs(sub,ses),rates.rec_h_srcCor_ro(sub,ses),rates.rec_h_srcCor_k(sub,ses),rates.rec_h_srcInc_rs(sub,ses),rates.rec_h_srcInc_ro(sub,ses),rates.rec_h_srcInc_k(sub,ses),rates.rec_cr_sure(sub,ses),rates.rec_cr_maybe(sub,ses),... % rates by confidence
          rates.irk_Fam(sub,ses),rates.irk_correct_Fam(sub,ses),rates.irk_incorrect_Fam(sub,ses),... % Independent Remember--Know Familiarity rates
          rates_wir.rec_h_srcCor_rs(sub,ses),rates_wir.rec_h_srcInc_rs(sub,ses),rates_wir.rec_h_srcCor_ro(sub,ses),rates_wir.rec_h_srcInc_ro(sub,ses),rates_wir.rec_h_srcCor_k(sub,ses),rates_wir.rec_h_srcInc_k(sub,ses),rates_wir.rec_cr_sure(sub,ses),rates_wir.rec_cr_maybe(sub,ses)... % within-response rates
          rt.rec_h(sub,ses),rt.rec_cr(sub,ses),... % rec reaction times
          rt.src_h(sub,ses),rt.src_m(sub,ses),rt.src_cr(sub,ses),rt.src_fa(sub,ses),... % src reaction times
          rt.rec_h_srcCor(sub,ses),rt.rec_h_srcInc(sub,ses),... % recognition hits reaction times
          rt.rec_h_srcCor_rs(sub,ses),rt.rec_h_srcCor_ro(sub,ses),rt.rec_h_srcCor_k(sub,ses),rt.rec_h_srcInc_rs(sub,ses),rt.rec_h_srcInc_ro(sub,ses),rt.rec_h_srcInc_k(sub,ses),rt.rec_cr_sure(sub,ses),rt.rec_cr_maybe(sub,ses),... % rt by confidence
          ];
      end
      
      if saveFiles == 1
        % print sub and ses
        if averageSes == 1
          allSes = sessions{1}(end);
          if length(sessions) > 1
            for i = 2:length(sessions)
              allSes = cat(2,sprintf('%s_%s',allSes,sessions{i}(end)));
            end
          end
          fprintf(outfile,'%s,%s',subjects{sub},allSes);
        else
          fprintf(outfile,'%s,%d',subjects{sub},ses-1);
        end
        % format for tableData and print
        tableDataStr = repmat(',%.4f',1,length(tableData));
        fprintf(outfile,tableDataStr,tableData);
        %for tabData = 1:length(tableData)
        %  fprintf(outfile,'%d\t',tableData(tabData));
        %end
        if rejectArt == 1
          artStr = sprintf(',%d,%d,%.4f',numEv.art(sub,ses),numEv.total(sub,ses),(numEv.art(sub,ses) / numEv.total(sub,ses)));
          fprintf(outfile,'%s\n',artStr);
        else
          fprintf(outfile,'\n');
        end
      end
      
    end % ses
  end % sub
  
  
  if saveFiles == 1
    fprintf('Saving %s...',eventsTable);
    fclose(outfile);
    fprintf('Done.\n');
  end
  
end % subsets
