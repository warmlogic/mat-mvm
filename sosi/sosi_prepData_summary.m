function sosi_prepData_summary(rejectArt,saveFiles,dataroot)
% function sosi_prepData_summary(rejectArt,saveFiles,dataroot)
%
% Save column-formatted data as a csv file
%
% Looks for each subject's events.mat and saves 3 summary files in dataroot
%

expName = 'SOSI';

subsets = {''};

if nargin < 3
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
  if nargin < 2
    saveFiles = 1;
    if nargin < 1
      rejectArt = 0;
    end
  end
end

% % behavioral
% subjects = {
%   'SOSI001';
%   'SOSI002';
%   'SOSI003';
%   'SOSI004';
%   'SOSI005';
%   'SOSI006';
%   'SOSI007';
%   'SOSI008';
%   'SOSI009';
%   'SOSI010';
%   'SOSI011';
%   'SOSI012';
%   'SOSI013';
%   'SOSI014';
%   'SOSI015';
%   'SOSI016';
%   'SOSI017';
%   'SOSI018';
%   'SOSI020';
%   };

subjects = {
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
  'SOSI020';
  'SOSI019';
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

sessions = {'session_0'};

%otherFilter = 'src_rt > 4000 | rkn_rt > 4000';
rejectManual = 1;

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
    if rejectArt == 1
      eventsTable = fullfile(dataroot,[expName,'_summary',subsets{s},'_rejArt.csv']);
    else
      eventsTable = fullfile(dataroot,[expName,'_summary',subsets{s},'.csv']);
    end
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
  
  %% for each subject
  for sub = 1:length(subjects)
    fprintf('Getting data for %s...',subjects{sub});
    
    %% for each session
    for ses = 1:length(sessions)
      fprintf('%s...\n',sessions{ses});
      % set the subject events directory
      eventsDir_sub = fullfile(dataroot,subjects{sub},sessions{ses},'events');
      events = loadEvents(fullfile(eventsDir_sub,'events.mat'));
      
      % include events with defined study color name and that was tested
      events = filterStruct(events,'src_correct ~= -1');
      
      if rejectManual == 1
        %events = filterStruct(events,otherFilter);
        manualArt = {'manual','badc,manual','eyem,manual','eyeb,manual','eyem,eyeb,manual'};
        events = filterStruct(events,'~ismember(nsBadReason,varargin{1})',manualArt);
      end
      
      if rejectArt == 1
        % remove artifacts
        subSesEv_all = events;
        events = filterStruct(subSesEv_all,'nsArt == 0');
      end
      
      % get only the events for this session
      %subSesEv = filterStruct(events,'ismember(subject,varargin{1}) & ismember(session,varargin{2})',subs_str{sub},sess(ses));
      
      % get only the events for this session
      fprintf('All\n');
      subSesEv = events;
      
      if rejectArt == 0
        %% RAW NUMBERS
        
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
        % sosi_createEvents.m
        src_targEv = filterStruct(rec_targEv,'src_isTarg == 1 & rec_correct == 1');
        src_lureEv = filterStruct(rec_targEv,'src_isTarg == 0 & rec_correct == 1');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEST-PRESENTATION EVENTS %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % recognition events only
        rec_h = filterStruct(rec_targEv,'rec_correct == 1');
        rec_m = filterStruct(rec_targEv,'rec_correct == 0');
        rec_cr = filterStruct(rec_lureEv,'rec_correct == 1');
        rec_fa = filterStruct(rec_lureEv,'rec_correct == 0');
        
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
        rec_h_rs = length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs);
        rec_h_ro = length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro);
        rec_h_k = length(rec_h_srcCor_k) + length(rec_h_srcInc_k);
        
        % recognition miss - old items called "new"
        rec_m_sure = filterStruct(rec_m,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'SURE'});
        rec_m_maybe = filterStruct(rec_m,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'MAYBE'});
        
        % recognition correct rejection - new items called "new"
        rec_cr_sure = filterStruct(rec_cr,'src_correct == 1 & ismember(rkn_resp,varargin{1})',{'SURE'});
        rec_cr_maybe = filterStruct(rec_cr,'src_correct == 1 & ismember(rkn_resp,varargin{1})',{'MAYBE'});
        
        % recognition false alarm - new items called "old"
        rec_fa_rs = filterStruct(rec_fa,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'REMEMBER_SOURCE'});
        rec_fa_ro = filterStruct(rec_fa,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'REMEMBER_OTHER'});
        rec_fa_k = filterStruct(rec_fa,'src_correct == 0 & ismember(rkn_resp,varargin{1})',{'KNOW'});
        
        %% RATES
        
        % accuracy rates
        rec_hr = length(rec_h) / length(rec_targEv);
        rec_mr = length(rec_m) / length(rec_targEv);
        rec_crr = length(rec_cr) / length(rec_lureEv);
        rec_far = length(rec_fa) / length(rec_lureEv);
        src_hr = length(src_h) / length(src_targEv);
        src_mr = length(src_m) / length(src_targEv);
        src_crr = length(src_cr) / length(src_lureEv);
        src_far = length(src_fa) / length(src_lureEv);
        
        % correct for rates = 1 or 0
        %
        % rec
        if rec_hr == 1
          rec_hr = 1 - (1/(2*length(rec_targEv)));
        elseif rec_hr == 0
          rec_hr = 1/(2*length(rec_targEv));
        end
        if rec_mr == 1
          rec_mr = 1 - (1/(2*length(rec_targEv)));
        elseif rec_mr == 0
          rec_mr = 1/(2*length(rec_targEv));
        end
        if rec_crr == 1
          rec_crr = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_crr == 0
          rec_crr = 1/(2*length(rec_lureEv));
        end
        if rec_far == 1
          rec_far = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_far == 0
          rec_far = 1/(2*length(rec_lureEv));
        end
        % src
        if src_hr == 1
          src_hr = 1 - (1/(2*length(src_targEv)));
        elseif src_hr == 0
          src_hr = 1/(2*length(src_targEv));
        end
        if src_mr == 1
          src_mr = 1 - (1/(2*length(src_targEv)));
        elseif src_mr == 0
          src_mr = 1/(2*length(src_targEv));
        end
        if src_crr == 1
          src_crr = 1 - (1/(2*length(src_lureEv)));
        elseif src_crr == 0
          src_crr = 1/(2*length(src_lureEv));
        end
        if src_far == 1
          src_far = 1 - (1/(2*length(src_lureEv)));
        elseif src_far == 0
          src_far = 1/(2*length(src_lureEv));
        end

        % d primes
        item_dp = norminv(rec_hr,0,1) - norminv(rec_far,0,1);
        source_dp = norminv(src_hr,0,1) - norminv(src_far,0,1);
        
        % response bias (criterion) from Macmillan & Creelman p. 29
        item_c = -0.5 * (norminv(rec_hr,0,1) + norminv(rec_far,0,1));
        source_c = -0.5 * (norminv(src_hr,0,1) + norminv(src_far,0,1));
        
        % From Mecklinger et al. (2007) and Corwin (1994)
        %
        % discrimination index
        Pr = rec_hr - rec_far;
        % response bias index
        Br = rec_far / (1 - Pr);
        
        % accuracy rates for recognition hits
        rec_h_srcCor_r = length(rec_h_srcCor) / length(rec_h);
        rec_h_srcInc_r = length(rec_h_srcInc) / length(rec_h);
        
        % rates by confidence
        rec_h_srcCor_rs_r = length(rec_h_srcCor_rs) / length(rec_targEv);
        rec_h_srcCor_ro_r = length(rec_h_srcCor_ro) / length(rec_targEv);
        rec_h_srcCor_k_r = length(rec_h_srcCor_k) / length(rec_targEv);
        rec_h_srcInc_rs_r = length(rec_h_srcInc_rs) / length(rec_targEv);
        rec_h_srcInc_ro_r = length(rec_h_srcInc_ro) / length(rec_targEv);
        rec_h_srcInc_k_r = length(rec_h_srcInc_k) / length(rec_targEv);
        rec_h_rs_r = rec_h_rs / length(rec_targEv);
        rec_h_ro_r = rec_h_ro / length(rec_targEv);
        rec_h_k_r = rec_h_k / length(rec_targEv);
        rec_m_sure_r = length(rec_m_sure) / length(rec_targEv);
        rec_m_maybe_r = length(rec_m_maybe) / length(rec_targEv);
        rec_cr_sure_r = length(rec_cr_sure) / length(rec_lureEv);
        rec_cr_maybe_r = length(rec_cr_maybe) / length(rec_lureEv);
        rec_fa_rs_r = length(rec_fa_rs) / length(rec_lureEv);
        rec_fa_ro_r = length(rec_fa_ro) / length(rec_lureEv);
        rec_fa_k_r = length(rec_fa_k) / length(rec_lureEv);
        
        %% fix for when rates are 1 or 0 (Macmillan and Creelman, 1991)
        %
        % srcCor
        if rec_h_srcCor_rs_r == 1
          rec_h_srcCor_rs_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcCor_rs_r == 0
          rec_h_srcCor_rs_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcCor_ro_r == 1
          rec_h_srcCor_ro_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcCor_ro_r == 0
          rec_h_srcCor_ro_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcCor_k_r == 1
          rec_h_srcCor_k_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcCor_k_r == 0
          rec_h_srcCor_k_r = 1/(2*length(rec_targEv));
        end
        % srcInc
        if rec_h_srcInc_rs_r == 1
          rec_h_srcInc_rs_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcInc_rs_r == 0
          rec_h_srcInc_rs_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcInc_ro_r == 1
          rec_h_srcInc_ro_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcInc_ro_r == 0
          rec_h_srcInc_ro_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcInc_k_r == 1
          rec_h_srcInc_k_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcInc_k_r == 0
          rec_h_srcInc_k_r = 1/(2*length(rec_targEv));
        end
        % miss
        if rec_m_sure_r == 1
          rec_m_sure_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_sure_r == 0
          rec_m_sure_r = 1/(2*length(rec_targEv));
        end
        if rec_m_maybe_r == 1
          rec_m_maybe_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_maybe_r == 0
          rec_m_maybe_r = 1/(2*length(rec_targEv));
        end
        % cr
        if rec_cr_sure_r == 1
          rec_cr_sure_r = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_cr_sure_r == 0
          rec_cr_sure_r = 1/(2*length(rec_lureEv));
        end
        if rec_cr_maybe_r == 1
          rec_cr_maybe_r = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_cr_maybe_r == 0
          rec_cr_maybe_r = 1/(2*length(rec_lureEv));
        end
        % fa
        if rec_fa_rs_r == 1
          rec_fa_rs_r = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_fa_rs_r == 0
          rec_fa_rs_r = 1/(2*length(rec_lureEv));
        end
        if rec_fa_ro_r == 1
          rec_fa_ro_r = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_fa_ro_r == 0
          rec_fa_ro_r = 1/(2*length(rec_lureEv));
        end
        if rec_fa_k_r == 1
          rec_fa_k_r = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_fa_k_r == 0
          rec_fa_k_r = 1/(2*length(rec_lureEv));
        end
        
        %% Independent Remember--Know Familiarity rate
        irk_Fam = rec_h_k_r / (1 - (rec_h_rs_r + rec_h_ro_r));
        
        %% Independent Remember--Know Familiarity rates - I think this is wrong
        irk_correct_Fam = rec_h_srcCor_k_r / (1 - (rec_h_srcCor_rs_r + rec_h_srcCor_ro_r));
        irk_incorrect_Fam = rec_h_srcInc_k_r / (1 - (rec_h_srcInc_rs_r + rec_h_srcInc_ro_r));
        
        % Independent Remember--Know Recollect Other rates
        %irk_correct_RcO = rec_h_srcCor_ro_r / (1 - (rec_h_srcCor_rs_r));
        %irk_incorrect_RcO = rec_h_srcInc_ro_r / (1 - (rec_h_srcInc_rs_r));

        %% within-response rates
        %
        % rs
        if isempty(rec_h_srcCor_rs)
          rec_h_srcCor_rs_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcCor_rs_wir = length(rec_h_srcCor_rs) / (length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs));
        end
        if rec_h_srcCor_rs_wir == 1
          rec_h_srcCor_rs_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_srcInc_rs)
          rec_h_srcInc_rs_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcInc_rs_wir = length(rec_h_srcInc_rs) / (length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs));
        end
        if rec_h_srcInc_rs_wir == 1
          rec_h_srcInc_rs_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % ro
        if isempty(rec_h_srcCor_ro)
          rec_h_srcCor_ro_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcCor_ro_wir = length(rec_h_srcCor_ro) / (length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro));
        end
        if rec_h_srcCor_ro_wir == 1
          rec_h_srcCor_ro_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_srcInc_ro)
          rec_h_srcInc_ro_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcInc_ro_wir = length(rec_h_srcInc_ro) / (length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro));
        end
        if rec_h_srcInc_ro_wir == 1
          rec_h_srcInc_ro_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % k
        if isempty(rec_h_srcCor_k)
          rec_h_srcCor_k_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcCor_k_wir = length(rec_h_srcCor_k) / (length(rec_h_srcCor_k) + length(rec_h_srcInc_k));
        end
        if rec_h_srcCor_k_wir == 1
          rec_h_srcCor_k_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_srcInc_k)
          rec_h_srcInc_k_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcInc_k_wir = length(rec_h_srcInc_k) / (length(rec_h_srcCor_k) + length(rec_h_srcInc_k));
        end
        if rec_h_srcInc_k_wir == 1
          rec_h_srcInc_k_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % miss
        if isempty(rec_m_sure)
          rec_m_sure_wir = 1/(2*length(rec_targEv));
        else
          rec_m_sure_wir = length(rec_m_sure) / (length(rec_m_sure) + length(rec_m_maybe));
        end
        if rec_m_sure_wir == 1
          rec_m_sure_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_m_maybe)
          rec_m_maybe_wir = 1/(2*length(rec_targEv));
        else
          rec_m_maybe_wir = length(rec_m_maybe) / (length(rec_m_sure) + length(rec_m_maybe));
        end
        if rec_m_maybe_wir == 1
          rec_m_maybe_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % cr
        if isempty(rec_cr_sure)
          rec_cr_sure_wir = 1/(2*length(rec_lureEv));
        else
          rec_cr_sure_wir = length(rec_cr_sure) / (length(rec_cr_sure) + length(rec_cr_maybe));
        end
        if rec_cr_sure_wir == 1
          rec_cr_sure_wir = 1 - (1/(2*length(rec_lureEv)));
        end
        if isempty(rec_cr_maybe)
          rec_cr_maybe_wir = 1/(2*length(rec_lureEv));
        else
          rec_cr_maybe_wir = length(rec_cr_maybe) / (length(rec_cr_sure) + length(rec_cr_maybe));
        end
        if rec_cr_maybe_wir == 1
          rec_cr_maybe_wir = 1 - (1/(2*length(rec_lureEv)));
        end
        % fa
        if isempty(rec_fa_rs)
          rec_fa_rs_wir = 1/(2*length(rec_lureEv));
        else
          rec_fa_rs_wir = length(rec_fa_rs) / (length(rec_fa_rs) + length(rec_fa_ro) + length(rec_fa_k));
        end
        if rec_fa_rs_wir == 1
          rec_fa_rs_wir = 1 - (1/(2*length(rec_lureEv)));
        end
        if isempty(rec_fa_ro)
          rec_fa_ro_wir = 1/(2*length(rec_lureEv));
        else
          rec_fa_ro_wir = length(rec_fa_ro) / (length(rec_fa_rs) + length(rec_fa_ro) + length(rec_fa_k));
        end
        if rec_fa_ro_wir == 1
          rec_fa_ro_wir = 1 - (1/(2*length(rec_lureEv)));
        end
        if isempty(rec_fa_k)
          rec_fa_k_wir = 1/(2*length(rec_lureEv));
        else
          rec_fa_k_wir = length(rec_fa_k) / (length(rec_fa_rs) + length(rec_fa_ro) + length(rec_fa_k));
        end
        if rec_fa_k_wir == 1
          rec_fa_k_wir = 1 - (1/(2*length(rec_lureEv)));
        end
        
        %% REACTION TIMES
        
        % test items
        rec_h_rt = mean(getStructField(rec_h,'src_rt'));
        rec_m_rt = mean(getStructField(rec_m,'src_rt'));
        rec_cr_rt = mean(getStructField(rec_cr,'src_rt'));
        rec_fa_rt = mean(getStructField(rec_fa,'src_rt'));
        
        src_h_rt = mean(getStructField(src_h,'src_rt'));
        src_m_rt = mean(getStructField(src_m,'src_rt'));
        src_cr_rt = mean(getStructField(src_cr,'src_rt'));
        src_fa_rt = mean(getStructField(src_fa,'src_rt'));
        
        rec_h_srcCor_rt = mean(getStructField(rec_h_srcCor,'src_rt'));
        rec_h_srcInc_rt = mean(getStructField(rec_h_srcInc,'src_rt'));
        
        rec_h_srcCor_rs_rt = mean(getStructField(rec_h_srcCor_rs,'src_rt'));
        rec_h_srcCor_ro_rt = mean(getStructField(rec_h_srcCor_ro,'src_rt'));
        rec_h_srcCor_k_rt = mean(getStructField(rec_h_srcCor_k,'src_rt'));
        rec_h_srcInc_rs_rt = mean(getStructField(rec_h_srcInc_rs,'src_rt'));
        rec_h_srcInc_ro_rt = mean(getStructField(rec_h_srcInc_ro,'src_rt'));
        rec_h_srcInc_k_rt = mean(getStructField(rec_h_srcInc_k,'src_rt'));
        
        rec_m_sure_rt = mean(getStructField(rec_m_sure,'src_rt'));
        rec_m_maybe_rt = mean(getStructField(rec_m_maybe,'src_rt'));
        rec_cr_sure_rt = mean(getStructField(rec_cr_sure,'src_rt'));
        rec_cr_maybe_rt = mean(getStructField(rec_cr_maybe,'src_rt'));
        rec_fa_rs_rt = mean(getStructField(rec_fa_rs,'src_rt'));
        rec_fa_ro_rt = mean(getStructField(rec_fa_ro,'src_rt'));
        rec_fa_k_rt = mean(getStructField(rec_fa_k,'src_rt'));
        
        %% print some info
        fprintf('\n\tRecognition accuracy info:\n');
        fprintf('\tHits: %.4f\n',rec_hr);
        fprintf('\tMisses: %.4f\n',rec_mr);
        fprintf('\tCorrect rejections: %.4f\n',rec_crr);
        fprintf('\tFalse alarms: %.4f\n',rec_far);
        
        fprintf('\n\tSource accuracy info:\n');
        fprintf('\tHits: %.4f\n',src_hr);
        fprintf('\tMisses: %.4f\n',src_mr);
        fprintf('\tCorrect rejections: %.4f\n',src_crr);
        fprintf('\tFalse alarms: %.4f\n\n',src_far);
        
      elseif rejectArt == 1
        %% Counts/rates rejecting artifacts
        
        % Don't currently have all event information if using rejectArt; only have
        % info for recognition hits and correct rejections. Source hits, CR,
        % missed, and FAs are ok because they come from recognition hits. This is
        % because we're not segmenting misses and false alarms in NS.
        
        %% RAW NUMBERS
        
        % TEST ITEMS: get all the targets (old study items presented at test)
        rec_targEv = filterStruct(subSesEv,'rec_isTarg == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
        % get all the lures (new test items)
        rec_lureEv = filterStruct(subSesEv,'rec_isTarg == 0 & ismember(type,varargin{1})',{'LURE_PRES'});
        
        % all source reponses (targets and lures) are conditionalized on
        % getting a recognition hit
        %
        % source targets (for Hs and Ms) and lures (for CRs and FAs) come
        % from the source coding of target items, as explained in
        % sosi_createEvents.m
        src_targEv = filterStruct(rec_targEv,'src_isTarg == 1 & rec_correct == 1');
        src_lureEv = filterStruct(rec_targEv,'src_isTarg == 0 & rec_correct == 1');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TEST-PRESENTATION EVENTS %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % recognition events only
        rec_h = filterStruct(rec_targEv,'rec_correct == 1');
        rec_cr = filterStruct(rec_lureEv,'rec_correct == 1');
        
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
        
        % recognition correct rejection - new items called "new"
        rec_cr_sure = filterStruct(rec_cr,'src_correct == 1 & ismember(rkn_resp,varargin{1})',{'SURE'});
        rec_cr_maybe = filterStruct(rec_cr,'src_correct == 1 & ismember(rkn_resp,varargin{1})',{'MAYBE'});
        
        % recognition numbers collapsed across source accuracy
        rec_h_rs = length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs);
        rec_h_ro = length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro);
        rec_h_k = length(rec_h_srcCor_k) + length(rec_h_srcInc_k);
        
        %% RATES
        
        % accuracy rates
        src_hr = length(src_h)/length(src_targEv);
        src_mr = length(src_m)/length(src_targEv);
        src_crr = length(src_cr)/length(src_lureEv);
        src_far = length(src_fa)/length(src_lureEv);
        
        % correct for rates = 1 or 0
        %
        % src
        if src_hr == 1
          src_hr = 1 - (1/(2*length(src_targEv)));
        elseif src_hr == 0
          src_hr = 1/(2*length(src_targEv));
        end
        if src_mr == 1
          src_mr = 1 - (1/(2*length(src_targEv)));
        elseif src_mr == 0
          src_mr = 1/(2*length(src_targEv));
        end
        if src_crr == 1
          src_crr = 1 - (1/(2*length(src_lureEv)));
        elseif src_crr == 0
          src_crr = 1/(2*length(src_lureEv));
        end
        if src_far == 1
          src_far = 1 - (1/(2*length(src_lureEv)));
        elseif src_far == 0
          src_far = 1/(2*length(src_lureEv));
        end
        
        % d primes
        source_dp = norminv(src_hr,0,1) - norminv(src_far,0,1);
        
        % response bias (criterion) from Macmillan & Creelman p. 29
        source_c = -0.5 * (norminv(src_hr,0,1) + norminv(src_far,0,1));
        
        % From Mecklinger et al. (2007) and Corwin (1994)
        %
        % discrimination index
        %Pr = rec_hr - rec_far;
        % response bias index
        %Br = rec_far / (1 - Pr);
        
        % accuracy rates for recognition hits
        rec_h_srcCor_r = length(rec_h_srcCor)/length(rec_h);
        rec_h_srcInc_r = length(rec_h_srcInc)/length(rec_h);
        
        % rates by confidence
        rec_h_srcCor_rs_r = length(rec_h_srcCor_rs) / length(rec_targEv);
        rec_h_srcCor_ro_r = length(rec_h_srcCor_ro) / length(rec_targEv);
        rec_h_srcCor_k_r = length(rec_h_srcCor_k) / length(rec_targEv);
        rec_h_srcInc_rs_r = length(rec_h_srcInc_rs) / length(rec_targEv);
        rec_h_srcInc_ro_r = length(rec_h_srcInc_ro) / length(rec_targEv);
        rec_h_srcInc_k_r = length(rec_h_srcInc_k) / length(rec_targEv);
        rec_cr_sure_r = length(rec_cr_sure) / length(rec_lureEv);
        rec_cr_maybe_r = length(rec_cr_maybe) / length(rec_lureEv);
        
        %% fix for when rates are 1 or 0 (Macmillan and Creelman, 1991)
        %
        % srcCor
        if rec_h_srcCor_rs_r == 1
          rec_h_srcCor_rs_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcCor_rs_r == 0
          rec_h_srcCor_rs_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcCor_ro_r == 1
          rec_h_srcCor_ro_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcCor_ro_r == 0
          rec_h_srcCor_ro_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcCor_k_r == 1
          rec_h_srcCor_k_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcCor_k_r == 0
          rec_h_srcCor_k_r = 1/(2*length(rec_targEv));
        end
        % srcInc
        if rec_h_srcInc_rs_r == 1
          rec_h_srcInc_rs_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcInc_rs_r == 0
          rec_h_srcInc_rs_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcInc_ro_r == 1
          rec_h_srcInc_ro_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcInc_ro_r == 0
          rec_h_srcInc_ro_r = 1/(2*length(rec_targEv));
        end
        if rec_h_srcInc_k_r == 1
          rec_h_srcInc_k_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_srcInc_k_r == 0
          rec_h_srcInc_k_r = 1/(2*length(rec_targEv));
        end
        % cr
        if rec_cr_sure_r == 1
          rec_cr_sure_r = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_cr_sure_r == 0
          rec_cr_sure_r = 1/(2*length(rec_lureEv));
        end
        if rec_cr_maybe_r == 1
          rec_cr_maybe_r = 1 - (1/(2*length(rec_lureEv)));
        elseif rec_cr_maybe_r == 0
          rec_cr_maybe_r = 1/(2*length(rec_lureEv));
        end
        
        %% Independent Remember--Know Familiarity rate
        irk_Fam = rec_h_k_r / (1 - (rec_h_rs_r + rec_h_ro_r));
        
        %% Independent Remember--Know Familiarity rates - I think this is wrong
        irk_correct_Fam = rec_h_srcCor_k_r / (1 - (rec_h_srcCor_rs_r + rec_h_srcCor_ro_r));
        irk_incorrect_Fam = rec_h_srcInc_k_r / (1 - (rec_h_srcInc_rs_r + rec_h_srcInc_ro_r));
        
        % Independent Remember--Know Recollect Other rates
        %irk_correct_RcO = rec_h_srcCor_ro_r / (1 - (rec_h_srcCor_rs_r));
        %irk_incorrect_RcO = rec_h_srcInc_ro_r / (1 - (rec_h_srcInc_rs_r));

        %% within-response rates
        %
        % rs
        if isempty(rec_h_srcCor_rs)
          rec_h_srcCor_rs_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcCor_rs_wir = length(rec_h_srcCor_rs) / (length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs));
        end
        if rec_h_srcCor_rs_wir == 1
          rec_h_srcCor_rs_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_srcInc_rs)
          rec_h_srcInc_rs_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcInc_rs_wir = length(rec_h_srcInc_rs) / (length(rec_h_srcCor_rs) + length(rec_h_srcInc_rs));
        end
        if rec_h_srcInc_rs_wir == 1
          rec_h_srcInc_rs_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % ro
        if isempty(rec_h_srcCor_ro)
          rec_h_srcCor_ro_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcCor_ro_wir = length(rec_h_srcCor_ro) / (length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro));
        end
        if rec_h_srcCor_ro_wir == 1
          rec_h_srcCor_ro_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_srcInc_ro)
          rec_h_srcInc_ro_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcInc_ro_wir = length(rec_h_srcInc_ro) / (length(rec_h_srcCor_ro) + length(rec_h_srcInc_ro));
        end
        if rec_h_srcInc_ro_wir == 1
          rec_h_srcInc_ro_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % k
        if isempty(rec_h_srcCor_k)
          rec_h_srcCor_k_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcCor_k_wir = length(rec_h_srcCor_k) / (length(rec_h_srcCor_k) + length(rec_h_srcInc_k));
        end
        if rec_h_srcCor_k_wir == 1
          rec_h_srcCor_k_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_srcInc_k)
          rec_h_srcInc_k_wir = 1/(2*length(rec_targEv));
        else
          rec_h_srcInc_k_wir = length(rec_h_srcInc_k) / (length(rec_h_srcCor_k) + length(rec_h_srcInc_k));
        end
        if rec_h_srcInc_k_wir == 1
          rec_h_srcInc_k_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % cr
        if isempty(rec_cr_sure)
          rec_cr_sure_wir = 1/(2*length(rec_lureEv));
        else
          rec_cr_sure_wir = length(rec_cr_sure) / (length(rec_cr_sure) + length(rec_cr_maybe));
        end
        if rec_cr_sure_wir == 1
          rec_cr_sure_wir = 1 - (1/(2*length(rec_lureEv)));
        end
        if isempty(rec_cr_maybe)
          rec_cr_maybe_wir = 1/(2*length(rec_lureEv));
        else
          rec_cr_maybe_wir = length(rec_cr_maybe) / (length(rec_cr_sure) + length(rec_cr_maybe));
        end
        if rec_cr_maybe_wir == 1
          rec_cr_maybe_wir = 1 - (1/(2*length(rec_lureEv)));
        end
        
        %% REACTION TIMES
        
        % test items
        rec_h_rt = mean(getStructField(rec_h,'src_rt'));
        rec_cr_rt = mean(getStructField(rec_cr,'src_rt'));
        
        src_h_rt = mean(getStructField(src_h,'src_rt'));
        src_m_rt = mean(getStructField(src_m,'src_rt'));
        src_cr_rt = mean(getStructField(src_cr,'src_rt'));
        src_fa_rt = mean(getStructField(src_fa,'src_rt'));
        
        rec_h_srcCor_rt = mean(getStructField(rec_h_srcCor,'src_rt'));
        rec_h_srcInc_rt = mean(getStructField(rec_h_srcInc,'src_rt'));
        
        rec_h_srcCor_rs_rt = mean(getStructField(rec_h_srcCor_rs,'src_rt'));
        rec_h_srcCor_ro_rt = mean(getStructField(rec_h_srcCor_ro,'src_rt'));
        rec_h_srcCor_k_rt = mean(getStructField(rec_h_srcCor_k,'src_rt'));
        rec_h_srcInc_rs_rt = mean(getStructField(rec_h_srcInc_rs,'src_rt'));
        rec_h_srcInc_ro_rt = mean(getStructField(rec_h_srcInc_ro,'src_rt'));
        rec_h_srcInc_k_rt = mean(getStructField(rec_h_srcInc_k,'src_rt'));
        
        rec_cr_sure_rt = mean(getStructField(rec_cr_sure,'src_rt'));
        rec_cr_maybe_rt = mean(getStructField(rec_cr_maybe,'src_rt'));
        
        %% print some info
        %fprintf('\n\tRecognition accuracy info:\n');
        %fprintf('\tHits: %.4f\n',rec_hr);
        %fprintf('\tMisses: %.4f\n',rec_mr);
        %fprintf('\tCorrect rejections: %.4f\n',rec_crr);
        %fprintf('\tFalse alarms: %.4f\n',rec_far);
        
        fprintf('\n\tSource accuracy info:\n');
        fprintf('\tHits: %.4f\n',src_hr);
        fprintf('\tMisses: %.4f\n',src_mr);
        fprintf('\tCorrect rejections: %.4f\n',src_crr);
        fprintf('\tFalse alarms: %.4f\n\n',src_far);
      end
      
      % ===================================================================
      
      if rejectArt == 1
        % count the artifacts
        
        % get all the targets (old study items presented at test)
        rec_targEv_art = filterStruct(subSesEv_all,'nsArt == 1 & rec_isTarg == 1 & rec_correct == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
        % get all the lures (new test items)
        rec_lureEv_art = filterStruct(subSesEv_all,'nsArt == 1 & rec_isTarg == 0 & rec_correct == 1 & ismember(type,varargin{1})',{'LURE_PRES'});
        
        % get all the targets (old study items presented at test)
        rec_targEv_all = filterStruct(subSesEv_all,'rec_isTarg == 1 & rec_correct == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
        % get all the lures (new test items)
        rec_lureEv_all = filterStruct(subSesEv_all,'rec_isTarg == 0 & rec_correct == 1 & ismember(type,varargin{1})',{'LURE_PRES'});
        
        numArt = length(rec_targEv_art) + length(rec_lureEv_art);
        
        numEv = length(rec_targEv_all) + length(rec_lureEv_all);
      end
      
      if rejectArt == 0
        tableData = [...
          length(rec_targEv),length(rec_lureEv),length(src_targEv),length(src_lureEv),... % raw numbers of targs, lures
          item_dp,source_dp,item_c,source_c,... % accuracy
          length(rec_h),length(rec_m),length(rec_cr),length(rec_fa),length(src_h),length(src_m),length(src_cr),length(src_fa),... % raw numbers for accuracy
          rec_hr,rec_mr,rec_crr,rec_far,src_hr,src_mr,src_crr,src_far,... % accuracy rates
          length(rec_h_srcCor),length(rec_h_srcInc),... % raw numbers for recognition hits
          rec_h_srcCor_r,rec_h_srcInc_r,... % accuracy rates for recognition hits
          rec_h_rs,rec_h_ro,rec_h_k,... % raw numbers for recognition hit response collapsed across source accuracy
          length(rec_h_srcCor_rs),length(rec_h_srcCor_ro),length(rec_h_srcCor_k),length(rec_h_srcInc_rs),length(rec_h_srcInc_ro),length(rec_h_srcInc_k),length(rec_m_sure),length(rec_m_maybe),length(rec_cr_sure),length(rec_cr_maybe),length(rec_fa_rs),length(rec_fa_ro),length(rec_fa_k),... % raw numbers by confidence
          rec_h_srcCor_rs_r,rec_h_srcCor_ro_r,rec_h_srcCor_k_r,rec_h_srcInc_rs_r,rec_h_srcInc_ro_r,rec_h_srcInc_k_r,rec_m_sure_r,rec_m_maybe_r,rec_cr_sure_r,rec_cr_maybe_r,rec_fa_rs_r,rec_fa_ro_r,rec_fa_k_r,... % rates by confidence
          irk_Fam,irk_correct_Fam,irk_incorrect_Fam,... % Independent Remember--Know Familiarity rates
          rec_h_srcCor_rs_wir,rec_h_srcInc_rs_wir,rec_h_srcCor_ro_wir,rec_h_srcInc_ro_wir,rec_h_srcCor_k_wir,rec_h_srcInc_k_wir,rec_m_sure_wir,rec_m_maybe_wir,rec_cr_sure_wir,rec_cr_maybe_wir,rec_fa_rs_wir,rec_fa_ro_wir,rec_fa_k_wir,... % within-response rates
          rec_h_rt,rec_m_rt,rec_cr_rt,rec_fa_rt,... % rec reaction times
          src_h_rt,src_m_rt,src_cr_rt,src_fa_rt,... % src reaction times
          rec_h_srcCor_rt,rec_h_srcInc_rt,... % recognition hits reaction times
          rec_h_srcCor_rs_rt,rec_h_srcCor_ro_rt,rec_h_srcCor_k_rt,rec_h_srcInc_rs_rt,rec_h_srcInc_ro_rt,rec_h_srcInc_k_rt,rec_m_sure_rt,rec_m_maybe_rt,rec_cr_sure_rt,rec_cr_maybe_rt,rec_fa_rs_rt,rec_fa_ro_rt,rec_fa_k_rt,... % rt by confidence
          Pr,Br,... % more accuracy
          ];
      elseif rejectArt == 1
        tableData = [...
          length(rec_targEv),length(rec_lureEv),length(src_targEv),length(src_lureEv),... % raw numbers of targs, lures
          source_dp,source_c,... % accuracy
          length(rec_h),length(rec_cr),length(src_h),length(src_m),length(src_cr),length(src_fa),... % raw numbers for accuracy
          src_hr,src_mr,src_crr,src_far,... % accuracy rates
          length(rec_h_srcCor),length(rec_h_srcInc),... % raw numbers for recognition hits
          rec_h_srcCor_r,rec_h_srcInc_r,... % accuracy rates for recognition hits
          rec_h_rs,rec_h_ro,rec_h_k,... % raw numbers for recognition hit response collapsed across source accuracy
          length(rec_h_srcCor_rs),length(rec_h_srcCor_ro),length(rec_h_srcCor_k),length(rec_h_srcInc_rs),length(rec_h_srcInc_ro),length(rec_h_srcInc_k),length(rec_cr_sure),length(rec_cr_maybe),... % raw numbers by confidence
          rec_h_srcCor_rs_r,rec_h_srcCor_ro_r,rec_h_srcCor_k_r,rec_h_srcInc_rs_r,rec_h_srcInc_ro_r,rec_h_srcInc_k_r,rec_cr_sure_r,rec_cr_maybe_r,... % rates by confidence
          irk_Fam,irk_correct_Fam,irk_incorrect_Fam,... % Independent Remember--Know Familiarity rates
          rec_h_srcCor_rs_wir,rec_h_srcInc_rs_wir,rec_h_srcCor_ro_wir,rec_h_srcInc_ro_wir,rec_h_srcCor_k_wir,rec_h_srcInc_k_wir,rec_cr_sure_wir,rec_cr_maybe_wir,... % within-response rates
          rec_h_rt,rec_cr_rt,... % rec reaction times
          src_h_rt,src_m_rt,src_cr_rt,src_fa_rt,... % src reaction times
          rec_h_srcCor_rt,rec_h_srcInc_rt,... % recognition hits reaction times
          rec_h_srcCor_rs_rt,rec_h_srcCor_ro_rt,rec_h_srcCor_k_rt,rec_h_srcInc_rs_rt,rec_h_srcInc_ro_rt,rec_h_srcInc_k_rt,rec_cr_sure_rt,rec_cr_maybe_rt,... % rt by confidence
          ];
      end
    end % ses
    
    if saveFiles == 1
      % print sub and ses
      fprintf(outfile,'%s,%d',subjects{sub},ses-1);
      % format for tableData and print
      tableDataStr = repmat(',%.4f',1,length(tableData));
      fprintf(outfile,tableDataStr,tableData);
      %for tabData = 1:length(tableData)
      %  fprintf(outfile,'%d\t',tableData(tabData));
      %end
      if rejectArt == 1
        artStr = sprintf(',%d,%d,%.4f',numArt,numEv,(numArt / numEv));
        fprintf(outfile,'%s\n',artStr);
      else
        fprintf(outfile,'\n');
      end
    end
    
  end % sub
  
  if saveFiles == 1
    fprintf('Saving %s...',eventsTable);
    fclose(outfile);
    fprintf('Done.\n');
  end
  
end % subsets
