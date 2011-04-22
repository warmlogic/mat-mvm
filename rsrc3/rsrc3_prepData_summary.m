function rsrc3_prepData_summary(session,rejectArt,saveFiles,dataroot)
% function rsrc3_prepData_summary(session,rejectArt,saveFiles,dataroot)
%
% Save column-formatted data as a csv file
%
% Looks for each subject's events.mat and saves a summary file in dataroot
%

expName = 'RSRC3';

if nargin < 4
  serverDir = fullfile('/Volumes/curranlab/Data',expName,'eeg/behavioral');
  serverLocalDir = fullfile('/Volumes/RAID/curranlab/Data',expName,'eeg/behavioral');
  if exist(serverDir,'dir')
    dataroot = serverDir;
  elseif exist(serverLocalDir,'dir')
    dataroot = serverLocalDir;
  else
    uroot = getenv('HOME');
    dataroot = fullfile(uroot,'data',expName,'eeg');
  end
  if nargin < 3
    saveFiles = 1;
    if nargin < 2
      rejectArt = 1;
    end
  end
end

subjects = {
    'RSRC3001'...
    'RSRC3002'...
    'RSRC3003'...
    'RSRC3004'...
    'RSRC3005'...
    'RSRC3006'...
    'RSRC3007'...
    'RSRC3008'...
    'RSRC3009'...
    'RSRC3010'...
    'RSRC3011'...
    'RSRC3012'...
    'RSRC3013'...
%     'RSRC3014'...
%     'RSRC3015'...
%     'RSRC3016'...
%     'RSRC3017'...
%     'RSRC3018'...
%     'RSRC3019'...
%     'RSRC3020'...
%     'RSRC3021'...
%     'RSRC3022'...
%     'RSRC3023'...
%     'RSRC3024'...
%     'RSRC3025'...
%     'RSRC3026'...
%     'RSRC3027'...
%     'RSRC3028'...
%     'RSRC3029'...
%     'RSRC3030'...
    };
% RSRC3006 - not sure if she understood the instructions

if session == 1
  %sessions = {'session_0','session_1'};
  sessions = {'session_1'};
elseif session == 2
  %sessions = {'session_2','session_3'};
  sessions = {'session_3'};
end

subsets = {''};


%% set up the header
if rejectArt == 0
  tableHeader = {'sub','ses',... % subject info
    'R targ','R lure','S targ','S lure',... % raw numbers of targs, lures
    'Item d''','Source d''','Item da','Item zROC slope','Item RespBias (c)','Source RespBias (c)','IRK Fam',... % accuracy
    'RH','RM','RCR','RFA','SH','SM','SCR','SFA',... % raw numbers for accuracy
    'RHR','RMR','RCRR','RFAR','SHR','SMR','SCRR','SFAR',... % accuracy rates
    'RH-SCor','RH-SInc','RM-SCor','RM-SInc',... % raw numbers for recognition hits/misses
    'RH-SCor R','RH-SInc R','RM-SCor R','RM-SInc R',... % rates for recognition hits/misses
    'RH-RS','RH-RO','RH-DF','RH-MF','RM-MU','RM-DU','RFA-RS','RFA-RO','RFA-DF','RFA-MF','RCR-MU','RCR-DU',... % raw numbers by confidence collapsed across source accuracy
    'RH-RS R','RH-RO R','RH-DF R','RH-MF R','RM-MU R','RM-DU R','RFA-RS R','RFA-RO R','RFA-DF R','RFA-MF R','RCR-MU R','RCR-DU R',... % rates by confidence collapsed across source accuracy
    'RH-RS SCor','RH-RO SCor','RH-DF SCor','RH-MF SCor','RM-MU SCor','RM-DU SCor',... % raw numbers for recognition hits/misses by confidence for source correct
    'RH-RS SInc','RH-RO SInc','RH-DF SInc','RH-MF SInc','RM-MU SInc','RM-DU SInc',... % raw numbers for recognition hits/misses by confidence for source incorrect
    'RH-RS SCor R','RH-RO SCor R','RH-DF SCor R','RH-MF SCor R','RM-MU SCor R','RM-DU SCor R',... % rates for recognition hits/misses by confidence for source correct
    'RH-RS SInc R','RH-RO SInc R','RH-DF SInc R','RH-MF SInc R','RM-MU SInc R','RM-DU SInc R',... % rates for recognition hits/misses by confidence for source incorrect
    'RH-RS SCor WIR','RH-RS SInc WIR','RH-RO SCor WIR','RH-RO SInc WIR','RH-DF SCor WIR','RH-DF SInc WIR','RH-MF SCor WIR','RH-MF SInc WIR','RM-MU SCor WIR','RM-MU SInc WIR','RM-DU SCor WIR','RM-DU SInc WIR',... % rates for recognition hits/misses by confidence for source correct
    'RH rt','RM rt','RCR rt','RFA rt',... % rec reaction times
    'SH rt','SM rt','SCR rt','SFA rt',... % src reaction times
    'RH-SCor rt','RH-SInc rt',... % recognition hits reaction times
    'Rec Discrimination (Pr)','Rec Response Bias (Br)',... % more accuracy
    };
elseif rejectArt == 1
  tableHeader = {'sub','ses',... % subject info
    'R targ','R lure','S targ','S lure',... % raw numbers of targs, lures
    'Source d''','Source da','Source zROC slope','Source RespBias (c)',... % accuracy
    'RH','RCR','SH','SM','SCR','SFA',... % raw numbers for accuracy
    'SHR','SMR','SCRR','SFAR',... % accuracy rates
    'RH-SCor','RH-SInc',... % raw numbers for recognition hits/misses
    'RH-SCor R','RH-SInc R',... % rates for recognition hits/misses
    'RH-RS','RH-RO','RH-DF','RH-MF','RCR-MU','RCR-DU',... % raw numbers by confidence collapsed across source accuracy
    'RH-RS R','RH-RO R','RH-DF R','RH-MF R','RCR-MU R','RCR-DU R',... % rates by confidence collapsed across source accuracy
    'RH-RS SCor','RH-RO SCor','RH-DF SCor','RH-MF SCor',... % raw numbers for recognition hits by confidence for source correct
    'RH-RS SInc','RH-RO SInc','RH-DF SInc','RH-MF SInc',... % raw numbers for recognition hits by confidence for source incorrect
    'RH-RS SCor R','RH-RO SCor R','RH-DF SCor R','RH-MF SCor R','RM-MU SCor R','RM-DU SCor R',... % rates for recognition hits/misses by confidence for source correct
    'RH-RS SInc R','RH-RO SInc R','RH-DF SInc R','RH-MF SInc R','RM-MU SInc R','RM-DU SInc R',... % rates for recognition hits/misses by confidence for source incorrect
    'RH rt','RCR rt',... % rec reaction times
    'SH rt','SM rt','SCR rt','SFA rt',... % src reaction times
    'RH-SCor rt','RH-SInc rt',... % recognition hits reaction times
    };
end

if rejectArt == 1
  artSegStr = sprintf('Segments marked bad for eye artifacts');
  tableHeader{end+1} = artSegStr;
end

% confidence numbers are [1 2 3 4 5 6]: RS, RO, DF, MF, MU, DU (old to new)
confRange_rec = [1 2 3 4 5 6];

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
      fprintf(outfile,'%s\n',expName);
      fprintf(outfile,'%s',tableHeader{1});
      for headCol = 2:length(tableHeader)
        fprintf(outfile,',%s',tableHeader{headCol});
      end
      fprintf(outfile,'\n');
    else
      fprintf('%s already exists, not overwriting!\n',eventsTable);
      continue
    end
  else
    fprintf('saveFiles = 0, not writing out any files\n')
  end
  
  %% for each subject
  for sub = 1:length(subjects)
    fprintf('Getting data for %s...',subjects{sub});
    
%     if str2double(subjects{sub}(end)) > 0 && str2double(subjects{sub}(end)) <= 5
%       % confidence numbers are [1 2 3 4 5 6]: RS, RO, DF, MF, MU, DU (old to new)
%       confRange_rec = [1 2 3 4 5 6];
%     elseif str2double(subjects{sub}(end)) == 0 || str2double(subjects{sub}(end)) > 5
%       % confidence numbers are [6 5 4 3 2 1]: DU, MU, MF, DF, RO, RS (new to old)
%       confRange_rec = [6 5 4 3 2 1];
%     end

    %% for each session
    for ses = 1:length(sessions)
      fprintf('%s...\n',sessions{ses});
      % set the subject events directory
      eventsDir_sub = fullfile(dataroot,subjects{sub},sessions{ses},'events');
      events = loadEvents(fullfile(eventsDir_sub,'test_events.mat'));
      
      if rejectArt == 1
        % remove artifacts
        subSesEv_all = events;
        subSesEv = filterStruct(subSesEv_all,'nsArt == 0');
      else
        subSesEv = events;
      end
      
%       if strcmp(subsets{s},'_task')
%         % get only the size study events
%         fprintf('Task\n');
%         subSesEv = filterStruct(events,'ismember(testType,varargin{1})',{'task'});
%       elseif strcmp(subsets{s},'_side')
%         % get only the life study events
%         fprintf('Side\n');
%         subSesEv = filterStruct(events,'ismember(testType,varargin{1})',{'side'});
%       else
%         fprintf('All\n');
%         subSesEv = events;
%       end
      
      if rejectArt == 0
        %% RAW NUMBERS
        
        % get all the targets (old study items presented at test)
        rec_targEv = filterStruct(subSesEv,'rec_isTarg == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
        % get all the lures (new test items)
        rec_lureEv = filterStruct(subSesEv,'rec_isTarg == 0 & ismember(type,varargin{1})',{'LURE_PRES'});
        
        % all source reponses (targets and lures) are conditionalized on
        % getting a recognition hit
        %
        % source targets (for Hs and Ms) and lures (for CRs and FAs) come
        % from the source coding of target items, as explained in
        % rsrc_createEvents.m
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
        % explained in rsrc_createEvents.m
        src_cr = filterStruct(src_lureEv,'src_correct == 1');
        src_fa = filterStruct(src_lureEv,'src_correct == 0');
        
        % recognition by confidence level
        rec_h_rs = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(1));
        rec_h_ro = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(2));
        rec_h_df = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(3));
        rec_h_mf = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(4));
        rec_m_mu = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(5));
        rec_m_du = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(6));
        rec_fa_rs = filterStruct(rec_lureEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(1));
        rec_fa_ro = filterStruct(rec_lureEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(2));
        rec_fa_df = filterStruct(rec_lureEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(3));
        rec_fa_mf = filterStruct(rec_lureEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(4));
        rec_cr_mu = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(5));
        rec_cr_du = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(6));
        
        % recognition hit and source correct/incorrect (correct = hit/CR,
        % incorrect= miss/FA)
        rec_h_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1');
        rec_h_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0');
        rec_m_srcCor = filterStruct(rec_targEv,'rec_correct == 0 & src_correct == 1');
        rec_m_srcInc = filterStruct(rec_targEv,'rec_correct == 0 & src_correct == 0');
        
        % recognition hit/miss, source correct
        rec_h_rs_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1',confRange_rec(1));
        rec_h_ro_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1',confRange_rec(2));
        rec_h_df_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1',confRange_rec(3));
        rec_h_mf_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1',confRange_rec(4));
        rec_m_mu_srcCor = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1',confRange_rec(5));
        rec_m_du_srcCor = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1',confRange_rec(6));
        % recognition hit/miss, source incorrect
        rec_h_rs_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0',confRange_rec(1));
        rec_h_ro_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0',confRange_rec(2));
        rec_h_df_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0',confRange_rec(3));
        rec_h_mf_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0',confRange_rec(4));
        rec_m_mu_srcInc = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0',confRange_rec(5));
        rec_m_du_srcInc = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0',confRange_rec(6));
        
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
        
        % rates for recognition hits/misses
        rec_h_srcCor_r = length(rec_h_srcCor) / length(rec_h);
        rec_h_srcInc_r = length(rec_h_srcInc) / length(rec_h);
        rec_m_srcCor_r = length(rec_m_srcCor) / length(rec_m);
        rec_m_srcInc_r = length(rec_m_srcInc) / length(rec_m);
        
        %% rates by confidence
        rec_h_rs_r = length(rec_h_rs) / length(rec_targEv);
        rec_h_ro_r = length(rec_h_ro) / length(rec_targEv);
        rec_h_df_r = length(rec_h_df) / length(rec_targEv);
        rec_h_mf_r = length(rec_h_mf) / length(rec_targEv);
        rec_m_mu_r = length(rec_m_mu) / length(rec_targEv);
        rec_m_du_r = length(rec_m_du) / length(rec_targEv);
        rec_fa_rs_r = length(rec_fa_rs) / length(rec_lureEv);
        rec_fa_ro_r = length(rec_fa_ro) / length(rec_lureEv);
        rec_fa_df_r = length(rec_fa_df) / length(rec_lureEv);
        rec_fa_mf_r = length(rec_fa_mf) / length(rec_lureEv);
        rec_cr_mu_r = length(rec_cr_mu) / length(rec_lureEv);
        rec_cr_du_r = length(rec_cr_du) / length(rec_lureEv);
        
        %% fix for when rates are 1 or 0 (Macmillan and Creelman, 1991)
        %
        % hit/miss
        if rec_h_rs_r == 1
          rec_h_rs_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_rs_r == 0
          rec_h_rs_r = 1/(2*length(rec_targEv));
        end
        if rec_h_ro_r == 1
          rec_h_ro_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_ro_r == 0
          rec_h_ro_r = 1/(2*length(rec_targEv));
        end
        if rec_h_df_r == 1
          rec_h_df_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_df_r == 0
          rec_h_df_r = 1/(2*length(rec_targEv));
        end
        if rec_h_mf_r == 1
          rec_h_mf_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_mf_r == 0
          rec_h_mf_r = 1/(2*length(rec_targEv));
        end
        if rec_m_mu_r == 1
          rec_m_mu_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_mu_r == 0
          rec_m_mu_r = 1/(2*length(rec_targEv));
        end
        if rec_m_du_r == 1
          rec_m_du_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_du_r == 0
          rec_m_du_r = 1/(2*length(rec_targEv));
        end
        % cr/fa
        if rec_fa_rs_r == 1
          rec_fa_rs_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_fa_rs_r == 0
          rec_fa_rs_r = 1/(2*length(rec_targEv));
        end
        if rec_fa_ro_r == 1
          rec_fa_ro_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_fa_ro_r == 0
          rec_fa_ro_r = 1/(2*length(rec_targEv));
        end
        if rec_fa_df_r == 1
          rec_fa_df_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_fa_df_r == 0
          rec_fa_df_r = 1/(2*length(rec_targEv));
        end
        if rec_fa_mf_r == 1
          rec_fa_mf_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_fa_mf_r == 0
          rec_fa_mf_r = 1/(2*length(rec_targEv));
        end
        if rec_cr_mu_r == 1
          rec_cr_mu_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_cr_mu_r == 0
          rec_cr_mu_r = 1/(2*length(rec_targEv));
        end
        if rec_cr_du_r == 1
          rec_cr_du_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_cr_du_r == 0
          rec_cr_du_r = 1/(2*length(rec_targEv));
        end
        
        %% Independent Remember--Know Familiarity rates
        irk_fam = (rec_h_df_r + rec_h_mf_r) / (1 - (rec_h_rs_r + rec_h_ro_r));
        
        %% RH-SCor: rates by confidence
        rec_h_rs_srcCor_r = length(rec_h_rs_srcCor) / length(rec_targEv);
        rec_h_ro_srcCor_r = length(rec_h_ro_srcCor) / length(rec_targEv);
        rec_h_df_srcCor_r = length(rec_h_df_srcCor) / length(rec_targEv);
        rec_h_mf_srcCor_r = length(rec_h_mf_srcCor) / length(rec_targEv);
        rec_m_mu_srcCor_r = length(rec_m_mu_srcCor) / length(rec_targEv);
        rec_m_du_srcCor_r = length(rec_m_du_srcCor) / length(rec_targEv);
        % RH-SInc: rates by confidence
        rec_h_rs_srcInc_r = length(rec_h_rs_srcInc) / length(rec_targEv);
        rec_h_ro_srcInc_r = length(rec_h_ro_srcInc) / length(rec_targEv);
        rec_h_df_srcInc_r = length(rec_h_df_srcInc) / length(rec_targEv);
        rec_h_mf_srcInc_r = length(rec_h_mf_srcInc) / length(rec_targEv);
        rec_m_mu_srcInc_r = length(rec_m_mu_srcInc) / length(rec_targEv);
        rec_m_du_srcInc_r = length(rec_m_du_srcInc) / length(rec_targEv);
        
        %% fix for when rates are 1 or 0 (Macmillan and Creelman, 1991)
        %
        % hit/miss srcCor
        if rec_h_rs_srcCor_r == 1
          rec_h_rs_srcCor_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_rs_srcCor_r == 0
          rec_h_rs_srcCor_r = 1/(2*length(rec_targEv));
        end
        if rec_h_ro_srcCor_r == 1
          rec_h_ro_srcCor_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_ro_srcCor_r == 0
          rec_h_ro_srcCor_r = 1/(2*length(rec_targEv));
        end
        if rec_h_df_srcCor_r == 1
          rec_h_df_srcCor_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_df_srcCor_r == 0
          rec_h_df_srcCor_r = 1/(2*length(rec_targEv));
        end
        if rec_h_mf_srcCor_r == 1
          rec_h_mf_srcCor_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_mf_srcCor_r == 0
          rec_h_mf_srcCor_r = 1/(2*length(rec_targEv));
        end
        if rec_m_mu_srcCor_r == 1
          rec_m_mu_srcCor_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_mu_srcCor_r == 0
          rec_m_mu_srcCor_r = 1/(2*length(rec_targEv));
        end
        if rec_m_du_srcCor_r == 1
          rec_m_du_srcCor_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_du_srcCor_r == 0
          rec_m_du_srcCor_r = 1/(2*length(rec_targEv));
        end
        % hit/miss srcInc
        if rec_h_rs_srcInc_r == 1
          rec_h_rs_srcInc_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_rs_srcInc_r == 0
          rec_h_rs_srcInc_r = 1/(2*length(rec_targEv));
        end
        if rec_h_ro_srcInc_r == 1
          rec_h_ro_srcInc_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_ro_srcInc_r == 0
          rec_h_ro_srcInc_r = 1/(2*length(rec_targEv));
        end
        if rec_h_df_srcInc_r == 1
          rec_h_df_srcInc_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_df_srcInc_r == 0
          rec_h_df_srcInc_r = 1/(2*length(rec_targEv));
        end
        if rec_h_mf_srcInc_r == 1
          rec_h_mf_srcInc_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_h_mf_srcInc_r == 0
          rec_h_mf_srcInc_r = 1/(2*length(rec_targEv));
        end
        if rec_m_mu_srcInc_r == 1
          rec_m_mu_srcInc_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_mu_srcInc_r == 0
          rec_m_mu_srcInc_r = 1/(2*length(rec_targEv));
        end
        if rec_m_du_srcInc_r == 1
          rec_m_du_srcInc_r = 1 - (1/(2*length(rec_targEv)));
        elseif rec_m_du_srcInc_r == 0
          rec_m_du_srcInc_r = 1/(2*length(rec_targEv));
        end
        
        %% within-response rates
        %
        % hit - rs
        if isempty(rec_h_rs_srcCor)
          rec_h_rs_srcCor_wir = 1/(2*length(rec_targEv));
        else
          rec_h_rs_srcCor_wir = length(rec_h_rs_srcCor) / (length(rec_h_rs_srcCor) + length(rec_h_rs_srcInc));
        end
        if rec_h_rs_srcCor_wir == 1
          rec_h_rs_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_rs_srcInc)
          rec_h_rs_srcInc_wir = 1/(2*length(rec_targEv));
        else
          rec_h_rs_srcInc_wir = length(rec_h_rs_srcInc) / (length(rec_h_rs_srcCor) + length(rec_h_rs_srcInc));
        end
        if rec_h_rs_srcInc_wir == 1
          rec_h_rs_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % hit - ro
        if isempty(rec_h_ro_srcCor)
          rec_h_ro_srcCor_wir = 1/(2*length(rec_targEv));
        else
          rec_h_ro_srcCor_wir = length(rec_h_ro_srcCor) / (length(rec_h_ro_srcCor) + length(rec_h_ro_srcInc));
        end
        if rec_h_ro_srcCor_wir == 1
          rec_h_ro_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_ro_srcInc)
          rec_h_ro_srcInc_wir = 1/(2*length(rec_targEv));
        else
          rec_h_ro_srcInc_wir = length(rec_h_ro_srcInc) / (length(rec_h_ro_srcCor) + length(rec_h_ro_srcInc));
        end
        if rec_h_ro_srcInc_wir == 1
          rec_h_ro_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % hit - df
        if isempty(rec_h_df_srcCor)
          rec_h_df_srcCor_wir = 1/(2*length(rec_targEv));
        else
          rec_h_df_srcCor_wir = length(rec_h_df_srcCor) / (length(rec_h_df_srcCor) + length(rec_h_df_srcInc));
        end
        if rec_h_df_srcCor_wir == 1
          rec_h_df_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_df_srcInc)
          rec_h_df_srcInc_wir = 1/(2*length(rec_targEv));
        else
          rec_h_df_srcInc_wir = length(rec_h_df_srcInc) / (length(rec_h_df_srcCor) + length(rec_h_df_srcInc));
        end
        if rec_h_df_srcInc_wir == 1
          rec_h_df_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % hit - mf
        if isempty(rec_h_mf_srcCor)
          rec_h_mf_srcCor_wir = 1/(2*length(rec_targEv));
        else
          rec_h_mf_srcCor_wir = length(rec_h_mf_srcCor) / (length(rec_h_mf_srcCor) + length(rec_h_mf_srcInc));
        end
        if rec_h_mf_srcCor_wir == 1
          rec_h_mf_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_h_mf_srcInc)
          rec_h_mf_srcInc_wir = 1/(2*length(rec_targEv));
        else
          rec_h_mf_srcInc_wir = length(rec_h_mf_srcInc) / (length(rec_h_mf_srcCor) + length(rec_h_mf_srcInc));
        end
        if rec_h_mf_srcInc_wir == 1
          rec_h_mf_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % miss - mu
        if isempty(rec_m_mu_srcCor)
          rec_m_mu_srcCor_wir = 1/(2*length(rec_targEv));
        else
          rec_m_mu_srcCor_wir = length(rec_m_mu_srcCor) / (length(rec_m_mu_srcCor) + length(rec_m_mu_srcInc));
        end
        if rec_m_mu_srcCor_wir == 1
          rec_m_mu_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_m_mu_srcInc)
          rec_m_mu_srcInc_wir = 1/(2*length(rec_targEv));
        else
          rec_m_mu_srcInc_wir = length(rec_m_mu_srcInc) / (length(rec_m_mu_srcCor) + length(rec_m_mu_srcInc));
        end
        if rec_m_mu_srcInc_wir == 1
          rec_m_mu_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
        end
        % miss - du
        if isempty(rec_m_du_srcCor)
          rec_m_du_srcCor_wir = 1/(2*length(rec_targEv));
        else
          rec_m_du_srcCor_wir = length(rec_m_du_srcCor) / (length(rec_m_du_srcCor) + length(rec_m_du_srcInc));
        end
        if rec_m_du_srcCor_wir == 1
          rec_m_du_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
        end
        if isempty(rec_m_du_srcInc)
          rec_m_du_srcInc_wir = 1/(2*length(rec_targEv));
        else
          rec_m_du_srcInc_wir = length(rec_m_du_srcInc) / (length(rec_m_du_srcCor) + length(rec_m_du_srcInc));
        end
        if rec_m_du_srcInc_wir == 1
          rec_m_du_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
        end
        
%%         % fa - rs
%         if isempty(rec_fa_rs_srcCor)
%           rec_fa_rs_srcCor_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_rs_srcCor_wir = length(rec_fa_rs_srcCor) / (length(rec_fa_rs_srcCor) + length(rec_fa_rs_srcInc));
%         end
%         if rec_fa_rs_srcCor_wir == 1
%           rec_fa_rs_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         if isempty(rec_fa_rs_srcInc)
%           rec_fa_rs_srcInc_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_rs_srcInc_wir = length(rec_fa_rs_srcInc) / (length(rec_fa_rs_srcCor) + length(rec_fa_rs_srcInc));
%         end
%         if rec_fa_rs_srcInc_wir == 1
%           rec_fa_rs_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         % fa - ro
%         if isempty(rec_fa_ro_srcCor)
%           rec_fa_ro_srcCor_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_ro_srcCor_wir = length(rec_fa_ro_srcCor) / (length(rec_fa_ro_srcCor) + length(rec_fa_ro_srcInc));
%         end
%         if rec_fa_ro_srcCor_wir == 1
%           rec_fa_ro_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         if isempty(rec_fa_ro_srcInc)
%           rec_fa_ro_srcInc_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_ro_srcInc_wir = length(rec_fa_ro_srcInc) / (length(rec_fa_ro_srcCor) + length(rec_fa_ro_srcInc));
%         end
%         if rec_fa_ro_srcInc_wir == 1
%           rec_fa_ro_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         % fa - df
%         if isempty(rec_fa_df_srcCor)
%           rec_fa_df_srcCor_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_df_srcCor_wir = length(rec_fa_df_srcCor) / (length(rec_fa_df_srcCor) + length(rec_fa_df_srcInc));
%         end
%         if rec_fa_df_srcCor_wir == 1
%           rec_fa_df_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         if isempty(rec_fa_df_srcInc)
%           rec_fa_df_srcInc_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_df_srcInc_wir = length(rec_fa_df_srcInc) / (length(rec_fa_df_srcCor) + length(rec_fa_df_srcInc));
%         end
%         if rec_fa_df_srcInc_wir == 1
%           rec_fa_df_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         % fa - mf
%         if isempty(rec_fa_mf_srcCor)
%           rec_fa_mf_srcCor_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_mf_srcCor_wir = length(rec_fa_mf_srcCor) / (length(rec_fa_mf_srcCor) + length(rec_fa_mf_srcInc));
%         end
%         if rec_fa_mf_srcCor_wir == 1
%           rec_fa_mf_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         if isempty(rec_fa_mf_srcInc)
%           rec_fa_mf_srcInc_wir = 1/(2*length(rec_targEv));
%         else
%           rec_fa_mf_srcInc_wir = length(rec_fa_mf_srcInc) / (length(rec_fa_mf_srcCor) + length(rec_fa_mf_srcInc));
%         end
%         if rec_fa_mf_srcInc_wir == 1
%           rec_fa_mf_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         % cr - mu
%         if isempty(rec_cr_mu_srcCor)
%           rec_cr_mu_srcCor_wir = 1/(2*length(rec_targEv));
%         else
%           rec_cr_mu_srcCor_wir = length(rec_cr_mu_srcCor) / (length(rec_cr_mu_srcCor) + length(rec_cr_mu_srcInc));
%         end
%         if rec_cr_mu_srcCor_wir == 1
%           rec_cr_mu_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         if isempty(rec_cr_mu_srcInc)
%           rec_cr_mu_srcInc_wir = 1/(2*length(rec_targEv));
%         else
%           rec_cr_mu_srcInc_wir = length(rec_cr_mu_srcInc) / (length(rec_cr_mu_srcCor) + length(rec_cr_mu_srcInc));
%         end
%         if rec_cr_mu_srcInc_wir == 1
%           rec_cr_mu_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         % cr - du
%         if isempty(rec_cr_du_srcCor)
%           rec_cr_du_srcCor_wir = 1/(2*length(rec_targEv));
%         else
%           rec_cr_du_srcCor_wir = length(rec_cr_du_srcCor) / (length(rec_cr_du_srcCor) + length(rec_cr_du_srcInc));
%         end
%         if rec_cr_du_srcCor_wir == 1
%           rec_cr_du_srcCor_wir = 1 - (1/(2*length(rec_targEv)));
%         end
%         if isempty(rec_cr_du_srcInc)
%           rec_cr_du_srcInc_wir = 1/(2*length(rec_targEv));
%         else
%           rec_cr_du_srcInc_wir = length(rec_cr_du_srcInc) / (length(rec_cr_du_srcCor) + length(rec_cr_du_srcInc));
%         end
%         if rec_cr_du_srcInc_wir == 1
%           rec_cr_du_srcInc_wir = 1 - (1/(2*length(rec_targEv)));
%         end
        
        %% Recognition HRs
        %
        % decreasing confidence that the stimulus was old
        %rec_roc_h = [length(rec_h_rs) length(rec_h_ro) length(rec_h_df) length(rec_h_mf) length(rec_m_mu) length(rec_m_du)];
        rec_roc_h = [length(rec_h_df) length(rec_h_mf) length(rec_m_mu) length(rec_m_du)];
        % fix for when a particular bin has no events (RscorePlus_doc p. 5)
        % option 1: 1/m
        rec_roc_h(rec_roc_h == 0) = 1/(length(rec_roc_h));
        % option 2: log-linear correction (Hautus & Lee, 1998)
        %if sum(rec_roc_h == 0) > 0
        %  rec_roc_h = rec_roc_h + 1/(length(confRange_rec));
        %end
        rec_roc_hr = cumsum(rec_roc_h(1:end-1) / sum(rec_roc_h));
        % fix for when rate == 1 or 0, because norminv can't handle that
        rec_roc_hr(rec_roc_hr >= 1) = 1 - 1/(2*length(rec_targEv));
        rec_roc_hr(rec_roc_hr == 0) = 1/(2*length(rec_targEv));
        rec_roc_zhr = norminv(rec_roc_hr,0,1);
        
        % Recognition FARs
        %
        % decreasing confidence that the stimulus was old
        %rec_roc_fa = [length(rec_fa_rs) length(rec_fa_ro) length(rec_fa_df) length(rec_fa_mf) length(rec_cr_mu) length(rec_cr_du)];
        rec_roc_fa = [length(rec_fa_df) length(rec_fa_mf) length(rec_cr_mu) length(rec_cr_du)];
        % fix for when a particular bin has no events (RscorePlus_doc p. 5)
        % option 1: 1/m
        rec_roc_fa(rec_roc_fa == 0) = 1/(length(rec_roc_fa));
        % option 2: log-linear correction (Hautus & Lee, 1998)
        %if sum(rec_roc_fa == 0) > 0
        %  rec_roc_fa = rec_roc_fa + 1/(length(confRange_rec));
        %end
        rec_roc_far = cumsum(rec_roc_fa(1:end-1) / sum(rec_roc_fa));
        % fix for when rate == 1 or 0, because norminv can't handle that
        rec_roc_far(rec_roc_far >= 1) = 1 - 1/(2*length(rec_lureEv));
        rec_roc_far(rec_roc_far == 0) = 1/(2*length(rec_lureEv));
        rec_roc_zfar = norminv(rec_roc_far,0,1);
        
        % Find da: Macmillan & Creelman, p. 61--62
        
        % find the slope and y-intercept of the zROC
        rec_MB = polyfit(rec_roc_zfar,rec_roc_zhr,1);
        item_z_slope = rec_MB(1);
        %rec_z_yint = rec_MB(2);
        
        % using the equation y = MB(1)x + MB(2):
        %
        % find x (dp1) when y = 0
        %rec_dp1 = polyval(rec_z_yint / item_z_slope,0); % shouldn't this be negative?
        % find y (dp2) when x = 0
        rec_dp2 = polyval(rec_MB,0);
        
        % alternative methods of calculating dp1 and dp2
        %rec_dp1 = ((1/item_z_slope) * rec_roc_zhr) - rec_roc_zfar;
        %rec_dp2 = rec_roc_zhr - (item_z_slope * rec_roc_zfar);
        %rec_dp1 = rec_dp2 / item_z_slope;
        
        % calculate da
        item_da = (2 / (1 + (item_z_slope^2)))^(1/2) * rec_dp2;
        %item_da = (2 / (1 + (item_z_slope^2)))^(1/2) * (rec_roc_zhr - (item_z_slope * rec_roc_zfar));
        
%         % plot zROC
%         figure
%         plot([0 0],[-3 3],'k--'); % vert
%         hold on
%         plot([-3 3],[0 0],'k--'); % horiz
%         plot([-3 3],[-3 3],'k-');% diag
%         X = [min(rec_roc_zfar) max(rec_roc_zfar)];
%         Ms = ones(1,2)*rec_MB(1);
%         Bs = ones(1,2)*rec_MB(2);
%         Y = Ms.*X + Bs;
%         plot(X,Y,'r-');
%         plot(rec_roc_zfar,rec_roc_zhr,'ko'); % zROC
%         %plot([0 rec_roc_far 1],[0 rec_roc_hr 1],'k*'); % ROC
%         title('Item zROC');
%         axis([-3 3 -3 3])
%         axis square
        
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
        
        %% REACTION TIMES
        rec_h_rt = mean(getStructField(rec_h,'rec_rt'));
        rec_m_rt = mean(getStructField(rec_m,'rec_rt'));
        rec_cr_rt = mean(getStructField(rec_cr,'rec_rt'));
        rec_fa_rt = mean(getStructField(rec_fa,'rec_rt'));
        
        src_h_rt = mean(getStructField(src_h,'src_rt'));
        src_m_rt = mean(getStructField(src_m,'src_rt'));
        src_cr_rt = mean(getStructField(src_cr,'src_rt'));
        src_fa_rt = mean(getStructField(src_fa,'src_rt'));
        
        rec_h_srcCor_rt = mean(getStructField(rec_h_srcCor,'rec_rt'));
        rec_h_srcInc_rt = mean(getStructField(rec_h_srcInc,'rec_rt'));
        
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
        % Don't currently have all event information if using rejectArt; only have
        % info for recognition hits and correct rejections. Source hits, CR,
        % missed, and FAs are ok because they come from recognition hits. This is
        % because we're not segmenting misses and false alarms in NS.
        
        %% RAW NUMBERS
        
        % get all the targets (old study items presented at test)
        rec_targEv = filterStruct(subSesEv,'rec_isTarg == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
        % get all the lures (new test items)
        rec_lureEv = filterStruct(subSesEv,'rec_isTarg == 0 & ismember(type,varargin{1})',{'LURE_PRES'});
        
        % all source reponses (targets and lures) are conditionalized on
        % getting a recognition hit
        %
        % source targets (for Hs and Ms) and lures (for CRs and FAs) come
        % from the source coding of target items, as explained in
        % rsrc_createEvents.m
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
        % explained in rsrc_createEvents.m
        src_cr = filterStruct(src_lureEv,'src_correct == 1');
        src_fa = filterStruct(src_lureEv,'src_correct == 0');
        
        % recognition by confidence level
        rec_h_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(1));
        rec_h_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(2));
        rec_h_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(3));
        rec_cr_hi = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(6));
        rec_cr_md = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(5));
        rec_cr_lo = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(4));
        
        % source by confidence level
        src_h_hi = filterStruct(src_targEv,'src_correct == 1 & src_conf == varargin{1}',confRange_src(6));
        src_h_md = filterStruct(src_targEv,'src_correct == 1 & src_conf == varargin{1}',confRange_src(5));
        src_h_lo = filterStruct(src_targEv,'src_correct == 1 & src_conf == varargin{1}',confRange_src(4));
        src_m_hi = filterStruct(src_targEv,'src_correct == 0 & src_conf == varargin{1}',confRange_src(1));
        src_m_md = filterStruct(src_targEv,'src_correct == 0 & src_conf == varargin{1}',confRange_src(2));
        src_m_lo = filterStruct(src_targEv,'src_correct == 0 & src_conf == varargin{1}',confRange_src(3));
        src_cr_hi = filterStruct(src_lureEv,'src_correct == 1 & src_conf == varargin{1}',confRange_src(1));
        src_cr_md = filterStruct(src_lureEv,'src_correct == 1 & src_conf == varargin{1}',confRange_src(2));
        src_cr_lo = filterStruct(src_lureEv,'src_correct == 1 & src_conf == varargin{1}',confRange_src(3));
        src_fa_hi = filterStruct(src_lureEv,'src_correct == 0 & src_conf == varargin{1}',confRange_src(6));
        src_fa_md = filterStruct(src_lureEv,'src_correct == 0 & src_conf == varargin{1}',confRange_src(5));
        src_fa_lo = filterStruct(src_lureEv,'src_correct == 0 & src_conf == varargin{1}',confRange_src(4));
        
        % recognition hit and source correct/incorrect (correct = hit/CR,
        % incorrect= miss/FA)
        rec_h_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1');
        rec_h_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0');
        
        % recognition hit/miss and source hit/miss/CR/FA combos
        %
        % rec hit hi
        %
        % recognition hit, source correct/incorrect
        rec_h_hi_srcCor_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(1),confRange_src(1),confRange_src(6));
        rec_h_hi_srcCor_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(1),confRange_src(2),confRange_src(5));
        rec_h_hi_srcCor_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(1),confRange_src(3),confRange_src(4));
        rec_h_hi_srcInc_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(1),confRange_src(1),confRange_src(6));
        rec_h_hi_srcInc_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(1),confRange_src(2),confRange_src(5));
        rec_h_hi_srcInc_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(1),confRange_src(3),confRange_src(4));
        % rec hit med
        %
        % recognition hit, source correct/incorrect
        rec_h_md_srcCor_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(2),confRange_src(1),confRange_src(6));
        rec_h_md_srcCor_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(2),confRange_src(2),confRange_src(5));
        rec_h_md_srcCor_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(2),confRange_src(3),confRange_src(4));
        rec_h_md_srcInc_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(2),confRange_src(1),confRange_src(6));
        rec_h_md_srcInc_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(2),confRange_src(2),confRange_src(5));
        rec_h_md_srcInc_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(2),confRange_src(3),confRange_src(4));
        % rec hit low
        %
        % recognition hit, source correct/incorrect
        rec_h_lo_srcCor_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(3),confRange_src(1),confRange_src(6));
        rec_h_lo_srcCor_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(3),confRange_src(2),confRange_src(5));
        rec_h_lo_srcCor_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(3),confRange_src(3),confRange_src(4));
        rec_h_lo_srcInc_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(3),confRange_src(1),confRange_src(6));
        rec_h_lo_srcInc_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(3),confRange_src(2),confRange_src(5));
        rec_h_lo_srcInc_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(3),confRange_src(3),confRange_src(4));
        
        %% RATES
        
        % accuracy rates
        src_hr = length(src_h) / length(src_targEv);
        src_mr = length(src_m) / length(src_targEv);
        src_crr = length(src_cr) / length(src_lureEv);
        src_far = length(src_fa) / length(src_lureEv);
        
        % rates for recognition hits/misses
        rec_h_srcCor_r = length(rec_h_srcCor) / length(rec_h);
        rec_h_srcInc_r = length(rec_h_srcInc) / length(rec_h);
        
        % rates by confidence
        rec_h_hi_r = length(rec_h_hi) / length(rec_targEv);
        rec_h_md_r = length(rec_h_md) / length(rec_targEv);
        rec_h_lo_r = length(rec_h_lo) / length(rec_targEv);
        rec_cr_hi_r = length(rec_cr_hi) / length(rec_lureEv);
        rec_cr_md_r = length(rec_cr_md) / length(rec_lureEv);
        rec_cr_lo_r = length(rec_cr_lo) / length(rec_lureEv);
        
%         % rates by confidence
%         src_h_hi_r = length(src_h_hi) / length(src_targEv);
%         src_h_md_r = length(src_h_md) / length(src_targEv);
%         src_h_lo_r = length(src_h_lo) / length(src_targEv);
%         src_m_hi_r = length(src_m_hi) / length(src_targEv);
%         src_m_md_r = length(src_m_md) / length(src_targEv);
%         src_m_lo_r = length(src_m_lo) / length(src_targEv);
%         src_cr_hi_r = length(src_cr_hi) / length(src_lureEv);
%         src_cr_md_r = length(src_cr_md) / length(src_lureEv);
%         src_cr_lo_r = length(src_cr_lo) / length(src_lureEv);
%         src_fa_hi_r = length(src_fa_hi) / length(src_lureEv);
%         src_fa_md_r = length(src_fa_md) / length(src_lureEv);
%         src_fa_lo_r = length(src_fa_lo) / length(src_lureEv);
        
        % RH-SCor: rates by confidence
        rec_h_hi_srcCor_hi_r = length(rec_h_hi_srcCor_hi) / length(rec_targEv);
        rec_h_hi_srcCor_md_r = length(rec_h_hi_srcCor_md) / length(rec_targEv);
        rec_h_hi_srcCor_lo_r = length(rec_h_hi_srcCor_lo) / length(rec_targEv);
        rec_h_md_srcCor_hi_r = length(rec_h_md_srcCor_hi) / length(rec_targEv);
        rec_h_md_srcCor_md_r = length(rec_h_md_srcCor_md) / length(rec_targEv);
        rec_h_md_srcCor_lo_r = length(rec_h_md_srcCor_lo) / length(rec_targEv);
        rec_h_lo_srcCor_hi_r = length(rec_h_lo_srcCor_hi) / length(rec_targEv);
        rec_h_lo_srcCor_md_r = length(rec_h_lo_srcCor_md) / length(rec_targEv);
        rec_h_lo_srcCor_lo_r = length(rec_h_lo_srcCor_lo) / length(rec_targEv);
        % RH-SInc: rates by confidence
        rec_h_hi_srcInc_hi_r = length(rec_h_hi_srcInc_hi) / length(rec_targEv);
        rec_h_hi_srcInc_md_r = length(rec_h_hi_srcInc_md) / length(rec_targEv);
        rec_h_hi_srcInc_lo_r = length(rec_h_hi_srcInc_lo) / length(rec_targEv);
        rec_h_md_srcInc_hi_r = length(rec_h_md_srcInc_hi) / length(rec_targEv);
        rec_h_md_srcInc_md_r = length(rec_h_md_srcInc_md) / length(rec_targEv);
        rec_h_md_srcInc_lo_r = length(rec_h_md_srcInc_lo) / length(rec_targEv);
        rec_h_lo_srcInc_hi_r = length(rec_h_lo_srcInc_hi) / length(rec_targEv);
        rec_h_lo_srcInc_md_r = length(rec_h_lo_srcInc_md) / length(rec_targEv);
        rec_h_lo_srcInc_lo_r = length(rec_h_lo_srcInc_lo) / length(rec_targEv);
        
        % Source HRs
        %
        % decreasing confidence that the stimulus was from target source
        src_roc_h = [length(src_h_hi) length(src_h_md) length(src_h_lo) length(src_m_lo) length(src_m_md) length(src_m_hi)];
        % fix for when a particular bin has no events (RscorePlus_doc p. 5)
        % option 1: 1/m
        src_roc_h(src_roc_h == 0) = 1/(length(confRange_rec));
        % option 2: log-linear correction (Hautus & Lee, 1998)
        %if sum(src_roc_h == 0) > 0
        %  src_roc_h = src_roc_h + 1/(length(confRange_rec));
        %end
        src_roc_hr = cumsum(src_roc_h(1:end-1) / length(src_targEv));
        % fix for when rate == 1 or 0, because norminv can't handle that
        src_roc_hr(src_roc_hr >= 1) = 1 - 1/(2*length(src_targEv));
        src_roc_hr(src_roc_hr == 0) = 1/(2*length(src_targEv));
        src_roc_zhr = norminv(src_roc_hr,0,1);
        
        % Source FARs
        %
        % decreasing confidence that the stimulus was from target source
        src_roc_fa = [length(src_fa_hi) length(src_fa_md) length(src_fa_lo) length(src_cr_lo) length(src_cr_md) length(src_cr_hi)];
        % fix for when a particular bin has no events (RscorePlus_doc p. 5)
        % option 1: 1/m
        src_roc_fa(src_roc_fa == 0) = 1/(length(confRange_rec));
        % option 2: log-linear correction (Hautus & Lee, 1998)
        %if sum(src_roc_fa == 0) > 0
        %  src_roc_fa = src_roc_fa + 1/(length(confRange_rec));
        %end
        src_roc_far = cumsum(src_roc_fa(1:end-1) / length(src_lureEv));
        % fix for when rate == 1 or 0, because norminv can't handle that
        src_roc_far(src_roc_far >= 1) = 1 - 1/(2*length(src_lureEv));
        src_roc_far(src_roc_far == 0) = 1/(2*length(src_lureEv));
        src_roc_zfar = norminv(src_roc_far,0,1);
        
        % Find da: Macmillan & Creelman, p. 61--62
        
        % find the slope and y-intercept of the zROC
        src_MB = polyfit(src_roc_zfar,src_roc_zhr,1);
        source_z_slope = src_MB(1);
        %src_z_yint = src_MB(2);
        
        % using the equation y = MB(1)x + MB(2):
        %
        % find x (dp1) when y = 0
        %src_dp1 = polyval(src_z_yint / source_z_slope,0); % shouldn't this be negative?
        % find y (dp2) when x = 0
        src_dp2 = polyval(src_MB,0);
        
        % alternative methods of calculating dp1 and dp2
        %src_dp1 = ((1/source_z_slope) * src_roc_zhr) - src_roc_zfar;
        %src_dp2 = src_roc_zhr - (source_z_slope * src_roc_zfar);
        %src_dp1 = src_dp2 / source_z_slope;
        
        % calculate da
        source_da = (2 / (1 + (source_z_slope^2)))^(1/2) * src_dp2;
        %source_da = (2 / (1 + (source_z_slope^2)))^(1/2) * (src_roc_zhr - (source_z_slope * src_roc_zfar));
        
%         % plot zROC
%         figure
%         plot([0 0],[-3 3],'k--'); % vert
%         hold on
%         plot([-3 3],[0 0],'k--'); % horiz
%         plot([-3 3],[-3 3],'k-');% diag
%         X = [min(src_roc_zfar) max(src_roc_zfar)];
%         Ms = ones(1,2)*src_MB(1);
%         Bs = ones(1,2)*src_MB(2);
%         Y = Ms.*X + Bs;
%         plot(X,Y,'r-');
%         plot(src_roc_zfar,src_roc_zhr,'ko'); % zROC
%         %plot([0 src_roc_far 1],[0 src_roc_hr 1],'k*'); % ROC
%         title('Item zROC');
%         axis([-3 3 -3 3])
%         axis square
        
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
        
        %% REACTION TIMES
        rec_h_rt = mean(getStructField(rec_h,'rec_rt'));
        rec_cr_rt = mean(getStructField(rec_cr,'rec_rt'));
        
        src_h_rt = mean(getStructField(src_h,'src_rt'));
        src_m_rt = mean(getStructField(src_m,'src_rt'));
        src_cr_rt = mean(getStructField(src_cr,'src_rt'));
        src_fa_rt = mean(getStructField(src_fa,'src_rt'));
        
        rec_h_srcCor_rt = mean(getStructField(rec_h_srcCor,'rec_rt'));
        rec_h_srcInc_rt = mean(getStructField(rec_h_srcInc,'rec_rt'));
        
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
          length(rec_targEv),length(rec_lureEv),length(src_targEv),length(src_lureEv),... % raw numbers of targets, lures
          item_dp,source_dp,item_da,item_z_slope,item_c,source_c,irk_fam,... % accuracy
          length(rec_h),length(rec_m),length(rec_cr),length(rec_fa),length(src_h),length(src_m),length(src_cr),length(src_fa),... % raw numbers for accuracy
          rec_hr,rec_mr,rec_crr,rec_far,src_hr,src_mr,src_crr,src_far,... % accuracy rates
          length(rec_h_srcCor),length(rec_h_srcInc),length(rec_m_srcCor),length(rec_m_srcInc),... % raw numbers for recognition hits/misses
          rec_h_srcCor_r,rec_h_srcInc_r,rec_m_srcCor_r,rec_m_srcInc_r,... % rates for recognition hits/misses
          length(rec_h_rs),length(rec_h_ro),length(rec_h_df),length(rec_h_mf),length(rec_m_mu),length(rec_m_du),length(rec_fa_rs),length(rec_fa_ro),length(rec_fa_df),length(rec_fa_mf),length(rec_cr_mu),length(rec_cr_du),... % raw numbers by confidence collapsed across source accuracy
          rec_h_rs_r,rec_h_ro_r,rec_h_df_r,rec_h_mf_r,rec_m_mu_r,rec_m_du_r,rec_fa_rs_r,rec_fa_ro_r,rec_fa_df_r,rec_fa_mf_r,rec_cr_mu_r,rec_cr_du_r,... % rates by confidence collapsed across source accuracy
          length(rec_h_rs_srcCor),length(rec_h_ro_srcCor),length(rec_h_df_srcCor),length(rec_h_mf_srcCor),length(rec_m_mu_srcCor),length(rec_m_du_srcCor),... % raw numbers for recognition hits/misses by confidence for source correct
          length(rec_h_rs_srcInc),length(rec_h_ro_srcInc),length(rec_h_df_srcInc),length(rec_h_mf_srcInc),length(rec_m_mu_srcInc),length(rec_m_du_srcInc),... % raw numbers for recognition hits/misses by confidence for source incorrect
          rec_h_rs_srcCor_r,rec_h_ro_srcCor_r,rec_h_df_srcCor_r,rec_h_mf_srcCor_r,rec_m_mu_srcCor_r,rec_m_du_srcCor_r,... % rates for recognition hits/misses by confidence for source correct
          rec_h_rs_srcInc_r,rec_h_ro_srcInc_r,rec_h_df_srcInc_r,rec_h_mf_srcInc_r,rec_m_mu_srcInc_r,rec_m_du_srcInc_r,... % rates for recognition hits/misses by confidence for source incorrect
          rec_h_rs_srcCor_wir,rec_h_rs_srcInc_wir,rec_h_ro_srcCor_wir,rec_h_ro_srcInc_wir,rec_h_df_srcCor_wir,rec_h_df_srcInc_wir,rec_h_mf_srcCor_wir,rec_h_mf_srcInc_wir,rec_m_mu_srcCor_wir,rec_m_mu_srcInc_wir,rec_m_du_srcCor_wir,rec_m_du_srcInc_wir,... % within-rates for recognition hits/misses by confidence
          rec_h_rt,rec_m_rt,rec_cr_rt,rec_fa_rt,... % rec reaction times
          src_h_rt,src_m_rt,src_cr_rt,src_fa_rt,... % src reaction times
          rec_h_srcCor_rt,rec_h_srcInc_rt,... % recognition hits reaction times
          Pr,Br,...
          ];
      elseif rejectArt == 1
        tableData = [...
          length(rec_targEv),length(rec_lureEv),length(src_targEv),length(src_lureEv),... % raw numbers of targets, lures
          source_dp,source_da,source_z_slope,source_c,... % accuracy
          length(rec_h),length(rec_cr),length(src_h),length(src_m),length(src_cr),length(src_fa),... % raw numbers for accuracy
          src_hr,src_mr,src_crr,src_far,... % accuracy rates
          length(rec_h_srcCor),length(rec_h_srcInc),... % raw numbers for recognition hits/misses
          rec_h_srcCor_r,rec_h_srcInc_r,... % rates for recognition hits/misses
          length(rec_h_hi),length(rec_h_md),length(rec_h_lo),length(rec_cr_hi),length(rec_cr_md),length(rec_cr_lo),... % raw numbers by confidence
          rec_h_hi_r,rec_h_md_r,rec_h_lo_r,rec_cr_hi_r,rec_cr_md_r,rec_cr_lo_r,... % rates by confidence
          length(rec_h_hi_srcCor_hi),length(rec_h_hi_srcCor_md),length(rec_h_hi_srcCor_lo),length(rec_h_md_srcCor_hi),length(rec_h_md_srcCor_md),length(rec_h_md_srcCor_lo),length(rec_h_lo_srcCor_hi),length(rec_h_lo_srcCor_md),length(rec_h_lo_srcCor_lo),... % RH-SCor: raw numbers by confidence
          length(rec_h_hi_srcInc_hi),length(rec_h_hi_srcInc_md),length(rec_h_hi_srcInc_lo),length(rec_h_md_srcInc_hi),length(rec_h_md_srcInc_md),length(rec_h_md_srcInc_lo),length(rec_h_lo_srcInc_hi),length(rec_h_lo_srcInc_md),length(rec_h_lo_srcInc_lo),... % RH-SInc: raw numbers by confidence
          rec_h_hi_srcCor_hi_r,rec_h_hi_srcCor_md_r,rec_h_hi_srcCor_lo_r,rec_h_md_srcCor_hi_r,rec_h_md_srcCor_md_r,rec_h_md_srcCor_lo_r,rec_h_lo_srcCor_hi_r,rec_h_lo_srcCor_md_r,rec_h_lo_srcCor_lo_r,... % RH-SCor: rates by confidence
          rec_h_hi_srcInc_hi_r,rec_h_hi_srcInc_md_r,rec_h_hi_srcInc_lo_r,rec_h_md_srcInc_hi_r,rec_h_md_srcInc_md_r,rec_h_md_srcInc_lo_r,rec_h_lo_srcInc_hi_r,rec_h_lo_srcInc_md_r,rec_h_lo_srcInc_lo_r,... % RH-SInc: rates by confidence
          rec_h_rt,rec_cr_rt,... % rec reaction times
          src_h_rt,src_m_rt,src_cr_rt,src_fa_rt,... % src reaction times
          rec_h_srcCor_rt,rec_h_srcInc_rt,... % recognition hits reaction times
          ];
      end
      
      %clear events subSesEv rec_* src_* study*
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
