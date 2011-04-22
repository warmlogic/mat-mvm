function rsrc2_prepData_summary(rejectArt,saveFiles,dataroot)
% function rsrc2_prepData_summary(dataroot)
%
% Save column-formatted data as a csv file
%
% Looks for each subject's events.mat and saves a summary file in dataroot
%

expName = 'RSRC2';

if nargin < 3
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
  if nargin < 2
    saveFiles = 1;
    if nargin < 1
      rejectArt = 1;
    end
  end
end

subjects = {
    'RSRC2001'...
    'RSRC2002'...
    'RSRC2003'...
    'RSRC2004'...
    'RSRC2005'...
    'RSRC2006'...
    'RSRC2007'...
    'RSRC2008'...
    'RSRC2009'...
    'RSRC2010'...
    'RSRC2011'...
    'RSRC2012'...
    'RSRC2014'...
    'RSRC2015'...
    'RSRC2016'...
    'RSRC2017'...
    'RSRC2018'...
    'RSRC2019'...
    'RSRC2020'...
    'RSRC2021'...
    'RSRC2022'...
    'RSRC2023'...
    'RSRC2024'...
    'RSRC2025'...
    'RSRC2026'...
    'RSRC2027'...
    'RSRC2028'...
    'RSRC2029'...
    'RSRC2030'...
    'RSRC2032'...
    'RSRC2033'...
    'RSRC2034'...
    'RSRC2041'...
    'RSRC2043'...
    'RSRC2047'...
    'RSRC2051'...
    };
%    'RSRC2013'... % fire alarm went off during blink period of final test block
%    'RSRC2031'... % did not finish session

sessions = {'session_0'};

subsets = {''};

% confidence numbers are [1 2 3 4 5 6], sure old to sure new
confRange_rec = [1 2 3 4 5 6];
% same for sure L to sure R
confRange_src = [1 2 3 4 5 6];

%events = loadEvents(fullfile(dataroot,'slork_events.mat'));
%
% make_report(dataroot,subjects,sessions,subsets,rejectArt);
%
% % make_report(dataroot,events,'',rejectArt);
% %
% % % get only the size study events
% % sizeEv = filterStruct(events,'ismember(study_resp,varargin{1})',{'BIGGER','SMALLER'});
% % % get only the life study events
% % lifeEv = filterStruct(events,'ismember(study_resp,varargin{1})',{'LIVING','NONLIVING'});
% %
% % make_report(dataroot,sizeEv,'_size',rejectArt);
% % make_report(dataroot,lifeEv,'_life',rejectArt);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_report(dataroot,subjects,sessions,subsets,rejectArt)

%% set up the header
if rejectArt == 0
  tableHeader = {'sub','ses',... % subject info
    'R targ','R lure','S targ','S lure',... % raw numbers of targs, lures
    'Item d''','Source d''','Item da','Source da','Item zROC slope','Source zROC slope','Item RespBias (c)','Source RespBias (c)',... % accuracy
    'RH','RM','RCR','RFA','SH','SM','SCR','SFA',... % raw numbers for accuracy
    'RHR','RMR','RCRR','RFAR','SHR','SMR','SCRR','SFAR',... % accuracy rates
    'RH-SCor','RH-SInc','RM-SCor','RM-SInc',... % raw numbers for recognition hits/misses
    'RH-SCor R','RH-SInc R','RM-SCor R','RM-SInc R',... % rates for recognition hits/misses
    'RH-hi','RH-md','RH-lo','RM-hi','RM-md','RM-lo','RCR-hi','RCR-md','RCR-lo','RFA-hi','RFA-md','RFA-lo',... % raw numbers by confidence
    'RH-hi R','RH-md R','RH-lo R','RM-hi R','RM-md R','RM-lo R','RCR-hi R','RCR-md R','RCR-lo R','RFA-hi R','RFA-md R','RFA-lo R',... % rates by confidence
    'RH-hi SCor-hi','RH-hi SCor-md','RH-hi SCor-lo','RH-md SCor-hi','RH-md SCor-md','RH-md SCor-lo','RH-lo SCor-hi','RH-lo SCor-md','RH-lo SCor-lo',... % RH-SCor: raw numbers by confidence
    'RH-hi SInc-hi','RH-hi SInc-md','RH-hi SInc-lo','RH-md SInc-hi','RH-md SInc-md','RH-md SInc-lo','RH-lo SInc-hi','RH-lo SInc-md','RH-lo SInc-lo',... % RH-SInc: raw numbers by confidence
    'RM-hi SCor-hi','RM-hi SCor-md','RM-hi SCor-lo','RM-md SCor-hi','RM-md SCor-md','RM-md SCor-lo','RM-lo SCor-hi','RM-lo SCor-md','RM-lo SCor-lo',... % RM-SCor: raw numbers by confidence
    'RM-hi SInc-hi','RM-hi SInc-md','RM-hi SInc-lo','RM-md SInc-hi','RM-md SInc-md','RM-md SInc-lo','RM-lo SInc-hi','RM-lo SInc-md','RM-lo SInc-lo',... % RM-SInc: raw numbers by confidence
    'RH-hi SCor-hi R','RH-hi SCor-md R','RH-hi SCor-lo R','RH-md SCor-hi R','RH-md SCor-md R','RH-md SCor-lo R','RH-lo SCor-hi R','RH-lo SCor-md R','RH-lo SCor-lo R',... % RH-SCor: rates by confidence
    'RH-hi SInc-hi R','RH-hi SInc-md R','RH-hi SInc-lo R','RH-md SInc-hi R','RH-md SInc-md R','RH-md SInc-lo R','RH-lo SInc-hi R','RH-lo SInc-md R','RH-lo SInc-lo R',... % RH-SInc: rates by confidence
    'RM-hi SCor-hi R','RM-hi SCor-md R','RM-hi SCor-lo R','RM-md SCor-hi R','RM-md SCor-md R','RM-md SCor-lo R','RM-lo SCor-hi R','RM-lo SCor-md R','RM-lo SCor-lo R',... % RM-SCor: rates by confidence
    'RM-hi SInc-hi R','RM-hi SInc-md R','RM-hi SInc-lo R','RM-md SInc-hi R','RM-md SInc-md R','RM-md SInc-lo R','RM-lo SInc-hi R','RM-lo SInc-md R','RM-lo SInc-lo R',... % RM-SInc: rates by confidence
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
    'RH-hi','RH-md','RH-lo','RCR-hi','RCR-md','RCR-lo',... % raw numbers by confidence
    'RH-hi R','RH-md R','RH-lo R','RCR-hi R','RCR-md R','RCR-lo R',... % rates by confidence
    'RH-hi SCor-hi','RH-hi SCor-md','RH-hi SCor-lo','RH-md SCor-hi','RH-md SCor-md','RH-md SCor-lo','RH-lo SCor-hi','RH-lo SCor-md','RH-lo SCor-lo',... % RH-SCor: raw numbers by confidence
    'RH-hi SInc-hi','RH-hi SInc-md','RH-hi SInc-lo','RH-md SInc-hi','RH-md SInc-md','RH-md SInc-lo','RH-lo SInc-hi','RH-lo SInc-md','RH-lo SInc-lo',... % RH-SInc: raw numbers by confidence
    'RH-hi SCor-hi R','RH-hi SCor-md R','RH-hi SCor-lo R','RH-md SCor-hi R','RH-md SCor-md R','RH-md SCor-lo R','RH-lo SCor-hi R','RH-lo SCor-md R','RH-lo SCor-lo R',... % RH-SCor: rates by confidence
    'RH-hi SInc-hi R','RH-hi SInc-md R','RH-hi SInc-lo R','RH-md SInc-hi R','RH-md SInc-md R','RH-md SInc-lo R','RH-lo SInc-hi R','RH-lo SInc-md R','RH-lo SInc-lo R',... % RH-SInc: rates by confidence
    'RH rt','RCR rt',... % rec reaction times
    'SH rt','SM rt','SCR rt','SFA rt',... % src reaction times
    'RH-SCor rt','RH-SInc rt',... % recognition hits reaction times
    };
end

if rejectArt == 1
  artSegStr = sprintf('Segments marked bad for eye artifacts');
  tableHeader{end+1} = artSegStr;
end

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
    
    %% for each session
    for ses = 1:length(sessions)
      fprintf('%s...\n',sessions{ses});
      % set the subject events directory
      eventsDir_sub = fullfile(dataroot,subjects{sub},sessions{ses},'events');
      events = loadEvents(fullfile(eventsDir_sub,'events.mat'));
      
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
        rec_h_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(1));
        rec_h_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(2));
        rec_h_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(3));
        rec_m_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(6));
        rec_m_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(5));
        rec_m_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(4));
        rec_cr_hi = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(6));
        rec_cr_md = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(5));
        rec_cr_lo = filterStruct(rec_lureEv,'rec_correct == 1 & rec_conf == varargin{1}',confRange_rec(4));
        rec_fa_hi = filterStruct(rec_lureEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(1));
        rec_fa_md = filterStruct(rec_lureEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(2));
        rec_fa_lo = filterStruct(rec_lureEv,'rec_correct == 0 & rec_conf == varargin{1}',confRange_rec(3));
        
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
        rec_m_srcCor = filterStruct(rec_targEv,'rec_correct == 0 & src_correct == 1');
        rec_m_srcInc = filterStruct(rec_targEv,'rec_correct == 0 & src_correct == 0');
        
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
        % rec miss hi
        %
        % recognition miss, source correct/incorrect
        rec_m_hi_srcCor_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(6),confRange_src(1),confRange_src(6));
        rec_m_hi_srcCor_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(6),confRange_src(2),confRange_src(5));
        rec_m_hi_srcCor_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(6),confRange_src(3),confRange_src(4));
        rec_m_hi_srcInc_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(6),confRange_src(1),confRange_src(6));
        rec_m_hi_srcInc_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(6),confRange_src(2),confRange_src(5));
        rec_m_hi_srcInc_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(6),confRange_src(3),confRange_src(4));
        % rec miss med,
        %
        % recognition miss, source correct/incorrect
        rec_m_md_srcCor_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(5),confRange_src(1),confRange_src(6));
        rec_m_md_srcCor_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(5),confRange_src(2),confRange_src(5));
        rec_m_md_srcCor_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(5),confRange_src(3),confRange_src(4));
        rec_m_md_srcInc_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(5),confRange_src(1),confRange_src(6));
        rec_m_md_srcInc_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(5),confRange_src(2),confRange_src(5));
        rec_m_md_srcInc_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(5),confRange_src(3),confRange_src(4));
        % rec miss low
        %
        % recognition miss, source correct/incorrect
        rec_m_lo_srcCor_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(4),confRange_src(1),confRange_src(6));
        rec_m_lo_srcCor_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(4),confRange_src(2),confRange_src(5));
        rec_m_lo_srcCor_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(4),confRange_src(3),confRange_src(4));
        rec_m_lo_srcInc_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(4),confRange_src(1),confRange_src(6));
        rec_m_lo_srcInc_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(4),confRange_src(2),confRange_src(5));
        rec_m_lo_srcInc_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{1} & src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_rec(4),confRange_src(3),confRange_src(4));
        
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
        
        % rates by confidence
        rec_h_hi_r = length(rec_h_hi) / length(rec_targEv);
        rec_h_md_r = length(rec_h_md) / length(rec_targEv);
        rec_h_lo_r = length(rec_h_lo) / length(rec_targEv);
        rec_m_hi_r = length(rec_m_hi) / length(rec_targEv);
        rec_m_md_r = length(rec_m_md) / length(rec_targEv);
        rec_m_lo_r = length(rec_m_lo) / length(rec_targEv);
        rec_cr_hi_r = length(rec_cr_hi) / length(rec_lureEv);
        rec_cr_md_r = length(rec_cr_md) / length(rec_lureEv);
        rec_cr_lo_r = length(rec_cr_lo) / length(rec_lureEv);
        rec_fa_hi_r = length(rec_fa_hi) / length(rec_lureEv);
        rec_fa_md_r = length(rec_fa_md) / length(rec_lureEv);
        rec_fa_lo_r = length(rec_fa_lo) / length(rec_lureEv);
        
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
        % RM-SCor: rates by confidence
        rec_m_hi_srcCor_hi_r = length(rec_m_hi_srcCor_hi) / length(rec_targEv);
        rec_m_hi_srcCor_md_r = length(rec_m_hi_srcCor_md) / length(rec_targEv);
        rec_m_hi_srcCor_lo_r = length(rec_m_hi_srcCor_lo) / length(rec_targEv);
        rec_m_md_srcCor_hi_r = length(rec_m_md_srcCor_hi) / length(rec_targEv);
        rec_m_md_srcCor_md_r = length(rec_m_md_srcCor_md) / length(rec_targEv);
        rec_m_md_srcCor_lo_r = length(rec_m_md_srcCor_lo) / length(rec_targEv);
        rec_m_lo_srcCor_hi_r = length(rec_m_lo_srcCor_hi) / length(rec_targEv);
        rec_m_lo_srcCor_md_r = length(rec_m_lo_srcCor_md) / length(rec_targEv);
        rec_m_lo_srcCor_lo_r = length(rec_m_lo_srcCor_lo) / length(rec_targEv);
        % RM-SInc: rates by confidence
        rec_m_hi_srcInc_hi_r = length(rec_m_hi_srcInc_hi) / length(rec_targEv);
        rec_m_hi_srcInc_md_r = length(rec_m_hi_srcInc_md) / length(rec_targEv);
        rec_m_hi_srcInc_lo_r = length(rec_m_hi_srcInc_lo) / length(rec_targEv);
        rec_m_md_srcInc_hi_r = length(rec_m_md_srcInc_hi) / length(rec_targEv);
        rec_m_md_srcInc_md_r = length(rec_m_md_srcInc_md) / length(rec_targEv);
        rec_m_md_srcInc_lo_r = length(rec_m_md_srcInc_lo) / length(rec_targEv);
        rec_m_lo_srcInc_hi_r = length(rec_m_lo_srcInc_hi) / length(rec_targEv);
        rec_m_lo_srcInc_md_r = length(rec_m_lo_srcInc_md) / length(rec_targEv);
        rec_m_lo_srcInc_lo_r = length(rec_m_lo_srcInc_lo) / length(rec_targEv);
        
        % Recognition HRs
        %
        % decreasing confidence that the stimulus was old
        rec_roc_h = [length(rec_h_hi) length(rec_h_md) length(rec_h_lo) length(rec_m_lo) length(rec_m_md) length(rec_m_hi)];
        % fix for when a particular bin has no events (RscorePlus_doc p. 5)
        % option 1: 1/m
        rec_roc_h(rec_roc_h == 0) = 1/(length(confRange_rec));
        % option 2: log-linear correction (Hautus & Lee, 1998)
        %if sum(rec_roc_h == 0) > 0
        %  rec_roc_h = rec_roc_h + 1/(length(confRange_rec));
        %end
        rec_roc_hr = cumsum(rec_roc_h(1:end-1) / length(rec_targEv));
        % fix for when rate == 1 or 0, because norminv can't handle that
        rec_roc_hr(rec_roc_hr >= 1) = 1 - 1/(2*length(rec_targEv));
        rec_roc_hr(rec_roc_hr == 0) = 1/(2*length(rec_targEv));
        rec_roc_zhr = norminv(rec_roc_hr,0,1);
        
        % Recognition FARs
        %
        % decreasing confidence that the stimulus was old
        rec_roc_fa = [length(rec_fa_hi) length(rec_fa_md) length(rec_fa_lo) length(rec_cr_lo) length(rec_cr_md) length(rec_cr_hi)];
        % fix for when a particular bin has no events (RscorePlus_doc p. 5)
        % option 1: 1/m
        rec_roc_fa(rec_roc_fa == 0) = 1/(length(confRange_rec));
        % option 2: log-linear correction (Hautus & Lee, 1998)
        %if sum(rec_roc_fa == 0) > 0
        %  rec_roc_fa = rec_roc_fa + 1/(length(confRange_rec));
        %end
        rec_roc_far = cumsum(rec_roc_fa(1:end-1) / length(rec_lureEv));
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
          item_dp,source_dp,item_da,source_da,item_z_slope,source_z_slope,item_c,source_c,... % accuracy
          length(rec_h),length(rec_m),length(rec_cr),length(rec_fa),length(src_h),length(src_m),length(src_cr),length(src_fa),... % raw numbers for accuracy
          rec_hr,rec_mr,rec_crr,rec_far,src_hr,src_mr,src_crr,src_far,... % accuracy rates
          length(rec_h_srcCor),length(rec_h_srcInc),length(rec_m_srcCor),length(rec_m_srcInc),... % raw numbers for recognition hits/misses
          rec_h_srcCor_r,rec_h_srcInc_r,rec_m_srcCor_r,rec_m_srcInc_r,... % rates for recognition hits/misses
          length(rec_h_hi),length(rec_h_md),length(rec_h_lo),length(rec_m_hi),length(rec_m_md),length(rec_m_lo),length(rec_cr_hi),length(rec_cr_md),length(rec_cr_lo),length(rec_fa_hi),length(rec_fa_md),length(rec_fa_lo),... % raw numbers by confidence
          rec_h_hi_r,rec_h_md_r,rec_h_lo_r,rec_m_hi_r,rec_m_md_r,rec_m_lo_r,rec_cr_hi_r,rec_cr_md_r,rec_cr_lo_r,rec_fa_hi_r,rec_fa_md_r,rec_fa_lo_r,... % rates by confidence
          length(rec_h_hi_srcCor_hi),length(rec_h_hi_srcCor_md),length(rec_h_hi_srcCor_lo),length(rec_h_md_srcCor_hi),length(rec_h_md_srcCor_md),length(rec_h_md_srcCor_lo),length(rec_h_lo_srcCor_hi),length(rec_h_lo_srcCor_md),length(rec_h_lo_srcCor_lo),... % RH-SCor: raw numbers by confidence
          length(rec_h_hi_srcInc_hi),length(rec_h_hi_srcInc_md),length(rec_h_hi_srcInc_lo),length(rec_h_md_srcInc_hi),length(rec_h_md_srcInc_md),length(rec_h_md_srcInc_lo),length(rec_h_lo_srcInc_hi),length(rec_h_lo_srcInc_md),length(rec_h_lo_srcInc_lo),... % RH-SInc: raw numbers by confidence
          length(rec_m_hi_srcCor_hi),length(rec_m_hi_srcCor_md),length(rec_m_hi_srcCor_lo),length(rec_m_md_srcCor_hi),length(rec_m_md_srcCor_md),length(rec_m_md_srcCor_lo),length(rec_m_lo_srcCor_hi),length(rec_m_lo_srcCor_md),length(rec_m_lo_srcCor_lo),... % RM-SCor: raw numbers by confidence
          length(rec_m_hi_srcInc_hi),length(rec_m_hi_srcInc_md),length(rec_m_hi_srcInc_lo),length(rec_m_md_srcInc_hi),length(rec_m_md_srcInc_md),length(rec_m_md_srcInc_lo),length(rec_m_lo_srcInc_hi),length(rec_m_lo_srcInc_md),length(rec_m_lo_srcInc_lo),... % RM-SInc: raw numbers by confidence
          rec_h_hi_srcCor_hi_r,rec_h_hi_srcCor_md_r,rec_h_hi_srcCor_lo_r,rec_h_md_srcCor_hi_r,rec_h_md_srcCor_md_r,rec_h_md_srcCor_lo_r,rec_h_lo_srcCor_hi_r,rec_h_lo_srcCor_md_r,rec_h_lo_srcCor_lo_r,... % RH-SCor: rates by confidence
          rec_h_hi_srcInc_hi_r,rec_h_hi_srcInc_md_r,rec_h_hi_srcInc_lo_r,rec_h_md_srcInc_hi_r,rec_h_md_srcInc_md_r,rec_h_md_srcInc_lo_r,rec_h_lo_srcInc_hi_r,rec_h_lo_srcInc_md_r,rec_h_lo_srcInc_lo_r,... % RH-SInc: rates by confidence
          rec_m_hi_srcCor_hi_r,rec_m_hi_srcCor_md_r,rec_m_hi_srcCor_lo_r,rec_m_md_srcCor_hi_r,rec_m_md_srcCor_md_r,rec_m_md_srcCor_lo_r,rec_m_lo_srcCor_hi_r,rec_m_lo_srcCor_md_r,rec_m_lo_srcCor_lo_r,... % RM-SCor: rates by confidence
          rec_m_hi_srcInc_hi_r,rec_m_hi_srcInc_md_r,rec_m_hi_srcInc_lo_r,rec_m_md_srcInc_hi_r,rec_m_md_srcInc_md_r,rec_m_md_srcInc_lo_r,rec_m_lo_srcInc_hi_r,rec_m_lo_srcInc_md_r,rec_m_lo_srcInc_lo_r,... % RM-SInc: rates by confidence
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

%srcCor = filterStruct(rec_targEv,'src_correct == 1');
%srcInc = filterStruct(rec_targEv,'src_correct == 0');
%srcCor_r = length(srcCor) / length(rec_targEv);
%srcInc_r = length(srcInc) / length(rec_targEv);
%         
%         srcCor_hi = filterStruct(rec_targEv,'src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_src(1),confRange_src(6));
%         srcCor_md = filterStruct(rec_targEv,'src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_src(2),confRange_src(5));
%         srcCor_lo = filterStruct(rec_targEv,'src_correct == 1 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_src(3),confRange_src(4));
%         srcInc_hi = filterStruct(rec_targEv,'src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_src(1),confRange_src(6));
%         srcInc_md = filterStruct(rec_targEv,'src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_src(2),confRange_src(5));
%         srcInc_lo = filterStruct(rec_targEv,'src_correct == 0 & (src_conf == varargin{2} | src_conf == varargin{3})',confRange_src(3),confRange_src(4));
%         % recognition hit/miss and source hit/miss combinations
%         rec_h_src_h = filterStruct(rec_targEv,'rec_correct == 1 & src_isTarg == 1 & src_correct == 1');
%         rec_h_src_m = filterStruct(rec_targEv,'rec_correct == 1 & src_isTarg == 1 & src_correct == 0');
%         rec_m_src_h = filterStruct(rec_targEv,'rec_correct == 0 & src_isTarg == 1 & src_correct == 1');
%         rec_m_src_m = filterStruct(rec_targEv,'rec_correct == 0 & src_isTarg == 1 & src_correct == 0');
%         % recognition hit/miss and source CR/FA combinations
%         rec_h_src_cr = filterStruct(rec_targEv,'rec_correct == 1 & src_isTarg == 0 & src_correct == 1');
%         rec_h_src_fa = filterStruct(rec_targEv,'rec_correct == 1 & src_isTarg == 0 & src_correct == 0');
%         rec_m_src_cr = filterStruct(rec_targEv,'rec_correct == 0 & src_isTarg == 0 & src_correct == 1');
%         rec_m_src_fa = filterStruct(rec_targEv,'rec_correct == 0 & src_isTarg == 0 & src_correct == 0');

%         % recognition hit, source hit/miss
%         rec_h_hi_src_h_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(6));
%         rec_h_hi_src_h_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(5));
%         rec_h_hi_src_h_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(4));
%         rec_h_hi_src_m_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(1));
%         rec_h_hi_src_m_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(2));
%         rec_h_hi_src_m_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(3));
%         % recognition hit, source CR/FA
%         rec_h_hi_src_cr_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(1));
%         rec_h_hi_src_cr_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(2));
%         rec_h_hi_src_cr_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(3));
%         rec_h_hi_src_fa_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(6));
%         rec_h_hi_src_fa_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(5));
%         rec_h_hi_src_fa_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{1} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(4));
%         % recognition hit, source hit/miss
%         rec_h_md_src_h_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(6));
%         rec_h_md_src_h_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(5));
%         rec_h_md_src_h_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(4));
%         rec_h_md_src_m_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(1));
%         rec_h_md_src_m_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(2));
%         rec_h_md_src_m_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(3));
%         % recognition hit, source CR/FA
%         rec_h_md_src_cr_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(1));
%         rec_h_md_src_cr_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(2));
%         rec_h_md_src_cr_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(3));
%         rec_h_md_src_fa_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(6));
%         rec_h_md_src_fa_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(5));
%         rec_h_md_src_fa_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{2} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(4));
%         % recognition hit, source hit/miss
%         rec_h_lo_src_h_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(6));
%         rec_h_lo_src_h_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(5));
%         rec_h_lo_src_h_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(4));
%         rec_h_lo_src_m_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(1));
%         rec_h_lo_src_m_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(2));
%         rec_h_lo_src_m_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(3));
%         % recognition hit, source CR/FA
%         rec_h_lo_src_cr_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(1));
%         rec_h_lo_src_cr_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(2));
%         rec_h_lo_src_cr_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(3));
%         rec_h_lo_src_fa_hi = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(6));
%         rec_h_lo_src_fa_md = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(5));
%         rec_h_lo_src_fa_lo = filterStruct(rec_targEv,'rec_correct == 1 & rec_conf == varargin{3} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(4));
%         % recognition miss, source hit/miss
%         rec_m_hi_src_h_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(6));
%         rec_m_hi_src_h_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(5));
%         rec_m_hi_src_h_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(4));
%         rec_m_hi_src_m_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(1));
%         rec_m_hi_src_m_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(2));
%         rec_m_hi_src_m_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(3));
%         % recognition miss, source CR/FA
%         rec_m_hi_src_cr_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(1));
%         rec_m_hi_src_cr_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(2));
%         rec_m_hi_src_cr_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(3));
%         rec_m_hi_src_fa_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(6));
%         rec_m_hi_src_fa_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(5));
%         rec_m_hi_src_fa_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{6} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(4));
%         % recognition miss, source hit/miss
%         rec_m_md_src_h_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(6));
%         rec_m_md_src_h_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(5));
%         rec_m_md_src_h_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(4));
%         rec_m_md_src_m_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(1));
%         rec_m_md_src_m_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(2));
%         rec_m_md_src_m_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(3));
%         % recognition miss, source CR/FA
%         rec_m_md_src_cr_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(1));
%         rec_m_md_src_cr_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(2));
%         rec_m_md_src_cr_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(3));
%         rec_m_md_src_fa_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(6));
%         rec_m_md_src_fa_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(5));
%         rec_m_md_src_fa_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{5} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(4));
%         % recognition miss, source hit/miss
%         rec_m_lo_src_h_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(6));
%         rec_m_lo_src_h_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(5));
%         rec_m_lo_src_h_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 1 & src_correct == 1 & src_conf == varargin{2}',confRange_src(4));
%         rec_m_lo_src_m_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(1));
%         rec_m_lo_src_m_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(2));
%         rec_m_lo_src_m_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 1 & src_correct == 0 & src_conf == varargin{2}',confRange_src(3));
%         % recognition miss, source CR/FA
%         rec_m_lo_src_cr_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(1));
%         rec_m_lo_src_cr_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(2));
%         rec_m_lo_src_cr_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 0 & src_correct == 1 & src_conf == varargin{2}',confRange_src(3));
%         rec_m_lo_src_fa_hi = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(6));
%         rec_m_lo_src_fa_md = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(5));
%         rec_m_lo_src_fa_lo = filterStruct(rec_targEv,'rec_correct == 0 & rec_conf == varargin{4} & src_isTarg == 0 & src_correct == 0 & src_conf == varargin{2}',confRange_src(4));

%           length(rec_h_src_h),length(rec_h_src_m),length(rec_h_src_cr),length(rec_h_src_fa),...
%           length(rec_m_src_h),length(rec_m_src_m),length(rec_m_src_cr),length(rec_m_src_fa),...
%           (length(rec_h_src_h)/length(src_targEv)),(length(rec_h_src_m)/length(src_targEv)),(length(rec_h_src_cr)/length(src_lureEv)),(length(rec_h_src_fa)/length(src_lureEv)),...
%           (length(rec_m_src_h)/length(src_targEv)),(length(rec_m_src_m)/length(src_targEv)),(length(rec_m_src_cr)/length(src_lureEv)),(length(rec_m_src_fa)/length(src_lureEv)),...
%           length(rec_h_hi_src_h_hi),length(rec_h_hi_src_h_md),length(rec_h_hi_src_h_lo),length(rec_h_hi_src_m_hi),length(rec_h_hi_src_m_md),length(rec_h_hi_src_m_lo),length(rec_h_hi_src_cr_hi),length(rec_h_hi_src_cr_md),length(rec_h_hi_src_cr_lo),length(rec_h_hi_src_fa_hi),length(rec_h_hi_src_fa_md),length(rec_h_hi_src_fa_lo),length(rec_h_md_src_h_hi),length(rec_h_md_src_h_md),length(rec_h_md_src_h_lo),length(rec_h_md_src_m_hi),length(rec_h_md_src_m_md),length(rec_h_md_src_m_lo),length(rec_h_md_src_cr_hi),length(rec_h_md_src_cr_md),length(rec_h_md_src_cr_lo),length(rec_h_md_src_fa_hi),length(rec_h_md_src_fa_md),length(rec_h_md_src_fa_lo),length(rec_h_lo_src_h_hi),length(rec_h_lo_src_h_md),length(rec_h_lo_src_h_lo),length(rec_h_lo_src_m_hi),length(rec_h_lo_src_m_md),length(rec_h_lo_src_m_lo),length(rec_h_lo_src_cr_hi),length(rec_h_lo_src_cr_md),length(rec_h_lo_src_cr_lo),length(rec_h_lo_src_fa_hi),length(rec_h_lo_src_fa_md),length(rec_h_lo_src_fa_lo),...
%           length(rec_m_hi_src_h_hi),length(rec_m_hi_src_h_md),length(rec_m_hi_src_h_lo),length(rec_m_hi_src_m_hi),length(rec_m_hi_src_m_md),length(rec_m_hi_src_m_lo),length(rec_m_hi_src_cr_hi),length(rec_m_hi_src_cr_md),length(rec_m_hi_src_cr_lo),length(rec_m_hi_src_fa_hi),length(rec_m_hi_src_fa_md),length(rec_m_hi_src_fa_lo),length(rec_m_md_src_h_hi),length(rec_m_md_src_h_md),length(rec_m_md_src_h_lo),length(rec_m_md_src_m_hi),length(rec_m_md_src_m_md),length(rec_m_md_src_m_lo),length(rec_m_md_src_cr_hi),length(rec_m_md_src_cr_md),length(rec_m_md_src_cr_lo),length(rec_m_md_src_fa_hi),length(rec_m_md_src_fa_md),length(rec_m_md_src_fa_lo),length(rec_m_lo_src_h_hi),length(rec_m_lo_src_h_md),length(rec_m_lo_src_h_lo),length(rec_m_lo_src_m_hi),length(rec_m_lo_src_m_md),length(rec_m_lo_src_m_lo),length(rec_m_lo_src_cr_hi),length(rec_m_lo_src_cr_md),length(rec_m_lo_src_cr_lo),length(rec_m_lo_src_fa_hi),length(rec_m_lo_src_fa_md),length(rec_m_lo_src_fa_lo),...
