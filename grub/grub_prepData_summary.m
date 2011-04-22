function grub_prepData_summary(rejectArt,saveFiles,dataroot)
% function grub_prepData_summary(dataroot)
%
% Save column-formatted data as a csv file
%
% Looks for each subject's events.mat and saves a summary file in dataroot
%

expName = 'GRUB';

%subsets = {'','_size','_life'};
subsets = {'_life'};

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
      rejectArt = 0;
    end
  end
end

subjects = {
  'GRUB002',...
  'GRUB003',...
  'GRUB004',...
  'GRUB005',...
  'GRUB006',...
  'GRUB007',...
  'GRUB008',...
  'GRUB009',...
  'GRUB010',...
  'GRUB011',...
  'GRUB012',...
  'GRUB014',...
  'GRUB015',...
  'GRUB016',...
  'GRUB017'...
  'GRUB018'...
  'GRUB019'...
  'GRUB020'...
  'GRUB021'...
  'GRUB022'...
  'GRUB023'...
  'GRUB024'...
  'GRUB025'...
  'GRUB026'...
  };
%    'GRUB001'... % exp was not set up correctly
%    'GRUB013'... % only half the session was recorded

sessions = {'session_0'};

%% set up the headers
if rejectArt == 0
  tableHeader = {'sub','ses',... % subject info
    'R targ','R lure','S targ','S lure',... % raw numbers of targs, lures
    'Item d''','Source d''','Item RespBias (c)','Source RespBias (c)',... % accuracy
    'RH','RM','RCR','RFA','SH','SM','SCR','SFA',... % raw numbers for accuracy
    'RHR','RMR','RCRR','RFAR','SHR','SMR','SCRR','SFAR',... % accuracy rates
    'RH-SCor','RH-SInc',... % raw numbers for recognition hits
    'RH-SCorR','RH-SIncR',... % accuracy rates for recognition hits
    'RH rt','RM rt','RCR rt','RFA rt',... % rec reaction times
    'SH rt','SM rt','SCR rt','SFA rt',... % src reaction times
    'RH-SCor rt','RH-SInc rt',... % recognition hits reaction times
    'Study rt','Study PoolCor rt','Study PoolInc rt',... % study reaction times
    'Study RH rt','Study RM rt','Study RH-SrcCor rt','Study RH-SrcInc rt',... % study rt by subsequent accuracy
    'Classification Accuracy','Rec Discrimination (Pr)','Rec Response Bias (Br)',... % more accuracy
    };
elseif rejectArt == 1
  tableHeader = {'sub','ses',... % subject info
    'R targ','R lure','S targ','S lure',... % raw numbers of targs, lures
    'Source d''','Source RespBias (c)',... % accuracy
    'RH','RCR','SH','SM','SCR','SFA',... % raw numbers for accuracy
    'SHR','SMR','SCRR','SFAR',... % accuracy rates
    'RH-SCor','RH-SInc',... % raw numbers for recognition hits
    'RH-SCorR','RH-SIncR',... % accuracy rates for recognition hits
    'RH rt','RCR rt',... % rec reaction times
    'SH rt','SM rt','SCR rt','SFA rt',... % src reaction times
    'RH-SCor rt','RH-SInc rt',... % recognition hits reaction times
    'Classification Accuracy',... % more accuracy
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
    
    % for each session
    for ses = 1:length(sessions)
      fprintf('%s...\n',sessions{ses});
      % set the subject events directory
      eventsDir_sub = fullfile(dataroot,subjects{sub},sessions{ses},'events');
      events = loadEvents(fullfile(eventsDir_sub,'events.mat'));
      
      if rejectArt == 1
        % remove artifacts
        subSesEv_all = events;
        events = filterStruct(subSesEv_all,'nsArt == 0');
      end
      
      % get overall FAR for use in calculating item d'; make sure to
      % calculate it after rejecting artifacts so that it is accurate
      rec_lureEv = filterStruct(events,'rec_isTarg == 0 & ismember(type,varargin{1})',{'LURE_PRES'});
      rec_fa = filterStruct(rec_lureEv,'rec_correct == 0');
      rec_far = length(rec_fa)/length(rec_lureEv);
      clear rec_lureEv rec_fa
      
      % get only the events for this subset
      if strcmp(subsets{s},'_size')
        % get only the size study events
        fprintf('Size\n');
        subSesEv = filterStruct(events,'ismember(study_resp,varargin{1})',{'BIGGER','SMALLER'});
      elseif strcmp(subsets{s},'_life')
        % get only the life study events
        fprintf('Life\n');
        subSesEv = filterStruct(events,'ismember(study_resp,varargin{1})',{'LIVING','NONLIVING'});
      else
        fprintf('All\n');
        subSesEv = events;
      end
      
      if rejectArt == 0
        %% RAW NUMBERS
        
        % get all the study events (excluding buffers)
        studyEv = filterStruct(subSesEv,'ismember(type,varargin{1})',{'STUDY_TARGET'});
        
        study_poolCor = filterStruct(studyEv,'pool_correct == 1');
        study_poolInc = filterStruct(studyEv,'pool_correct == 0');
        
        study_rh = filterStruct(studyEv,'rec_correct == 1');
        study_rm = filterStruct(studyEv,'rec_correct == 0');
        
        study_rh_srcCor = filterStruct(study_rh,'src_correct == 1');
        study_rh_srcInc = filterStruct(study_rh,'src_correct == 0');
        
        % get all the targets (old study items presented at test)
        rec_targEv = filterStruct(subSesEv,'rec_isTarg == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
        % get all the lures (new test items)
        rec_lureEv = filterStruct(subSesEv,'rec_isTarg == 0 & ismember(type,varargin{1})',{'LURE_PRES'});
        
        % all source reponses (targets and lures) are conditionalized on
        % getting a recognition hit
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
        %
        % source targets (for Hs and Ms) and lures (for CRs and FAs) come
        % from the source coding of target items, as explained in
        % grub_createEvents.m
        src_h = filterStruct(src_targEv,'src_correct == 1');
        src_m = filterStruct(src_targEv,'src_correct == 0');
        src_cr = filterStruct(src_lureEv,'src_correct == 1');
        src_fa = filterStruct(src_lureEv,'src_correct == 0');
        
        % recognition hit - old items called "old"
        %
        % collapsed across source hits and CRs
        rec_h_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1');
        % collapsed across source misses and FAs
        rec_h_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0');
        % can't have rec_m_srcCor/Inc because missed don't get source
        % judgments
        %
        % can't do rec_cr_srcCor/Inc and rec_fa_srcCor/Inc because lures do
        % not have sources assocaited with them
        
        %% RATES
        
        % accuracy rates
        rec_hr = length(rec_h)/length(rec_targEv);
        rec_mr = length(rec_m)/length(rec_targEv);
        rec_crr = length(rec_cr)/length(rec_lureEv);
        % instead of calculating FAR here, use the one calculated up above
        % using all items so that we can have a FAR even when we've split
        % out the two different sources
        %rec_far = length(rec_fa)/length(rec_lureEv);
        src_hr = length(src_h)/length(src_targEv);
        src_mr = length(src_m)/length(src_targEv);
        src_crr = length(src_cr)/length(src_lureEv);
        src_far = length(src_fa)/length(src_lureEv);
        
        % d primes
        item_dp = norminv(rec_hr,0,1) - norminv(rec_far,0,1);
        source_dp = norminv(src_hr,0,1) - norminv(src_far,0,1);
        
        % response bias (criterion) from Macmillan & Creelman p. 29
        item_c = -0.5 * (norminv(rec_hr,0,1) + norminv(rec_far,0,1));
        source_c = -0.5 * (norminv(src_hr,0,1) + norminv(src_far,0,1));
        
        % Classification accuracy
        pool_accuracy = mean(getStructField(rec_targEv,'pool_correct'));
        
        % From Mecklinger et al. (2007) and Corwin (1994)
        %
        % discrimination index
        Pr = rec_hr - rec_far;
        % response bias index
        Br = rec_far / (1 - Pr);
        
        % accuracy rates for recognition hits
        rec_h_srcCor_r = length(rec_h_srcCor)/length(rec_h);
        rec_h_srcInc_r = length(rec_h_srcInc)/length(rec_h);
        
        %% REACTION TIMES
        
        % study items
        studyEv_rt = mean(getStructField(studyEv,'study_rt'));
        
        study_poolCor_rt = mean(getStructField(study_poolCor,'study_rt'));
        study_poolInc_rt = mean(getStructField(study_poolInc,'study_rt'));
        
        study_rh_rt = mean(getStructField(study_rh,'study_rt'));
        study_rm_rt = mean(getStructField(study_rm,'study_rt'));
        
        study_rh_srcCor_rt = mean(getStructField(study_rh_srcCor,'study_rt'));
        study_rh_srcInc_rt = mean(getStructField(study_rh_srcInc,'study_rt'));
        
        % test items
        rec_h_rt = mean(getStructField(rec_h,'rec_rt'));
        rec_m_rt = mean(getStructField(rec_m,'rec_rt'));
        rec_cr_rt = mean(getStructField(rec_cr,'rec_rt'));
        rec_fa_rt = mean(getStructField(rec_fa,'rec_rt'));
        
        rec_h_srcCor_rt = mean(getStructField(rec_h_srcCor,'rec_rt'));
        rec_h_srcInc_rt = mean(getStructField(rec_h_srcInc,'rec_rt'));
        
        src_h_rt = mean(getStructField(src_h,'src_rt'));
        src_m_rt = mean(getStructField(src_m,'src_rt'));
        src_cr_rt = mean(getStructField(src_cr,'src_rt'));
        src_fa_rt = mean(getStructField(src_fa,'src_rt'));
        
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
        % Don't currently have all event information if using rejectArt;
        % only have info for recognition hits and correct rejections.
        % Source hits, source misses, source CRs, and source FAs are ok
        % because they come from recognition hits. This is because we're
        % not segmenting misses and false alarms in NS.
        
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
        % grub_createEvents.m
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
        % explained in grub_createEvents.m
        src_cr = filterStruct(src_lureEv,'src_correct == 1');
        src_fa = filterStruct(src_lureEv,'src_correct == 0');
        
        % recognition hit - old items called "old"
        %
        % collapsed across source hits and CRs
        rec_h_srcCor = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 1');
        % collapsed across source misses and FAs
        rec_h_srcInc = filterStruct(rec_targEv,'rec_correct == 1 & src_correct == 0');
        % can't have rec_m_srcCor/Inc because missed don't get source judgments
        %
        % can't do rec_cr_srcCor/Inc and rec_fa_srcCor/Inc because lures do not
        % have sources assocaited with them
        
        %% RATES
        
        % accuracy rates
        src_hr = length(src_h)/length(src_targEv);
        src_mr = length(src_m)/length(src_targEv);
        src_crr = length(src_cr)/length(src_lureEv);
        src_far = length(src_fa)/length(src_lureEv);
        
        % d primes
        source_dp = norminv(src_hr,0,1) - norminv(src_far,0,1);
        
        % response bias (criterion) from Macmillan & Creelman p. 29
        source_c = -0.5 * (norminv(src_hr,0,1) + norminv(src_far,0,1));
        
        % Classification accuracy
        pool_accuracy = mean(getStructField(rec_targEv,'pool_correct'));
        
        % accuracy rates for recognition hits
        rec_h_srcCor_r = length(rec_h_srcCor)/length(rec_h);
        rec_h_srcInc_r = length(rec_h_srcInc)/length(rec_h);
        
        %% REACTION TIMES
        
        % test items
        rec_h_rt = mean(getStructField(rec_h,'rec_rt'));
        rec_cr_rt = mean(getStructField(rec_cr,'rec_rt'));
        
        rec_h_srcCor_rt = mean(getStructField(rec_h_srcCor,'rec_rt'));
        rec_h_srcInc_rt = mean(getStructField(rec_h_srcInc,'rec_rt'));
        
        src_h_rt = mean(getStructField(src_h,'src_rt'));
        src_m_rt = mean(getStructField(src_m,'src_rt'));
        src_cr_rt = mean(getStructField(src_cr,'src_rt'));
        src_fa_rt = mean(getStructField(src_fa,'src_rt'));
        
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
          rec_h_rt,rec_m_rt,rec_cr_rt,rec_fa_rt,... % rec reaction times
          src_h_rt,src_m_rt,src_cr_rt,src_fa_rt,... % src reaction times
          rec_h_srcCor_rt,rec_h_srcInc_rt,... % recognition hits reaction times
          studyEv_rt,study_poolCor_rt,study_poolInc_rt,... % study reaction times
          study_rh_rt,study_rm_rt,study_rh_srcCor_rt,study_rh_srcInc_rt,... % study rt by subsequent accuracy
          pool_accuracy,Pr,Br,... % more accuracy
          ];
      elseif rejectArt == 1
        tableData = [...
          length(rec_targEv),length(rec_lureEv),length(src_targEv),length(src_lureEv),... % raw numbers of targs, lures
          source_dp,source_c,... % accuracy
          length(rec_h),length(rec_cr),length(src_h),length(src_m),length(src_cr),length(src_fa),... % raw numbers for accuracy
          src_hr,src_mr,src_crr,src_far,... % accuracy rates
          length(rec_h_srcCor),length(rec_h_srcInc),... % raw numbers for recognition hits
          rec_h_srcCor_r,rec_h_srcInc_r,... % accuracy rates for recognition hits
          rec_h_rt,rec_cr_rt,... % rec reaction times
          src_h_rt,src_m_rt,src_cr_rt,src_fa_rt,... % src reaction times
          rec_h_srcCor_rt,rec_h_srcInc_rt,... % recognition hits reaction times
          pool_accuracy,... % more accuracy
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

%matlabpool close