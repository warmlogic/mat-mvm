function slork_prepData_events(subjects,prep_eeg)
% slork_prepData_events(subjects,prep_eeg)
%
% Purpose
%   Create behavioral events; if prep_eeg == 1: export Net Station events
%
% Inputs
%   subjects: a cell of subject numbers
%
%   prep_eeg: boolean, whether or not to prepare the EEG data
%
% Outputs
%   Events (struct and NetStation) will be saved in:
%     /Users/username/data/SLORK/subject/session/events
%
% Assumptions
%   Each subject only ran one session (session_0)
%
%   The behavioral data is located in:
%     /Users/username/data/SLORK/subject/session
%

serverDir = '/Volumes/curranlab/Data/SLORK/eeg/behavioral';
serverLocalDir = '/Volumes/RAID/curranlab/Data/SLORK/eeg/behavioral';
if exist(serverDir,'dir')
  dataroot = serverDir;
elseif exist(serverLocalDir,'dir')
  dataroot = serverLocalDir;
else
  %uname = getenv('USER');
  uroot = getenv('HOME');
  dataroot = fullfile(uroot,'data/SLORK/eeg');
end
saveDir = dataroot;

if nargin == 0
  subjects = {...
      'SLORK002'... % 1
      'SLORK003'...
      'SLORK004'...
      'SLORK005'...
      'SLORK006'...
      'SLORK007'...
      'SLORK008'...
      'SLORK009'...
      'SLORK010'...
      'SLORK011'...
      'SLORK012'...
      'SLORK013'...
      'SLORK014'...
      'SLORK015'...
      'SLORK016'...
      'SLORK017'...
      'SLORK018'...
      'SLORK019'...
      'SLORK020'...
      'SLORK022'... % 21 had problems
      'SLORK023'...
      'SLORK024'...
      'SLORK025'...
      'SLORK026'...
      'SLORK027'... % 28 had problems
      'SLORK029'...
      'SLORK030'...
      'SLORK031'...
      'SLORK032'...
      'SLORK033'... % master's analyses included up to 33
      'SLORK034'...
      'SLORK035'...
      'SLORK036'...
      'SLORK038'...
      'SLORK039'...
      'SLORK040'...
      'SLORK041'...
      'SLORK042'...
      'SLORK043'...
             };
  
  prep_eeg = 1;
end

% behavioral pilot
% %    'SLORK001'... % behavioral pilot didn't do last trial

% eeg
% %    'SLORK001'... % eeg crashed in the middle
% %    'SLORK021'... % no eeg recorded
% %    'SLORK028'... % no eeg recorded

% 37 does not exist?
% 41 was rerun using a different person because didn't hit record

% if prep_eeg == 1
%   nsroot = fullfile(dataroot,'ns_files');
% end

sessions = {'session_0'};

% dr = dir(dataroot);
% s = filterStruct(dr,{'strncmp(name,''SLORK'',4)'});
% subjects = getStructField(s,'name');

% if prep_eeg == 1
%   % get bad channel info
%   [subBC,sesBC,badChan] = textread(fullfile(dataroot,'slork_badChan.txt'),'%s%s%s','delimiter','\t');
%   for i = 1:length(badChan)
%     badChan{i} = str2num(badChan{i});
%   end
%   infoStruct = struct('subject',subBC,'session',sesBC,'badChan',badChan);
% end

% do_accuracy = 0;
% 
% if do_accuracy
%   headerlines = {'Sub','Hit','Hit','Hit','Hit','Hit','Miss','Miss','Miss','Miss','FA','FA','FA','FA','FA','CR','CR','CR','CR','Hit SC','Hit SC','Hit SC','Hit SC','Hit SC','Hit SI','Hit SI','Hit SI','Hit SI','Hit SI','Discrimination','Response Bias';...
%                  'Sub','Overall','rt','RS','RO','Fam','Overall','rt','Sure','Maybe','Overall','rt','RS','RO','Fam','Overall','rt','Sure','Maybe','Overall','rt','RS','RO','Fam','Overall','rt','RS','RO','Fam','Pr','Br'};
%   resp_data = nan(length(subjects),size(headerlines,2));
% end

%matlabpool open

for sub = 1:length(subjects)
  fprintf('Getting data for %s...',subjects{sub});
  for ses = 1:length(sessions)
    fprintf('%s...\n',sessions{ses});
    
    if prep_eeg == 1
      % find the bad channels for this subject and session
%       sesStruct = filterStruct(infoStruct,'ismember(subject,varargin{1}) & ismember(session,varargin{2})',subjects{sub},sessions{ses});
%       subSesBadChan = sesStruct.badChan;
%       if isempty(sesStruct)
%         error('no subject listing found for this session');
%       end
      subSesBadChan = [];
    end
    
    % set the subject events directory
    eventsOutdir_sub = fullfile(saveDir,subjects{sub},sessions{ses},'events');
    if ~exist(eventsOutdir_sub,'dir')
      mkdir(eventsOutdir_sub);
    end
    
    % set the subject events file
    eventsOutfile_sub = fullfile(eventsOutdir_sub,'events.mat');
    if exist(eventsOutfile_sub,'file')
      %fprintf('%s already exists! Skipping this subject!\n',eventsOutfile_sub);
      %continue
      fprintf('%s already exists! Next subject...\n',eventsOutfile_sub);
      continue
    else
      %if ~lockFile(eventsOutfile_sub)
      fprintf('Creating events...\n');
      % create the events
      events = slork_createEvents(dataroot,subjects{sub},sessions{ses});
      fprintf('Saving %s...\n',eventsOutfile_sub);
      % save each subject's events
      saveEvents(events,eventsOutfile_sub);
      
      % release the lockFile
      %releaseFile(eventsOutfile_sub);
    end
    
    %% old accuracy section
%     if do_accuracy
%       %events = filterStruct(events,'ismember(testType,varargin{1})',{'task'});
%       %events = filterStruct(events,'ismember(testType,varargin{1})',{'side'});
%       
%       subNum = str2num(subjects{sub}(end-2:end));
%       % collect data in sub_data vector
%       sub_data = subNum;
%       
%       rec_targEv = filterStruct(events,'rec_isTarg == 1 & ismember(type,varargin{1})',{'TARG_PRES'});
%       rec_lureEv = filterStruct(events,'rec_isTarg == 0 & ismember(type,varargin{1})',{'LURE_PRES'});
%       
%       src_correct_targ = getStructField(rec_targEv,'src_correct'); % h and m
%       src_correct_lure = getStructField(rec_lureEv,'src_correct'); % cr and fa
%       src_rt_targ = getStructField(rec_targEv,'src_rt');
%       src_rt_lure = getStructField(rec_lureEv,'src_rt');
%       
%       rkn_resp_targ = getStructField(rec_targEv,'rkn_resp');
%       rkn_resp_lure = getStructField(rec_lureEv,'rkn_resp');
%       rkn_rt_targ = getStructField(rec_targEv,'rkn_rt');
%       rkn_rt_lure = getStructField(rec_lureEv,'rkn_rt');
%       
%       fprintf('\n\nSource accuracy info:\n');
%       
%       % get the counts for hit source correct and hit source incorrect
%       hit_sc = sum(src_correct_targ == 1);
%       hit_si = sum(src_correct_targ == 0 & (strcmp(rkn_resp_targ,'REMEMBER_SOURCE') | strcmp(rkn_resp_targ,'REMEMBER_OTHER') | strcmp(rkn_resp_targ,'KNOW')));
%       
%       % Hit
%       hit = hit_sc + hit_si;
%       rt_h = mean(src_rt_targ(src_correct_targ == 1 | (src_correct_targ == 0 & (strcmp(rkn_resp_targ,'REMEMBER_SOURCE') | strcmp(rkn_resp_targ,'REMEMBER_OTHER') | strcmp(rkn_resp_targ,'KNOW')))));
%       rs_h = sum(src_correct_targ == 1 & strcmp(rkn_resp_targ,'REMEMBER_SOURCE'));
%       ro_h = sum(src_correct_targ == 1 & strcmp(rkn_resp_targ,'REMEMBER_OTHER'));
%       fam_h = sum(src_correct_targ == 1 & strcmp(rkn_resp_targ,'KNOW'));
%       fprintf('Hits (%d): %.4f\n',hit,hit / length(rec_targEv));
%       fprintf('\tRS (%d): %.4f\n',rs_h,rs_h / sum(src_correct_targ == 1));
%       fprintf('\tRO (%d): %.4f\n',ro_h,ro_h / sum(src_correct_targ == 1));
%       fprintf('\t F (%d): %.4f\n',fam_h,fam_h / sum(src_correct_targ == 1));
%       
%       sub_data = [sub_data,hit / length(rec_targEv),rt_h,rs_h,ro_h,fam_h];
%       
%       % Miss
%       miss = sum(src_correct_targ == 0 & (strcmp(rkn_resp_targ,'SURE') | strcmp(rkn_resp_targ,'MAYBE')));
%       rt_m = mean(src_rt_targ(src_correct_targ == 0 & (strcmp(rkn_resp_targ,'SURE') | strcmp(rkn_resp_targ,'MAYBE'))));
%       sure_m = sum(src_correct_targ == 0 & strcmp(rkn_resp_targ,'SURE'));
%       maybe_m = sum(src_correct_targ == 0 & strcmp(rkn_resp_targ,'MAYBE'));
%       fprintf('Misses (%d): %.4f\n',miss,miss / length(rec_targEv));
%       fprintf('\t S (%d): %.4f\n',sure_m,sure_m / sum(src_correct_targ == 0));
%       fprintf('\t M (%d): %.4f\n',maybe_m,maybe_m / sum(src_correct_targ == 0));
%       
%       sub_data = [sub_data,miss / length(rec_targEv),rt_m,sure_m,maybe_m];
%       
%       % False Alarm
%       fa = sum(src_correct_lure == 0);
%       rt_fa = mean(src_rt_lure(src_correct_lure == 0));
%       rs_fa = sum(src_correct_lure == 0 & strcmp(rkn_resp_lure,'REMEMBER_SOURCE'));
%       ro_fa = sum(src_correct_lure == 0 & strcmp(rkn_resp_lure,'REMEMBER_OTHER'));
%       fam_fa = sum(src_correct_lure == 0 & strcmp(rkn_resp_lure,'KNOW'));
%       fprintf('False alarms (%d): %.4f\n',fa,fa / length(rec_lureEv));
%       fprintf('\tRS (%d): %.4f\n',rs_fa,rs_fa / sum(src_correct_lure == 0));
%       fprintf('\tRO (%d): %.4f\n',ro_fa,ro_fa / sum(src_correct_lure == 0));
%       fprintf('\t F (%d): %.4f\n',fam_fa,fam_fa / sum(src_correct_lure == 0));
%       
%       sub_data = [sub_data,fa / length(rec_lureEv),rt_fa,rs_fa,ro_fa,fam_fa];
%       
%       % Correct Rejection
%       cr = sum(src_correct_lure == 1);
%       rt_cr = mean(src_rt_lure(src_correct_lure == 1));
%       sure_cr = sum(src_correct_lure == 1 & strcmp(rkn_resp_lure,'SURE'));
%       maybe_cr = sum(src_correct_lure == 1 & strcmp(rkn_resp_lure,'MAYBE'));
%       fprintf('Correct rejections (%d): %.4f\n',cr,cr / length(rec_lureEv));
%       fprintf('\t S (%d): %.4f\n',sure_cr,sure_cr / sum(src_correct_lure == 1));
%       fprintf('\t M (%d): %.4f\n',maybe_cr,maybe_cr / sum(src_correct_lure == 1));
%       
%       sub_data = [sub_data,cr / length(rec_lureEv),rt_cr,sure_cr,maybe_cr];
%       
%       % Hit, Source Correct
%       rt_hsi = mean(src_rt_targ(src_correct_targ == 1));
%       rs_hsc = sum(src_correct_targ == 1 & strcmp(rkn_resp_targ,'REMEMBER_SOURCE'));
%       ro_hsc = sum(src_correct_targ == 1 & strcmp(rkn_resp_targ,'REMEMBER_OTHER'));
%       fam_hsc = sum(src_correct_targ == 1 & strcmp(rkn_resp_targ,'KNOW'));
%       fprintf('Hits, Source Correct (%d): %.4f\n',hit_sc,hit_sc / (hit_sc + hit_si));
%       fprintf('\tRS (%d): %.4f\n',rs_hsc,rs_hsc / sum(src_correct_targ == 1));
%       fprintf('\tRO (%d): %.4f\n',ro_hsc,ro_hsc / sum(src_correct_targ == 1));
%       fprintf('\t F (%d): %.4f\n',fam_hsc,fam_hsc / sum(src_correct_targ == 1));
%       
%       sub_data = [sub_data,hit_sc / (hit_sc + hit_si),rt_hsi,rs_hsc,ro_hsc,fam_hsc];
%       
%       % Hit, Source Incorrect
%       rt_hsi = mean(src_rt_targ(src_correct_targ == 0 & (strcmp(rkn_resp_targ,'REMEMBER_SOURCE') | strcmp(rkn_resp_targ,'REMEMBER_OTHER') | strcmp(rkn_resp_targ,'KNOW'))));
%       rs_hsi = sum(src_correct_targ == 0 & strcmp(rkn_resp_targ,'REMEMBER_SOURCE'));
%       ro_hsi = sum(src_correct_targ == 0 & strcmp(rkn_resp_targ,'REMEMBER_OTHER'));
%       fam_hsi = sum(src_correct_targ == 0 & strcmp(rkn_resp_targ,'KNOW'));
%       fprintf('Hits, Source Incorrect (%d): %.4f\n',hit_si,hit_si / (hit_sc + hit_si));
%       fprintf('\tRS (%d): %.4f\n',rs_hsi,rs_hsi / sum(src_correct_targ == 0));
%       fprintf('\tRO (%d): %.4f\n',ro_hsi,ro_hsi / sum(src_correct_targ == 0));
%       fprintf('\t F (%d): %.4f\n',fam_hsi,fam_hsi / sum(src_correct_targ == 0));
%       
%       sub_data = [sub_data,hit_si / (hit_sc + hit_si),rt_hsi,rs_hsi,ro_hsi,fam_hsi];
%       
%       % From Mecklinger et al. (2007) and Corwin (1994)
%       %
%       % discrimination index
%       Pr = (hit / length(rec_targEv)) - (fa / length(rec_lureEv));
%       % response bias index
%       Br = (fa / length(rec_lureEv)) / (1 - Pr);
%       
%       sub_data = [sub_data,Pr,Br];
%       
%       resp_data(sub,:) = sub_data;
%       
%       %     src_targEv = filterStruct(rec_targEv,'rec_correct == 1');
%       %     src_lureEv = filterStruct(rec_lureEv,'rec_correct == 0');
%       %     fprintf('\n\nSource accuracy info:\n');
%       %     src_correct_targ = getStructField(src_targEv,'src_correct');
%       %     %src_correct_lure = getStructField(src_lureEv,'src_correct');
%       %     fprintf('Hits: %.4f\n',sum(src_correct_targ == 1) / length(src_correct_targ(src_correct_targ ~= -1)));
%       %     fprintf('Misses: %.4f\n',sum(src_correct_targ == 0) / length(src_correct_targ(src_correct_targ ~= -1)));
%       %     %fprintf('Correct rejections: %.4f\n',sum(src_correct_lure == 1) / length(src_correct_lure(src_correct_lure ~= -1)));
%       %     %fprintf('False alarms: %.4f\n',sum(src_correct_lure == 0) / length(src_correct_lure(src_correct_lure ~= -1)));
%     end
    
    %% prep the EEG data
    if prep_eeg == 1
      fprintf('Prepping EEG data...\n');
      % get this subject's session dir
      sesDir = fullfile(dataroot,subjects{sub},sessions{ses});
%       % this is where the RAW file is currently stored
%       oldEEGDir = fullfile(sesDir,'eeg');
%       % this is where we want to move the EEG data to
%       newEEGRoot = fullfile(dataroot,subjects{sub},sessions{ses});
%       newEEGDir = fullfile(newEEGRoot,'eeg');
%       if exist(newEEGDir,'dir')
%         sprintf('%s already exists, skipping!\n',newEEGDir);
%         continue
%       else
%         mkdir(newEEGRoot);
%       end
      
      subEegDir = fullfile(sesDir,'eeg','eeg.noreref');
      pfile = dir(fullfile(subEegDir,[subjects{sub},'*params.txt']));
      
      if ~exist(fullfile(subEegDir,pfile.name),'file')
        curDir = pwd;
        % cd to the session directory since prep_egi_data needs to be there
        cd(sesDir);
        % align the behavioral and EEG data
        prep_egi_data_CU(subjects{sub},sesDir,{fullfile(sesDir,'events/events.mat')},subSesBadChan,'mstime','HCGSN');
        % go back to the previous working directory
        cd(curDir);
      end
      
      % export the events for netstation; saves to the session's events dir
      slork_events2ns(dataroot,subjects{sub},sessions{ses});
      
%       % change the location of EEG files
%       if ~exist(newEEGDir,'dir')
%         % do the actual move
%         unix(sprintf('mv %s %s',oldEEGDir,newEEGRoot));
%         % resave the events
%         events = loadEvents(fullfile(sesDir,'events/events.mat'),{oldEEGDir,newEEGDir});
%         saveEvents(events,fullfile(sesDir,'events/events.mat'));
%       else
%         error('Error: %s already exists, not moving any directories!\n',newEEGDir);
%       end

    end % prep_eeg
  end % ses
  fprintf('Done.\n');
end % sub

%matlabpool close

% if do_accuracy
%   if ~exist(fullfile(dataroot,'resp_distrib_head.csv'),'file')
%     fid = fopen(fullfile(dataroot,'resp_distrib_head.csv'),'w');
%     for i = 1:size(headerlines,1)
%       fprintf(fid,'%s',headerlines{i,1});
%       for j = 2:size(headerlines,2)
%         fprintf(fid,',%s',headerlines{i,j});
%       end
%       fprintf(fid,'\n');
%     end
%     fclose(fid);
%   end
%   %xlswrite(fullfile(dataroot,'resp_distrib_head.xls'),headerlines);
%   
%   csvwrite(fullfile(dataroot,'resp_distrib.csv'),resp_data);
% end
