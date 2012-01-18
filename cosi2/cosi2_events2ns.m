function cosi2_events2ns(dataroot,subject,session)
% cosi_events2ns(dataroot,subject,session)
%2
% Create a NetStation .evt file so we can import events from a
% PyEPL-based experiment into NetStation using the Markup File tool
%
% NB: when writing to the .evt file, this script uses a carriage
% return (\r) instead of a linefeed return (\n) because NetStation
% expects a CR
%

% set the sample rate that the recording should be
recordingSR = 500;

overwrite_nsevents = 1;

% dataroot = '/Users/matt/data/COSI2/eeg/behavioral';
% subject = 'COSI2001';
% session = 'session_0';

eventsDir = fullfile(dataroot,subject,session,'events');

% load in the events
fprintf('Loading %s...\n',fullfile(eventsDir,'events.mat'));
events = loadEvents(fullfile(eventsDir,'events.mat'));

% predefined strings for the netstation file
stimEvStr = 'Stimulus Event';
trackStr = 'Events';
durStr = '_00:00:00.001';

% get the name of the netstation .raw file
eegdir = fullfile(dataroot,subject,session,'eeg');
rawfile = dir(fullfile(eegdir,'*.raw'));
% make sure that we actually found the netstation .raw file
if isempty(rawfile)
  % if we don't find it there, look in the eeg directory
  eegdir = fullfile(dataroot,subject,session,'eeg');
  rawfile = dir(fullfile(eegdir,'*.raw'));
  if ~isempty(rawfile)
    rawfile_name = rawfile.name;
    % take off the .raw extension
    ns_filename = rawfile_name(1:end-4);
    fprintf('Using info from %s\n',fullfile(eegdir,rawfile_name));
  else
    fprintf('Cannot find .raw file! Using subject name instead. You must rename .evt file AND put the proper name in the file before importing events!\n');
    rawfile_name = [subject,'_YYYYMMDD_HHMM_rename_me'];
    ns_filename = rawfile_name;
  end
else
  rawfile_name = rawfile.name;
  % take off the .raw extension
  ns_filename = rawfile_name(1:end-4);
  fprintf('Using info from %s\n',fullfile(eegdir,rawfile_name));
end
% get the samplerate
[samplerate] = GetRateAndFormat(fullfile(eegdir,'eeg.noreref'));
%[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fullfile(eegdir,'eeg.noreref'));

% make sure it's the right sample rate, otherwise throw an error because
% we're not reading params.txt correctly
if samplerate ~= recordingSR
  error('NetStation file samplerate is not %d Hz. Something is wrong.',recordingSR);
end

% set the name of the file that we will import back into netstation
netstation_events = fullfile(eventsDir,sprintf('%s.evt',ns_filename));
if overwrite_nsevents == 1
  if exist(netstation_events,'file')
    fprintf('%s already exists! Overwriting this file!\n',netstation_events);
  else
    fprintf('Saving NetStation event data to %s\n',netstation_events);
  end
else
  if exist(netstation_events,'file')
    fprintf('%s already exists! Skipping this file and not overwriting!\n',netstation_events);
    return
  else
    fprintf('Saving NetStation event data to %s\n',netstation_events);
  end
end

% open the file that we'll save
outfile = fopen(netstation_events,'wt');

% put in the basic info at the top
fprintf(outfile,'%s\t\t\t\t\t\r',ns_filename);
fprintf(outfile,'Time Mode: Relative Time\t\t\t\t\t\r');
fprintf(outfile,'Code\tLabel\tType\tTrack\tOnset\tDuration\r');

for i = 1:length(events)
  code = '';
  label = '';
  
  % convert the sample offset to ms time
  rel_mstime = (events(i).eegoffset / samplerate) * 1000;
  
  % create the timestamp for this event from ms time
  hrs = sprintf('%02.0f',floor(rel_mstime / (1000 * 60 * 60)));
  mins = sprintf('%02.0f',floor(mod(rel_mstime,(1000 * 60 * 60)) / (1000 * 60)));
  secs = sprintf('%06.3f',mod(mod(rel_mstime,(1000 * 60 * 60)),(1000 * 60)) / 1000);
  timeStr = sprintf('_%s:%s:%s',hrs,mins,secs);
  
  % session info
  subStr = 'subj';
  sesStr = 'sess';
  % list num
  listStr = 'list';
  % trial num
  trialStr = 'tril';
  % item info
  ncolStr = 'ncol'; % number of colors
  typeStr = 'type'; % type of response (STUDY_RESP, SOURCE_RESP, etc.)
  condStr = 'cond'; % condition (color or side)
  itemStr = 'item'; % item name
  srpsStr = 'srps'; % serial position
  scolStr = 'scol'; % study color
  lcolStr = 'lcol'; % lure color
  xcrdStr = 'xcrd'; % x-coordinate during study
  ycrdStr = 'ycrd'; % y-coordinate during study
  % recognition response (implicit)
  rec_targStr = 'rtrg'; % recognition target? 1 or 0
  rcorStr = 'rcor'; % was it correct?
  % source response
  src_targStr = 'strg'; % source target? 1 or 0
  srcrStr = 'srcr'; % the response
  srrtStr = 'srrt'; % reaction time
  scorStr = 'scor'; % was it correct?
  % remember/know (familiar), new response
  rknrStr = 'rknr'; % the response
  rkrtStr = 'rkrt'; % reaction time
  rkncStr = 'rknc'; % was it correct? (always = 1)
  redoStr = 'redo'; % was this a redo trial? (the "corrected" answer, not the accidental one)
  % eeg_toolbox artifact info
  %artiStr = 'arti'; % did eeg_toolbox find an artifact?
  
  % session info
  subject = num2str(str2double(events(i).subject(end-2:end)));
  session = num2str(events(i).session);
  % list num
  list = num2str(events(i).list);
  % trial num
  trial = num2str(events(i).trial);
  % item info
  ncol = num2str(events(i).numColors);
  type = events(i).type;
  cond = events(i).cond;
  item = events(i).item;
  serpos = num2str(events(i).serialpos);
  scol = events(i).study_color;
  lcol = events(i).lure_color;
  xcrd = num2str(events(i).study_loc_x);
  ycrd = num2str(events(i).study_loc_y);
  
  % study and lure colors
  if isempty(scol)
    scol = '-1';
  end
  if isempty(lcol)
    lcol = '-1';
  end
  
  % implicit item recognition response
  rec_isTarg = num2str(events(i).rec_isTarg);
  rec_correct = num2str(events(i).rec_correct);
  % put in -1 if this field doesn't apply
  if isempty(rec_isTarg)
    rec_isTarg = '-1';
  end
  if isempty(rec_correct)
    rec_correct = '-1';
  end
  
  % source response
  src_isTarg = num2str(events(i).src_isTarg);
  % put in -1 if this field doesn't apply
  if isempty(src_isTarg)
    src_isTarg = '-1';
  end
  src_resp = events(i).src_resp;
  src_rt = num2str(events(i).src_rt);
  src_correct = num2str(events(i).src_correct);
  % put in -1 if this field doesn't apply
  if isempty(src_resp)
    src_resp = '-1';
  elseif isnumeric(src_resp)
    src_resp = num2str(src_resp);
  end
  if isempty(src_rt)
    src_rt = '-1';
  end
  if isempty(src_correct)
    src_correct = '-1';
  end
  
  % remember/know/new response
  rkn_resp = events(i).rkn_resp;
  rkn_rt = num2str(events(i).rkn_rt);
  rkn_correct = num2str(events(i).rkn_correct);
  % put in -1 if this field doesn't apply
  if isempty(rkn_resp)
    rkn_resp = '-1';
  elseif isnumeric(rkn_resp)
    rkn_resp = num2str(rkn_resp);
  end
  if isempty(rkn_rt)
    rkn_rt = '-1';
  end
  if isempty(rkn_correct)
    rkn_correct = '-1';
  end
  
  % was this a "corrected" redo trial?
  redo = num2str(events(i).redo);
  
  % % eeg_toolbox artifact information
  % artifactMS = num2str(events(i).artifactMS);
  
  % find each type of event and properly create the code and label;
  % NOTE: code is limited to 4 characters
  switch events(i).type
    
    case 'STUDY_BUFFER'
      code = 'sbuf';
      label = 'study_buff';
      
    case 'STUDY_TARGET'
      code = 'strg';
      label = 'study_targ';
      
    case 'STUDY_RESP'
      code = 'rstu';
      
    case 'TEST_LURE'
      code = 'tlur';
      label = 'test_lure';
      
    case 'TEST_TARGET'
      code = 'ttrg';
      label = 'test_targ';
      
    case 'SOURCE_RESP'
      code = 'rsrc';
      label = 'resp_src';
      
    case 'RK_RESP'
      code = 'rrk_';
      label = 'resp_rk';
      
    case 'NEW_RESP'
      code = 'rnew';
      label = 'resp_new';
      
    otherwise
      error('Event type ''%s'' was not defined.',events(i).type);
  end
  
  eventInfo = {code,label,stimEvStr,trackStr,timeStr,durStr};
  
  %cellLabels = {subStr,sesStr,listStr,trialStr,ncolStr,typeStr,itemStr,poolStr,pcorStr,scolStr,lcolStr,srpsStr,rec_targStr,scndStr,sturStr,strtStr,src_targStr,srcrStr,srrtStr,scorStr,rknrStr,rkrtStr,rcorStr,artiStr};
  %cellData = {subject,session,list,trial,ncol,type,item,pool,pool_correct,scol,lcol,serpos,rec_isTarg,study_cond,study_resp,study_rt,src_isTarg,src_resp,src_rt,src_correct,rkn_resp,rkn_rt,rkn_correct,artifactMS};
  cellLabels = {subStr,sesStr,listStr,trialStr,ncolStr,typeStr,condStr,itemStr,srpsStr,scolStr,lcolStr,xcrdStr,ycrdStr,rec_targStr,rcorStr,src_targStr,srcrStr,srrtStr,scorStr,rknrStr,rkrtStr,rkncStr,redoStr};
  cellData = {subject,session,list,trial,ncol,type,cond,item,serpos,scol,lcol,xcrd,ycrd,rec_isTarg,rec_correct,src_isTarg,src_resp,src_rt,src_correct,rkn_resp,rkn_rt,rkn_correct,redo};
  
  % print the event info
  for c = 1:length(eventInfo)
    fprintf(outfile,'%s\t',eventInfo{c});
  end
  % print the meta data for this event
  for c = 1:length(cellLabels)
    fprintf(outfile,'%s\t%s\t',cellLabels{c},cellData{c});
  end
  
  % put in a carriage return if we're not on the last event
  if i < length(events)
    fprintf(outfile,'\r');
  end
  
end
fclose(outfile);
