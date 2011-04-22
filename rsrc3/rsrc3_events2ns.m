function rsrc3_events2ns(dataroot,subject,session)
% rsrc3_events2ns(dataroot,subject,session)
%
% Create a NetStation .evt file so we can import events from a
% PyEPL-based experiment into NetStation using the Markup File tool
%
% Note: when writing to the .evt file, uses a carriage return (\r)
% instead of a linefeed return (\n) because NetStation expects a CR
%

% dataroot = '/Volumes/curranlab/Data/RSRC3/Behavioral/Sessions';
% subject = 'RSRC3001';
% session = 1;

eventsDir = fullfile(dataroot,subject,session,'events');

% load in the events
fprintf('Loading %s...\n',fullfile(eventsDir,'test_events.mat'));
events = loadEvents(fullfile(eventsDir,'test_events.mat'));

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
  eegdir = fullfile(strrep(dataroot,'behavioral','eeg'),subject,session,'eeg');
  rawfile = dir(fullfile(eegdir,'*.raw'));
  if ~isempty(rawfile)
    rawfile_name = rawfile.name;
    % take off the .raw extension
    ns_filename = rawfile_name(1:end-4);
    fprintf('Using info from %s\n',fullfile(eegdir,rawfile_name));
  else
    fprintf('Cannot find .raw file! Using subject name instead. You must rename .evt file AND put the proper name in the file before importing events!\n');
    rawfile_name = [subject,'_YYYYMMDD_HHMM_rename_me'];
%     if isfield(events,'eegfile')
%       %rawfile_name_file = events(1).eegfile(end-19:end);
%       %rawfile_date = strrep(rawfile_name_file(1:end-5),[subject,'_'],'');
%       %rawfile_name = [subject,'_',datestr(rawfile_date,'yyyymmdd'),rawfile_name_file(end-4:end)];
%       if strfind(events(1).eegfile,'eeg.reref')
%         rawfile_name_file = strrep(events(1).eegfile,[fullfile(eegdir,'eeg.reref'),filesep],'');
%         rawfile_date = strrep(rawfile_name_file(1:end-5),[subject,'_'],'');
%         rawfile_name = [subject,'_',datestr(rawfile_date,'yyyymmdd'),rawfile_name_file(end-4:end)];
%       elseif strfind(events(1).eegfile,'eeg.noreref')
%         rawfile_namefile = strrep(events(1).eegfile,[fullfile(eegdir,'eeg.noreref'),filesep],'');
%         rawfile_date = strrep(rawfile_name_file(1:end-5),[subject,'_'],'');
%         rawfile_name = [subject,'_',datestr(rawfile_date,'yyyymmdd'),rawfile_name_file(end-4:end)];
%       else
%         fprintf('Cannot find .raw file or eegfile field with eeg.reref or eeg.noreref directory! Using subject name instead. You must rename .evt file AND put the proper name in the file before importing events!\n');
%         rawfile_name = [subject,'_YYYYMMDD_HHMM_rename_me.raw'];
%       end
%     else
%       fprintf('Cannot find .raw file or any eegfile field! Using subject name instead. You must rename .evt file AND put the proper name in the file before importing events!\n');
%       rawfile_name = [subject,'_YYYYMMDD_HHMM_rename_me.raw'];
%     end
    ns_filename = rawfile_name;
  end
else
  rawfile_name = rawfile.name;
  % take off the .raw extension
  ns_filename = rawfile_name(1:end-4);
  fprintf('Using info from %s\n',fullfile(eegdir,rawfile_name));
end
% get the samplerate
% [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fullfile(eegdir,'eeg.reref'));
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fullfile(eegdir,'eeg.noreref'));

% make sure it's 250, otherwise throw an error because we're not
% reading params.txt correctly
if samplerate ~= 250
  error('NetStation file samplerate is not 250 Hz. Something is wrong.');
end

% set the name of the file that we will import back into netstation
netstation_events = fullfile(eventsDir,sprintf('%s.evt',ns_filename));
overwrite_nsevents = 1;
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
fprintf(outfile,'%s\r',ns_filename);
fprintf(outfile,'Time Mode: Relative Time\r');
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
  typeStr = 'type';
  itemStr = 'item';
  xcrdStr = 'xcrd';
  ycrdStr = 'ycrd';
  srpsStr = 'srps';
  % recog response
  rtrgStr = 'rtrg';
  rrspStr = 'rrsp';
  rcnfStr = 'rcnf';
  rcorStr = 'rcor';
  rrtmStr = 'rrtm';
  rslwStr = 'rslw';
  % source response
  strgStr = 'strg';
  srspStr = 'srsp';
  %scnfStr = 'scnf';
  scorStr = 'scor';
  srtmStr = 'srtm';
  sslwStr = 'sslw';
  % eeg_toolbox artifact info
  %artiStr = 'arti';
  
  % session info
  subject = num2str(str2double(events(i).subject(end-2:end)));
  session = num2str(events(i).session);
  % list num
  list = num2str(events(i).list);
  % trial num
  trial = num2str(events(i).trial);
  % item info
  type = events(i).type;
  item = events(i).item;
  xcoord = num2str(events(i).xcoord);
  ycoord = num2str(events(i).ycoord);
  serpos = num2str(events(i).serialpos);
  % recog response
  rec_isTarg = num2str(events(i).rec_isTarg);
  rec_resp = num2str(events(i).rec_resp);
  rec_conf = num2str(events(i).rec_conf);
  rec_correct = num2str(events(i).rec_correct);
  rec_rt = num2str(events(i).rec_rt);
  rec_tooslow = num2str(events(i).rec_tooslow);
  % put in -1 if this field doesn't apply
  if isempty(rec_isTarg)
    rec_isTarg = '-1';
  end
  if isempty(rec_resp)
    rec_resp = '-1';
  end
  if isempty(rec_conf)
    rec_conf = '-1';
  end
  if isempty(rec_correct)
    rec_correct = '-1';
  end
  if isempty(rec_rt)
    rec_rt = '-1';
  end
  if isempty(rec_tooslow)
    rec_tooslow = '-1';
  end
  % source response
  src_isTarg = num2str(events(i).src_isTarg);
  src_correct = num2str(events(i).src_correct);
  src_resp = num2str(events(i).src_resp);
  src_rt = num2str(events(i).src_rt);
  src_tooslow = num2str(events(i).src_tooslow);
  % put in -1 if this field doesn't apply
  if isempty(src_isTarg)
    src_isTarg = '-1';
  end
  if isempty(src_resp)
    src_resp = '-1';
  end
  if isempty(src_correct)
    src_correct = '-1';
  end
  if isempty(src_rt)
    src_rt = '-1';
  end
  if isempty(src_tooslow)
    src_tooslow = '-1';
  end
  %artifactMS = num2str(events(i).artifactMS);
  
  % find each type of event and properly create the code and label;
  % NOTE: code is limited to 4 characters
  switch events(i).type
   case 'STUDY_BUFFER'
    code = 'sbuf';
    
    % only see buffers once, so there is no info to include in the labal
    label = 'study_buff';
    
   case 'STUDY_TARGET'
    code = 'strg';
    label = 'study_targ';
    
   case 'LURE_PRES'
    code = 'tlur';
    label = 'test_lure';
    
   case 'TARG_PRES'
    code = 'ttrg';
    label = 'test_targ';
    
   case 'RK_RESP'
    code = 'rrec';
    label = 'resp_rec';
    
   case 'SOURCE_RESP'
    code = 'rsrc';
    label = 'resp_src';
    
  end
  
  eventInfo = {code,label,stimEvStr,trackStr,timeStr,durStr};
  
  cellLabels = {subStr,sesStr,listStr,trialStr,typeStr,itemStr,xcrdStr,ycrdStr,srpsStr,rtrgStr,rrspStr,rcnfStr,rcorStr,rrtmStr,rslwStr,strgStr,srspStr,scorStr,srtmStr,sslwStr};
  cellData = {subject,session,list,trial,type,item,xcoord,ycoord,serpos,rec_isTarg,rec_resp,rec_conf,rec_correct,rec_rt,rec_tooslow,src_isTarg,src_resp,src_correct,src_rt,src_tooslow};
  
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
