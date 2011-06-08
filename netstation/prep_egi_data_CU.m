function prep_egi_data_CU(subject,session,eventfiles,badchan,ms_field,captype)
%PREP_EGI_DATA - Prepare data collected with EGI system and pyepl.
%
% FUNCTION:
%   prep_egi_data(subject,session,eventfiles,badchan,ms_field)
%
% INPUT ARGS:
%   subject = 'FR062';
%   session = '.';  % set that when called from session dir
%   eventfiles = {'events.mat','revents.mat'};
%   badchan = [1 4 8];
%   ms_field = 'mstime';
%   captype = 'GSN200'; % or 'HCGSN' for HydroCel caps
%
% This function expects that the events structures are in the
% session_? directory and the raw eeg data are in the session_?/eeg
% directory.  There may be multiple EEG files.

if ~exist('captype','var')
	captype = 'GSN200';
end
if ~exist('ms_field','var')
  ms_field = 'mstime';
end
if ~exist('badchan','var')
  badchan = [];
end

if isstr(session)
  sessdir = session;
else
  sessdir = fullfile('../data',subject,['session_' num2str(session-1)]);
end

% set up directory structure
eegdir = fullfile(sessdir,'eeg');
norerefdir = fullfile(eegdir,'eeg.noreref');
rerefdir = fullfile(eegdir,'eeg.reref');

% process the eeg.eeglog file
fixEEGLog(fullfile(sessdir,'eeg.eeglog'), fullfile(sessdir,'eeg.eeglog.up'));

% extract the raw data
% find file to extract
rawfile = dir(fullfile(eegdir,'*.raw'));
if length(rawfile) == 0
  error('No raw eeg file found.');
  return
end

for i=1:length(rawfile)
  % split this file into channels
  basename = egi_split_CU(fullfile(eegdir,rawfile(i).name),subject,norerefdir);

  % get the filename of one of the channels
  chan_file{i} = fullfile(norerefdir, [basename '.001']);
  
  % get the EEG sync pulse file
  %d = dir(fullfile(norerefdir, [basename '.DIN*']));
  d = dir(fullfile(norerefdir, [basename '.D255']));
  %d = dir(fullfile(norerefdir, [basename '.D136']));
  eeg_file{i} = fullfile(norerefdir, d(1).name);
end

% get the rate info (assume same for all eeg files)
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(norerefdir);
beh_file = {'eeg.eeglog.up'};

% align data
isfrei = 0;
isegi = 1;
moreAccurateAlign = 1;
runAlign_CU(samplerate, beh_file, eeg_file, chan_file, eventfiles, ms_field, isfrei, isegi, moreAccurateAlign);

% % set variables based on which cap was used
% [eog,perif] = getcapinfo(captype);
% 
% % add artifact info to the events structs
% for i = 1:length(eventfiles)
%   addArtifacts(eventfiles{i},eog,100,0);
% end
% 
% % deal with bad channels
% % see if badchan is a file to open, or to save
% if ~isempty(badchan)
%   if isstr(badchan)
%     % open from file
%     badchanfile = badchan;
%     badchan = textread(badchanfile)';
%   else
%     % save to file
%     badchanfile = fullfile(eegdir,[basename '.bad_chan']);
%     chan2file(badchanfile,badchan);
%   end
% end
% % set to bad plus perif chans
% badchan = unique([badchan perif]);
% 
% % rereference the data
% % read in and combine all events structures
% fileroots = {};
% channels = 1:129;
% for i = 1:length(eventfiles)
%   ev = loadEvents(eventfiles{i});
%   newroots = unique(getStructField(ev,'eegfile','~strcmp(eegfile,'''')'));
%   fileroots(length(fileroots)+1:length(fileroots)+length(newroots)) = newroots;
% end
% fileroots = unique(fileroots);
% reref(fileroots,channels,rerefdir,{1:129,setdiff(channels,badchan)});
% 
% % redo the events to point to the rereferenced data
% for i = 1:length(eventfiles)
%   % save backup of old eventfile
%   copyfile(eventfiles{i},[eventfiles{i} '.noreref'],'f');
%   
%   % load in the new file, replacing the eegfile
%   ev = loadEvents(eventfiles{i},{rerefdir});
%   
%   % save out again
%   saveEvents(ev,eventfiles{i});
% end




function fixEEGLog(infile,outfile)
%FIXEEGLOG - Fix pyepl EEG Logfile leaving only UP pulses.
%
%
% FUNCTION:
%   fixEEGLog(infile,outfile)
%
% INPUT ARGS:
%   infile = 'eeg.eeglog';
%   outfile = 'eeg.eeglog.up';
%

% read in the logfile
[mstime, maxoffset, type] = textread(infile,'%s%n%s%*[^\n]','emptyvalue',NaN);

% write out new file
fid = fopen(outfile,'w');
for i = 1:length(type)
  if strcmp(type{i},'UP')
    % save it to file
    fprintf(fid,'%s\t%d\t%s\n',mstime{i},maxoffset(i),type{i});
  end
end
fclose(fid);


function chan2file(filename,chans)
%CHAN2FILE - Write channels to a file
%
%
% FUNCTION:
%   chan2file(filename,chans)
%
% INPUT ARGS:
%   filename
%   chans
%
%


fid = fopen(filename,'w');
fprintf(fid,'%d\n',chans);
fclose(fid);


