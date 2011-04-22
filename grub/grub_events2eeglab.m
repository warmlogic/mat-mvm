function grub_events2eeglab(dataroot,subject,session)
% grub_events2eeglab(dataroot,subject,session)
%
% Create a .txt file so we can import a NetStation RAW file into
% eeglab
%

% dataroot = '/Users/matt/data/GRUB';
% subject = 'GRUB002';
% session = 'session_0';

eventsDir = fullfile(dataroot,subject,session,'events');

% load in the events
fprintf('Loading %s...\n',fullfile(eventsDir,'events.mat'));
events = loadEvents(fullfile(eventsDir,'events.mat'));

testEvents = filterStruct(events,'ismember(type,varargin{1})',{'TARG_PRES','LURE_PRES'});

eeglab_filename = sprintf('%s_%s.txt',subject,session(end));
eeglab_events = fullfile(eventsDir,sprintf('%s',eeglab_filename));

% start the file
outfile = fopen(eeglab_events,'wt');
% write the header
fprintf(outfile,'Latency\tType\tCorrect\n');

for ev = 1:length(testEvents)
  type = [];
  % RHSC
  if testEvents(ev).isTarg == 1 & testEvents(ev).rec_correct == 1 &  testEvents(ev).src_correct == 1
    type = 'RHSC';
  % RHSI
  elseif testEvents(ev).isTarg == 1 & testEvents(ev).rec_correct == 1 &  testEvents(ev).src_correct == 0
    type = 'RHSI';
  % RCR
  elseif testEvents(ev).isTarg == 0 & testEvents(ev).rec_correct == 1
    type = 'RCR';
  else
    continue
  end
  
  fprintf(outfile,'%d\t%s\t%d\n',testEvents(ev).eegoffset,type,sum([testEvents(ev).rec_correct,testEvents(ev).src_correct]));
end

fclose(outfile);
