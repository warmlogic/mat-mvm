function frce_ns_addTRSP2evt(dataroot)

% This function reads in an .evt file with multiple event tags denoting
% trial-specific (TRSP) information, backs up the original .evt, and writes
% out a new .evt file (with the same name as the original) with TRSP
% information collapsed into a single event tag. The new .evt file is
% useful for doing segmentation in NS.
%
% NB: The original .evt file needs to be exported using relative time.
%
% TODO: don't separate into targets & buffers, let separation happen in NS
%
%
%
% USAGE STEPS:
%
% This function is specific to the FRCE experiment. To properly process the
% original continuous EEG files with multiple tags per stimulus, you must
% 1: markup the file with subsequent memory performance information. These
%    .evt files are in /curranlab/Data/FRCE/EEG/evt files/.
% 2: export .evt files again
% 3: run this function
% 4: markup the continuous files again with the new .evt files.
%
% NB: Keep track of the different .evt files. Don't mix them up!
%
% FUNCTION SPECIFICS:
%
% Event codes in the new .evt file are changed to STIM. The folliwng old
% TRSP tag data will be transfered to the key codes mentioned below:
%
% cue type
% Old tag: C001 (ear/auditory) or C002 = eye/visual
% New key code: CUET (ear or eye)
%
% start of prestimulus period
% Old tag: PRES
% New key code: Don't need to transfer this info
%
% length of the prestimulus period
% Old tag: 1250, 1500, or 1750
% New key code: PSTM (1250, 1500, or 1750)
%
% stimulus number
% Old tag: P001-P320
% New key code: STMN
%
% subsequent memory perofrmance
% Old tag: reca or forg
% New key code: SMEM (0 or 1)
%
% stimulus modality (confounded with cue type)
% Old tag: M001 (auditory) or M002 (visual)
% New key code: SMOD (auditory or visual)
%
% block number
% Old tag: B001-B020
% New key code: BLOC (1-20)
%
% serial position
% Old tag: S001-S016
% New key code: SPOS (1-16)
%
% Also, adds a stimulus type key code, where it's a buffer if serial
% position is 1, 2, 15, or 16, or it's a target if serial position is 3-14.
% New key code: TYPE (STUDY_TARGET or STUDY_BUFFER)


if nargin < 1
  % where the evt files are located
  dataroot = pwd;
end

%% Read in the old evt and write out the new evt

% Find the evt files
fprintf('Searching for .evt files in %s...',dataroot);
evt = dir(fullfile(dataroot,'*.evt'));
if ~isempty(evt)
  fprintf('found %d.\n',length(evt));
else
  error('No .evt files found!');
end

numHeaderlines = 5;

for i = 1:length(evt)
  fprintf('Processing %s...',evt(i).name);
  
  %% Read in the evt file
  [path,name,ext] = fileparts(evt(i).name);
  
  % figure out how many columns there are
  fid = fopen(fullfile(dataroot,evt(i).name),'r');
  if fid == -1
    error('Could not open the file. Make sure you do not have it open in another application.');
  end
  % get the file name line
  fgetl(fid);
  % get the time mode line
  fgetl(fid);
  % get the header line
  headerline = fgetl(fid);
  % headers
  if strcmp(headerline(end),sprintf('\t'))
    headerline = headerline(1:end-1);
  end
  hdr = regexp(headerline,'\t','split');
  cols.code = find(strcmp(hdr,'Code'));
  %cols.label = find(strcmp(hdr,'Label'));
  %cols.type = find(strcmp(hdr,'Type'));
  %cols.track = find(strcmp(hdr,'Track'));
  %cols.onset = find(strcmp(hdr,'Onset'));
  %cols.duration = find(strcmp(hdr,'Duration'));
  
  % get these lines because we want to put them back in the file
  sessline = fgetl(fid);
  cellline = fgetl(fid);
  
  % get the first data line so we can count the number of columns
  dataline = fgetl(fid);
  
  % close the file
  fclose(fid);
  
  % if the last character is a tab, remove it
  if strcmp(sessline(end),sprintf('\t'))
    sessline = sessline(1:end-1);
  end
  % if the last character is a tab, remove it
  if strcmp(cellline(end),sprintf('\t'))
    cellline = cellline(1:end-1);
  end
  % if the last character is a tab, remove it
  if strcmp(dataline(end),sprintf('\t'))
    dataline = dataline(1:end-1);
  end
  % since it's tab delimited, split it and count the length
  numCols = length(regexp(dataline,'\t','split'));
  
  fid = fopen(fullfile(dataroot,evt(i).name),'r');
  data = textscan(fid,repmat('%s',1,numCols),'Delimiter','\t');
  fclose(fid);
  
  % rename the old evt file
  feval(str2func('unix'),sprintf('mv %s %s',strrep(fullfile(dataroot,evt(i).name),' ','\ '),strrep(fullfile(dataroot,[name,'_backup',ext]),' ','\ ')));
  
  % start the new file to write out
  new_evt = fullfile(dataroot,evt(i).name);
  outfile = fopen(new_evt,'wt');
  
  % add the first five lines
  fprintf(outfile,'%s\t\t\t\t\t\r',name);
  fprintf(outfile,'Time Mode: Relative Time\t\t\t\t\t\r');
  fprintf(outfile,'Code\tLabel\tType\tTrack\tOnset\tDuration\r');
  fprintf(outfile,'%s\r',sessline);
  fprintf(outfile,'%s\r',cellline);
  
  
  %% start looking for the data we need
  % first stim number
  stimNum = 1;
  blockNum = 1;
  serPos = 1;
  
  % initialize the struct to hold the stimulus information
  allStim = struct;
  
  % go through each line and collect the info for each stimulus
  for j = (numHeaderlines + 1):length(data{cols.code})
    
    % set the stimulus number tag so we know what to look for
    stimTag = sprintf('P%03d',stimNum);
    
    if strcmp(data{cols.code}{j},stimTag)
      allStim(stimNum).prevInfo = [];
      for k = 2:numCols
        allStim(stimNum).prevInfo = sprintf('%s\t%s',allStim(stimNum).prevInfo,data{k}{j});
      end
      allStim(stimNum).prevInfo = allStim(stimNum).prevInfo(2:end);
      
      % set the cue type (-3)
      if strcmp(data{cols.code}{j-3},'C001')
        allStim(stimNum).CUET = 'ear';
      elseif strcmp(data{cols.code}{j-3},'C002')
        allStim(stimNum).CUET = 'eye';
      end
      
      % skip the PRES tag (-2)
      
      % set the prestimulus period (-1)
      allStim(stimNum).PSTM = data{cols.code}{j-1};
      
      % set the stimulus number
      thisStimNum = data{cols.code}{j};
      thisStimNum = str2double(thisStimNum(2:end));
      if thisStimNum ~= stimNum
        error('Stimulus number counter (%d) does not match the data (%d)',stimNum,thisStimNum);
      end
      allStim(stimNum).STMN = num2str(thisStimNum);
      %allStim(stimNum).STMN = num2str(stimNum);
      
      % set the subsequent memory performance (+1)
      if strcmp(data{cols.code}{j+1},'reca')
        allStim(stimNum).SMEM = '1';
      elseif strcmp(data{cols.code}{j+1},'forg')
        allStim(stimNum).SMEM = '0';
      else
        % break out of this file if subsequent memory performance is not
        % known
        continue
      end
      
      % set the stimulus modality (+2)
      if strcmp(data{cols.code}{j+2},'M001')
        allStim(stimNum).SMOD = 'auditory';
      elseif strcmp(data{cols.code}{j+2},'M002')
        allStim(stimNum).SMOD = 'visual';
      end
      
      % set the block number (+3)
      thisBlockNum = data{cols.code}{j+3};
      thisBlockNum = str2double(thisBlockNum(2:end));
      if thisBlockNum ~= blockNum
        error('Block number counter (%d) does not match the data (%d)',blockNum,thisBlockNum);
      end
      allStim(stimNum).BLOC = num2str(thisBlockNum);
      %allStim(stimNum).BLOC = num2str(blockNum);
      
      % set the serial position for this block (+4)
      thisSerPos = data{cols.code}{j+4};
      thisSerPos = str2double(thisSerPos(2:end));
      if thisSerPos ~= serPos
        error('Serial position counter (%d) does not match the data (%d)',serPos,thisSerPos);
      end
      allStim(stimNum).SPOS = num2str(thisSerPos);
      %allStim(stimNum).SPOS = num2str(serPos);
      
      % set the stimulus type (buffers are the first and last 2 items)
      %
      % TODO: don't separate into buffers, let separation happen in NS
      if thisSerPos == 1 || thisSerPos == 2 || thisSerPos == 15 || thisSerPos == 16
        allStim(stimNum).TYPE = 'STUDY_BUFFER';
      else
        allStim(stimNum).TYPE = 'STUDY_TARGET';
      end
      
      % now assemble the full line for writing to file
      fn = fieldnames(allStim);
      for k = 1:length(fn)
        if strcmp(fn{k},'prevInfo')
          continue
        else
          allStim(stimNum).prevInfo = sprintf('%s\t%s\t%s',allStim(stimNum).prevInfo,fn{k},allStim(stimNum).(fn{k}));
        end
      end
      
      % write the full line to file with the STIM code
      fprintf(outfile,'%s\t%s\r','STIM',allStim(stimNum).prevInfo);
      
      % advance the counter to look for the next stimulus
      stimNum = stimNum + 1;
      if serPos == 16
        serPos = 1;
        % only advance the block number after going through all positions
        blockNum = blockNum + 1;
      else
        serPos = serPos + 1;
      end
    end % if this is a stimulus presentation (P###)
    
  end % for j (each line of the original .evt)
  
  % close the output file
  fclose(outfile);
  
  fprintf('Done.\n');
end % for i (each .evt)
