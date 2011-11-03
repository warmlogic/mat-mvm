function ns_addTRSP2evt(dataroot)

% script to read in an .evt file with multiple event tags denoting
% trial-specific (TRSP) information and write out .evt files with TRSP
% information in a single event tag. Useful for doing segmentation in NS.

% I need to collapse all codes into a single Code, based on information in
% the P### Code (stimulus presentation). The new code is STIM.

% STIM event will have these keys:

% cue type
% C001 = ear/auditory
% C002 = eye/visual
% key: CUET (string: ear/eye)

% start of prestimulus period
% PRES
% key: Don't need to transfer this info

% {1250, 1500, 1750}
% length of the prestimulus period
% key: PSTM (double: 1250, 1500, 1750)

% stimulus number
% P001-P320
% key: STMN

% subsequent memory perofrmance
% reca
% forg
% key: SMEM (double: 0, 1)

% target modality (same as cue type)
% M001 = auditory
% M002 = visual
% key: TMOD (string: auditory/visual)

% block
% B001-B020
% key: BLOC (double: 1-20)

% serial position
% S001-S016
% key: SPOS (double: 1-16)

% also add target type, either target buffer if S001, S002,
% S015, S016, or target stim if anything else
% key: STMT (string: STUDY_TARGET, STUDY_BUFFER)


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
  
  % initialize the timeMode
  timeModeStr = [];
  
  %% Read in the evt file
  [path,name,ext] = fileparts(evt(i).name);
  
  % figure out how many columns there are
  fid = fopen(fullfile(dataroot,evt(i).name),'r');
  % get the file name line
  fgetl(fid);
  % get the time mode line
  fgetl(fid);
  % get the header line
  tline = fgetl(fid);
  % headers
  if strcmp(tline(end),sprintf('\t'))
    tline = tline(1:end-1);
  end
  hdr = regexp(tline,'\t','split');
  cols.code = find(strcmp(hdr,'Code'));
  cols.label = find(strcmp(hdr,'Label'));
  cols.type = find(strcmp(hdr,'Track'));
  cols.track = find(strcmp(hdr,'Track'));
  cols.onset = find(strcmp(hdr,'Onset'));
  cols.duration = find(strcmp(hdr,'Duration'));
  
  % get these lines because we want to put them back in the file
  sessline = fgetl(fid);
  cellline = fgetl(fid);
  
  % get the first data line so we can count the number of columns
  tline = fgetl(fid);
  
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
  if strcmp(tline(end),sprintf('\t'))
    tline = tline(1:end-1);
  end
  % since it's tab delimited, split it and count the length
  numCols = length(regexp(tline,'\t','split'));
  
  fid = fopen(fullfile(dataroot,evt(i).name),'r');
  data = textscan(fid,repmat('%s',1,numCols),'Delimiter','\t');
  fclose(fid);
  
  % rename the old evt file
  feval(str2func('unix'),sprintf('mv %s %s',strrep(fullfile(dataroot,evt(i).name),' ','\ '),strrep(fullfile(dataroot,[name,'_separateTags',ext]),' ','\ ')));
  
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
    
    stimTag = sprintf('P%03d',stimNum);
    blocTag = sprintf('B%03d',blockNum);
    serpTag = sprintf('S%03d',serPos);
    
    switch cell2mat(data{cols.code}(j))
      case {'C001'}
        allStim(stimNum).CUET = 'ear';
      case {'C002'}
        allStim(stimNum).CUET = 'eye';
      case {'1250', '1500', '1750'}
        allStim(stimNum).PSTM = data{cols.code}(j);
      case {stimTag}
        allStim(stimNum).STMN = num2str(stimNum);
        prevInfo = [];
        for k = 2:numCols
          prevInfo = sprintf('%s\t%s',prevInfo,cell2mat(data{k}(j)));
        end
        allStim(stimNum).prevInfo = prevInfo(2:end);
        stimNum = stimNum + 1;
      case {'reca'}
        allStim(stimNum).SMEM = '1';
      case {'forg'}
        allStim(stimNum).SMEM = '0';
      case {'M001'}
        allStim(stimNum).TMOD = 'auditory';
      case {'M002'}
        allStim(stimNum).TMOD = 'visual';
      case {blocTag}
        allStim(stimNum).BLOC = num2str(blockNum);
        blockNum = blockNum + 1;
      case {serpTag}
        allStim(stimNum).SPOS = num2str(serPos);
        if strcmp(serpTag,'S001') || strcmp(serpTag,'S002') || strcmp(serpTag,'S015') || strcmp(serpTag,'S015')
          allStim(stimNum).STMT = 'STUDY_BUFFER';
        else
          allStim(stimNum).STMT = 'STUDY_TARGET';
        end
        
        % advance after finding the serial position tag because it's last
        serPos = serPos + 1;
    end % switch
    
  end
  
end


