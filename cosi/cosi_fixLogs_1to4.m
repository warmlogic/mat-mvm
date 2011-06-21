% study
%
% side - remove 'none'
%
% color - add color string before RGB
%
% both - swap columns 3 and 4

% test
%
% side - remove color info
% side - add study coordinates
%
% both - change TARG_PRES and LURE_PRES to TEST_TARG and TEST_LURE

expName = 'COSI';

% both sessions
subjects = {'COSI001','COSI002','COSI003','COSI004'};
sessions = {'session_0','session_1'};

numCols = 12;

dataroot = fullfile(getenv('HOME'),'data',expName,'eeg','behavioral');

for sub = 1:length(subjects)
  for ses = 1:length(sessions)
    fprintf('Working on %s, %s...\n',subjects{sub},sessions{ses});
    
    sesDir = fullfile(dataroot,subjects{sub},sessions{ses});
    
    logfilename = 'session.log';
    logfile = fullfile(sesDir,logfilename);
    
    % copy the log file to the old file
    unix(sprintf('cp %s %s',fullfile(sesDir,logfilename),fullfile(sesDir,sprintf('%s.old',logfilename))));
    
    fid = fopen(logfile);
    logdata = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
    fclose(fid);
    
    % start the new log file
    logfilenew = fullfile(sesDir,sprintf('%s',logfilename));
    outfile = fopen(logfilenew,'wt');
    
    for i = 1:length(logdata{1})
      %%%%%%%%%%%%%%%%%
      % study
      %%%%%%%%%%%%%%%%%
      
      % swap columns 3 and 4
      if strcmp(logdata{3}(i),'side') || strcmp(logdata{3}(i),'color')
        cond = logdata{3}(i);
        type = logdata{4}(i);
        logdata{3}(i) = type;
        logdata{4}(i) = cond;
      end
    end
      
    for i = 1:length(logdata{1})
      %%%%%%%%%%%%%%%%%
      % study
      %%%%%%%%%%%%%%%%%
      
      % remove 'none'
      if strcmp(logdata{4}(i),'side') && strcmp(logdata{11}(i),'none')
        logdata{11}(i) = {''};
      end
      
      % add color name before RGB
      if strcmp(logdata{4}(i),'color') && (strcmp(logdata{3}(i),'STUDY_BUFFER') || strcmp(logdata{3}(i),'STUDY_TARGET'))
        if strcmp(logdata{11}(i),'(3, 67, 225, 255)')
          rgb = logdata{11}(i);
          col = {'blue'};
          logdata{11}(i) = col;
          logdata{12}(i) = rgb;
        elseif strcmp(logdata{11}(i),'(229, 0, 0, 255)')
          rgb = logdata{11}(i);
          col = {'red'};
          logdata{11}(i) = col;
          logdata{12}(i) = rgb;
        end
      end
      
      %%%%%%%%%%%%%%%%%
      % test
      %%%%%%%%%%%%%%%%%
      
      if strcmp(logdata{3}(i),'TARG_PRES')
        logdata{3}(i) = {'TEST_TARGET'};
        if strcmp(logdata{4}(i),'side')
          name = logdata{6}(i);
          for j = 1:length(logdata{1})
            if strcmp(logdata{6}(j),name) && (strcmp(logdata{3}(j),'STUDY_BUFFER') || strcmp(logdata{3}(j),'STUDY_TARGET'))
              xcoord = logdata{8}(j);
              ycoord = logdata{9}(j);
              break
            end
          end
          logdata{9}(i) = xcoord;
          logdata{10}(i) = ycoord;
          logdata{11}(i) = {''};
          logdata{12}(i) = {''};
        end
      elseif strcmp(logdata{3}(i),'LURE_PRES')
        logdata{3}(i) = {'TEST_LURE'};
        if strcmp(logdata{4}(i),'side')
          logdata{9}(i) = {'-1'};
          logdata{10}(i) = {'-1'};
          logdata{11}(i) = {''};
          logdata{12}(i) = {''};
        end
      end
      
      %%%%%%%%%%%%%%%%%
      % write it out
      %%%%%%%%%%%%%%%%%
      
      lineStr = [];
      for j = 1:numCols
        if ~isempty(cell2mat(logdata{j}(i)))
          lineStr = sprintf('%s%s',lineStr,cell2mat(logdata{j}(i)));
          if j < numCols
            if ~isempty(cell2mat(logdata{j+1}(i)))
              lineStr = sprintf('%s\t',lineStr);
            end
          end
        end
      end
      fprintf(outfile,'%s\n',lineStr);
      
    end % logdata
    fclose(outfile);
    fprintf('Done.\n');
    
  end
end
