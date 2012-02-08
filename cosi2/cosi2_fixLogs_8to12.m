% COSI2, subjects 008-012, test period
%
% 009: The keys that were recorded (and the responses) for Y and B source
% responses were reversed. So subject 9 thought X was for B and / was for Y
% (because that's what they saw on the screen), but it was recorded as X=Y
% and /=B. The fix is to swap the recorded color responses and the accuracy
% (and not the keys).
%
% the others just had color responses for false alarms recorded incorrectly

expName = 'COSI2';

logfilename = 'session.log';

% both sessions
subjects = {'COSI2008','COSI2009','COSI2010','COSI2011','COSI2012'};
sessions = {'session_0','session_1'};

numCols = 13;

%dataroot = fullfile(getenv('HOME'),'data',expName,'eeg','behavioral');
dataroot = fullfile(filesep,'Volumes','curranlab','Data',expName,'eeg','behavioral');

for sub = 1:length(subjects)
  for ses = 1:length(sessions)
    fprintf('Working on %s, %s...\n',subjects{sub},sessions{ses});
    
    sesDir = fullfile(dataroot,subjects{sub},sessions{ses});
    
    logfile = fullfile(sesDir,logfilename);
    
    % copy the log file to the old file
    if exist(logfile,'file')
      logfileOld = fullfile(sesDir,sprintf('%s.old',logfilename));
      if ~exist(logfileOld,'file')
        unix(sprintf('cp %s %s',fullfile(sesDir,logfilename),logfileOld));
      else
        error('%s already exists!',logfileOld);
      end
    else
      fprintf('%s does not exist. Moving on.\n',logfile);
      continue
    end
    
    fid = fopen(logfile,'r');
    logdata = textscan(fid,repmat('%s',1,numCols),'delimiter','\t','emptyvalue',NaN);
    fclose(fid);
    
    % start the new log file
    outfile = fopen(logfile,'wt');
    
    fprintf('Fixing log files');
    for i = 1:length(logdata{1})
      %%%%%%%%%%%%%%%%%
      % test
      %%%%%%%%%%%%%%%%%
      
      % set the false alarm response to the right color
      if strcmp(subjects{sub},'COSI2008') || strcmp(subjects{sub},'COSI2009') || strcmp(subjects{sub},'COSI2010')
        if strcmp(logdata{3}(i),'SOURCE_RESP') || strcmp(logdata{3}(i),'PRACTICE_SOURCE_RESP')
          if strcmp(logdata{4}(i-1),'color')
            if strcmp(logdata{4}{i},'/') && strcmp(logdata{5}{i},'Yellow')
              logdata{5}{i} = 'Blue';
              fprintf('.');
            end
          end
        end
      elseif strcmp(subjects{sub},'COSI2011')
        if strcmp(logdata{3}(i),'SOURCE_RESP') || strcmp(logdata{3}(i),'PRACTICE_SOURCE_RESP')
          if strcmp(logdata{4}(i-1),'color')
            if strcmp(logdata{4}{i},'.') && strcmp(logdata{5}{i},'Blue')
              logdata{5}{i} = 'Yellow';
              fprintf('.');
            end
          end
        end
      elseif strcmp(subjects{sub},'COSI2012')
        if strcmp(logdata{3}(i),'SOURCE_RESP') || strcmp(logdata{3}(i),'PRACTICE_SOURCE_RESP')
          if strcmp(logdata{4}(i-1),'color')
            if strcmp(logdata{4}{i},'.') && strcmp(logdata{5}{i},'Yellow')
              logdata{5}{i} = 'Blue';
              fprintf('.');
            end
          end
        end
      end
      
      if strcmp(subjects{sub},'COSI2009')
        % swap blue and yellow responses (not the keys)
        if strcmp(logdata{3}(i),'SOURCE_RESP') || strcmp(logdata{3}(i),'PRACTICE_SOURCE_RESP')
          if strcmp(logdata{4}(i-1),'color')
            if strcmp(logdata{4}{i},'/') && strcmp(logdata{5}{i},'Blue')
              % change the color response
              logdata{5}{i} = 'Yellow';
              % and fix the accuracy
              if strcmp(logdata{9}{i-1},logdata{5}{i}) && strcmp(logdata{8}{i},'0')
                logdata{8}{i} = '1';
              end
              fprintf('.');
            elseif strcmp(logdata{4}{i},'X') && strcmp(logdata{5}{i},'Yellow')
              % change the color response
              logdata{5}{i} = 'Blue';
              % and fix the accuracy
              if strcmp(logdata{9}{i-1},logdata{5}{i}) && strcmp(logdata{8}{i},'0')
                logdata{8}{i} = '1';
              end
              fprintf('.');
            end
          end
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
    
  end % ses
end % sub
