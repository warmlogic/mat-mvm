
% the first 4 subjects had study presentation time logged incorrectly
fixStudyPresSubs = {'RSRC3001','RSRC3002','RSRC3003','RSRC3004'};

% subjects 6-10 had the wrong RK_RESP recorded
% key  right_coding  wrong_coding
%  /        DU           RS
%  .        MU           RO
%  ,        MF           DF
%  M        DF           MF
%  X        RO           MU
%  Z        RS           DU
fixRKRespSubs = {'RSRC3006','RSRC3007','RSRC3008','RSRC3009','RSRC3010',};

confRange_rec_inc = {'REMEMBER_SOURCE','REMEMBER_OTHER','DEFINITELY_FAMILIAR','MAYBE_FAMILIAR','MAYBE_UNFAMILIAR','DEFINITELY_UNFAMILIAR'};
confRange_rec_cor = {'DEFINITELY_UNFAMILIAR','MAYBE_UNFAMILIAR','MAYBE_FAMILIAR','DEFINITELY_FAMILIAR','REMEMBER_OTHER','REMEMBER_SOURCE'};

%corKeys = {'/','.',',','M','X','Z'}; % new to old
newKeys = {'/','.'};
oldKeys = {',','M','X','Z'};

newStrs = {'DEFINITELY_UNFAMILIAR','MAYBE_UNFAMILIAR'};
oldStrs = {'MAYBE_FAMILIAR','DEFINITELY_FAMILIAR','REMEMBER_OTHER','REMEMBER_SOURCE'};

dataroot = '/Volumes/curranlab/Data/RSRC3/eeg/behavioral';
subjects = unique(cat(2,fixStudyPresSubs,fixRKRespSubs));


for sub = 1:length(subjects)
  subject = subjects{sub};
  
  if ismember(subject,fixStudyPresSubs)
    sessions = {'session_0'};
    for ses = 1:length(sessions)
      fprintf('Working on %s, %s...\n',subjects{sub},sessions{ses});
      
      sesDir = fullfile(dataroot,subject,sessions{ses});
      
      logfilename = 'session.log';
      logfile = fullfile(sesDir,logfilename);
      
      fid = fopen(logfile,'r');
      logdata = textscan(fid,'%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
      fclose(fid);
      
      % copy the log file to the old file
      unix(sprintf('cp %s %s',fullfile(sesDir,logfilename),fullfile(sesDir,sprintf('%s.old',logfilename))));
      
      % start the new log file
      logfilenew = fullfile(sesDir,sprintf('%s',logfilename));
      outfile = fopen(logfilenew,'wt');
      
      for i = 1:length(logdata{1})
        if strcmp(logdata{3}(i),'STUDY_BUFFER') || strcmp(logdata{3}(i),'STUDY_TARGET')
          presTime = str2double(logdata{1}(i+1)) - str2double(logdata{6}(i+1));
          logdata{1}(i) = {num2str(presTime)};
        end
        
        lineStr = [];
        for j = 1:9
          if ~isempty(cell2mat(logdata{j}(i)))
            lineStr = sprintf('%s%s',lineStr,cell2mat(logdata{j}(i)));
            if j < 9
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
    end % sessions
    
  elseif ismember(subject,fixRKRespSubs)
    sessions = {'session_1'};
  
    for ses = 1:length(sessions)
      fprintf('Working on %s, %s...\n',subjects{sub},sessions{ses});
      
      sesDir = fullfile(dataroot,subject,sessions{ses});
      
      logfilename = 'session.log';
      logfile = fullfile(sesDir,logfilename);
      
      fid = fopen(logfile,'r');
      logdata = textscan(fid,'%s%s%s%s%s%s%s%s%s','delimiter','\t','emptyvalue',NaN);
      fclose(fid);
      
      % copy the log file to the old file
      unix(sprintf('cp %s %s',fullfile(sesDir,logfilename),fullfile(sesDir,sprintf('%s.old',logfilename))));
      
      % start the new log file
      logfilenew = fullfile(sesDir,sprintf('%s',logfilename));
      outfile = fopen(logfilenew,'wt');
      
      for i = 1:length(logdata{3})
        if strcmp(logdata{3}(i),'RK_RESP')
          
          %if str2double(logdata{6}(i-1)) == 1
          %  fprintf('This is actually a TARGET\n');
          %elseif str2double(logdata{6}(i-1)) == 0
          %  fprintf('This is actually a LURE\n');
          %end
          
          respIndex = find(ismember(confRange_rec_inc,logdata{5}(i)));
          %fprintf('Changing %s to %s (%s)\n',confRange_rec_inc{respIndex},confRange_rec_cor{respIndex},cell2mat(logdata{4}(i)));
          logdata{5}(i) = confRange_rec_cor(respIndex);
          
          %         if strcmp(logdata{5}(i),confRange_rec_inc{1})
          %           %fprintf('Changing %s to %s (%s)\n',cell2mat(logdata{5}(i)),confRange_rec_cor{1},cell2mat(logdata{4}(i)));
          %           logdata{5}(i) = confRange_rec_cor(1);
          %         elseif strcmp(logdata{5}(i),confRange_rec_inc{2})
          %           %fprintf('Changing %s to %s (%s)\n',cell2mat(logdata{5}(i)),confRange_rec_cor{2},cell2mat(logdata{4}(i)));
          %           logdata{5}(i) = confRange_rec_cor(2);
          %         elseif strcmp(logdata{5}(i),confRange_rec_inc{3})
          %           %fprintf('Changing %s to %s (%s)\n',cell2mat(logdata{5}(i)),confRange_rec_cor{3},cell2mat(logdata{4}(i)));
          %           logdata{5}(i) = confRange_rec_cor(3);
          %         elseif strcmp(logdata{5}(i),confRange_rec_inc{4})
          %           %fprintf('Changing %s to %s (%s)\n',cell2mat(logdata{5}(i)),confRange_rec_cor{4},cell2mat(logdata{4}(i)));
          %           logdata{5}(i) = confRange_rec_cor(4);
          %         elseif strcmp(logdata{5}(i),confRange_rec_inc{5})
          %           %fprintf('Changing %s to %s (%s)\n',cell2mat(logdata{5}(i)),confRange_rec_cor{5},cell2mat(logdata{4}(i)));
          %           logdata{5}(i) = confRange_rec_cor(5);
          %         elseif strcmp(logdata{5}(i),confRange_rec_inc{6})
          %           %fprintf('Changing %s to %s (%s)\n',cell2mat(logdata{5}(i)),confRange_rec_cor{6},cell2mat(logdata{4}(i)));
          %           logdata{5}(i) = confRange_rec_cor(6);
          %         end
          
          if ismember(logdata{4}(i),oldKeys) && str2double(logdata{6}(i-1)) == 1
            %fprintf('Correct\n');
            logdata{8}(i) = {'1'};
          elseif ismember(logdata{4}(i),oldKeys) && str2double(logdata{6}(i-1)) == 0
            %fprintf('Incorrect\n');
            logdata{8}(i) = {'0'};
          elseif ismember(logdata{4}(i),newKeys) && str2double(logdata{6}(i-1)) == 0
            %fprintf('Correct\n');
            logdata{8}(i) = {'1'};
          elseif ismember(logdata{4}(i),newKeys) && str2double(logdata{6}(i-1)) == 1
            %fprintf('Incorrect\n');
            logdata{8}(i) = {'0'};
          end
          
        end % if RK_RESP
        
        lineStr = [];
        for j = 1:9
          if ~isempty(cell2mat(logdata{j}(i)))
            lineStr = sprintf('%s%s',lineStr,cell2mat(logdata{j}(i)));
            if j < 9
              if ~isempty(cell2mat(logdata{j+1}(i)))
                lineStr = sprintf('%s\t',lineStr);
              end
            end
          end
        end
        
        %fprintf(outfile,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',cell2mat(logdata{1}(i)),cell2mat(logdata{2}(i)),cell2mat(logdata{3}(i)),cell2mat(logdata{4}(i)),cell2mat(logdata{5}(i)),cell2mat(logdata{6}(i)),cell2mat(logdata{7}(i)),cell2mat(logdata{8}(i)),cell2mat(logdata{9}(i)));
        fprintf(outfile,'%s\n',lineStr);
        
      end % logdata
      
      fclose(outfile);
      fprintf('Done.\n');
    end % sessions
  end % ismember
end % subjects

