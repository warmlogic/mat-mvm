function [exper] = mm_getBadChan(cfg,exper,dirs)
%MM_GETBADCHAN - Finds bad channel listings
%
% [exper] = mm_getBadChan(cfg,exper,dirs)
%
% Input:
%
% cfg.badChanManual=true requires a tab-delimited file titled
% [exper.name,'_badChan.txt'] to reside in
% fullfile(dirs.dataroot,dirs.dataDir). The three tab sections are subject
% name (e.g., EXPER001), session name (e.g., session_0), and bad channel
% numbers listed as integers in brackets (e.g., [56 93]).
%
% cfg.badChanEP=true requires the Artifact_Correction_Log output from EP
% Toolkit artifact processing, and must reside in a directory labeled with
% the session name (from exper.sessions) which is in a directory called
% 'ep_art' in fullfile(dirs.dataroot,dirs.dataDir). This will only look for
% channels listed as being globally bad.
%
% Output:
%
% Adds badChan field to exper struct. This will overwrite the bad channel
% field that gets put there using mm_ft_artifact.
%
% See also: MM_FT_ARTIFACT
%

% TODO: Deal with shorted channels
% "Warning: shorted channels: E1-E2; E3-E4"

if ~isfield(cfg,'badChanManual')
  cfg.badChanManual = false;
end
if ~isfield(cfg,'badChanEP')
  cfg.badChanEP = false;
end

if cfg.badChanManual == false && cfg.badChanEP == false
  error('One or both of cfg.badChanManual and cfg.badChanEP must be set to true.');
end

% initialize
if isfield(exper,'badChan')
  error('Bad channel information already exists in this exper struct. Remove the field and run %s again if you want new bad channel information.',mfilename);
else
  exper.badChan = cell(length(exper.subjects),length(exper.sessions));
end

fprintf('Finding ');
if cfg.badChanManual
  fprintf('manually entered ');
  if cfg.badChanEP
    fprintf('and ');
  end
end
if cfg.badChanEP
  fprintf('EP Toolkit ');
end
fprintf('bad channel information.\n');

for sub = 1:length(exper.subjects)
  subject = exper.subjects{sub};
  
  % initialize
  badChanAllSes = [];
  
  for ses = 1:length(exper.sessions)
    session = exper.sessions(ses);
    for thisSes = 1:length(session)
      sesName = session{thisSes};
      
      if cfg.badChanManual
        %clear badChanManual
        foundBadChan = false;
        
        fid = fopen(fullfile(dirs.dataroot,dirs.dataDir,[exper.name,'_badChan.txt']));
        
        if fid == -1
          error('Could not open %s. Make sure you do not have it open in another application.',fullfile(dirs.dataroot,dirs.dataDir,[exper.name,'_badChan.txt']));
        else
          [badChanInfo] = textscan(fid,'%s%s%s','delimiter','\t');
          fclose(fid);
          for i = 1:length(badChanInfo{1})
            if strcmp(badChanInfo{1}{i},subject) && strcmp(badChanInfo{2}{i},sesName)
              badChanInfo{3}{i} = strrep(strrep(badChanInfo{3}{i},'[',''),']','');
              badChanManual = eval(sprintf('[%s];',badChanInfo{3}{i}));
              foundBadChan = true;
              break
            else
              badChanManual = [];
            end
          end
          if ~foundBadChan
            %error('No manual bad channel information found for %s %s.',subject,sesName);
            warning([mfilename,':noManualBadChan'],'No manual bad channel information found for %s %s.\n',subject,sesName);
          end
        end
      end
      
      if cfg.badChanEP
        foundBadChan = false;
        foundNoBadChan = false;
        
        epGBC_str = 'Global bad Channels: ';
        ep_prefix = 'Artifact_Correction_Log';
        
        %   datasetFile = ft_findcfg(data.cfg,'dataset');
        %   if ~isempty(datasetFile)
        %     datasetSep = strfind(datasetFile,filesep);
        %     if ~isempty(datasetSep)
        %       datasetFile = datasetFile(datasetSep(end)+1:end);
        %     end
        %     datasetSep = strfind(datasetFile,'.');
        %     datasetFile = datasetFile(1:datasetSep(end)-1);
        %     datasetFile = strrep(datasetFile,'_e','');
        %   else
        %     datasetFile = subject;
        %   end
        %   epArtFile = dir(fullfile(dirs.dataroot,dirs.dataDir,'ep_art',sesName,sprintf('%s %s*.txt',ep_prefix,datasetFile)));
        epArtFile = dir(fullfile(dirs.dataroot,dirs.dataDir,'ep_art',sesName,sprintf('%s %s*.txt',ep_prefix,subject)));
        
        if ~isempty(epArtFile)
          if length(epArtFile) == 1
            epArtFile = epArtFile.name;
            fid = fopen(fullfile(dirs.dataroot,dirs.dataDir,'ep_art',sesName,epArtFile),'r');
            if fid == -1
              error('Could not open the file for %s %s. Make sure you do not have it open in another application.',subject,sesName);
            else
              
              while ~foundBadChan
                % get the next line
                tline = fgetl(fid);
                if strncmpi(tline,epGBC_str,length(epGBC_str))
                  foundBadChan = true;
                  break
                end
              end
              fclose(fid);
              
              badChanInfo = strrep(tline,epGBC_str,'');
              if strcmp(badChanInfo(end),sprintf('\t'))
                badChanInfo = badChanInfo(1:end-1);
              end
              if strcmp(badChanInfo,'None')
                foundNoBadChan = true;
              end
              if foundBadChan
                if ~foundNoBadChan
                  %badChanEP = regexp(badChanEP,'\t','split');
                  badChanEP = eval(sprintf('[%s];',badChanInfo));
                elseif foundNoBadChan
                  badChanEP = [];
                end
              elseif ~foundBadChan
                error('No bad channel information found.');
              end
            end
          elseif length(epArtFile) > 1
            error('More than one EP Articifact Correction Log found for %s %s, probably multiple sessions in the same directory. Make sure each subject has one Artifact Log file in individual session folders, or just put the info in a text file.',subject,sesName);
          end
        elseif isempty(epArtFile)
          %error('No EP Articifact Correction Log found for %s %s.',subject,sesName);
          warning([mfilename,':noEPBadChan'],'No EP Toolkit bad channel information found for %s %s.\n',subject,sesName);
        end
      end
      
      % collect the bad channels
      if cfg.badChanManual || cfg.badChanEP
        if ~exist('badChanManual','var')
          badChanManual = [];
        end
        if ~exist('badChanEP','var')
          badChanEP = [];
        end
        badChan = unique(cat(2,badChanManual,badChanEP));
      else
        badChan = [];
      end
      badChanAllSes = cat(2,badChanAllSes,badChan);
      
    end % thisSes
    
    exper.badChan{sub,ses} = badChanAllSes;
  end % ses
end % sub

end
