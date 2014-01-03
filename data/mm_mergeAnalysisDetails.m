function [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails(orig_anaDetails_file,add_anaDetails_file,backup_orig_AD,sortBySubj,replaceOrig)

% concatenate the additional analysis details onto the end of the original
% analysis details, and save in the location of the original one

if ~exist('backup_orig_AD','var') || isempty(backup_orig_AD)
  backup_orig_AD = true;
end

if ~exist('sortBySubj','var') || isempty(sortBySubj)
  sortBySubj = false;
end

if ~exist('replaceOrig','var') || isempty(replaceOrig)
  replaceOrig = false;
end

if exist(orig_anaDetails_file,'file')
  orig_anaDetails = load(orig_anaDetails_file);
else
  error('original AD file (%s) does not exist!',orig_anaDetails_file);
end

if exist(add_anaDetails_file,'file')
  add_anaDetails = load(add_anaDetails_file);
else
  error('additional AD file (%s) does not exist!',add_anaDetails_file);
end

if all(ismember(orig_anaDetails.exper.eventValues,add_anaDetails.exper.eventValues))
  subjects_new = ~ismember(add_anaDetails.exper.subjects,orig_anaDetails.exper.subjects);
  subjects_exist = ismember(orig_anaDetails.exper.subjects,add_anaDetails.exper.subjects);
  
  if any(subjects_new) || (replaceOrig && any(subjects_exist))
    % don't backup until we're sure we're going to re-save the file
    if backup_orig_AD
      [pathstr,name,ext] = fileparts(orig_anaDetails_file);
      orig_AD_backup_file = fullfile(pathstr,sprintf('%s_backup%s',name,ext));
      fprintf('Backing up original analysisDetails file to %s... ',orig_AD_backup_file);
      s = unix(sprintf('cp %s %s',orig_anaDetails_file,orig_AD_backup_file));
      if s ~= 0
        error('Could not back up original AD file using the unix command!');
      else
        fprintf('Done.\n');
      end
    end
    
    if any(subjects_new)
      addSub_str = add_anaDetails.exper.subjects(subjects_new);
      if sortBySubj
        fprintf('Sorting (alphabetically)%s together with%s...\n',sprintf(repmat(' %s',1,sum(subjects_new)),addSub_str{:}),sprintf(repmat(' %s',1,sum(subjects_new)),orig_anaDetails.exper.subjects{:}));
      else
        fprintf('Concatenating%s after%s...\n',sprintf(repmat(' %s',1,sum(subjects_new)),addSub_str{:}),sprintf(repmat(' %s',1,sum(subjects_new)),orig_anaDetails.exper.subjects{:}));
      end
      
      subjects_all = cat(1,orig_anaDetails.exper.subjects,add_anaDetails.exper.subjects(subjects_new));
      
      newSubInd = find(subjects_new);
      
      if sortBySubj
        [~, subjInd] = sort(subjects_all);
      end
      
      nTrials_all = orig_anaDetails.exper.nTrials;
      nTr_fn = fieldnames(nTrials_all);
      for i = 1:length(nTr_fn)
        for j = 1:length(newSubInd)
          nTr_thisSub = add_anaDetails.exper.nTrials.(nTr_fn{i});
          nTrials_all.(nTr_fn{i}) = cat(1,nTrials_all.(nTr_fn{i}),nTr_thisSub(newSubInd(j),:));
        end
        if sortBySubj
          nTrials_all.(nTr_fn{i}) = nTrials_all.(nTr_fn{i})(subjInd);
        end
      end
      
      badChan_all = cat(1,orig_anaDetails.exper.badChan,add_anaDetails.exper.badChan(newSubInd,:));
      badEv_all = cat(1,orig_anaDetails.exper.badEv,add_anaDetails.exper.badEv(newSubInd,:));
      
      if sortBySubj
        subjects_all = subjects_all(subjInd);
        badChan_all = badChan_all(subjInd);
        badEv_all = badEv_all(subjInd);
      end
      
      % add the combined info into the struct we want to save
      orig_anaDetails.exper.subjects = subjects_all;
      orig_anaDetails.exper.nTrials = nTrials_all;
      orig_anaDetails.exper.badChan = badChan_all;
      orig_anaDetails.exper.badEv = badEv_all;
      
      % fill in the variables to save
      exper = orig_anaDetails.exper;
      ana = orig_anaDetails.ana;
      dirs = orig_anaDetails.dirs;
      files = orig_anaDetails.files;
      cfg_proc = orig_anaDetails.cfg_proc;
      cfg_pp = orig_anaDetails.cfg_pp;
      
      fprintf('\nResaving analysis details with additional data: %s...',orig_anaDetails_file);
      save(orig_anaDetails_file,'exper','ana','dirs','files','cfg_proc','cfg_pp');
      fprintf('Done.\n');
    else
      fprintf('There are no new subjects in exper.subjects (everyone is already in the original analysisDetails.mat).  Not saving anything with additional subjects.\n');
    end
    
    if replaceOrig
      if any(subjects_exist)
        % get the original data
        nTrials_orig = orig_anaDetails.exper.nTrials;
        nTr_fn = fieldnames(nTrials_orig);
        badChan_orig = orig_anaDetails.exper.badChan;
        badEv_orig = orig_anaDetails.exper.badEv;
        
        for os = 1:length(orig_anaDetails.exper.subjects)
          if subjects_exist(os)
            fprintf('Replacing %s in original data...\n',orig_anaDetails.exper.subjects{os});
            replSubInd = find(ismember(add_anaDetails.exper.subjects,orig_anaDetails.exper.subjects(os)));
            % move nTrials
            for fn = 1:length(nTr_fn)
              nTrials_orig.(nTr_fn{fn})(os) = add_anaDetails.exper.nTrials.(nTr_fn{fn})(replSubInd);
            end
            % move badChan
            badChan_orig{os} = add_anaDetails.exper.badChan{replSubInd,:};
            % move badEv
            badEv_orig{os} = add_anaDetails.exper.badEv{replSubInd,:};
          end
        end
        
        % add the combined info into the struct we want to save
        %exper.subjects = orig_anaDetails.exper.subjects;
        orig_anaDetails.exper.nTrials = nTrials_orig;
        orig_anaDetails.exper.badChan = badChan_orig;
        orig_anaDetails.exper.badEv = badEv_orig;
        
        % fill in the variables to save
        exper = orig_anaDetails.exper;
        ana = orig_anaDetails.ana;
        dirs = orig_anaDetails.dirs;
        files = orig_anaDetails.files;
        cfg_proc = orig_anaDetails.cfg_proc;
        cfg_pp = orig_anaDetails.cfg_pp;
        
        fprintf('\nResaving analysis details with replacement data: %s...',orig_anaDetails_file);
        save(orig_anaDetails_file,'exper','ana','dirs','files','cfg_proc','cfg_pp');
        fprintf('Done.\n');
      else
        fprintf('\nThere are no subjects to replace in the original analysis details file. Not saving anything with replacements.\n');
      end
    end
  end
  
else
  fprintf('exper.eventValues does not match up between original and additional subjects! Not re-saving analysisDetails.mat.\n');
end
