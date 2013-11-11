function mm_mergeAnalysisDetails(orig_anaDetails_file,add_anaDetails_file,backup_orig_AD)

% concatenate the additional analysis details onto the end of the original
% analysis details, and save in the location of the original one

if ~exist('backup_orig_AD','var') || isempty(backup_orig_AD)
  backup_orig_AD = true;
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
  
  if any(subjects_new)
    % don't backup until we're sure we're going to re-save the file
    if backup_orig_AD
      [pathstr,name,ext] = fileparts(orig_anaDetails_file);
      orig_AD_backup_file = fullfile(pathstr,sprintf('%s_backup',name),ext);
      fprintf('Backing up original analysisDetails file to %s...\n',orig_AD_backup_file);
      s = unix(sprintf('cp %s %s',orig_anaDetails_file,orig_AD_backup_file));
      if s ~= 0
        error('Could not back up original AD file using the unix command!');
      else
        fprintf('Done.\n');
      end
    end
    
    subjects_all = cat(1,orig_anaDetails.exper.subjects,add_anaDetails.exper.subjects(subjects_new));
    
    newSubInd = find(subjects_new);
    
    nTrials_all = orig_anaDetails.exper.nTrials;
    nTr_fn = fieldnames(nTrials_all);
    for i = 1:length(nTr_fn)
      for j = 1:length(newSubInd)
        nTr_thisSub = add_anaDetails.exper.nTrials.(nTr_fn{i});
        nTrials_all.(nTr_fn{i}) = cat(1,nTrials_all.(nTr_fn{i}),nTr_thisSub(newSubInd(j),:));
      end
    end
    
    badChan_all = cat(1,orig_anaDetails.exper.badChan,add_anaDetails.exper.badChan(newSubInd,:));
    badEv_all = cat(1,orig_anaDetails.exper.badEv,add_anaDetails.exper.badEv(newSubInd,:));
    
    % add the combined info into the struct we want to save
    orig_anaDetails.exper.subjects = subjects_all;
    orig_anaDetails.exper.nTrials = nTrials_all;
    orig_anaDetails.exper.badChan = badChan_all;
    orig_anaDetails.exper.badEv = badEv_all;
    
    % enumerate the variables to save
    exper = orig_anaDetails.exper;
    ana = orig_anaDetails.ana;
    dirs = orig_anaDetails.dirs;
    files = orig_anaDetails.files;
    cfg_proc = orig_anaDetails.cfg_proc;
    cfg_pp = orig_anaDetails.cfg_pp;
    
    save(orig_anaDetails_file,'exper','ana','dirs','files','cfg_proc','cfg_pp');
    fprintf('Done.\n');
  else
    fprintf('There are no new subjects in exper.subjects (everyone is already in the old analysisDetails.mat).  Not re-saving analysisDetails.mat.\n');
  end
else
  fprintf('exper.eventValues does not match up between old and new subjects! Not re-saving analysisDetails.mat.\n');
end
