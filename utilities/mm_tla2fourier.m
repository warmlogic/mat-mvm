function mm_tla2fourier(exper,dirs,cfg_ana,cfg_ft,cfg_alt)

% function mm_tla2fourier(exper,dirs,cfg_ana,cfg_ft,cfg_alt)
%
% Take timelock data and convert it into fourier data. Can also
% subsequently turn it into power data.
%

if nargin ~= 5
  error('Not the correct number of input arguments!');
end

if ~isfield(cfg_ana,'orig_ftype')
  cfg_ana.orig_ftype = 'tla';
end
if ~isfield(cfg_ana,'orig_param')
  cfg_ana.orig_param = 'trial';
end

if ~isfield(cfg_ana,'param')
  cfg_ana.out_param = 'fourierspctrm';
end

if ~isfield(cfg_ana,'output2alt')
  cfg_ana.output2alt = true;
end
if ~isfield(cfg_ana,'alt_ftype')
  cfg_ana.alt_ftype = 'pow';
end
if ~isfield(cfg_ana,'alt_param')
  cfg_ana.alt_param = 'powspctrm';
end

if ~isfield(cfg_ana,'useLockFiles')
  cfg_ana.useLockFiles = false;
end

if ~isfield(cfg_ana,'splitTrials')
  cfg_ana.splitTrials = false;
end
if cfg_ana.splitTrials
  if ~isfield(cfg_ana,'splitSize')
    cfg_ana.splitSize = 100;
  end
  % number of trials to lump in with the last trial split. Must be >= 1.
  if ~isfield(cfg_ana,'splitRemainderLump')
    cfg_ana.splitRemainderLump = 10;
  end
  
  % whether to see if files are too big to combine (RAM limit issue)
  if ~isfield(cfg_ana,'checkSplitFileSizeForSaving')
    cfg_ana.checkSplitFileSizeForSaving = false;
  end
  % approximate limit in MB that the combined files can reach; otherwise will
  % keep the split files
  if cfg_ana.checkSplitFileSizeForSaving
    if ~isfield(cfg_ana,'combineSavingLimitMB')
      cfg_ana.combineSavingLimitMB = 3000;
    end
  end
end

if ~isfield(cfg_ana,'rmPreviousCfg')
  %warning('Upon loading of data, the data.cfg.previous field will NOT be removed! This may cause increased memory usage!');
  %cfg_ana.rmPreviousCfg = false;
  warning('Upon loading of data, the data.cfg.previous field WILL be removed to save memory!');
  cfg_ana.rmPreviousCfg = true;
end
if ~cfg_ana.rmPreviousCfg
  warning('The data.cfg.previous field will NOT be removed! This may cause increased memory usage!');
end

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    fprintf('%s %s...\n',exper.subjects{sub},exper.sesStr{ses});
    
    sesDir = fullfile(dirs.saveDirProc,exper.subjects{sub},exper.sesStr{ses});
    origData = dir(fullfile(sesDir,sprintf('data_%s*.mat',cfg_ana.orig_ftype)));
    
    saveDir = fullfile(cfg_ana.saveroot,exper.subjects{sub},exper.sesStr{ses});
    if ~exist(saveDir,'dir')
      [s] = mkdir(saveDir);
      if ~s
        error('Could not make %s',saveDir);
      end
    end
    
    saveDir_alt = fullfile(cfg_ana.saveroot_alt,exper.subjects{sub},exper.sesStr{ses});
    if ~exist(saveDir_alt,'dir')
      [s] = mkdir(saveDir_alt);
      if ~s
        error('Could not make %s',saveDir_alt);
      end
    end
    
    for od = 1:length(origData)
      % set up file names
      origFile = origData(od).name;
      origFile_full = fullfile(sesDir,origFile);
      
      outputFile = strrep(origFile,cfg_ana.orig_ftype,cfg_ft.output);
      outputFile_full_orig = fullfile(saveDir,outputFile);
      
      if cfg_ana.output2alt
        outputFile_alt = strrep(origFile,cfg_ana.orig_ftype,cfg_ana.alt_ftype);
        outputFile_alt_full_orig = fullfile(saveDir_alt,outputFile_alt);
      end
      
      if cfg_ana.useLockFiles
        % check the lock on the final file that we should compute
        if cfg_ana.output2alt
          lockedFile = outputFile_alt_full_orig;
        else
          lockedFile = outputFile_full_orig;
        end
        
        % file could be in limbo; see if it's already been started
        if ~lockFile(lockedFile)
          fprintf('Moving on, final file already exists: %s\n',lockedFile);
          continue
        end
      end
      
      if ~exist(outputFile_full_orig,'file') || (cfg_ana.output2alt && ~exist(outputFile_alt_full_orig,'file'))
        % load original regardless of what output we have
        fprintf('\tLoading %s: %s...\n',cfg_ana.orig_ftype,origFile);
        orig = load(origFile_full);
        origType = fieldnames(orig);
        fprintf('Done.\n');
        
        if length(origType) == 1
          data_fn = cell2mat(origType);
        else
          error('More than one field in data struct! There should only be one.\n');
        end
        
        % save space by removing a potentially large 'previous' field
        if cfg_ana.rmPreviousCfg
          if isfield(orig.(data_fn).cfg,'previous')
            fprintf('Removing ''previous'' field from data.cfg\n');
            orig.(data_fn).cfg = rmfield(orig.(data_fn).cfg,'previous');
          end
        end
        
        if cfg_ana.splitTrials
          % split up trials if there are more than the RAM can handle
          nTrials = size(orig.(data_fn).(cfg_ana.orig_param),1);
          if nTrials > (cfg_ana.splitSize + cfg_ana.splitRemainderLump)
            fprintf('Splitting trials...\n');
            
            tooBigToCombine = false;
            
            outputFile_full = outputFile_full_orig;
            if cfg_ana.output2alt
              outputFile_alt_full = outputFile_alt_full_orig;
            end
            
            nSplits = floor(nTrials / cfg_ana.splitSize);
            remainSize = nTrials - (cfg_ana.splitSize * nSplits);
            splitSizes = repmat(cfg_ana.splitSize,1,nSplits);
            if remainSize > 0
              % if it's only a few trials remaining, just lump them in with
              % the last split; also, ft_freqdescriptives will error if
              % only 1 trial is processed by itself because it expects some
              % other type of data
              if remainSize > cfg_ana.splitRemainderLump
                splitSizes = cat(2,splitSizes,remainSize);
              else
                splitSizes(end) = splitSizes(end) + remainSize;
              end
            end
            
            % trial counter
            count = 0;
            for sp = 1:length(splitSizes)
              cfg_ft.trials = false(1,nTrials);
              
              cfg_ft.trials(count+1:count+splitSizes(sp)) = true;
              
              if sp == 1
                origExt = '.mat';
                newExt = sprintf('_%d.mat',sp);
                outputFile_full = strrep(outputFile_full,origExt,newExt);
              elseif sp > 1
                origExt = sprintf('_%d.mat',sp-1);
                newExt = sprintf('_%d.mat',sp);
                outputFile_full = strrep(outputFile_full,origExt,newExt);
              end
              
              % do the first conversion if it is not already calculated
              if ~exist(outputFile_full,'file') && ~exist(outputFile_full_orig,'file')
                fprintf('\tSplit %d: Calculating %s for trials %d to %d...\n',sp,cfg_ana.out_param,count+1,count+splitSizes(sp));
                if isfield(cfg_ft,'outputfile')
                  cfg_ft = rmfield(cfg_ft,'outputfile');
                end
                cfg_ft.outputfile = outputFile_full;
                freq = ft_freqanalysis(cfg_ft,orig.(data_fn));
                fprintf('Done.\n');
                
                if sp == 1 && cfg_ana.checkSplitFileSizeForSaving
                  varinfo = whos('freq');
                  % convert to Megabytes
                  if (varinfo.bytes / (1024^2)) * length(splitSizes) > cfg_ana.combineSavingLimitMB
                    tooBigToCombine = true;
                  end
                end
                
              else
                if exist(outputFile_full,'file')
                  fprintf('%s file (split) already exist, not recalculating: %s\n',cfg_ana.out_param,outputFile_full);
                end
                if exist(outputFile_full_orig,'file')
                  fprintf('%s file (combined) already exist, not recalculating: %s\n',cfg_ana.out_param,outputFile_full_orig);
                end
              end
              
              % calculate power for this split
              if cfg_ana.output2alt
                outputFile_alt_full = strrep(outputFile_alt_full,origExt,newExt);
                
                % if it is not already calculated
                if ~exist(outputFile_alt_full,'file') && ~exist(outputFile_alt_full_orig,'file')
                  
                  % load the output file so we have data to work with
                  if ~exist('freq','var')
                    if exist(outputFile_full,'file')
                      load(outputFile_full);
                    else
                      if exist(outputFile_full_orig,'file')
                        fprintf('Combined %s file already exists. Need split files, so loading full file and splitting...\n',cfg_ana.out_param);
                        % reloading the file for each split may not be
                        % efficient, but it may be better than keeping a
                        % backup of the full struct in memory (untested)
                        load(outputFile_full_orig);
                        
                        cfg_sel = [];
                        cfg_sel.trials = false(1,nTrials);
                        cfg_sel.trials(count+1:count+splitSizes(sp)) = true;
                        freq = ft_selectdata_new(cfg_sel,freq);
                        
                        % % delete the file and re-calculate output; since
                        % % the files are going to get re-combined anyway,
                        % % seems easier and more efficient to start over
                        % fprintf('Deleting the combined file...\n');
                        % %delete(outputFile_full);
                        % [s] = unix(sprintf('rm %s',outputFile_full_orig));
                        % if s ~= 0
                        %   warning('Something went wrong when deleting %s',outputFile_full_orig);
                        % end
                        % fprintf('Done.\n');
                        %
                        % fprintf('\tSplit %d: Re-calculating %s for trials %d to %d...\n',sp,cfg_ana.out_param,count+1,count+splitSizes(sp));
                        % if isfield(cfg_ft,'outputfile')
                        %   cfg_ft = rmfield(cfg_ft,'outputfile');
                        % end
                        % cfg_ft.outputfile = outputFile_full;
                        % freq = ft_freqanalysis(cfg_ft,orig.(data_fn));
                      else
                        fprintf('\tSplit %d: Re-calculating %s for trials %d to %d...\n',sp,cfg_ana.out_param,count+1,count+splitSizes(sp));
                        if isfield(cfg_ft,'outputfile')
                          cfg_ft = rmfield(cfg_ft,'outputfile');
                        end
                        cfg_ft.outputfile = outputFile_full;
                        freq = ft_freqanalysis(cfg_ft,orig.(data_fn));
                      end
                      fprintf('Done.\n');
                    end
                  end
                  
                  fprintf('\tSplit %d: Calculating %s for trials %d to %d...\n',sp,cfg_ana.alt_param,count+1,count+splitSizes(sp));
                  if isfield(cfg_alt,'outputfile')
                    cfg_alt = rmfield(cfg_alt,'outputfile');
                  end
                  cfg_alt.outputfile = outputFile_alt_full;
                  %pow.(cfg_ana.alt_param) = (abs(freq.(cfg_ana.out_param))).^2;
                  %pow = ft_freqdescriptives(cfg_alt,freq);
                  ft_freqdescriptives(cfg_alt,freq);
                  fprintf('Done.\n');
                else
                  if exist(outputFile_alt_full,'file')
                    fprintf('%s file (split) already exist, not recalculating: %s\n',cfg_ana.alt_param,outputFile_alt_full);
                  end
                  if exist(outputFile_alt_full_orig,'file')
                    fprintf('%s file (combined) already exist, not recalculating: %s\n',cfg_ana.alt_param,outputFile_alt_full_orig);
                  end
                end
              end
              
              % increase the count for next time
              count = count + splitSizes(sp);
              
              % clear some memory
              clear freq
            end % for sp
            
            % combine the split files
            if ~tooBigToCombine
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % Main output (fourier)
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
              outputFile_full = outputFile_full_orig;
              
              if ~exist(outputFile_full,'file')
                % initialize to hold all files
                freq_all = struct;
                
                for sp = 1:length(splitSizes)
                  
                  if sp == 1
                    origExt = '.mat';
                    newExt = sprintf('_%d.mat',sp);
                    outputFile_full = strrep(outputFile_full,origExt,newExt);
                  elseif sp > 1
                    origExt = sprintf('_%d.mat',sp-1);
                    newExt = sprintf('_%d.mat',sp);
                    outputFile_full = strrep(outputFile_full,origExt,newExt);
                  end
                  
                  freq_all.(sprintf('freq_%d',sp)) = load(outputFile_full);
                  
                end % for sp
                
                % append the data together
                fn = fieldnames(freq_all);
                data_str = [];
                for sp = 1:length(fn)
                  data_str = cat(2,data_str,sprintf(',freq_all.%s.freq',fn{sp}));
                end
                
                cfg_af = [];
                cfg_af.parameter = cfg_ana.out_param;
                freq = eval(sprintf('ft_appendfreq(cfg_af%s)',data_str));
                
                % see if freq struct is too big to save in v7 format
                varinfo = whos('freq');
                if (varinfo.bytes / (1024^2)) > 2000
                  saveVersion = '-v7.3';
                else
                  saveVersion = '-v7';
                end
                
                % save it
                outputFile_full = outputFile_full_orig;
                save(outputFile_full,'freq',saveVersion);
                
                % clear some memory
                clear freq freq_all
              end
              
              % remove the split files
              for sp = 1:length(splitSizes)
                if sp == 1
                  origExt = '.mat';
                  newExt = sprintf('_%d.mat',sp);
                  outputFile_full = strrep(outputFile_full,origExt,newExt);
                elseif sp > 1
                  origExt = sprintf('_%d.mat',sp-1);
                  newExt = sprintf('_%d.mat',sp);
                  outputFile_full = strrep(outputFile_full,origExt,newExt);
                end
                
                if exist(outputFile_full,'file')
                  %delete(outputFile_full);
                  [s] = unix(sprintf('rm %s',outputFile_full));
                  if s ~= 0
                    warning('Something went wrong when deleting %s',outputFile_full);
                  end
                end
              end % remove sp
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % Alternative output (Power)
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
              % combine power files
              if cfg_ana.output2alt
                outputFile_alt_full = outputFile_alt_full_orig;
                
                if ~exist(outputFile_alt_full,'file')
                  % initialize to hold all files
                  freq_all = struct;
                  
                  for sp = 1:length(splitSizes)
                    
                    if sp == 1
                      origExt = '.mat';
                      newExt = sprintf('_%d.mat',sp);
                      outputFile_alt_full = strrep(outputFile_alt_full,origExt,newExt);
                    elseif sp > 1
                      origExt = sprintf('_%d.mat',sp-1);
                      newExt = sprintf('_%d.mat',sp);
                      outputFile_alt_full = strrep(outputFile_alt_full,origExt,newExt);
                    end
                    
                    freq_all.(sprintf('freq_%d',sp)) = load(outputFile_alt_full);
                    
                  end % for sp
                  
                  % append the data together
                  fn = fieldnames(freq_all);
                  data_str = [];
                  for sp = 1:length(fn)
                    data_str = cat(2,data_str,sprintf(',freq_all.%s.freq',fn{sp}));
                  end
                  
                  cfg_af = [];
                  cfg_af.parameter = cfg_ana.alt_param;
                  freq = eval(sprintf('ft_appendfreq(cfg_af%s)',data_str));
                  
                  % see if freq struct is too big to save in v7 format
                  varinfo = whos('freq');
                  if (varinfo.bytes / (1024^2)) > 2000
                    saveVersion = '-v7.3';
                  else
                    saveVersion = '-v7';
                  end
                  
                  % save it
                  outputFile_alt_full = outputFile_alt_full_orig;
                  save(outputFile_alt_full,'freq',saveVersion);
                  
                  % clear some memory
                  clear freq freq_all
                end
                
                % remove the split files
                for sp = 1:length(splitSizes)
                  if sp == 1
                    origExt = '.mat';
                    newExt = sprintf('_%d.mat',sp);
                    outputFile_alt_full = strrep(outputFile_alt_full,origExt,newExt);
                  elseif sp > 1
                    origExt = sprintf('_%d.mat',sp-1);
                    newExt = sprintf('_%d.mat',sp);
                    outputFile_alt_full = strrep(outputFile_alt_full,origExt,newExt);
                  end
                  
                  if exist(outputFile_alt_full,'file')
                    %delete(outputFile_alt_full);
                    [s] = unix(sprintf('rm %s',outputFile_alt_full));
                    if s ~= 0
                      warning('Something went wrong when deleting %s',outputFile_alt_full);
                    end
                  end
                end % remove sp
              end % if output2alt
            end % tooBigToCombine
            
          else
            fprintf('No need to split trials...\n');
            
            outputFile_full = outputFile_full_orig;
            
            % operate on all trials
            cfg_ft.trials = 'all';
            
            % do first conversion if it is not already calculated
            if ~exist(outputFile_full,'file')
              fprintf('\tCalculating %s...\n',cfg_ana.out_param);
              if isfield(cfg_ft,'outputfile')
                cfg_ft = rmfield(cfg_ft,'outputfile');
              end
              cfg_ft.outputfile = outputFile_full;
              freq = ft_freqanalysis(cfg_ft,orig.(data_fn));
              fprintf('Done.\n');
            else
              fprintf('%s file already exist, not recalculating: %s\n',cfg_ana.out_param,outputFile_full);
            end
            
            % calculate power
            if cfg_ana.output2alt
              outputFile_alt_full = outputFile_alt_full_orig;
              
              % if it is not already calculated
              if ~exist(outputFile_alt_full,'file')
                
                % load the output file so we have data to work with
                if ~exist('freq','var')
                  load(outputFile_full);
                end
                
                fprintf('\tCalculating %s...\n',cfg_ana.alt_param);
                if isfield(cfg_alt,'outputfile')
                  cfg_alt = rmfield(cfg_alt,'outputfile');
                end
                cfg_alt.outputfile = outputFile_alt_full;
                %pow.(cfg_ana.alt_param) = (abs(freq.(cfg_ana.out_param))).^2;
                %pow = ft_freqdescriptives(cfg_alt,freq);
                ft_freqdescriptives(cfg_alt,freq);
                fprintf('Done.\n');
              else
                fprintf('%s file already exist, not recalculating: %s\n',cfg_ana.alt_param,outputFile_alt_full);
              end
            end
            
            % clear some memory
            clear freq
          end % nTrials > splitSize
        else
          outputFile_full = outputFile_full_orig;
          
          % first conversion: operate on all trials
          cfg_ft.trials = 'all';
          
          % do the main output if it is not already calculated
          if ~exist(outputFile_full,'file')
            fprintf('\tCalculating %s...\n',cfg_ana.out_param);
            if isfield(cfg_ft,'outputfile')
              cfg_ft = rmfield(cfg_ft,'outputfile');
            end
            cfg_ft.outputfile = outputFile_full;
            freq = ft_freqanalysis(cfg_ft,orig.(data_fn));
            fprintf('Done.\n');
          else
            fprintf('%s file already exist, not recalculating: %s\n',cfg_ana.out_param,outputFile_full);
          end
          
          % calculate power
          if cfg_ana.output2alt
            outputFile_alt_full = outputFile_alt_full_orig;
            
            % if it is not already calculated
            if ~exist(outputFile_alt_full,'file')
              
              % TODO: check on this exist variable name
              
              % load the output file so we have data to work with
              if ~exist('freq','var')
                load(outputFile_full);
              end
              
              fprintf('\tCalculating %s...\n',cfg_ana.alt_param);
              if isfield(cfg_alt,'outputfile')
                cfg_alt = rmfield(cfg_alt,'outputfile');
              end
              cfg_alt.outputfile = outputFile_alt_full;
              %pow.(cfg_ana.alt_param) = (abs(freq.(cfg_ana.out_param))).^2;
              %pow = ft_freqdescriptives(cfg_alt,freq);
              ft_freqdescriptives(cfg_alt,freq);
              fprintf('Done.\n');
            else
              fprintf('%s file already exist, not recalculating: %s\n',cfg_ana.alt_param,outputFile_alt_full);
            end
          end
          
          % clear some memory
          clear freq
        end % split trials
        
        % clear some memory
        clear orig
      else
        fprintf('\tAll necessary data for %s already exists! Moving on...\n',origFile);
      end
      
      if cfg_ana.useLockFiles
        releaseFile(lockedFile);
      end
    end % for od
    
  end % ses
end % sub

end % function
