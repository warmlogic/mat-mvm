function mm_tla2fourier(exper,dirs,cfg_ana,cfg_ft,cfg_fd)

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
      
      if cfg_ana.fourier2pow
        outputFile_alt = strrep(origFile,cfg_ana.orig_ftype,cfg_ana.alt_ftype);
        outputFile_alt_full_orig = fullfile(saveDir_alt,outputFile_alt);
      end
      
      if ~exist(outputFile_full_orig,'file') || (cfg_ana.fourier2pow && ~exist(outputFile_alt_full_orig,'file'))
        % load timelock regardless of what output we have
        fprintf('\tLoading timelock: %s...\n',origFile);
        load(origFile_full);
        fprintf('Done.\n');
        
        if cfg_ana.splitTrials
          % split up trials if there are more than the RAM can handle
          nTrials = size(timelock.trial,1);
          if nTrials > cfg_ana.splitSize
            fprintf('Splitting trials...\n');
            
            outputFile_full = outputFile_full_orig;
            if cfg_ana.fourier2pow
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
              if remainSize > 10
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
              
              % calculate fourier if it is not already calculated
              if ~exist(outputFile_full,'file') && ~exist(outputFile_full_orig,'file')
                fprintf('\tSplit %d: Calculating fourier for trials %d to %d...\n',sp,count+1,count+splitSizes(sp));
                if isfield(cfg_ft,'outputfile')
                  cfg_ft = rmfield(cfg_ft,'outputfile');
                end
                cfg_ft.outputfile = outputFile_full;
                freq = ft_freqanalysis(cfg_ft,timelock);
                fprintf('Done.\n');
                
                %freq = runFreqanalysis(cfg_ft,timelock,outputFile_full);
              else
                if exist(outputFile_full,'file')
                  fprintf('fourier file already exist, not recalculating: %s\n',outputFile_full);
                end
                if exist(outputFile_full_orig,'file')
                  fprintf('fourier file already exist, not recalculating: %s\n',outputFile_full_orig);
                end
              end
              
              % calculate power for this split
              if cfg_ana.fourier2pow
                outputFile_alt_full = strrep(outputFile_alt_full,origExt,newExt);
                
                % if it is not already calculated
                if ~exist(outputFile_alt_full,'file') && ~exist(outputFile_alt_full_orig,'file')
                  
                  % load the fourier file so we have data to work with
                  if ~exist('freq','var')
                    if exist(outputFile_full,'file')
                      load(outputFile_full);
                    else
                      if exist(outputFile_full_orig,'file')
                        % fprintf('Combined fourier file already exists. Need split files, so loading full file and splitting...\n');
                        % % reloading the file for each split may not be
                        % % efficient, but it may be better than keeping a
                        % % backup of the full struct in memory (untested)
                        % load(outputFile_full_orig);
                        
                        % cfg_sel = [];
                        % cfg_sel.trials = false(1,nTrials);
                        % cfg_sel.trials(count+1:count+splitSizes(sp)) = true;
                        % freq = ft_selectdata_new(cfg_sel,freq);
                        
                        % delete the file and re-calculate fourier
                        fprintf('Deleting the combined file...\n');
                        %delete(outputFile_full);
                        [s] = unix(sprintf('rm %s',outputFile_full_orig));
                        if s ~= 0
                          warning('Something went wrong when deleting %s',outputFile_full_orig);
                        end
                        fprintf('Done.\n');
                        
                        fprintf('\tSplit %d: Re-calculating fourier for trials %d to %d...\n',sp,count+1,count+splitSizes(sp));
                        if isfield(cfg_ft,'outputfile')
                          cfg_ft = rmfield(cfg_ft,'outputfile');
                        end
                        cfg_ft.outputfile = outputFile_full;
                        freq = ft_freqanalysis(cfg_ft,timelock);
                        
                      else
                        fprintf('\tSplit %d: Re-calculating fourier for trials %d to %d...\n',sp,count+1,count+splitSizes(sp));
                        if isfield(cfg_ft,'outputfile')
                          cfg_ft = rmfield(cfg_ft,'outputfile');
                        end
                        cfg_ft.outputfile = outputFile_full;
                        freq = ft_freqanalysis(cfg_ft,timelock);
                      end
                      fprintf('Done.\n');
                    end
                  end
                  
                  fprintf('\tSplit %d: Calculating power for trials %d to %d...\n',sp,count+1,count+splitSizes(sp));
                  if isfield(cfg_fd,'outputfile')
                    cfg_fd = rmfield(cfg_fd,'outputfile');
                  end
                  cfg_fd.outputfile = outputFile_alt_full;
                  %pow.(cfg_ana.alt_param) = (abs(freq.(cfg_ana.param))).^2;
                  %pow = ft_freqdescriptives(cfg_fd,freq);
                  ft_freqdescriptives(cfg_fd,freq);
                  fprintf('Done.\n');
                else
                  if exist(outputFile_alt_full,'file')
                    fprintf('power file already exist, not recalculating: %s\n',outputFile_alt_full);
                  end
                  if exist(outputFile_alt_full_orig,'file')
                    fprintf('power file already exist, not recalculating: %s\n',outputFile_alt_full_orig);
                  end
                end
              end
              
              % increase the count for next time
              count = count + splitSizes(sp);
              
              % clear some memory
              clear freq
            end % for sp
            
            % combine the split files
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fourier
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
              cfg_af.parameter = 'fourierspctrm';
              freq = eval(sprintf('ft_appendfreq(cfg_af%s)',data_str));
              
              % save it
              outputFile_full = outputFile_full_orig;
              save(outputFile_full,'freq');
              
              % clear some memory
              clear freq freq_all
              
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
                
                %delete(outputFile_full);
                [s] = unix(sprintf('rm %s',outputFile_full));
                if s ~= 0
                  warning('Something went wrong when deleting %s',outputFile_full);
                end
              end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Power
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % combine power files
            if cfg_ana.fourier2pow
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
                
                % save it
                outputFile_alt_full = outputFile_alt_full_orig;
                save(outputFile_alt_full,'freq');
                
                % clear some memory
                clear freq freq_all
                
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
                  
                  %delete(outputFile_alt_full);
                  [s] = unix(sprintf('rm %s',outputFile_alt_full));
                  if s ~= 0
                    warning('Something went wrong when deleting %s',outputFile_alt_full);
                  end
                end
              end
            end
            
          else
            fprintf('No need to split trials...\n');
            
            outputFile_full = outputFile_full_orig;
            
            % fourier: operate on all trials
            cfg_ft.trials = 'all';
            
            % calculate fourier if it is not already calculated
            if ~exist(outputFile_full,'file')
              fprintf('\tCalculating fourier...\n');
              if isfield(cfg_ft,'outputfile')
                cfg_ft = rmfield(cfg_ft,'outputfile');
              end
              cfg_ft.outputfile = outputFile_full;
              freq = ft_freqanalysis(cfg_ft,timelock);
              fprintf('Done.\n');
            else
              fprintf('fourier file already exist, not recalculating: %s\n',outputFile_full);
            end
            
            % calculate power
            if cfg_ana.fourier2pow
              outputFile_alt_full = outputFile_alt_full_orig;
              
              % if it is not already calculated
              if ~exist(outputFile_alt_full,'file')
                
                % load the fourier file so we have data to work with
                if ~exist('freq','var')
                  load(outputFile_full);
                end
                
                fprintf('\tCalculating power...\n');
                if isfield(cfg_fd,'outputfile')
                  cfg_fd = rmfield(cfg_fd,'outputfile');
                end
                cfg_fd.outputfile = outputFile_alt_full;
                %pow.(cfg_ana.alt_param) = (abs(freq.(cfg_ana.param))).^2;
                %pow = ft_freqdescriptives(cfg_fd,freq);
                ft_freqdescriptives(cfg_fd,freq);
                fprintf('Done.\n');
              else
                fprintf('power file already exist, not recalculating: %s\n',outputFile_alt_full);
              end
            end
            
            % clear some memory
            clear freq
          end % nTrials > splitSize
        else
          outputFile_full = outputFile_full_orig;
          
          % fourier: operate on all trials
          cfg_ft.trials = 'all';
          
          % calculate fourier if it is not already calculated
          if ~exist(outputFile_full,'file')
            fprintf('\tCalculating fourier...\n');
            if isfield(cfg_ft,'outputfile')
              cfg_ft = rmfield(cfg_ft,'outputfile');
            end
            cfg_ft.outputfile = outputFile_full;
            freq = ft_freqanalysis(cfg_ft,timelock);
            fprintf('Done.\n');
          else
            fprintf('fourier file already exist, not recalculating: %s\n',outputFile_full);
          end
          
          % calculate power
          if cfg_ana.fourier2pow
            outputFile_alt_full = outputFile_alt_full_orig;
            
            % if it is not already calculated
            if ~exist(outputFile_alt_full,'file')
              
              % load the fourier file so we have data to work with
              if ~exist('fourier','var')
                load(outputFile_full);
              end
              
              fprintf('\tCalculating power...\n');
              if isfield(cfg_fd,'outputfile')
                cfg_fd = rmfield(cfg_fd,'outputfile');
              end
              cfg_fd.outputfile = outputFile_alt_full;
              %pow.(cfg_ana.alt_param) = (abs(freq.(cfg_ana.param))).^2;
              %pow = ft_freqdescriptives(cfg_fd,freq);
              ft_freqdescriptives(cfg_fd,freq);
              fprintf('Done.\n');
            else
              fprintf('power file already exist, not recalculating: %s\n',outputFile_alt_full);
            end
          end
          
          % clear some memory
          clear freq
        end % split trials
        
        % clear some memory
        clear timelock
      else
        fprintf('\tAll necessary data for %s already exists! Moving on...\n',origFile);
      end
      
    end % for od
  end % ses
end % sub

end % function mm_tla2fourier
