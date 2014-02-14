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
      
      % load timelock regardless of what output we have
      fprintf('\tLoading timelock...\n');
      load(origFile_full);
      fprintf('Done.\n');
        
      if cfg_ana.splitTrials
        % split the trials up if there are more trials than RAM can handle
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
            splitSizes = cat(2,splitSizes,remainSize);
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
            if ~exist(outputFile_full,'file')
              fprintf('\tSplit %d: Calculating fourier for trials %d to %d...\n',sp,count+1,count+splitSizes(sp));
              
              if isfield(cfg_ft,'outputfile')
                cfg_ft = rmfield(cfg_ft,'outputfile');
              end
              cfg_ft.outputfile = outputFile_full;
              freq = ft_freqanalysis(cfg_ft,timelock);
              fprintf('Done.\n');
            else
              fprintf('fourier file already exist, not recalculating: %s\n',outputFile_full);
            end
            
            % calculate power for this split
            if cfg_ana.fourier2pow
              outputFile_alt_full = strrep(outputFile_alt_full,origExt,newExt);
              
              % if it is not already calculated
              if ~exist(outputFile_alt_full,'file')
                
                % load the fourier file so we have data to work with
                if ~exist('freq','var')
                  load(outputFile_full);
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
                fprintf('power file already exist, not recalculating: %s\n',outputFile_alt_full);
              end
            end
            
            % increase the count for next time
            count = count + splitSizes(sp);
            
            % clear some memory
            clear freq
          end
        else
          fprintf('No need to split trials...\n');
          
          outputFile_full = outputFile_full_orig;
          if cfg_ana.fourier2pow
            outputFile_alt_full = outputFile_alt_full_orig;
          end
          
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
        if cfg_ana.fourier2pow
          outputFile_alt_full = outputFile_alt_full_orig;
        end
        
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
      
    end % for od
  end % ses
end % sub

end % function
