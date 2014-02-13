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
      outputFile_full = fullfile(saveDir,outputFile);
      
      outputFile_alt = strrep(origFile,cfg_ana.orig_ftype,cfg_ana.alt_ftype);
      outputFile_alt_full = fullfile(saveDir_alt,outputFile_alt);
      
      % load timelock
      %fprintf('\tLoading timelock...\n');
      %load(origFile_full);
      %fprintf('Done.\n');
      
      % calculate fourier
      fprintf('\tCalculating fourier...\n');
      cfg_ft.inputfile = origFile_full;
      cfg_ft.outputfile = outputFile_full;
      fourier = ft_freqanalysis(cfg_ft);
      %fourier = ft_freqanalysis(cfg_ft,timelock);
      fprintf('Done.\n');
      
      % save fourier
      %fprintf('\tSaving fourier...');
      %save(outputFile_full,'fourier');
      %fprintf('Done.\n');
      
      if cfg_ana.fourier2pow
        if cfg_ana.alt_splitTrials
          % split it in two if there are too many trials
          nTrials = size(fourier.fourierspctrm,1);
          if nTrials > 300
            fprintf('Splitting trials...\n');
            cfg_fd.trials = false(1,nTrials);
            firstHalf = floor(nTrials / 2);
            cfg_fd.trials(1:firstHalf) = true;
            outputFile_alt_full = strrep(outputFile_alt_full,'.mat','_1.mat');
            
            % calculate pow
            fprintf('\tCalculating power for first half of trials...\n');
            if isfield(cfg_fd,'outputfile')
              cfg_fd = rmfield(cfg_fd,'outputfile');
            end
            cfg_fd.outputfile = outputFile_alt_full;
            %pow.(cfg_ana.alt_param) = (abs(fourier.(cfg_ana.param))).^2;
            %pow = ft_freqdescriptives(cfg_fd,fourier);
            ft_freqdescriptives(cfg_fd,fourier);
            fprintf('Done.\n');
            
            
            cfg_fd.trials = false(1,nTrials);
            cfg_fd.trials(firstHalf+1:nTrials) = true;
            outputFile_alt_full = strrep(outputFile_alt_full,'_1.mat','_2.mat');
            
            % calculate pow
            fprintf('\tCalculating power for second half of trials...\n');
            if isfield(cfg_fd,'outputfile')
              cfg_fd = rmfield(cfg_fd,'outputfile');
            end
            cfg_fd.outputfile = outputFile_alt_full;
            %pow.(cfg_ana.alt_param) = (abs(fourier.(cfg_ana.param))).^2;
            %pow = ft_freqdescriptives(cfg_fd,fourier);
            ft_freqdescriptives(cfg_fd,fourier);
            fprintf('Done.\n');
          else
            cfg_fd.trials = 'all';
            
            % calculate pow
            fprintf('\tCalculating power...\n');
            if isfield(cfg_fd,'outputfile')
              cfg_fd = rmfield(cfg_fd,'outputfile');
            end
            cfg_fd.outputfile = outputFile_alt_full;
            %pow.(cfg_ana.alt_param) = (abs(fourier.(cfg_ana.param))).^2;
            %pow = ft_freqdescriptives(cfg_fd,fourier);
            ft_freqdescriptives(cfg_fd,fourier);
            fprintf('Done.\n');
          end
        else
          cfg_fd.trials = 'all';
          
          % calculate pow
          fprintf('\tCalculating power...\n');
          if isfield(cfg_fd,'outputfile')
            cfg_fd = rmfield(cfg_fd,'outputfile');
          end
          cfg_fd.outputfile = outputFile_alt_full;
          %pow.(cfg_ana.alt_param) = (abs(fourier.(cfg_ana.param))).^2;
          %pow = ft_freqdescriptives(cfg_fd,fourier);
          ft_freqdescriptives(cfg_fd,fourier);
          fprintf('Done.\n');
        end
        
        % save pow
        %fprintf('\tSaving power...\n');
        %save(outputFile_alt_full,'pow');
        %fprintf('Done.\n');
      end
      
      clear fourier
      
    end
  end
end

end