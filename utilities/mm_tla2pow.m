function mm_tla2pow(exper,dirs,cfg_ana)

% function mm_tla2pow(exper,dirs,cfg_ana)
%
% Take timelock data and convert it into power data.
%

if nargin ~= 3
  error('Not the correct number of input arguments!');
end

if ~isfield(cfg_ana,'method') 
  error('Must set cfg_ana.method to ''wavelet_ndtools'' or ''wavelet_ants''.');
elseif isfield(cfg_ana,'method')
  if isempty(cfg_ana.method)
    error('Must set cfg_ana.method to ''wavelet_ndtools'' or ''wavelet_ants''.');
  elseif ~ismember(cfg_ana.method,{'wavelet_ndtools','wavelet_ants'})
    error('Must set cfg_ana.method to ''wavelet_ndtools'' or ''wavelet_ants''. You set it to: ''%s''.',cfg_ana.method);
  elseif strcmp(cfg_ana.method,'wavelet_ants')
    if isfield(cfg_ana,'wavelet_cycles') && ~isempty(cfg_ana.wavelet_cycles)
      warning('cfg_ana.wavelet_cycles is only used with cfg_ana.method=''wavelet_ants''. Removing this field.');
      cfg_ana = rmfield(cfg_ana,'wavelet_cycles');
    end
  elseif strcmp(cfg_ana.method,'wavelet_ndtools')
    warning('Automatically using wavelet cycles = 6 because it is hardcoded in NDTools.');
  end
end

if ~isfield(cfg_ana,'orig_ftype')
  cfg_ana.orig_ftype = 'tla';
end
if ~isfield(cfg_ana,'orig_param')
  cfg_ana.orig_param = 'trial';
end

if ~isfield(cfg_ana,'sampleRate')
  error('need to set cfg_ana.sampleRate for the sample rate of original data');
end

% not recommended to resample before calculating power because I don't know
% how it will change the data
if ~isfield(cfg_ana,'resample_tla')
  cfg_ana.resample_tla = false;
else
  if cfg_ana.resample_tla && ~isfield(cfg_ana,'resampleRate_tla')
    error('cfg_ana.resample_tla=true but cfg_ana.resampleRate_tla is not set.');
  end
end

% it might be a good idea to make your resampleRate_pow a multiple of your
% original sampleRate (e.g., resampleRate_pow = sampleRate/5)
if ~isfield(cfg_ana,'resample_pow')
  cfg_ana.resample_pow = false;
else
  if cfg_ana.resample_pow && ~isfield(cfg_ana,'resampleRate_pow')
    error('cfg_ana.resample_pow=true but cfg_ana.resampleRate_pow is not set.');
  elseif cfg_ana.resample_pow && isfield(cfg_ana,'resampleRate_pow') && cfg_ana.resampleRate_pow == cfg_ana.sampleRate
    warning('cfg_ana.resampleRate_pow and cfg_ana.sampleRate are the same. Not resampling.');
    cfg_ana.resample_pow = false;
  elseif cfg_ana.resample_pow && isfield(cfg_ana,'resampleRate_pow') && cfg_ana.resampleRate_pow > cfg_ana.sampleRate
    error('Cannot set cfg_ana.resampleRate_pow higher than cfg_ana.sampleRate.');
  end
end

if ~isfield(cfg_ana,'param')
  cfg_ana.out_param = 'powspctrm';
end

if ~isfield(cfg_ana,'output_ftype')
  cfg_ana.output_ftype = 'pow';
end

if ~isfield(cfg_ana,'keepTimeSec')
  cfg_ana.keepTimeSec = [];
else
  if length(cfg_ana.keepTimeSec) ~= 2
    error('cfg_ana.keepTimeSec not set correctly. Set start and end times in seconds');
  end
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
    cfg_ana.splitRemainderLump = 50;
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
    
    for od = 1:length(origData)
      % set up file names
      origFile = origData(od).name;
      origFile_full = fullfile(sesDir,origFile);
      
      outputFile = strrep(origFile,cfg_ana.orig_ftype,cfg_ana.output_ftype);
      outputFile_full_orig = fullfile(saveDir,outputFile);
      
      if cfg_ana.useLockFiles
        lockedFile = outputFile_full_orig;
        
        % file could be in limbo; see if it's already been started
        if ~lockFile(lockedFile)
          fprintf('Moving on, final file already exists: %s\n',lockedFile);
          continue
        end
      end
      
      if ~exist(outputFile_full_orig,'file')
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
        
        if cfg_ana.resample_tla
          cfg_rs = [];
          cfg_rs.resamplefs = cfg_ana.resampleRate_tla;
          % cfg_rs.time = repmat({-0.5:0.04:1.0},size(orig.(data_fn).trial,1),1);
          % cfg_rs.method = 'pchip'; % see help interp1
          cfg_rs.detrend = 'no';
          orig.(data_fn) = ft_resampledata(cfg_rs,orig.(data_fn));
          
          sampleRate = cfg_ana.resampleRate_tla;
        else
          sampleRate = cfg_ana.sampleRate;
        end
        
        % other wavelet parameters
        if strcmp(cfg_ana.freq_spacing,'log')
          frequencies = logspace(log10(cfg_ana.min_freq),log10(cfg_ana.max_freq),cfg_ana.num_freq);
        elseif strcmp(cfg_ana.freq_spacing,'lin')
          frequencies = linspace(cfg_ana.min_freq,cfg_ana.max_freq,cfg_ana.num_freq);
        else
          error('Need to set frequency spacing (''log'' or ''lin'').');
        end
        
        % set the number of samples
        n_samples = size(orig.(data_fn).trial,3);
        
        if strcmp(cfg_ana.method,'wavelet_ndtools')
          wobj = mkwobj('morl', n_samples, sampleRate, 1./frequencies);
        elseif strcmp(cfg_ana.method,'wavelet_ants')
          time = -1:1/sampleRate:1;
          half_of_wavelet_size = (length(time)-1)/2;
          if mod(half_of_wavelet_size,1) ~= 0
            error('half_of_wavelet_size=%.2f. Need to make it an integer. (This has to do with cfg_ana.sampleRate.)',half_of_wavelet_size);
          end
          
          % FFT parameters (use next-power-of-2)
          n_wavelet     = length(time);
          n_convolution = n_wavelet+n_samples-1;
          n_conv_pow2   = pow2(nextpow2(n_convolution));
          
          % create wavelet bank
          fft_wavelet = zeros(length(frequencies),n_conv_pow2);
          for fi=1:length(frequencies)
            % create wavelet and get its FFT
            wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( cfg_ana.wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
            fft_wavelet(fi,:) = fft(wavelet,n_conv_pow2);
          end
        end
        
        % find out about resampling power
        keepTime = true(1,n_samples);
        if cfg_ana.resample_pow
          keepTime = false(1,n_samples);
          keepTimeInd = [];
          if ~isempty(cfg_ana.keepTimeSec)
            timeLimits = cfg_ana.keepTimeSec;
          else
            timeLimits = [orig.(data_fn).time(1) orig.(data_fn).time(end)];
          end
          
          findTheseTimes = timeLimits(1):1/cfg_ana.resampleRate_pow:timeLimits(end);
          
          % find the nearest samples
          for i = 1:length(findTheseTimes)
            keepTimeInd = cat(2,keepTimeInd,nearest(orig.(data_fn).time,findTheseTimes(i)));
          end
          keepTime(keepTimeInd) = true;
        elseif ~cfg_ana.resample_pow && ~isempty(cfg_ana.keepTimeSec)
          keepTime = false(1,n_samples);
          tbeg = nearest(orig.(data_fn).time,cfg_ana.keepTimeSec(1));
          tend = nearest(orig.(data_fn).time,cfg_ana.keepTimeSec(2));
          keepTime(tbeg:tend) = true;
        end
        
        if cfg_ana.splitTrials
          % split up trials if there are more than the RAM can handle
          nTrials = size(orig.(data_fn).(cfg_ana.orig_param),1);
          if nTrials > (cfg_ana.splitSize + cfg_ana.splitRemainderLump)
            fprintf('%d trials found. Splitting trials...\n',nTrials);
            
            tooBigToCombine = false;
            
            outputFile_full = outputFile_full_orig;
            
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
              theseTrials = false(1,nTrials);
              
              theseTrials(count+1:count+splitSizes(sp)) = true;
              
              thisData = orig.(data_fn).trial(theseTrials,:,:);

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
                
                % initialize output time-frequency data: trials, channels, freq, time
                tf_data = zeros(size(thisData,1),size(orig.(data_fn).trial,2),length(frequencies),sum(keepTime));
                
                for t = 1:size(tf_data,1)
                  fprintf('.');
                  
                  if strcmp(cfg_ana.method,'wavelet_ndtools')
                    % can process all channels inside qwave
                    q = qwave(squeeze(thisData(t,:,:)), wobj);
                    tf_data(t,:,:,:) = abs(q.cfs(:,:,keepTime)).^2;
                  elseif strcmp(cfg_ana.method,'wavelet_ants')
                    for c = 1:size(tf_data,2)
                      % get each channel for each trial
                      
                      % result of squeeze(thisData(t,c,:)) needs to be a row vector
                      
                      % fft_data = repmat(fft(squeeze(thisData(t,c,:))',n_conv_pow2),size(fft_wavelet,1),1);
                      % convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2,2);
                      
                      % consolidate commands for speed
                      convolution_result_fft = ifft(fft_wavelet.*repmat(fft(squeeze(thisData(t,c,:))',n_conv_pow2),size(fft_wavelet,1),1),n_conv_pow2,2);
                      
                      convolution_result_fft = convolution_result_fft(:,1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
                      convolution_result_fft = convolution_result_fft(:,half_of_wavelet_size+1:end-half_of_wavelet_size);
                      tf_data(t,c,:,:) = abs(convolution_result_fft(keepTime)).^2;
                      
%                       % old slow for loop over frequencies
%                       % FFT of data (note: this doesn't change on frequency iteration)
%                       % result of squeeze(thisData(t,c,:)) needs to be a row vector
%                       fft_data = fft(squeeze(thisData(t,c,:))',n_conv_pow2);
%                       
%                       for fi=1:length(frequencies)
%                         % run convolution
%                         convolution_result_fft = ifft(fft_wavelet(fi,:).*fft_data,n_conv_pow2);
%                         convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
%                         convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
%                         
%                         % put power data into time-frequency matrix
%                         tf_data(t,c,fi,:) = abs(convolution_result_fft(keepTime)).^2;
%                       end
                    end
                  end
                end
                fprintf('\n');
                
                freq = orig.(data_fn);
                if isfield(freq,'avg')
                  freq = rmfield(freq,'avg');
                end
                if isfield(freq,'var')
                  freq = rmfield(freq,'var');
                end
                if isfield(freq,'dof')
                  freq = rmfield(freq,'dof');
                end
                if isfield(freq,'elec')
                  freq = rmfield(freq,'elec');
                end
                if isfield(freq,'trial')
                  freq = rmfield(freq,'trial');
                end
                if ~isfield(freq,'fsample')
                  if cfg_ana.resample_pow
                    freq.fsample = cfg_ana.resampleRate_pow;
                  else
                    freq.fsample = sampleRate;
                  end
                end
                freq.dimord = 'rpt_chan_freq_time';
                if strcmp(cfg_ana.method,'wavelet_ndtools')
                  % store the center frequencies
                  freq.freq = wobj.cF;
                elseif strcmp(cfg_ana.method,'wavelet_ants')
                  freq.freq = frequencies;
                end
                freq.time = freq.time(keepTime);
                freq.(cfg_ana.out_param) = tf_data;
                if isfield(freq,'trialinfo')
                  freq.trialinfo = freq.trialinfo(theseTrials,:);
                end
                if isfield(freq,'sampleinfo')
                  freq.sampleinfo = freq.sampleinfo(theseTrials,:);
                end
                
                varinfo = whos('freq');
                if (varinfo.bytes / (1024^2)) > 2000
                  saveVersion = '-v7.3';
                else
                  saveVersion = '-v7';
                end
                
                save(outputFile_full,'freq',saveVersion);
                
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
              
              % increase the count for next time
              count = count + splitSizes(sp);
              
              % clear some memory
              clear freq
            end % for sp
            
            % combine the split files
            if ~tooBigToCombine
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % Output power
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
              
            end % tooBigToCombine
            
          else
            fprintf('No need to split trials...\n');
            
            outputFile_full = outputFile_full_orig;
            
            % do first conversion if it is not already calculated
            if ~exist(outputFile_full,'file')
              fprintf('\tCalculating %s...\n',cfg_ana.out_param);
              
              thisData = orig.(data_fn).trial;
              
              % initialize output time-frequency data: trials, channels, freq, time
              tf_data = zeros(size(orig.(data_fn).trial,1),size(orig.(data_fn).trial,2),length(frequencies),sum(keepTime));
              
              for t = 1:size(tf_data,1)
                fprintf('.');
                
                if strcmp(cfg_ana.method,'wavelet_ndtools')
                  % can process all channels inside qwave
                  q = qwave(squeeze(thisData(t,:,:)), wobj);
                  tf_data(t,:,:,:) = abs(q.cfs(:,:,keepTime)).^2;
                elseif strcmp(cfg_ana.method,'wavelet_ants')
                  for c = 1:size(tf_data,2)
                    % get each channel for each trial
                    
                    % result of squeeze(thisData(t,c,:)) needs to be a row vector
                    
                    % fft_data = repmat(fft(squeeze(thisData(t,c,:))',n_conv_pow2),size(fft_wavelet,1),1);
                    % convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2,2);
                    
                    % consolidate commands for speed
                    convolution_result_fft = ifft(fft_wavelet.*repmat(fft(squeeze(thisData(t,c,:))',n_conv_pow2),size(fft_wavelet,1),1),n_conv_pow2,2);
                    
                    convolution_result_fft = convolution_result_fft(:,1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
                    convolution_result_fft = convolution_result_fft(:,half_of_wavelet_size+1:end-half_of_wavelet_size);
                    tf_data(t,c,:,:) = abs(convolution_result_fft(keepTime)).^2;
                    
%                     % old slow for loop over frequencies
%                     % FFT of data (note: this doesn't change on frequency iteration)
%                     % result of squeeze needs to be a row vector
%                     fft_data = fft(squeeze(orig.(data_fn).trial(t,c,:))',n_conv_pow2);
%                     
%                     for fi=1:length(frequencies)
%                       % run convolution
%                       convolution_result_fft = ifft(fft_wavelet(fi,:).*fft_data,n_conv_pow2);
%                       convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
%                       convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
%                       
%                       % put power data into time-frequency matrix
%                       tf_data(t,c,fi,:) = abs(convolution_result_fft(keepTime)).^2;
%                     end
                  end
                end
              end
              fprintf('\n');
              
              freq = orig.(data_fn);
              if isfield(freq,'avg')
                freq = rmfield(freq,'avg');
              end
              if isfield(freq,'var')
                freq = rmfield(freq,'var');
              end
              if isfield(freq,'dof')
                freq = rmfield(freq,'dof');
              end
              if isfield(freq,'elec')
                freq = rmfield(freq,'elec');
              end
              if isfield(freq,'trial')
                freq = rmfield(freq,'trial');
              end
              if ~isfield(freq,'fsample')
                if cfg_ana.resample_pow
                  freq.fsample = cfg_ana.resampleRate_pow;
                else
                  freq.fsample = sampleRate;
                end
              end
              freq.dimord = 'rpt_chan_freq_time';
              if strcmp(cfg_ana.method,'wavelet_ndtools')
                % store the center frequencies
                freq.freq = wobj.cF;
              elseif strcmp(cfg_ana.method,'wavelet_ants')
                freq.freq = frequencies;
              end
              freq.time = freq.time(keepTime);
              freq.(cfg_ana.out_param) = tf_data;
              
              varinfo = whos('freq');
              if (varinfo.bytes / (1024^2)) > 2000
                saveVersion = '-v7.3';
              else
                saveVersion = '-v7';
              end
              
              save(outputFile_full,'freq',saveVersion);

              fprintf('Done.\n');
            else
              fprintf('%s file already exist, not recalculating: %s\n',cfg_ana.out_param,outputFile_full);
            end
            
            % clear some memory
            clear freq
          end % nTrials > splitSize
        else
          outputFile_full = outputFile_full_orig;
          
          % do the main output if it is not already calculated
          if ~exist(outputFile_full,'file')
            fprintf('\tCalculating %s...\n',cfg_ana.out_param);
            
            thisData = orig.(data_fn).trial;
            
            % initialize output time-frequency data: trials, channels, freq, time
            tf_data = zeros(size(orig.(data_fn).trial,1),size(orig.(data_fn).trial,2),length(frequencies),sum(keepTime));
            
            for t = 1:size(tf_data,1)
              fprintf('.');
              
              if strcmp(cfg_ana.method,'wavelet_ndtools')
                % can process all channels inside qwave
                q = qwave(squeeze(thisData(t,:,:)), wobj);
                tf_data(t,:,:,:) = abs(q.cfs(:,:,keepTime)).^2;
              elseif strcmp(cfg_ana.method,'wavelet_ants')
                for c = 1:size(tf_data,2)
                  % get each channel for each trial
                  
                  % result of squeeze(thisData(t,c,:)) needs to be a row vector
                  
                  % fft_data = repmat(fft(squeeze(thisData(t,c,:))',n_conv_pow2),size(fft_wavelet,1),1);
                  % convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2,2);
                  
                  % consolidate commands for speed
                  convolution_result_fft = ifft(fft_wavelet.*repmat(fft(squeeze(thisData(t,c,:))',n_conv_pow2),size(fft_wavelet,1),1),n_conv_pow2,2);
                  
                  convolution_result_fft = convolution_result_fft(:,1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
                  convolution_result_fft = convolution_result_fft(:,half_of_wavelet_size+1:end-half_of_wavelet_size);
                  tf_data(t,c,:,:) = abs(convolution_result_fft(keepTime)).^2;
                  
%                   % old slow for loop over frequencies
%                   % FFT of data (note: this doesn't change on frequency iteration)
%                   % result of squeeze needs to be a row vector
%                   fft_data = fft(squeeze(orig.(data_fn).trial(t,c,:))',n_conv_pow2);
%                   
%                   for fi=1:length(frequencies)
%                     % run convolution
%                     convolution_result_fft = ifft(fft_wavelet(fi,:).*fft_data,n_conv_pow2);
%                     convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
%                     convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
%                     
%                     % put power data into time-frequency matrix
%                     tf_data(t,c,fi,:) = abs(convolution_result_fft(keepTime)).^2;
%                   end
                end
              end
            end
            fprintf('\n');
            
            freq = orig.(data_fn);
            if isfield(freq,'avg')
              freq = rmfield(freq,'avg');
            end
            if isfield(freq,'var')
              freq = rmfield(freq,'var');
            end
            if isfield(freq,'dof')
              freq = rmfield(freq,'dof');
            end
            if isfield(freq,'elec')
              freq = rmfield(freq,'elec');
            end
            if isfield(freq,'trial')
              freq = rmfield(freq,'trial');
            end
            if ~isfield(freq,'fsample')
              if cfg_ana.resample_pow
                freq.fsample = cfg_ana.resampleRate_pow;
              else
                freq.fsample = sampleRate;
              end
            end
            freq.dimord = 'rpt_chan_freq_time';
            if strcmp(cfg_ana.method,'wavelet_ndtools')
              % store the center frequencies
              freq.freq = wobj.cF;
            elseif strcmp(cfg_ana.method,'wavelet_ants')
              freq.freq = frequencies;
            end
            freq.time = freq.time(keepTime);
            freq.(cfg_ana.out_param) = tf_data;
            
            varinfo = whos('freq');
            if (varinfo.bytes / (1024^2)) > 2000
              saveVersion = '-v7.3';
            else
              saveVersion = '-v7';
            end
            
            save(outputFile_full,'freq',saveVersion);
            
            fprintf('Done.\n');
          else
            fprintf('%s file already exist, not recalculating: %s\n',cfg_ana.out_param,outputFile_full);
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
