function cfg = mm_freqcheck(exper,cfg,baseline,plotit)
%MM_FREQCHECK - Checks your freq cfg for ft_freqanalysis
%
% cfg = mm_freqcheck(exper,cfg,baseline,plotit)
%
% Input:
%   exper      = exper struct; only needs prepost and sampleRate fields
%   cfg        = cfg that will go into FT_FREQANALYSIS; only needs fields:
%                method, foi (or foilim), toi, and output fields, plus
%                fields relevant to the method (e.g., width for wavelet)
%   baseline   = e.g., [-1.0 -0.75]
%   plotit     = whether to plot some info (default: true)
%
% Output:
%   cfg        = same cfg with some parameters/defaults set if you didn't
%                set them (optional)
%
% NB: I think this is correct, but the FieldTrip tutorials do not pay such
%     close attention to baseline periods. This is at least a conservative
%     approach.
%
% NB: if you just want to figure out how much time a given wavelet will
%     occupy, the equation is
%     (sample rate * width) / 2;
%     so: ((1 second / frequency of interest) * width) / 2;
%
% See also: FT_FREQANALYSIS
%

% check on the input
if ~exist('exper','var') || isempty(exper)
  error('Must set exper input.');
end
if ~exist('cfg','var') || isempty(cfg)
  error('Must set cfg input.');
end
if ~exist('baseline','var') || isempty(baseline)
  error('Must set baseline input.');
end
if ~exist('plotit','var') || isempty(plotit)
  plotit = true;
end

if ~isfield(cfg,'method')
  error('Must set cfg.methd (e.g., ''wavelet'' or ''mtmconvol'').');
end

% check on trial padding
if ~isfield(cfg,'pad')
  cfg.pad = 'maxperlen';
end
if strcmp(cfg.pad, 'maxperlen')
  padding = diff(exper.prepost) * exper.sampleRate;
  cfg.pad = padding / exper.sampleRate;
else
  padding = cfg.pad * exper.sampleRate;
  if padding < diff(exper.prepost)
    error('the specified padding is too short');
  end
end

% check on the frequencies of interest (foi or foilim)
if ~isfield(cfg,'foi') && ~isfield(cfg,'foilim')
  error('Must set cfg.foi (e.g., 4:1:8) or cfg.foilim (e.g., [4 8]).');
elseif isfield(cfg,'foi') && isfield(cfg,'foilim') && ~isempty(cfg.foi) && ~isempty(cfg.foilim)
  error('Do not set both cfg.foi and cfg.foilim.');
elseif isfield(cfg,'foi') && ~isfield(cfg,'foilim')
  cfg.foilim = [];
elseif ~isfield(cfg,'foi') && isfield(cfg,'foilim')
  cfg.foi = [];
end

% check on the times of interest (toi)
if ~isfield(cfg,'toi')
  error('Must set cfg.toi (e.g., -1.0:0.04:2.0). Every 40 or 50 ms is good.');
end

% modify foi, if necessary
if isempty(cfg.foi) && ~isempty(cfg.foilim)
  % % fieldtrip default step size
  % freqstep = (exper.sampleRate/(cfg.pad*exper.sampleRate));
  
  % double the step size
  freqstep = (exper.sampleRate/(cfg.pad*exper.sampleRate)) * 2;
  
  fprintf('You set cfg.foilim=[%.1f %.1f]. Removing this field and setting cfg.foi=%.1f:%.2f:%.1f;\n',cfg.foilim(1),cfg.foilim(2),cfg.foilim(1),freqstep,cfg.foilim(2));
  cfg.foi = cfg.foilim(1):freqstep:cfg.foilim(2);
  cfg = rmfield(cfg,'foilim');
elseif ~isempty(cfg.foi) && isempty(cfg.foilim)
  oldfoi = cfg.foi;
  fboi = round(cfg.foi ./ (exper.sampleRate ./ (cfg.pad * exper.sampleRate))) + 1;
  cfg.foi = (fboi-1) ./ cfg.pad; % boi - 1 because 0 Hz is included in fourier output
  if isfield(cfg,'t_ftimwin') && ~isempty(cfg.t_ftimwin) && strcmp(cfg.correctt_ftimwin,'yes')
    cyclenum = oldfoi .* cfg.t_ftimwin;
    cfg.t_ftimwin = cyclenum ./ cfg.foi;
  end
end

if ~isfield(cfg,'output')
  error('Must set set cfg.output (''pow'', ''powandcsd'', or ''fourier'').');
  %cfg.output = 'pow';
  %cfg.output = 'powandcsd';
  %cfg.output = 'fourier';
end

if strcmp(cfg.output,'pow') && ~isfield(cfg,'channel')
  %fprintf('You set cfg.output=''%s'', setting cfg.channel=''all''.\n',cfg.output);
  cfg.channel = 'all';
elseif strcmp(cfg.output,'powandcsd') && ~isfield(cfg,'channelcmb')
  %fprintf('You set cfg.output=''%s'', setting cfg.channelcmb={''all'',''all''}.\n',cfg.output);
  cfg.channelcmb = {'all','all'};
end

% % some test parameters for debugging

% exper = [];
% exper.prepost = [-1 1.25];
% exper.sampleRate = 500;

% % wavelet
% cfg = [];
% % cfg.foi = 4:1:8;
% cfg.foilim = [4 8];
% cfg.toi = -1.0:0.04:1.25;
% baseline = [-0.3 -0.1];
% cfg.method = 'wavelet';
% cfg.width = 6;
% cfg.output = 'pow';

% % multitaper - hanning taper
% cfg = [];
% % cfg.foi = 2:2:30;
% cfg.foilim = [2 30];
% cfg.toi = -1.0:0.04:2.0;
% baseline = [-0.3 -0.1];
% cfg.method = 'mtmconvol';
% cfg.taper = 'hanning';
% cfg.t_ftimwin = 6 ./ cfg.foi;
% cfg.output = 'pow';

% % multitaper - slepian tapers
% cfg = [];
% % cfg.foi = 2:2:30;
% cfg.foilim = [2 30];
% cfg.toi = -1.0:0.04:2.0;
% baseline = [-0.3 -0.1];
% cfg.method = 'mtmconvol';
% cfg.taper = 'dpss';
% cfg.t_ftimwin = 6 ./ cfg.foi;
% cfg.tapsmofrq = 0.4 * cfg.foi;
% cfg.output = 'pow';

fprintf('\nUsing the ''%s'' method with ''%s'' output.\n',cfg.method,cfg.output);

fprintf('ft_freqanalysis will process:\n');
fprintf('Time: %.2f to %.2f s at every %.1f ms\n',cfg.toi(1),cfg.toi(end),diff([cfg.toi(1) cfg.toi(2)]) * 1000);
fprintf('Freq: %.2f to %.2f Hz\n',cfg.foi(1),cfg.foi(end));

if strcmp(cfg.pad,'maxperlen')
  pfacs = factor(length(cfg.toi));
  pfacs = pfacs(isprime(pfacs));
  fprintf('Number of samples = %d. Prime factor sum = %d (factors =%s).\n',length(cfg.toi),sum(pfacs),sprintf(repmat(' %d',1,length(pfacs)),pfacs));
  fprintf('The FFTs could potentially run faster if your number of samples was a power of 2 (e.g., 64, 128).\n');
  % TODO
  fprintf('I need to actually check on how this works in FieldTrip, so don''t put too much stock in it.\n');
end

if strcmp(cfg.method,'wavelet')
  % wavelet
  if ~isfield(cfg,'width')
    fprintf('Setting cfg.width=6.\n');
    cfg.width = 6;
  end
  fprintf('Using wavelets with %d cycles.\n',cfg.width);
  wavelet_widths_s = (1 ./ cfg.foi * cfg.width);
  min_baseline_s = wavelet_widths_s ./ 2;
  
  %fprintf('For these frequencies, these are your wavelet widths, and the baseline period must end at least this many seconds before stimulus onset:\n');
  %fprintf('Freqs%s\n',sprintf(repmat('\t%.1f',1,length(cfg.foi)),cfg.foi));
  %fprintf('Width%s\n',sprintf(repmat('\t%.2f',1,length(wavelet_widths_s)),wavelet_widths_s));
  %fprintf('BL-end%s\n',sprintf(repmat('\t%.2f',1,length(min_baseline_s)),min_baseline_s));
  
  if plotit
    figure;
    subplot(2,2,1)
    plot(cfg.foi,wavelet_widths_s,'b*');
    xlabel('Freq (Hz)');
    ylabel('Wavelet widths (s)');
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    ylim([0 wavelet_widths_s(1)]);
    
    subplot(2,2,2)
    plot(cfg.foi,wavelet_widths_s ./ 2,'b*');
    xlabel('Freq (Hz)');
    ylabel('BL end before stim (s): wavelet\_widths ./ 2');
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    ylim([0 wavelet_widths_s(1)]);
    
    subplot(2,2,3)
    plot(cfg.foi,cfg.width,'b*');
    xlabel('Freq (Hz)');
    ylabel('Cycles per wavelet: cfg.width');
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    
    subplot(2,2,4)
    plot(cfg.foi,cfg.toi(end) - min_baseline_s,'b*');
    xlabel('Freq (Hz)');
    ylabel(sprintf('%s calculated out to (s)',cfg.output));
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    
    publishfig(gcf,0,12,16);
  end
  
elseif strcmp(cfg.method,'mtmconvol')
  % multitaper method
  
  if ~isfield(cfg,'t_ftimwin')
    fprintf('Setting frequency-dependent time window (cfg.t_ftimwin; temporal smoothing) to 6./cfg.foi (6 cycles per time window).\n');
    cfg.t_ftimwin = 6 ./ cfg.foi;
  end
  
  min_baseline_s = cfg.t_ftimwin ./ 2;
  
  %fprintf('For these frequencies, these are your time-window lengths, and the baseline period must end at least this many seconds before stimulus onset:\n');
  %fprintf('Freqs%s\n',sprintf(repmat('\t%.1f',1,length(cfg.foi)),cfg.foi));
  %fprintf('Length%s\n',sprintf(repmat('\t%.2f',1,length(cfg.t_ftimwin)),cfg.t_ftimwin));
  %fprintf('BL-end%s\n',sprintf(repmat('\t%.2f',1,length(min_baseline_s)),min_baseline_s));
  
  % check on the taper
  if ~isfield(cfg,'taper')
    fprintf('Using default taper: ''dpss''');
    cfg.taper = 'dpss';
  end
  
  if strcmp(cfg.taper,'hanning')
    % hanning taper
    fprintf('Using the ''%s'' taper. This is typically good for processing below 30 Hz.\n',cfg.taper);
    
    if min(cfg.foi) > 30
      fprintf('You''re only processing frequencies above 30 Hz. Consider using cfg.taper=''dpss'' (multitaper) instead of ''%s''.\n',cfg.taper);
    end
    
    if min(cfg.foi) < 30 && max(cfg.foi) > 30
      fprintf('You''re processing frequencies both below and above 30 Hz. Consider splitting the analyses into two parts.\n');
    end
    
    if plotit
      figure;
      subplot(2,2,1)
      plot(cfg.foi,cfg.t_ftimwin,'b*');
      xlabel('Freq (Hz)');
      ylabel('Time window length (s): t\_ftimwin');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      
      subplot(2,2,2)
      plot(cfg.foi,cfg.t_ftimwin ./ 2,'b*');
      xlabel('Freq (Hz)');
      ylabel('BL end before stim (s): t\_ftimwin ./ 2');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      ylim([0 cfg.t_ftimwin(1)]);
      
      subplot(2,2,3)
      plot(cfg.foi,cfg.t_ftimwin .* cfg.foi,'b*');
      xlabel('Freq (Hz)');
      ylabel('Cycles per time window: t\_ftimwin .* foi');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      
      subplot(2,2,4)
      plot(cfg.foi,cfg.toi(end) - min_baseline_s,'b*');
      xlabel('Freq (Hz)');
      ylabel(sprintf('%s calculated out to (s)',cfg.output));
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      
      publishfig(gcf,0,12,16);
    end
    
  elseif strcmp(cfg.taper,'dpss')
    % dpss taper, using multitapers
    fprintf('Using the ''%s'' taper, which is based on Slepian sequences as tapers. This is typically good for processing above 30 Hz.\n',cfg.taper);
    
    if max(cfg.foi) < 30
      fprintf('You''re only processing frequencies below 30 Hz. Consider using cfg.taper=''hanning'' instead of ''%s''.\n',cfg.taper);
    end
    
    if min(cfg.foi) < 30 && max(cfg.foi) > 30
      fprintf('You''re processing frequencies both below and above 30 Hz. Consider splitting the analyses into two parts.\n');
    end
    
    if ~isfield(cfg,'tapsmofrq')
      fprintf('Setting spectral smoothing (cfg.tapsmofrq) to 0.4*cfg.foi.\n');
      cfg.tapsmofrq = 0.4 * cfg.foi;
    end
    
    % K is the number of multitapers that will be applied
    K = diag(2 * cfg.t_ftimwin' * cfg.tapsmofrq - 1);
    if sum(K < 0) > 0
      fprintf('One frequency has fewer than zero tapers: That is bad! Set the multiplier in cfg.tapsmofrq to something larger than %.4f.\n',cfg.tapsmofrq/cfg.foi);
    end
    
    if plotit
      subplot(2,3,1)
      plot(cfg.foi,cfg.t_ftimwin,'b*');
      xlabel('Freq (Hz)');
      ylabel('Time window length (s): t\_ftimwin');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      
      subplot(2,3,2)
      plot(cfg.foi,cfg.t_ftimwin ./ 2,'b*');
      xlabel('Freq (Hz)');
      ylabel('BL end before stim (s): t\_ftimwin ./ 2');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      ylim([0 cfg.t_ftimwin(1)]);
      
      subplot(2,3,3)
      plot(cfg.foi,cfg.t_ftimwin .* cfg.foi,'b*');
      xlabel('Freq (Hz)');
      ylabel('Cycles per time window: t\_ftimwin .* foi');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      
      subplot(2,3,4)
      plot(cfg.foi,cfg.tapsmofrq,'b*');
      xlabel('Freq (Hz)');
      ylabel('Frequency smoothing: tapsmofrq');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      
      subplot(2,3,5)
      plot(cfg.foi,K,'b*');
      xlabel('Freq (Hz)');
      ylabel('Number of tapers');
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      ylim([(min(K) - 1) (max(K) + 1)]);
      
      subplot(2,3,6)
      plot(cfg.foi,cfg.toi(end) - min_baseline_s,'b*');
      xlabel('Freq (Hz)');
      ylabel(sprintf('%s calculated out to (s)',cfg.output));
      xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
      
      publishfig(gcf,0,12,16);
    end
    
  else
    fprintf('cfg.taper=%s is not currently supported in %s.\n',cfg.taper,mfilename);
  end
  
else
  fprintf('cfg.method=%s is not currently supported in %s.\n',cfg.method,mfilename);
end

% tell us about the baseline
% With your settings, you need to process
fprintf('\nThe average baseline period that you want to use is [%.2f %.2f] (length: %.2f s). ',baseline(1),baseline(end),diff(baseline));
if baseline(end) <= -1 * min_baseline_s(1)
  fprintf('Technically, that''s fine.\n');
elseif baseline(end) > -1 * min_baseline_s(1)
  fprintf('\n*** YOU MUST CHANGE IT TO END AT %.2f SECONDS OR EARLIER ***\n',-1*min_baseline_s(1));
  if strcmp(cfg.method,'wavelet')
    cycleStr = 'wavelet';
  elseif strcmp(cfg.method,'mtmconvol')
    cycleStr = 'time window';
  end
  fprintf('At the current setting of ending at %.2f s, the baseline period %s overlaps with ''%s'' calculation for the stimulus period.\n',baseline(end),cycleStr,cfg.output);
  fprintf('(Yeah, it takes a long chunk of data to process low frequencies properly.)\n');
end

fprintf('\nTo summarize:\n');
fprintf('The lowest frequency is %.2f Hz.\n',cfg.foi(1));
fprintf('To properly baseline correct at %.1f Hz, the baseline period must end at least %.4f s prior to stimulus onset.\n',cfg.foi(1),min_baseline_s(1));
fprintf('Thus, if you want to use data at %.1f Hz with this [%.2f %.2f] baseline correction (averaged across %.2f s):\n',cfg.foi(1),baseline(1),baseline(end),diff(baseline));
fprintf('\t- start processing at %.4f seconds (or earlier).\n',-2*min_baseline_s(1) - diff(baseline));
fprintf('\t- set the baseline to [%.4f %.4f].\n',-1*min_baseline_s(1) - diff(baseline),-1*min_baseline_s(1));

fprintf('When processing the data out to %.2f s, ''%s'' will be calculated out to %.4f s at %.2f Hz.\n',cfg.toi(end),cfg.output,cfg.toi(end) - min_baseline_s(1),cfg.foi(1));
fprintf('\t- If you want more usable ''%s'' data at %.2f Hz, take the last timepoint you want and add %.4f.\n',cfg.output,cfg.foi(1),min_baseline_s(1));

fprintf('\nNB: I think this is correct, but the FieldTrip tutorials do not pay such close attention to baseline periods.\nNeither do many published papers, and ultimately it probably doesn''t matter so much.\n');
end

