function cfg = mm_freqcheck(cfg,baseline)
%MM_FREQCHECK - Checks your freq cfg for ft_freqanalysis
%
% cfg = mm_freqcheck(cfg,baseline)
%
% Input:
%   cfg        = cfg that will go into FT_FREQANALYSIS
%   baseline   = e.g., [-1.0 -0.75]
%   samplerate = sample rate of data
%
% Output:
%   cfg        = same cfg with some defaults set if you didn't set them
%
% See also: FT_FREQANALYSIS
%

if ~isfield(cfg,'method')
  error('Must set cfg.methd (e.g., ''wavelet'' or ''mtmconvol'')');
end

if ~isfield(cfg,'foi') && ~isfield(cfg,'foilim')
  error('Must set cfg.foi (e.g., 4:1:8) or cfg.foilim (e.g., [4 8])');
end
% TODO: how does FieldTrip determine the spacing of foilim?
if ~isfield(cfg,'foi') && isfield(cfg,'foilim')
  fprintf('You set cfg.foilim=[%.1f %.1f]. Removing this field and setting cfg.foi=%.1f:1:%.1f\n',cfg.foilim(1),cfg.foilim(2),cfg.foilim(1),cfg.foilim(2));
  cfg.foi = cfg.foilim(1):1:cfg.foilim(2);
  cfg = rmfield(cfg,'foilim');
end

if ~isfield(cfg,'toi')
  error('Must set cfg.toi (e.g., -1.0:0.04:2.0)');
end

if ~isfield(cfg,'output')
  cfg.output = 'pow';
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

if ~isfield(cfg,'pad')
  cfg.pad = 'maxperlen';
end


% % test parameters
%
% % wavelet
% cfg = [];
% cfg.foi = 4:1:8;
% cfg.toi = -1.0:0.04:2.0;
% baseline = [-1.0 -0.75];
% cfg.method = 'wavelet';
% cfg.width = 6;
%
% % multitaper - hanning taper
% cfg = [];
% cfg.foi = 1:2:30;
% cfg.toi = -1.0:0.04:2.0;
% baseline = [-0.3 -0.1];
% cfg.method = 'mtmconvol';
% cfg.taper = 'hanning';
% cfg.t_ftimwin = 6 ./ cfg.foi;
%
% % multitaper - slepian tapers
% cfg = [];
% cfg.foi = 1:2:30;
% cfg.toi = -1.0:0.04:2.0;
% baseline = [-0.3 -0.1];
% cfg.method = 'mtmconvol';
% cfg.taper = 'dpss';
% cfg.t_ftimwin = 6 ./ cfg.foi;
% cfg.tapsmofrq = 0.4 * cfg.foi;

fprintf('\nUsing the ''%s'' method with ''%s'' output.\n',cfg.method,cfg.output);

fprintf('ft_freqanalysis will process:\n');
fprintf('Time: %.2f to %.2f s at every %.1f ms\n',cfg.toi(1),cfg.toi(end),diff([cfg.toi(1) cfg.toi(2)]) * 1000);
fprintf('Freq: %.2f to %.2f Hz\n',cfg.foi(1),cfg.foi(end));

if strcmp(cfg.method,'wavelet')
  % wavelet
  if ~isfield(cfg,'width')
    fprintf('Setting cfg.width=6, a good default.\n');
    cfg.width = 6;
  end
  fprintf('Using wavelets with %d cycles.\n',cfg.width);
  wavelet_widths_s = (1 ./ cfg.foi * cfg.width);
  min_baseline_s = wavelet_widths_s ./ 2;
  
  %fprintf('For these frequencies, these are your wavelet widths, and the baseline period must end at least this many seconds before stimulus onset:\n');
  %fprintf('Freqs%s\n',sprintf(repmat('\t%.1f',1,length(cfg.foi)),cfg.foi));
  %fprintf('Width%s\n',sprintf(repmat('\t%.2f',1,length(wavelet_widths_s)),wavelet_widths_s));
  %fprintf('BL-end%s\n',sprintf(repmat('\t%.2f',1,length(min_baseline_s)),min_baseline_s));
  
  subplot(1,3,1)
  plot(cfg.foi,wavelet_widths_s,'b*');
  xlabel('Freq (Hz)');
  ylabel('Wavelet widths (s)');
  xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
  ylim([0 wavelet_widths_s(1)]);
  
  subplot(1,3,2)
  plot(cfg.foi,wavelet_widths_s ./ 2,'b*');
  xlabel('Freq (Hz)');
  ylabel('BL end before stimulus (s): wavelet\_widths ./ 2');
  xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
  ylim([0 wavelet_widths_s(1)]);
  
  subplot(1,3,3)
  plot(cfg.foi,cfg.width,'b*');
  xlabel('Freq (Hz)');
  ylabel('Cycles per wavelet: cfg.width');
  xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
  
elseif strcmp(cfg.method,'mtmconvol')
  % multitaper method
  
  if ~isfield(cfg,'t_ftimwin')
    fprintf('Setting frequency-dependent time window (cfg.t_ftimwin; temporal smoothing) to 6./cfg.foi (6 cycles per time window), a good default.\n');
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
    fprintf('Using the ''%s'' taper. This is good for processing below 30 Hz.\n',cfg.taper);
    
    if min(cfg.foi) > 30
      fprintf('You''re only processing frequencies above 30 Hz. You should consider using the dpss taper/multitaper instead.\n');
    end
    
    subplot(1,3,1)
    plot(cfg.foi,cfg.t_ftimwin,'b*');
    xlabel('Freq (Hz)');
    ylabel('Time window length (s): t\_ftimwin');
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    
    subplot(1,3,2)
    plot(cfg.foi,cfg.t_ftimwin ./ 2,'b*');
    xlabel('Freq (Hz)');
    ylabel('BL end before stimulus (s): t\_ftimwin ./ 2');
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    ylim([0 cfg.t_ftimwin(1)]);
    
    subplot(1,3,3)
    plot(cfg.foi,cfg.t_ftimwin .* cfg.foi,'b*');
    xlabel('Freq (Hz)');
    ylabel('Cycles per time window: t\_ftimwin .* foi');
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    
  elseif strcmp(cfg.taper,'dpss')
    % dpss taper, using multitapers
    fprintf('Using the ''%s'' taper, which is based on Slepian sequences as tapers. This is good for processing above 30 Hz.\n',cfg.taper);
    
    if max(cfg.foi) < 30
      fprintf('You''re only processing frequencies below 30 Hz. You should consider using the hanning taper instead.\n');
    end
    
    if ~isfield(cfg,'tapsmofrq')
      fprintf('Setting spectral smoothing (cfg.tapsmofrq) to 0.4*cfg.foi, a good default.\n');
      cfg.tapsmofrq = 0.4 * cfg.foi;
    end
    
    % K is the number of multitapers that will be applied
    K = diag(2 * cfg.t_ftimwin' * cfg.tapsmofrq - 1);
    if sum(K < 0) > 0
      fprintf('One frequency has fewer than zero tapers: That is bad! Set the multiplier in cfg.tapsmofrq to something larger than %.4f.\n',cfg.tapsmofrq/cfg.foi);
    end
    
    subplot(2,3,1)
    plot(cfg.foi,cfg.t_ftimwin,'b*');
    xlabel('Freq (Hz)');
    ylabel('Time window length (s): t\_ftimwin');
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    
    subplot(2,3,2)
    plot(cfg.foi,cfg.t_ftimwin ./ 2,'b*');
    xlabel('Freq (Hz)');
    ylabel('BL end before stimulus (s): t\_ftimwin ./ 2');
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
    ylabel('Number of tapers (K)');
    ylim([(min(K) - 1) (max(K) + 1)]);
    xlim([(cfg.foi(1) - 1) (cfg.foi(end) + 1)]);
    
  else
    fprintf('cfg.taper=%s is not supported.\n',cfg.taper);
  end
  
else
  fprintf('cfg.method=%s is not supported in %s.\n',cfg.method,mfilename);
end

% tell us about the baseline
% With your settings, you need to process
fprintf('The average baseline period that you want to use is [%.2f %.2f] (length: %.2f s). ',baseline(1),baseline(end),diff(baseline));
if baseline(end) <= -1 * min_baseline_s(1)
  fprintf('That''s fine.\n');
elseif baseline(end) > -1 * min_baseline_s(1)
  fprintf('\n*** YOU MUST CHANGE IT TO END AT %.4f OR EARLIER ***\nAt the current setting of ending at %.2f s, it overlaps with power calculation for the stimulus period.\n', -1*min_baseline_s(1),baseline(end));
end

fprintf('\nTo summarize:\n');
fprintf('To properly baseline correct at the lowest frequency (%.1f Hz), your baseline must end at least %.4f s prior to the stimulus.\n',cfg.foi(1),min_baseline_s(1));
fprintf('Thus, if you want to use data at the lowest frequency (%.1f Hz) with %.2f s avg baseline correction, you need to start processing at %.4f seconds or earlier.\n',cfg.foi(1),diff(baseline),-2*min_baseline_s(1) - diff(baseline));

fprintf('If you processing the data out to %.2f s, the data will be useable out to %.4f s at the lowest frequency (%.1f Hz).\n',cfg.toi(end),cfg.toi(end) - min_baseline_s(1),cfg.foi(1));

fprintf('Note that I think this is correct, but the FieldTrip tutorials do not pay such close attention to baseline periods.\n');
end

