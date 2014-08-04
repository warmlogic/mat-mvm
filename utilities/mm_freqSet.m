function ana = mm_freqSet(ana,freqFxn)

% function ana = mm_freqSet(ana,freqFxn)
%
% Hard-coded frequency ranges from 3 to 80 Hz. Assumes there are 30
% log-spaced frequencies.
%
% ana:        ana struct
%
% freqFxn:    'ndtools' or 'ants'. Method used in mm_tla2pow.
%             default: 'ndtools' because it recalculates the actual center
%             frequency of each wavelet.
%
% Automatically splits alpha, beta, and gamma frequencies.
%
% default frequency range based on ndtools: wobj = mkwobj('morl', n_samples, sampleRate, 1./frequencies);
% ana.freq.theta = [4.1 7.7];
% ana.freq.alpha = [8.4 12];
% ana.freq.alpha_lower = [8.4 10.1];
% ana.freq.alpha_upper = [11 12];
% ana.freq.beta = [13.1 29.2];
% ana.freq.beta_lower = [13.1 20.5];
% ana.freq.beta_upper = [22.4 29.2];
% ana.freq.gamma = [31.9 77.4];
% ana.freq.gamma_lower = [31.9 49.7];
% ana.freq.gamma_upper = [54.3 77.4];
%
% 'ants' will group similarly using: logspace(log10(3),log10(80),30)
% ana.freq.theta = [4.2 7.4];
% ana.freq.alpha = [8.3 11.7];
% ana.freq.alpha_lower = [8.3 9.3];
% ana.freq.alpha_upper = [10.4 11.7];
% ana.freq.beta = [13.1 28.8];
% ana.freq.beta_lower = [13.1 20.6];
% ana.freq.beta_upper = [23 28.8];
% ana.freq.gamma = [32.3 80];
% ana.freq.gamma_lower = [32.3 50.9];
% ana.freq.gamma_upper = [57 80];

if nargin < 2
  freqFxn = 'ndtools';
end

if strcmp(freqFxn,'ndtools')
  ana.freq.theta = [4.1 7.7];
  
  ana.freq.alpha = [8.4 12];
  ana.freq.alpha_lower = [8.4 10.1];
  ana.freq.alpha_upper = [11 12];
  
  ana.freq.beta = [13.1 29.2];
  ana.freq.beta_lower = [13.1 20.5];
  ana.freq.beta_upper = [22.4 29.2];
  
  ana.freq.gamma = [31.9 77.4];
  ana.freq.gamma_lower = [31.9 49.7];
  ana.freq.gamma_upper = [54.3 77.4];
elseif strcmp(freqFxn,'ants')
  ana.freq.theta = [4.2 7.4];
  
  ana.freq.alpha = [8.3 11.7];
  ana.freq.alpha_lower = [8.3 9.3];
  ana.freq.alpha_upper = [10.4 11.7];
  
  ana.freq.beta = [13.1 28.8];
  ana.freq.beta_lower = [13.1 20.6];
  ana.freq.beta_upper = [23 28.8];
  
  ana.freq.gamma = [32.3 80];
  ana.freq.gamma_lower = [32.3 50.9];
  ana.freq.gamma_upper = [57 80];
end
