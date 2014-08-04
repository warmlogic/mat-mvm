function freq = mm_freqSet(freqFxn)

% function freq = mm_freqSet(freqFxn)
%
% Hard-coded frequency ranges from 3 to 80 Hz. Assumes there are 30
% log-spaced frequencies.
%
% freqFxn:    string: 'ndtools' or 'ants'. Method used in mm_tla2pow.
%             default: 'ndtools' because it recalculates the actual center
%             frequency of each wavelet using:
%             wobj = mkwobj('morl', n_samples, sampleRate, 1./frequencies);
%             NB: n_samples and sampleRate do not affect the frequencies
%
% Automatically splits alpha, beta, and gamma frequencies.
%
% default frequency range based on ndtools:
% freq.fullRange = [4.1 77.4];
% freq.theta = [4.1 7.7];
% freq.alpha = [8.4 12];
% freq.alpha_lower = [8.4 10.1];
% freq.alpha_upper = [11 12];
% freq.beta = [13.1 29.2];
% freq.beta_lower = [13.1 20.5];
% freq.beta_upper = [22.4 29.2];
% freq.gamma = [31.9 77.4];
% freq.gamma_lower = [31.9 45.5];
% freq.gamma_upper = [49.7 77.4];
%
% 'ants' will group similarly using: logspace(log10(3),log10(80),30)
% freq.fullRange = [4.2 80];
% freq.theta = [4.2 7.4];
% freq.alpha = [8.3 11.7];
% freq.alpha_lower = [8.3 9.3];
% freq.alpha_upper = [10.4 11.7];
% freq.beta = [13.1 28.8];
% freq.beta_lower = [13.1 20.6];
% freq.beta_upper = [23 28.8];
% freq.gamma = [32.3 80];
% freq.gamma_lower = [32.3 45.4];
% freq.gamma_upper = [50.9 80];

if nargin < 1
  freqFxn = 'ndtools';
end

freq = struct;

if strcmp(freqFxn,'ndtools')
  freq.fullRange = [4.1 77.4];
  
  freq.theta = [4.1 7.7];
  
  freq.alpha = [8.4 12];
  freq.alpha_lower = [8.4 10.1];
  freq.alpha_upper = [11 12];
  
  freq.beta = [13.1 29.2];
  freq.beta_lower = [13.1 20.5];
  freq.beta_upper = [22.4 29.2];
  
  freq.gamma = [31.9 77.4];
  freq.gamma_lower = [31.9 45.5];
  freq.gamma_upper = [49.7 77.4];
  
elseif strcmp(freqFxn,'ants')
  freq.fullRange = [4.2 80];
  
  freq.theta = [4.2 7.4];
  
  freq.alpha = [8.3 11.7];
  freq.alpha_lower = [8.3 9.3];
  freq.alpha_upper = [10.4 11.7];
  
  freq.beta = [13.1 28.8];
  freq.beta_lower = [13.1 20.6];
  freq.beta_upper = [23 28.8];
  
  freq.gamma = [32.3 80];
  freq.gamma_lower = [32.3 45.4];
  freq.gamma_upper = [50.9 80];
end
