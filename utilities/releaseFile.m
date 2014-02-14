function releaseFile(filename)
%RELEASEFILE - Release a file locked with lockFile.
%
% This function helps to keep multiple procs from working on the same
% file when on a cluster.  When used along with lockFile it will
% release the lock on a file once you are done working on it.
%
% FUNCTION:
%   releaseFile(filename)
%
% INPUT ARGS:
%   filename- The file to release.
%
%
% Attribution: Borrowed from the Kahana Lab's eeg_toolbox

if strcmp(computer,'MACI64')
  % test name
  lockname = [filename '.lock'];
  
  if exist(lockname,'file')
    % release it
    system(['rm -f ' lockname]);
  else
    % tell them
    fprintf('No lock found for %s.\n',filename);
  end
elseif strcmp(computer,'GLNXA64') || strcmp(computer,'GLNX86')
  % for the dream cluster
  system(['lockfile-remove' filename]);
end
