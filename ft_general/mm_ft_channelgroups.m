function ana = mm_ft_channelgroups(varargin)
%
% This is a backwards compatibility wrapper.
%
% This function has been renamed mm_ft_elecGroups. Please update your code.
%
% See also: MM_FT_ELECGROUPS
%

fprintf('mm_ft_channelgroups has been renamed mm_ft_elecGroups. Please update your code.\n');

ana = mm_ft_elecGroups(varargin{1});
