function ana = mm_ft_channelgroups(ana,elecGroupsExtra,elecGroupsStrExtra)
%MM_FT_CHANNELGROUPS Set up channel groups
%
%   ana = mm_ft_channelgroups(ana,elecGroupsExtra,elecGroupsStrExtra)
%
%  Input:
%    elecGroupsExtra    = Optional. Muse be a cell array containing cell
%                         arrays of strings of electrode labels, even if
%                         there is only one group.
%                         e.g., {{'E1','E2'},{'E3','E4'}}.
%
%    elecGroupsStrExtra = Optional. Identification strings for the groups.
%                         A cell array of strings, one for each group.
%
%  Output:
%    adds the elecGroups and elecGroupsStr fields to the ana struct
%
% NB: See the function code to view the pre-defined electrode groups
%

if nargin == 3
  if (isempty(elecGroupsExtra) && ~isempty(elecGroupsStrExtra)) || (~isempty(elecGroupsExtra) && isempty(elecGroupsStrExtra))
    error('elecGroupsExtra and elecGroupsStrExtra both need to contain something.');
  elseif ~isempty(elecGroupsExtra) && ~isempty(elecGroupsStrExtra)
    if length(elecGroupsExtra) ~= length(elecGroupsStrExtra)
    error('elecGroupsExtra and elecGroupsStrExtra must be the same length.');
    end
  end
elseif nargin == 2
  error('Cannot define elecGroupsExtra without elecGroupsStrExtra');
elseif nargin == 1
  elecGroupsExtra = {};
  elecGroupsStrExtra = {};
end

ana.elecGroups = {...
  {'E32','E33','E38','E39','E43','E44','E128'},... % LAI
  {'E1','E114','E115','E120','E121','E122','E125'},... % RAI
  {'E12','E13','E19','E20','E24','E28','E29'},... % LAS
  {'E4','E5','E111','E112','E117','E118','E124'},... % RAS
  {'E37','E42','E52','E53','E54','E60','E61'},... % LPS
  {'E78','E79','E85','E86','E87','E92','E93'},... % RPS
  {'E57','E58','E63','E64','E65','E68','E69'},... %LPI
  {'E89','E90','E94','E95','E96','E99','E100'},... % RPI
  {'E52','E53','E60','E61','E59','E66','E67'},... % LPS2
  {'E77','E78','E84','E85','E86','E91','E92'},... % RPS2
  {'E4','E10','E11','E15','E16','E18','E19'},... % Frontal Inferior
  {'E4','E5','E6','E11','E12','E13','E19','E112'},... % Frontal Superior
  {'E7','E31','E55','E80','E106','E129'},... % Central (central)
  {'E54','E61','E62','E67','E72','E77','E78','E79'},... % Posterior Superior
  {'E70','E71','E74','E75','E76','E82','E83'},... % Posterior Inferior
  {'E11'},... % Fz (frontocentral)
  {'Cz'},... % Cz (central)
  {'E75'},... % Oz - P1 effect (occipitocentral)
  {'E62'},... % Pz (posteriocentral)
  {'E64'},... % P9 - N1 effect? (left lateral)
  {'E95'},... % P10 - N1 effect? (right lateral)
  };

ana.elecGroupsStr = {'LAI','RAI','LAS','RAS','LPS','RPS','LPI','RPI','LPS2','RPS2','FI','FS','C','PS','PI','Fz','Cz','Oz','Pz','P9','P10'};

if ~isempty(elecGroupsExtra)
  if sum(ismember(elecGroupsStrExtra,ana.elecGroupsStr)) > 0
    repeats = elecGroupsStrExtra(ismember(elecGroupsStrExtra,ana.elecGroupsStr));
    error('elecGroupsStrExtra group(s) %salready exists in the predefined elecGroupsStr field. Unfortunately you need to rename the offending electrode group.',sprintf(repmat('''%s'' ',1,sum(ismember(elecGroupsStrExtra,ana.elecGroupsStr))),repeats{:}));
  end
  
  ana.elecGroups = cat(2,ana.elecGroups,elecGroupsExtra);
  ana.elecGroupsStr = cat(2,ana.elecGroupsStr,elecGroupsStrExtra);
end

if sum(size(ana.elecGroups) ~= size(ana.elecGroupsStr)) ~= 0
  error('The number of groups in ana.elecGroups is not the same as in ana.elecGroupsStr.')
end

end

% {'E62','E67','E71','E72','E75','E76','E77'},... % Posterior Inferior
