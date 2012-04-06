function [data] = mm_rmBadChan(exper,ana,data)
%MM_RMBADCHAN - Remove bad channels from individual subject data
%
% [data] = mm_rmBadChan(exper,ana,data)
%
% Input:
%   exper   = must contain badChan field
%   ana     = used for event values
%   data    = ERP or TF data
%
% Output
%   data    = data structs with bad channels removed (deleted)
%
% NB: Currently this runs FieldTrip's ft_selectdata_old because there is a
%     bug in ft_selectdata_new (dimord field gets removed)
%     http://bugzilla.fcdonders.nl/show_bug.cgi?id=1409
%
%
% See also: MM_GETBADCHAN

if ~isfield(exper,'badChan')
  error('No bad channel information in exper struct. Run mm_getBadChan.');
end

for sub = 1:length(exper.subjects)
  for ses = 1:length(exper.sesStr)
    if ~isempty(exper.badChan{sub,ses})
      for typ = 1:length(ana.eventValues)
        for evVal = 1:length(ana.eventValues{typ})
          cfg_sd = [];
          cfg_sd.channel = ft_channelselection(eval(sprintf('{''all''%s};',sprintf(repmat(' ''-E%d''',1,length(exper.badChan{sub,ses})),exper.badChan{sub,ses}))),data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.label);
          % remove the bad channels
          fprintf('Removing %d bad channel(s) from %s, %s, %s\n',length(exper.badChan{sub,ses}),exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
          
          % % ft_selectdata_new
          % data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_selectdata(cfg_sd,data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data);
          
          % ft_selectdata_old
          data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data = ft_selectdata(data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data,'channel',cfg_sd.channel);
          
          % % setting bad channels to nan - doesn't work with GA
          % fprintf('Setting bad channels to NaNs for %s, %s, %s\n',exper.subjects{sub},exper.sesStr{ses},ana.eventValues{typ}{evVal});
          % data.(ana.eventValues{typ}{evVal}).sub(sub).ses(ses).data.(param)(exper.badChan{sub,ses},:) = NaN;
        end % evVal
      end % typ
    else
      fprintf('%s %s has no bad channels.\n',exper.subjects{sub},exper.sesStr{ses});
      continue
    end % if
  end % ses
end % sub

end

