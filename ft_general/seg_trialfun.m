function trl = seg_trialfun(cfg)

if ~iscell(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% offset should be negative
offsetSamp = round(-cfg.trialdef.prestim*hdr.Fs);
% duration should be 1 sample less than the whole length of an event
durationSamp = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

trl = [];

for i = 1:length(event)
  if strcmp(event(i).type,cfg.trialdef.eventtype)
    if ~isempty(intersect(event(i).value, cfg.trialdef.eventvalue))
      %keyboard
      %if strcmp(strrep(event(i).value,' ',''),cfg.trialdef.eventvalue)
      
      % get the type of event for the trialinfo field
      if length(cfg.trialdef.eventvalue) > 1
        eventNumber = find(ismember(cfg.trialdef.eventvalue,event(i).value));
      else
        eventNumber = [];
      end
      
      % add this trial
      begsample = event(i).sample;
      endsample = begsample + durationSamp;
      
      trl = [trl; [begsample endsample offsetSamp eventNumber]];
    end
  end
end

