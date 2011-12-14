function trl = seg_trialfun(cfg)

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% offset should be negative
offsetSamp = round(-cfg.trialdef.prestim*hdr.Fs);
% duration should be 1 sample less than the whole length of an event
durationSamp = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
% TODO: should this be ceil instead of round?

% initialize the trl matrix
trl = [];
%trl = nan(length(ismember({event(:).value},cfg.trialdef.eventvalue)),4);

for i = 1:length(event)
  if strcmp(event(i).type,cfg.trialdef.eventtype)
    %if ~isempty(intersect(event(i).value, cfg.trialdef.eventvalue))
    %keyboard
    %if strcmp(deblank(event(i).value),cfg.trialdef.eventvalue)
    if ismember(deblank(event(i).value),cfg.trialdef.eventvalue)
      
      % get the type of event for the trialinfo field
      if length(cfg.trialdef.eventvalue) > 1
        eventNumber = find(ismember(cfg.trialdef.eventvalue,deblank(event(i).value)));
      else
        eventNumber = [];
      end
      
      % add this trial [beginning sample, ending sample, offset, evNum]
      trl = cat(1,trl,[event(i).sample, (event(i).sample + durationSamp), offsetSamp, eventNumber]);
      
%       % add this trial
%       %
%       % beginning sample
%       trl(i,1) = event(i).sample;
%       % ending sample
%       trl(i,2) = event(i).sample + durationSamp;
%       % offset
%       trl(i,3) = offsetSamp;
%       % type of trial
%       trl(i,4) = eventNumber;
      
    end
  end
end
