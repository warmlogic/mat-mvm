function trl = seg_trialfun(cfg)

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
  cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% IMPORTANT: the way offsetSamp, durationSamp, and event(i).sample get
% combined in this function are a special case for processing pre-segmented
% data. It differs from the normal continuous data segmentation described
% here:
% http://fieldtrip.fcdonders.nl/example/making_your_own_trialfun_for_conditional_trial_definition
%
% This is because when using cfg.trialdef.eventtype='trial' (as we do here)
% the event tag in the pre-segmented file is at the beginning of the
% segment that we want to demarcate with a particular sample, meaning it
% likely already has some pre-stimulus period built in (if you segmented it
% this way). This is in contrast to the typical FT setup using
% cfg.trialdef.eventtype='trigger', where the event we want to segment
% around is actually at that trigger, and so the prestim perior needs to be
% subtracted from the event's sample location.

% offset should be negative
offsetSamp = -round(cfg.trialdef.prestim * hdr.Fs);
% duration should be 1 sample less than the whole length of an event
durationSamp = round((cfg.trialdef.poststim + cfg.trialdef.prestim) * hdr.Fs) - 1;
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
      
      if ~isempty(eventNumber) && length(eventNumber) > 1
        warning('More than one event value matched the current value (%s)',deblank(event(i).value));
        keyboard
      end
      
      % TODO: add in a check comparing
      % range(event(i).sample,event(i).sample + durationSamp), maybe use
      % ft_read_data; or maybe durationSamp vs event(i).sample - event(i-1).sample
      
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
