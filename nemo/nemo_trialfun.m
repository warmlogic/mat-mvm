function trl = nemo_trialfun(cfg)

% operates using Net Station evt files and event structs

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
    cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
fprintf('Reading flags from EEG file using FieldTrip...');
ft_hdr = ft_read_header(cfg.dataset);
ft_event = ft_read_event(cfg.dataset);
fprintf('Done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in external data, if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg.eventinfo.useMetadata
    md = cfg.eventinfo.metadata;
    
    if ismember('nsEvt', md.types)
        % read the file with information about all events
        evtDir = 'ns_evt';
        % find the evt file
        evtfile = dir(fullfile(md.dataroot,md.sesName,evtDir,[md.subject,'*.evt']));
        if isempty(evtfile)
            error('Cannot find %s*.evt file in %s',md.subject,fullfile(md.dataroot,md.sesName,evtDir));
        elseif length(evtfile) > 1
            error('More than one %s*.evt file found in %s',md.subject,fullfile(md.dataroot,md.sesName,evtDir));
        elseif length(evtfile) == 1
            evtfile = fullfile(md.dataroot,md.sesName,evtDir,evtfile.name);
        end
        fprintf('Reading evt file: %s...',evtfile);
        % figure out how many columns there are
        fid = fopen(evtfile,'r');
        maxNumCols = -Inf;
        while 1
            % get each line of the file
            tline = fgetl(fid);
            if ischar(tline)
                if strcmp(tline(end),sprintf('\t'))
                    tline = tline(1:end-1);
                end
                % since it's tab delimited, split it and count the length
                numCols = length(regexp(tline,'\t','split'));
                % change the max number
                if numCols > maxNumCols
                    maxNumCols = numCols;
                end
            else
                break
            end
        end
        fclose(fid);
        % read the evt file
        fid = fopen(evtfile,'r');
        ns_evt = textscan(fid,repmat('%s',1,maxNumCols),'Headerlines',3,'Delimiter','\t');
        fclose(fid);
        fprintf('Done.\n');
    end
end

% read the types of triggers (flags) in the Net Station evt file
if ismember('nsEvt', md.types)
    triggers = unique(ns_evt{1});
else
    warning('Need to set up %s to find trigger types!',mfilename);
    keyboard
    triggers = {};
    %triggers = {'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND', 'DIN '};
end

% initialize the trl matrix
trl = [];

% all trls need to have the same length
maxTrlCols = -Inf;
fn_trl_ord = fieldnames(cfg.eventinfo.trl_order);
for fn = 1:length(fn_trl_ord)
    if ismember(fn_trl_ord{fn},cfg.eventinfo.eventValues)
        if length(cfg.eventinfo.trl_order.(fn_trl_ord{fn})) > maxTrlCols
            maxTrlCols = length(cfg.eventinfo.trl_order.(fn_trl_ord{fn}));
        end
    end
end
if maxTrlCols == -Inf
    fprintf('Did not set maximum number of trialinfo columns!\n');
    keyboard
end
timeCols = 3;
trl_ini = -1 * ones(1, timeCols + maxTrlCols);

% only keep the ft events with triggers
ft_event = ft_event(ismember({ft_event.value},triggers));

%% Alignment 1/3: read evt and put events with the same sample in alphabetical order

% The header, as read by FieldTrip, is (usually) sorted alphabetically when
% events with different values occurred at the same sample (but not
% necessarily the same millisecond). Therefore, events in the Net Station
% evt file need to be rearranged, especially because the next section of
% this function will remove "duplicate" events (defined as events with the
% same value that occurred at the same sample). There is a third section of
% this function to catch any events that should not be in alphabetical
% order.

% backup so we can access original data
ns_evt_orig = ns_evt;

% % wait to change the values
% ecInd = 1:length(ns_evt{1});

for ec = 1:length(ns_evt{1})
    if ec > 1 && ~strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
        this_time_ms = (str2double(ns_evt{5}{ec}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec}(8:9)) * 1000) + (str2double(ns_evt{5}{ec}(11:13)));
        this_time_samp = fix((this_time_ms / 1000) * ft_hdr.Fs);
        prev_time_ms = (str2double(ns_evt{5}{ec-1}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(8:9)) * 1000) + (str2double(ns_evt{5}{ec-1}(11:13)));
        prev_time_samp = fix((prev_time_ms / 1000) * ft_hdr.Fs);
        
        if this_time_samp == prev_time_samp || abs(this_time_ms - prev_time_ms) < (1000 / ft_hdr.Fs)
            if ~issorted({ns_evt{1}{ec-1},ns_evt{1}{ec}})
                % change the values on the fly
                for i = 1:length(ns_evt_orig)
                    ns_evt{i}(ec - 1) = ns_evt_orig{i}(ec);
                    ns_evt{i}(ec) = ns_evt_orig{i}(ec - 1);
                end
                
                % % wait to change the values
                % prevInd = ec - 1;
                % ecInd(ec) = prevInd;
                % ecInd(prevInd) = ec;
            end
        end
    end
    
end

% % wait to change the values
% for i = 1:length(ns_evt)
%   ns_evt{i} = ns_evt{i}(ecInd);
% end

%% Alignment 2/3: remove duplicate events (from FieldTrip's sample-based perspective)

% If two (or more) events with the same value string occurred at the same
% sample, FieldTrip will only keep one of them. Therefore, we need to
% remove the duplicate events from the Net Station evt file.

ecInd = true(length(ns_evt{1}),1);

% keep track of how many real evt events we have counted
ec = 0;

for i = 1:length(ft_event)
    if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
        % found a trigger in the EEG file events; increment index if
        % value is correct.
        
        if ismember(ft_event(i).value,triggers)
            ec = ec + 1;
        else
            continue
        end
        
        if ec > 1 && strcmp(ns_evt{1}(ec),ns_evt{1}(ec-1))
            this_time_ms = (str2double(ns_evt{5}{ec}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec}(8:9)) * 1000) + (str2double(ns_evt{5}{ec}(11:13)));
            this_time_samp = fix((this_time_ms / 1000) * ft_hdr.Fs);
            prev_time_ms = (str2double(ns_evt{5}{ec-1}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec-1}(8:9)) * 1000) + (str2double(ns_evt{5}{ec-1}(11:13)));
            prev_time_samp = fix((prev_time_ms / 1000) * ft_hdr.Fs);
            
            if (this_time_samp == prev_time_samp || abs(this_time_ms - prev_time_ms) < (1000 / ft_hdr.Fs))
                % need to check whether these are the same events; can't rely on
                % comparing current NS and FT events with ~strcmp because they
                % might have the same value even though they're not the same exact
                % event. So, check the event value for the next NS event vs the
                % current FT event, and make it more robust by checking the sample
                % offset
                next_time_ms = (str2double(ns_evt{5}{ec+1}(2:3)) * 60 * 60 * 1000) + (str2double(ns_evt{5}{ec+1}(5:6)) * 60 * 1000) + (str2double(ns_evt{5}{ec+1}(8:9)) * 1000) + (str2double(ns_evt{5}{ec+1}(11:13)));
                next_time_samp = fix((next_time_ms / 1000) * ft_hdr.Fs);
                
                if strcmp(ns_evt{1}(ec+1),ft_event(i).value) && abs(next_time_samp - prev_time_samp) == abs(ft_event(i).sample - ft_event(i-1).sample)
                    ecInd(ec) = false;
                    ec = ec + 1;
                end
            end
            
        end
    end
end

% only keep the ns_evt indices that align with ft_event
for i = 1:length(ns_evt)
    ns_evt{i} = ns_evt{i}(ecInd);
end

%% Alignment 3/3: final comparison of Net Station and FieldTrip data

% make sure the Net Station evt Code column and ft_event value field are in
% the same order. If they are not (e.g., if the events should not have been
% rearranged alphabetically), then put the Net Station events in the same
% order as the FieldTrip events, but only if the neighboring events are
% aligned.

ns_evt_backup = ns_evt;

% verify that our lists have the same event values

% if length(ns_evt{1}) == length(ft_event(~ismember({ft_event.value},'epoc')))
if length(ns_evt{1}) == length(ft_event)
    for i = 1:length(ft_event)
        if ~strcmp(ns_evt{1}{i},ft_event(i).value)
            if strcmp(ns_evt{1}{i-1},ft_event(i-1).value) && strcmp(ns_evt{1}{i},ft_event(i+1).value) && strcmp(ns_evt{1}{i+1},ft_event(i).value) && strcmp(ns_evt{1}{i+2},ft_event(i+2).value)
                % they're just out of order for some reason (probably due to
                % alphabetical sorting), so put them back
                for j = 1:length(ns_evt_orig)
                    ns_evt{j}(i) = ns_evt_backup{j}(i + 1);
                    ns_evt{j}(i + 1) = ns_evt_backup{j}(i);
                end
                
            else
                warning('Index %d does not have the same value! ns_evt: %s, ft_event: %s',i,ns_evt{1}{i},ft_event(i).value);
                keyboard
            end
        end
    end
else
    warning('ns_evt and ft_event are not the same length! Need to fix alignment code.');
    keyboard
end

%% go through the events

ses = cfg.eventinfo.sessionNum;
% sesName = cfg.eventinfo.sessionNames{ses};
sesType = find(ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses}));

for pha = 1:length(cfg.eventinfo.phaseNames{sesType})
    phaseName = cfg.eventinfo.phaseNames{sesType}{pha};
    %phaseType = find(ismember(cfg.eventinfo.phaseNames{sesType},phaseName));
    phaseType = pha;
    
    %% process the study phase
    
    fprintf('Processing %s...\n',phaseName);
    
    % keep track of how many real evt events we have counted
    ec = 0;
    
    fprintf('%s NS event flag count: %s',phaseName,repmat(' ',1,length(num2str(length(ft_event)))));
    
    for i = 1:length(ft_event)
        fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
        
        if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
            % found a trigger in the EEG file events; increment index if
            % value is correct.
            
            %if ~ismember(event(i).value,{'epoc'})
            if ismember(ft_event(i).value,triggers)
                ec = ec + 1;
            else
                continue
            end
            %%%Deleted Hack%%%
            switch ft_event(i).value
                %case 'TRSP'
                case {'prm+','trg+','stm+'}
                    
                    if strcmp(ns_evt{1}(ec),ft_event(i).value)
                        
                        % find the following TRSP because that row contains
                        % the metadata that we want to store
                        trspInd = ec+1;
                        while 1
                            if strcmp(ns_evt{1}(trspInd),'TRSP')
                                break
                            else
                                trspInd = trspInd + 1;
                            end
                        end
                        
                        % set column types because Net Station evt files can vary
                        ns_evt_cols = {};
                        for ns = 1:size(ns_evt,2)
                            ns_evt_cols = cat(1,ns_evt_cols,ns_evt{ns}(trspInd));
                        end
                        cols.(phaseName).phase = find(strcmp(ns_evt_cols,'exid'));
                        if isempty(cols.(phaseName).phase)
                            keyboard
                        end
                        
                        if strcmp(ns_evt{cols.(phaseName).phase+1}(trspInd),phaseName)
                            
                            cols.(phaseName).trial = find(strcmp(ns_evt_cols,'trl#'));
                            cols.(phaseName).cell_label = find(strcmp(ns_evt_cols,'cel#'));
                            cols.(phaseName).times_studied = find(strcmp(ns_evt_cols,'obs#'));
                            cols.(phaseName).resp_value = find(strcmp(ns_evt_cols,'rsp#'));
                            cols.(phaseName).accuracy = find(strcmp(ns_evt_cols,'eval'));
                            cols.(phaseName).reaction_time = find(strcmp(ns_evt_cols,'rtim'));
                            cols.(phaseName).random_number = find(strcmp(ns_evt_cols,'rnum'));
                            
                            if isempty(cols.(phaseName).trial)
                                keyboard
                            end
                            
                            switch phaseName
                                
                                case {'TC_NEMO_fN400study','TC_NEMO_fN400test'}
                                    % set column types because Net Station evt files can vary
                                    %cols.(phaseName).prime_string = find(strcmp(ns_evt_cols,'prim')); %string
                                    %cols.(phaseName).target_string = find(strcmp(ns_evt_cols,'targ')); %string
                                    cols.(phaseName).FSG_target = find(strcmp(ns_evt_cols,'fstg'));
                                    if isempty(cols.(phaseName).FSG_target)
                                        FSG_target = -1;
                                    end
                                    cols.(phaseName).BSG_target = find(strcmp(ns_evt_cols,'bstg'));
                                    if isempty(cols.(phaseName).BSG_target)
                                        BSG_target = -1;
                                    end
                                    %cols.(phaseName).related_target = find(strcmp(ns_evt_cols,'relt')); %string
                                    cols.(phaseName).CNC_target = find(strcmp(ns_evt_cols,'cnct'));
                                    if isempty(cols.(phaseName).CNC_target)
                                        CNC_target = -1;
                                    end
                                    cols.(phaseName).KFFRQ_target = find(strcmp(ns_evt_cols,'frqt'));
                                    if isempty(cols.(phaseName).KFFRQ_target)
                                        KFFRQ_target = -1;
                                    end
                                    cols.(phaseName).NLET_target = find(strcmp(ns_evt_cols,'lent'));
                                    if isempty(cols.(phaseName).NLET_target)
                                        NLET_target = -1;
                                    end
                                    cols.(phaseName).NPHON_target = find(strcmp(ns_evt_cols,'phot'));
                                    if isempty(cols.(phaseName).NPHON_target)
                                        NPHON_target = -1;
                                    end
                                    cols.(phaseName).NSYLL_target = find(strcmp(ns_evt_cols,'sylt'));
                                    if isempty(cols.(phaseName).NSYLL_target)
                                        NSYLL_target = -1;
                                    end
                                    cols.(phaseName).orthoN_target = find(strcmp(ns_evt_cols,'ornt'));
                                    if isempty(cols.(phaseName).orthoN_target)
                                        orthoN_target = -1;
                                    end
                                    cols.(phaseName).phonoN_target = find(strcmp(ns_evt_cols,'phnt'));
                                    if isempty(cols.(phaseName).phonoN_target)
                                        phonoN_target = -1;
                                    end
                                    %cols.(phaseName).con_cat_target = find(strcmp(ns_evt_cols,'ccat')); %string
                                    cols.(phaseName).AOA_target = find(strcmp(ns_evt_cols,'aoat'));
                                    if isempty(cols.(phaseName).AOA_target)
                                        AOA_target = -1;
                                    end
                                    cols.(phaseName).BFRQ_target = find(strcmp(ns_evt_cols,'bfqt'));
                                    if isempty(cols.(phaseName).BFRQ_target)
                                        BFRQ_target = -1;
                                    end
                                    cols.(phaseName).FAM_target = find(strcmp(ns_evt_cols,'famt'));
                                    if isempty(cols.(phaseName).FAM_target)
                                        FAM_target = -1;
                                    end
                                    cols.(phaseName).IMG_target = find(strcmp(ns_evt_cols,'imgt'));
                                    if isempty(cols.(phaseName).IMG_target)
                                        IMG_target = -1;
                                    end
                                    cols.(phaseName).CMEAN_target = find(strcmp(ns_evt_cols,'cmnt'));
                                    if isempty(cols.(phaseName).CMEAN_target)
                                        CMEAN_target = -1;
                                    end
                                    cols.(phaseName).PMEAN_target = find(strcmp(ns_evt_cols,'pmnt'));
                                    if isempty(cols.(phaseName).PMEAN_target)
                                        PMEAN_target = -1;
                                    end
                                    cols.(phaseName).TLFRQ_target = find(strcmp(ns_evt_cols,'tfrt'));
                                    if isempty(cols.(phaseName).TLFRQ_target)
                                        TLFRQ_target = -1;
                                    end
                                    %cols.(phaseName).valence_target = find(strcmp(ns_evt_cols,'valt')); %string
                                    %cols.(phaseName).response_target = find(strcmp(ns_evt_cols,'rspt')); %string
                                    cols.(phaseName).wordid_target = find(strcmp(ns_evt_cols,'widt'));
                                    if isempty(cols.(phaseName).wordid_target)
                                        wordid_target = -1;
                                    end
                                    cols.(phaseName).CNC_prime = find(strcmp(ns_evt_cols,'cncp'));
                                    if isempty(cols.(phaseName).CNC_prime)
                                        CNC_prime = -1;
                                    end
                                    cols.(phaseName).KFFRQ_prime = find(strcmp(ns_evt_cols,'frqp'));
                                    if isempty(cols.(phaseName).KFFRQ_prime)
                                        KFFRQ_prime = -1;
                                    end
                                    cols.(phaseName).NLET_prime = find(strcmp(ns_evt_cols,'lenp'));
                                    if isempty(cols.(phaseName).NLET_prime)
                                        NLET_prime = -1;
                                    end
                                    cols.(phaseName).NPHON_prime = find(strcmp(ns_evt_cols,'phop'));
                                    if isempty(cols.(phaseName).NPHON_prime)
                                        NPHON_prime = -1;
                                    end
                                    cols.(phaseName).NSYLL_prime = find(strcmp(ns_evt_cols,'sylp'));
                                    if isempty(cols.(phaseName).NSYLL_prime)
                                        NSYLL_prime = -1;
                                    end
                                    cols.(phaseName).orthoN_prime = find(strcmp(ns_evt_cols,'ornp'));
                                    if isempty(cols.(phaseName).orthoN_prime)
                                        orthoN_prime = -1;
                                    end
                                    cols.(phaseName).phonoN_prime = find(strcmp(ns_evt_cols,'phnp'));
                                    if isempty(cols.(phaseName).phonoN_prime)
                                        phonoN_prime = -1;
                                    end
                                    %cols.(phaseName).con_cat_prime = find(strcmp(ns_evt_cols,'ccap')); %string
                                    cols.(phaseName).AOA_prime = find(strcmp(ns_evt_cols,'aoap'));
                                    if isempty(cols.(phaseName).AOA_prime)
                                        AOA_prime = -1;
                                    end
                                    cols.(phaseName).BFRQ_prime = find(strcmp(ns_evt_cols,'bfqp'));
                                    if isempty(cols.(phaseName).BFRQ_prime)
                                        BFRQ_prime = -1;
                                    end
                                    cols.(phaseName).FAM_prime = find(strcmp(ns_evt_cols,'famp'));
                                    if isempty(cols.(phaseName).FAM_prime)
                                        FAM_prime = -1;
                                    end
                                    cols.(phaseName).IMG_prime = find(strcmp(ns_evt_cols,'imgp'));
                                    if isempty(cols.(phaseName).IMG_prime)
                                        IMG_prime = -1;
                                    end
                                    cols.(phaseName).CMEAN_prime = find(strcmp(ns_evt_cols,'cmnp'));
                                    if isempty(cols.(phaseName).CMEAN_prime)
                                        CMEAN_prime = -1;
                                    end
                                    cols.(phaseName).PMEAN_prime = find(strcmp(ns_evt_cols,'pmnp'));
                                    if isempty(cols.(phaseName).PMEAN_prime)
                                        PMEAN_prime = -1;
                                    end
                                    cols.(phaseName).TLFRQ_prime = find(strcmp(ns_evt_cols,'tfrp'));
                                    if isempty(cols.(phaseName).TLFRQ_prime)
                                        TLFRQ_prime = -1;
                                    end
                                    %cols.(phaseName).valence_prime = find(strcmp(ns_evt_cols,'valp')); %string
                                    %cols.(phaseName).response_prime = find(strcmp(ns_evt_cols,'rspp')); %string
                                    cols.(phaseName).wordid_prime = find(strcmp(ns_evt_cols,'widp'));
                                    if isempty(cols.(phaseName).wordid_prime)
                                        wordid_prime = -1;
                                    end
                                    cols.(phaseName).word_pair_id = find(strcmp(ns_evt_cols,'pair'));
                                    if isempty(cols.(phaseName).word_pair_id)
                                        word_pair_id = -1;
                                    end
                                    %cols.(phaseName).modality = find(strcmp(ns_evt_cols,'mody')); %string
                                    %cols.(phaseName).study_task = find(strcmp(ns_evt_cols,'task')); %string
                                    
                            end
                            
                            % Critical: set up the stimulus type, as well as the
                            % event string to match eventValues
                            if strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'11') && strcmp(phaseName,'TC_NEMO_AO') && strcmp(ns_evt{cols.(phaseName).accuracy+1}(trspInd),'1') && (str2double(ns_evt{cols.(phaseName).random_number+1}(trspInd)) <= 150)
                                evVal = 'ao_standard_corr';
                            elseif strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'12') && strcmp(phaseName,'TC_NEMO_AO') && strcmp(ns_evt{cols.(phaseName).accuracy+1}(trspInd),'1')
                                evVal = 'ao_target_corr';
                            elseif strcmp(ns_evt{1}(ec),'prm+') && strcmp(phaseName,'TC_NEMO_fN400study')
                                evVal = 'study_prime';
                            elseif strcmp(ns_evt{1}(ec),'trg+') && strcmp(phaseName,'TC_NEMO_fN400study') && strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'6')
                                evVal = 'study_targ_AA_rel';
                            elseif strcmp(ns_evt{1}(ec),'trg+') && strcmp(phaseName,'TC_NEMO_fN400study') && strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'8')
                                evVal = 'study_targ_AA_unrel';
                            elseif strcmp(ns_evt{1}(ec),'trg+') && strcmp(phaseName,'TC_NEMO_fN400study') && strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'1')
                                evVal = 'study_targ_CA_rel';
                            elseif strcmp(ns_evt{1}(ec),'trg+') && strcmp(phaseName,'TC_NEMO_fN400study') && strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'3')
                                evVal = 'study_targ_CA_unrel';
                            elseif strcmp(ns_evt{1}(ec),'trg+') && strcmp(phaseName,'TC_NEMO_fN400test') && strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'10') && strcmp(ns_evt{cols.(phaseName).accuracy+1}(trspInd),'1')
                                evVal = 'test_targ_A_new_corr';
                            elseif strcmp(ns_evt{1}(ec),'trg+') && strcmp(phaseName,'TC_NEMO_fN400test') && strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'8') && strcmp(ns_evt{cols.(phaseName).accuracy+1}(trspInd),'1')
                                evVal = 'test_targ_AA_unrel_corr';
                            elseif strcmp(ns_evt{1}(ec),'trg+') && strcmp(phaseName,'TC_NEMO_fN400test') && strcmp(ns_evt{cols.(phaseName).cell_label+1}(trspInd),'3') && strcmp(ns_evt{cols.(phaseName).accuracy+1}(trspInd),'1')
                                evVal = 'test_targ_CA_unrel_corr';
                            %else
                                %keyboard
                            end
                            
                            trl_order = cfg.eventinfo.trl_order.(evVal);
                            
                            % find where this event type occurs in the list
                            eventNumber = find(ismember(cfg.trialdef.eventvalue,evVal));
                            if isempty(eventNumber)
                                eventNumber = -1;
                            end
                            
                            if length(eventNumber) == 1 && eventNumber ~= -1
                                % set the times we need to segment before and after the
                                % trigger
                                prestimSec = abs(cfg.eventinfo.prepost(eventNumber,1));
                                poststimSec = cfg.eventinfo.prepost(eventNumber,2);
                                
                                % prestimulus period should be negative
                                prestimSamp = -round(prestimSec * ft_hdr.Fs);
                                poststimSamp = round(poststimSec * ft_hdr.Fs);
                            else
                                fprintf('event number not found for %s!\n',evVal);
                                keyboard
                            end
                            
                            % dynamically assign values to variables, as
                            % named by the field names under
                            % cols.(phaseName)
                            colFn = fieldnames(cols.(phaseName));
                            for fn = 1:length(colFn)
                                if ~exist(colFn{fn},'var')
                                    if isstrprop(ns_evt{cols.(phaseName).(colFn{fn})+1}{trspInd},'digit')
                                        eval(sprintf('%s = str2double(ns_evt{cols.%s.%s+1}{%d});',colFn{fn},phaseName,colFn{fn},trspInd));
                                    else
                                        eval(sprintf('%s = -1;',colFn{fn}))
                                    end
                                end
                            end
                            
                            % add it to the trial definition
                            this_trl = trl_ini;
                            
                            % get the time of this event
                            this_sample = ft_event(i).sample;
                            
                            % prestimulus sample
                            this_trl(1) = this_sample + prestimSamp;
                            % poststimulus sample
                            this_trl(2) = this_sample + poststimSamp;
                            % offset in samples
                            this_trl(3) = prestimSamp;
                            
                            for to = 1:length(trl_order)
                                thisInd = find(ismember(trl_order,trl_order{to}));
                                if ~isempty(thisInd)
                                    if exist(trl_order{to},'var')
                                        this_trl(timeCols + thisInd) = eval(trl_order{to});
                                    else
                                        fprintf('variable %s does not exist!\n',trl_order{to});
                                        keyboard
                                    end
                                end
                            end
                            
                            % put all the trials together
                            trl = cat(1,trl,double(this_trl));
                            
                        end % check the evt event
                        
                    else
                        % the count is off? EEG and evt events don't line up
                        keyboard
                    end
            end
            
        end
    end % for
    fprintf('\n');
end % pha
