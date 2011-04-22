function [exper] = mm_ft_concatTrialCounts_cluster(job,exper,allSubjects)
%MM_FT_CONCATTRIALCOUNTS_CLUSTER Concatenate the trial counts for subjects
%and event values on a computer cluster when using distributed jobs
%
% mm_ft_concatTrialCounts_cluster(job,exper,allSubjects)

%% concat exper.eventValues and exper.eventValuesExtra.newValue

eventValuesAll = exper.eventValues;
if ~isempty(exper.eventValuesExtra.newValue)
  for nVal = 1:length(exper.eventValuesExtra.newValue)
    eventValuesAll = cat(2,eventValuesAll,exper.eventValuesExtra.newValue{nVal});
  end
end

%% overwrite exper with the first subject's create_ft_struct output

% option 1
tasks = get(job,'Tasks');
exper = tasks(1).OutputArguments{1};

% % option 2
% % need to use: t(i) = createTask(job,@create_ft_struct,1,inArg);
% taskOutput = get(t(1),'OutputArguments');
% exper = taskOutput{1};

% % option 3
% exper = job.Tasks(1).OutputArguments{1};

% put all subjects in the new exper struct
exper.subjects = allSubjects;

%% get the trial counts for subjects and event values

if length(exper.subjects) > 1
  for i = 2:length(exper.subjects)
    % option 1
    % get exper struct for this subject
    exper_sub = tasks(i).OutputArguments{1};
    
    % % option 2
    % % get the task outputs for this subject as cell array
    % taskOutput = get(t(i),'OutputArguments');
    % % get exper struct for this subject
    % exper_sub = taskOutput{1};
    
    % % option 3
    % % get exper struct for this subject
    % exper_sub = job.Tasks(i).OutputArguments{1};
    
    for evVal = 1:length(eventValuesAll)
      % concatenate the number of trials for each event value across all
      % subjects
      exper.nTrials.(eventValuesAll{evVal}) = cat(1,exper.nTrials.(eventValuesAll{evVal}),exper_sub.nTrials.(eventValuesAll{evVal}));
    end
  end
end

end

