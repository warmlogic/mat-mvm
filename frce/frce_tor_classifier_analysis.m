%save matt_data_flat data_flat

dataroot = fullfile(getenv('HOME'),'Downloads','FRCE_data_classification');

%load matt_data_flat

% time within freq within channel
% 73 x 9 x 20

% note which rows in the colIndex represent which type of data
rows.chan = 1;
rows.freq = 2;
rows.time = 3;

% just to see what it looks like
%data_flat.subjects(1).colIndex; % channel, freq, time

%%

% spider parameters
algorithm_name = 'cv_svm';
nfolds = 5;

% initialize
accuracy = nan(length(cfg_prep.subjs),length(cfg_prep.freqs));

for sub = 1:length(cfg_prep.subjs)
  for frq = 1:size(cfg_prep.freqs,1)
    
    X = double(data_flat.subjects(sub).data);
    Y = data_flat.subjects(sub).outcome;  % 1 is recalled, -1 is forgotten
    
    % select features/freqs
    
    wh = data_flat.subjects(sub).colIndex(rows.freq, :) == frq; % alpha power
    X = X(:, wh);
    
    %
    
    dat = fmri_data;
    dat.dat = X';
    dat.Y = Y;
    
    %[cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_lassopcr','lasso_num', 5, 'nfolds', 5, 'error_type', 'mse')
    
    %STATS = xval_regression_multisubject('lasso', {Y}, {X}, 'pca', 'ndims', 'variable', 'holdout_method', 'balanced4');
    
    % [cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_pcr', 'nfolds', 5, 'error_type', 'mse');
    % pred_outcome_r = stats.pred_outcome_r
    
    [cverr, stats, optout] = predict(dat, 'algorithm_name', algorithm_name, 'nfolds', nfolds, 'error_type', 'mse');
    
    accuracy(sub,frq) = sum(stats.Y == stats.yfit) / length(stats.Y);
    fprintf('\n%s, %d-%d Hz, accuracy = %.2f',cfg_prep.subjs{sub},cfg_prep.freqs(frq,1),cfg_prep.freqs(frq,2),accuracy(sub,frq));
    if accuracy(sub,frq) > 0.5
      fprintf(' ***\n\n\n');
    else
      fprintf('\n\n\n');
    end
    
    fprintf('\t\tPredicted\n');
    fprintf('\t\t%s\t%s\n','VisForg','VisReca');
    fprintf('True\t%s\t%d\t%d\n','VisForg',sum((stats.Y == -1) & (stats.yfit == -1)),sum((stats.Y == -1) & (stats.yfit ~= -1)));
    fprintf('\t%s\t%d\t%d\n','VisReca',sum((stats.Y == 1) & (stats.yfit == 1)),sum((stats.Y == 1) & (stats.yfit ~= 1)));
    
  end
end

%% print out the results

fid = fopen(fullfile(dataroot,sprintf('spider_results_%s_%dfolds_%dsubjs_%dfreqs_%dchans.txt',algorithm_name,nfolds,length(cfg_prep.subjs),length(cfg_prep.freqs),length(cfg_prep.chans))),'w+');

fprintf(fid,'%s',sprintf(repmat('\t%d-%d Hz',1,size(cfg_prep.freqs,1)),cfg_prep.freqs'));
fprintf(fid,'\tAverage\n');
for sub = 1:length(cfg_prep.subjs)
  fprintf(fid,'%s%s',cfg_prep.subjs{sub},sprintf(repmat('\t%.2f',1,size(accuracy,2)),accuracy(sub,:)));
  fprintf(fid,'\t%.2f\n',mean(accuracy(sub,:),2));
end
fprintf(fid,'Average%s',sprintf(repmat('\t%.2f',1,size(accuracy,2)),mean(accuracy,1)));
fprintf(fid,'\t%.2f\n',mean(mean(accuracy,1),2));

fclose(fid);


%% can't run because the files are on Tor's computer
% STATS = xval_regression_multisubject('lasso', {Y}, {X}, 'pca', 'ndims', 20, 'holdout_method', 'balanced4');
% 
% yfit = contrast_code(double(STATS.subjfit{1} > 0) - .5);
% yobs = double(round(STATS.INPUTS.Y{1}));
% accuracy = sum(yobs == yfit) ./ length(yobs)

%% fishing

n = size(dat.dat, 1);
z = ones(size(dat.dat, 2), 1);
Xr = [dat.Y];

[b, p] = deal(zeros(n, 1));

for i = 1:n
    
    [bb, dev, stats] = glmfit(Xr, dat.dat(i, :)');
    b(i, 1) = b(2);
    p(i, 1) = stats.p(2);
    
end

%%
whsig = p < .001;

%
% 
% X = double(data_flat.subjects(1).data);
% Y = data_flat.subjects(1).outcome;  % 1 is correct, 0 is wrong
% 
% % select features/freqs
% 
% wh = data_flat.subjects(1).colIndex(1, :) == 3; % alpha power
% X = X(:, wh);
% 
% X = X(:, whsig);
% 
% Y(Y == 0) = -1;

%

dat = fmri_data;
dat.dat = X';
dat.dat = dat.dat(whsig, :);
dat.Y = Y;

[cverr, stats, optout] = predict(dat, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'error_type', 'mse');
accuracy = sum(stats.Y == stats.yfit) ./ length(stats.Y)

xx = 1:size(dat.dat, 1);
create_figure; plot(xx-.1, dat.dat(:, dat.Y > 0), 'ko');
hold on; plot(xx+.1, dat.dat(:, dat.Y < 0), 'rs');

sc1 = dat.dat(1, dat.Y > 0)';
sc2 = dat.dat(1, dat.Y < 0)';
[H,p,ci,stats] = ttest2_printout(sc1, sc2, 1);

