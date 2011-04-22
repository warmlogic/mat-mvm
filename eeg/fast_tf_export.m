uname = getenv('USER');
dataroot = fullfile('/Users',uname,'data/RSRC/behavioral');
eegroot = strrep(dataroot,'behavioral','eeg');
if ~strcmp(uname,'admin')
  events = loadEvents(fullfile(dataroot,'rsrc_events.mat'),{'admin',uname});
else
  events = loadEvents(fullfile(dataroot,'rsrc_events.mat'));
end
% exclude the blinking subjects
blinkSub = {'RSRC032','RSRC034','RSRC044','RSRC051'};
badSub = {'RSRC047','RSRC030','RSRC040','RSRC049'};
events = filterStruct(events,'~ismember(subject,varargin{1}) & ~ismember(subject,varargin{2})',blinkSub,badSub);
subjects = unique(getStructField(events,'subject'));
fprintf('%d subjects included in this analysis\n',length(subjects));

checkForArt_ms = 1500;

plotit = 0;

saverooteeg = fullfile(eegroot,'data_fast_tf');
if ~exist(saverooteeg,'dir')
  mkdir(saverooteeg);
end

%% get the events
testEvents = filterStruct(events,['ismember(type,varargin{1}) & (artifactMS == -1 | artifactMS >= ',num2str(checkForArt_ms),')'],{'TARG_PRES','LURE_PRES'});

targ_testEvents = filterStruct(testEvents,'ismember(type,varargin{1})',{'TARG_PRES'});
lure_testEvents = filterStruct(testEvents,'ismember(type,varargin{1})',{'LURE_PRES'});

recHit_srcCor = filterStruct(targ_testEvents,'rec_correct == 1 & src_correct == 1');
recHit_srcInc = filterStruct(targ_testEvents,'rec_correct == 1 & src_correct == 0');

%rec_hit = filterStruct(targ_testEvents,'rec_isTarg == 1 & rec_correct == 1');
%rec_miss = filterStruct(targ_testEvents,'rec_isTarg == 1 & rec_correct == 0');
rec_cr = filterStruct(lure_testEvents,'rec_isTarg == 0 & rec_correct == 1');
%rec_fa = filterStruct(lure_testEvents,'rec_isTarg == 0 & rec_correct == 0');


%% set EEG params
sampleRate = GetRateAndFormat(events(1));
resampledRate = sampleRate;
%resampledRate = 200;

durationMS = 2000;
offsetMS = -500;
bufferMS = 1000;
filtfreq = [58 62];
filttype = 'stop';
filtorder = 4;

% plotting with gete_ms
timeStepMS = 1000/resampledRate;
timeMS = [offsetMS:timeStepMS:(durationMS+offsetMS-1)];
timeS = timeMS / 1000;
%timeS = [(offsetMS/1000):(1000/sampleRate)/1000:((durationMS+offsetMS-1)/1000)];

%elecNum = [1:129];

gruber_parietal = [53,54,55,60,61,62,66,67,70,71,72,75,76,77,78,79,83,84,85,86];
gruber_frontal = [4,5,6,9,10,11,12,13,14,15,16,17,18,19,21,22,112];
elecNum = unique([gruber_parietal,gruber_frontal]);

elecNum_str = {};
for el = 1:length(elecNum)
  if ismember(elecNum(el),gruber_parietal)
    elecNum_str{el} = sprintf('PAR%d',elecNum(el));
  elseif ismember(elecNum(el),gruber_frontal)
    elecNum_str{el} = sprintf('FR%d',elecNum(el));
  else
    elecNum_str{el} = sprintf('%d',elecNum(el));
  end
end

%% Set parameters for the fast_tf executable
outfile_fast_tf = fullfile(saverooteeg,'run_fast_tf.sh');
% set the general parameters
fast_tf_exe = 'fast_tf.4.5.4';
ana_type = '--z_score';
blackman_window = 0.05;
wavelet_m = 6;
begin_baseline = -0.4;
end_baseline = -0.1;
begin_analysis = -0.5;
end_analysis = 1.5;
other_params = '--verbose --rewrite';
% set the parameters for gamma
first_freq_gam = 25;
last_freq_gam = 95;
freq_step_gam = 1;
% set the parameters for theta
first_freq_th = 2;
last_freq_th = 12;
freq_step_th = 0.5;

%% Open the fast_tf executable
fid_tf = fopen(outfile_fast_tf,'w');
fprintf(fid_tf,'#!/bin/bash\n');

%% Grab the EEG and save it out
for sub = 1:length(subjects)
    fprintf('Subject: %s\n',subjects{sub});
    
    outfile_sub_RHSC = fullfile(saverooteeg,[subjects{sub},'_RHSC.txt']);
    outfile_sub_RHSI = fullfile(saverooteeg,[subjects{sub},'_RHSI.txt']);
    outfile_sub_RCR = fullfile(saverooteeg,[subjects{sub},'_RCR.txt']);
    
    if ~lockFile(outfile_sub_RHSC)
        continue
    end
    
    recHit_srcCor_sub = filterStruct(recHit_srcCor,'ismember(subject,varargin{1})',subjects{sub});
    recHit_srcInc_sub = filterStruct(recHit_srcInc,'ismember(subject,varargin{1})',subjects{sub});
    recCR_sub = filterStruct(rec_cr,'ismember(subject,varargin{1})',subjects{sub});
    
    % preallocate space to collect EEG
    eeg_RHSC = nan(length(elecNum),length(recHit_srcCor_sub),size(timeMS,2));
    eeg_RHSI = nan(length(elecNum),length(recHit_srcInc_sub),size(timeMS,2));
    eeg_RCR = nan(length(elecNum),length(recCR_sub),size(timeMS,2));
    fprintf('Getting raw EEG\n');
    for el = 1:length(elecNum)
        fprintf('%d',elecNum(el));
        eeg_RHSC(el,:,:) = gete_ms(elecNum(el),recHit_srcCor_sub,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate);
        eeg_RHSI(el,:,:) = gete_ms(elecNum(el),recHit_srcInc_sub,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate);
        fprintf('.');
        eeg_RCR(el,:,:) = gete_ms(elecNum(el),recCR_sub,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate);
    end
    fprintf('\n');
    
    %% Hit, Source Correct
    fprintf('Writing %s...',outfile_sub_RHSC);
    fid = fopen(outfile_sub_RHSC,'w');
    fprintf(fid,'ascii\n');
    %timeSstr = sprintf(repmat('%s',1,length(timeS)),num2str(timeS));
    timeSstr = strrep(strrep(mat2str(timeS,5),'[',''),']','');
    fprintf(fid,'Time %d %s\n',length(timeS),timeSstr);
    fprintf(fid,'Trials %d\n',length(recHit_srcCor_sub));
    fprintf(fid,'Channels %d',length(elecNum));
    for el = 1:length(elecNum)
        fprintf(fid,' %s',elecNum_str{el});
    end
    fprintf(fid,'\n');
    for ev = 1:length(recHit_srcCor_sub)
        for el = 1:length(elecNum)
            eeg_str = strrep(strrep(mat2str(squeeze(eeg_RHSC(el,ev,:))',5),'[',''),']','');
            fprintf(fid,'%s\n',eeg_str);
        end
    end
    fclose(fid);
    fprintf('done.\n');
    
    %% Hit, Source Incorrect
    fprintf('Writing %s...',outfile_sub_RHSI);
    fid = fopen(outfile_sub_RHSI,'w');
    fprintf(fid,'ascii\n');
    %timeSstr = sprintf(repmat('%s',1,length(timeS)),num2str(timeS));
    timeSstr = strrep(strrep(mat2str(timeS,5),'[',''),']','');
    fprintf(fid,'Time %d %s\n',length(timeS),timeSstr);
    fprintf(fid,'Trials %d\n',length(recHit_srcInc_sub));
    fprintf(fid,'Channels %d',length(elecNum));
    for el = 1:length(elecNum)
        fprintf(fid,' %s',elecNum_str{el});
    end
    fprintf(fid,'\n');
    for ev = 1:length(recHit_srcInc_sub)
        for el = 1:length(elecNum)
            eeg_str = strrep(strrep(mat2str(squeeze(eeg_RHSI(el,ev,:))',5),'[',''),']','');
            fprintf(fid,'%s\n',eeg_str);
        end
    end
    fclose(fid);
    fprintf('done.\n');
    
    %% Correct rejections
    fprintf('Writing %s...',outfile_sub_RCR);
    fid = fopen(outfile_sub_RCR,'w');
    fprintf(fid,'ascii\n');
    %timeSstr = sprintf(repmat('%s',1,length(timeS)),num2str(timeS));
    timeSstr = strrep(strrep(mat2str(timeS,5),'[',''),']','');
    fprintf(fid,'Time %d %s\n',length(timeS),timeSstr);
    fprintf(fid,'Trials %d\n',length(recCR_sub));
    fprintf(fid,'Channels %d',length(elecNum));
    for el = 1:length(elecNum)
        fprintf(fid,' %s',elecNum_str{el});
    end
    fprintf(fid,'\n');
    for ev = 1:length(recCR_sub)
        for el = 1:length(elecNum)
            eeg_str = strrep(strrep(mat2str(squeeze(eeg_RCR(el,ev,:))',5),'[',''),']','');
            fprintf(fid,'%s\n',eeg_str);
        end
    end
    fclose(fid);
    fprintf('done.\n');
    
    releaseFile(outfile_sub_RHSC);
    
    %% Write to the fast_tf executable
    % fast_tf inputs and outputs
    fast_tf_input_file_RHSC = outfile_sub_RHSC;
    fast_tf_input_file_RHSI = outfile_sub_RHSI;
    fast_tf_input_file_RCR = outfile_sub_RCR;
    fast_tf_output_file_RHSC_gam = fullfile(saverooteeg,[subjects{sub},'_RHSC_gam']);
    fast_tf_output_file_RHSI_gam = fullfile(saverooteeg,[subjects{sub},'_RHSI_gam']);
    fast_tf_output_file_RCR_gam = fullfile(saverooteeg,[subjects{sub},'_RCR_gam']);
    fast_tf_output_file_RHSC_th = fullfile(saverooteeg,[subjects{sub},'_RHSC_th']);
    fast_tf_output_file_RHSI_th = fullfile(saverooteeg,[subjects{sub},'_RHSI_th']);
    fast_tf_output_file_RCR_th = fullfile(saverooteeg,[subjects{sub},'_RCR_th']);
    % fast_tf: RHSC gam
    fprintf(fid_tf,'%s %s --first_frequency %.2f --last_frequency %.2f --frequency_step %.2f --blackman_win %.2f --wavelet_m %d --begin_baseline %.2f --end_baseline %.2f --begin_analysis %.2f --end_analysis %.2f %s --output_file %s --stdin < %s\n\n',fast_tf_exe,ana_type,first_freq_gam,last_freq_gam,freq_step_gam,blackman_window,wavelet_m,begin_baseline,end_baseline,begin_analysis,end_analysis,other_params,fast_tf_output_file_RHSC_gam,fast_tf_input_file_RHSC);
    % fast_tf: RHSI gam
    fprintf(fid_tf,'%s %s --first_frequency %.2f --last_frequency %.2f --frequency_step %.2f --blackman_win %.2f --wavelet_m %d --begin_baseline %.2f --end_baseline %.2f --begin_analysis %.2f --end_analysis %.2f %s --output_file %s --stdin < %s\n\n',fast_tf_exe,ana_type,first_freq_gam,last_freq_gam,freq_step_gam,blackman_window,wavelet_m,begin_baseline,end_baseline,begin_analysis,end_analysis,other_params,fast_tf_output_file_RHSI_gam,fast_tf_input_file_RHSI);
    % fast_tf: RCR gam
    fprintf(fid_tf,'%s %s --first_frequency %.2f --last_frequency %.2f --frequency_step %.2f --blackman_win %.2f --wavelet_m %d --begin_baseline %.2f --end_baseline %.2f --begin_analysis %.2f --end_analysis %.2f %s --output_file %s --stdin < %s\n\n',fast_tf_exe,ana_type,first_freq_gam,last_freq_gam,freq_step_gam,blackman_window,wavelet_m,begin_baseline,end_baseline,begin_analysis,end_analysis,other_params,fast_tf_output_file_RCR_gam,fast_tf_input_file_RCR);
    % fast_tf: RHSC th
    fprintf(fid_tf,'%s %s --first_frequency %.2f --last_frequency %.2f --frequency_step %.2f --blackman_win %.2f --wavelet_m %d --begin_baseline %.2f --end_baseline %.2f --begin_analysis %.2f --end_analysis %.2f %s --output_file %s --stdin < %s\n\n',fast_tf_exe,ana_type,first_freq_th,last_freq_th,freq_step_th,blackman_window,wavelet_m,begin_baseline,end_baseline,begin_analysis,end_analysis,other_params,fast_tf_output_file_RHSC_th,fast_tf_input_file_RHSC);
    % fast_tf: RHSI th
    fprintf(fid_tf,'%s %s --first_frequency %.2f --last_frequency %.2f --frequency_step %.2f --blackman_win %.2f --wavelet_m %d --begin_baseline %.2f --end_baseline %.2f --begin_analysis %.2f --end_analysis %.2f %s --output_file %s --stdin < %s\n\n',fast_tf_exe,ana_type,first_freq_th,last_freq_th,freq_step_th,blackman_window,wavelet_m,begin_baseline,end_baseline,begin_analysis,end_analysis,other_params,fast_tf_output_file_RHSI_th,fast_tf_input_file_RHSI);
    % fast_tf: RCR th
    fprintf(fid_tf,'%s %s --first_frequency %.2f --last_frequency %.2f --frequency_step %.2f --blackman_win %.2f --wavelet_m %d --begin_baseline %.2f --end_baseline %.2f --begin_analysis %.2f --end_analysis %.2f %s --output_file %s --stdin < %s\n\n',fast_tf_exe,ana_type,first_freq_th,last_freq_th,freq_step_th,blackman_window,wavelet_m,begin_baseline,end_baseline,begin_analysis,end_analysis,other_params,fast_tf_output_file_RCR_th,fast_tf_input_file_RCR);
    
end

% close the fast_tf executable
fclose(fid_tf);

if plotit == 1

  %% Load subjects
  for sub = 1:length(subjects)

    %% Plotting
    
    % RCR
    [tf_RCR_gam,v_t_RCR_gam,v_f_RCR_gam,channels_RCR_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RCR_gam_zscore.lena']));
    figure(1)
    title_str = sprintf('RCR, Parietal (Gamma), %s',subjects{sub});
    h = pcolor(v_t_RCR_gam,v_f_RCR_gam,squeeze(mean(tf_RCR_gam(:,:,ismember(elecNum,gruber_parietal)),3))');
    shading interp
    % imagesc(v_t_RCR_gam,v_f_RCR_gam,squeeze(mean(tf_RCR_gam(:,:,ismember(elecNum,gruber_parietal)),3))');
    % axis xy
    colorbar
    title(title_str);
    [tf_RCR_th,v_t_RCR_th,v_f_RCR_th,channels_RCR_th] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RCR_th_zscore.lena']));
    figure(2)
    title_str = sprintf('RCR, Frontal (Theta), %s',subjects{sub});
    h = pcolor(v_t_RCR_th,v_f_RCR_th,squeeze(mean(tf_RCR_th(:,:,ismember(elecNum,gruber_frontal)),3))');
    shading interp
    % imagesc(v_t_RCR_th,v_f_RCR_th,squeeze(mean(tf_RCR_th(:,:,ismember(elecNum,gruber_parietal)),3))');
    % axis xy
    colorbar
    title(title_str);
    
    % RHSI
    [tf_RHSI_gam,v_t_RHSI_gam,v_f_RHSI_gam,channels_RHSI_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSI_gam_zscore.lena']));
    figure(3)
    title_str = sprintf('RHSI, Parietal (Gamma), %s',subjects{sub});
    h = pcolor(v_t_RHSI_gam,v_f_RHSI_gam,squeeze(mean(tf_RHSI_gam(:,:,ismember(elecNum,gruber_parietal)),3))');
    shading interp
    % imagesc(v_t_RHSI_gam,v_f_RHSI_gam,squeeze(mean(tf_RHSI_gam(:,:,ismember(elecNum,gruber_parietal)),3))');
    % axis xy
    colorbar
    title(title_str);
    [tf_RHSI_th,v_t_RHSI_th,v_f_RHSI_th,channels_RHSI_th] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSI_th_zscore.lena']));
    figure(4)
    title_str = sprintf('RHSI, Frontal (Theta), %s',subjects{sub});
    h = pcolor(v_t_RHSI_th,v_f_RHSI_th,squeeze(mean(tf_RHSI_th(:,:,ismember(elecNum,gruber_frontal)),3))');
    shading interp
    % imagesc(v_t_RHSI_th,v_f_RHSI_th,squeeze(mean(tf_RHSI_th(:,:,ismember(elecNum,gruber_frontal)),3))');
    % axis xy
    colorbar
    title(title_str);
    
    % RHSC
    [tf_RHSC_gam,v_t_RHSC_gam,v_f_RHSC_gam,channels_RHSC_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSC_gam_zscore.lena']));
    figure(5)
    title_str = sprintf('RHSC, Parietal (Gamma), %s',subjects{sub});
    h = pcolor(v_t_RHSC_gam,v_f_RHSC_gam,squeeze(mean(tf_RHSC_gam(:,:,ismember(elecNum,gruber_parietal)),3))');
    shading interp
    % imagesc(v_t_RHSC_gam,v_f_RHSC_gam,squeeze(mean(tf_RHSC_gam(:,:,ismember(elecNum,gruber_parietal)),3))');
    % axis xy
    colorbar
    title(title_str);
    [tf_RHSC_th,v_t_RHSC_th,v_f_RHSC_th,channels_RHSC_th] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSC_th_zscore.lena']));
    figure(6)
    title_str = sprintf('RHSC, Frontal (Theta), %s',subjects{sub});
    h = pcolor(v_t_RHSC_th,v_f_RHSC_th,squeeze(mean(tf_RHSC_th(:,:,ismember(elecNum,gruber_frontal)),3))');
    shading interp
    % imagesc(v_t_RHSC_th,v_f_RHSC_th,squeeze(mean(tf_RHSC_th(:,:,ismember(elecNum,gruber_frontal)),3))');
    % axis xy
    colorbar
    title(title_str);
    
    keyboard
  end


  %% Load subjects
  % get an example
  [tf_RCR_gam,v_t_RCR_gam,v_f_RCR_gam,channels_RCR_gam] = read_lena(fullfile(saverooteeg,[subjects{1},'_RCR_gam_zscore.lena']));
  [tf_RCR_th,v_t_RCR_th,v_f_RCR_th,channels_RCR_th] = read_lena(fullfile(saverooteeg,[subjects{1},'_RCR_th_zscore.lena']));
  % preallocate using the size of the example file
  tf_RCR_gam_all = nan(length(subjects),size(tf_RCR_gam,3),size(tf_RCR_gam,2),size(tf_RCR_gam,1));
  tf_RHSI_gam_all = nan(length(subjects),size(tf_RCR_gam,3),size(tf_RCR_gam,2),size(tf_RCR_gam,1));
  tf_RHSC_gam_all = nan(length(subjects),size(tf_RCR_gam,3),size(tf_RCR_gam,2),size(tf_RCR_gam,1));
  tf_RCR_th_all = nan(length(subjects),size(tf_RCR_th,3),size(tf_RCR_th,2),size(tf_RCR_th,1));
  tf_RHSI_th_all = nan(length(subjects),size(tf_RCR_th,3),size(tf_RCR_th,2),size(tf_RCR_th,1));
  tf_RHSC_th_all = nan(length(subjects),size(tf_RCR_th,3),size(tf_RCR_th,2),size(tf_RCR_th,1));
%   tf_RCR_gam_all = nan(length(subjects),size(tf_RCR_gam,2),size(tf_RCR_gam,1));
%   tf_RHSI_gam_all = nan(length(subjects),size(tf_RCR_gam,2),size(tf_RCR_gam,1));
%   tf_RHSC_gam_all = nan(length(subjects),size(tf_RCR_gam,2),size(tf_RCR_gam,1));
%   tf_RCR_th_all = nan(length(subjects),size(tf_RCR_th,2),size(tf_RCR_th,1));
%   tf_RHSI_th_all = nan(length(subjects),size(tf_RCR_th,2),size(tf_RCR_th,1));
%   tf_RHSC_th_all = nan(length(subjects),size(tf_RCR_th,2),size(tf_RCR_th,1));

  % read the files and average within ROIs
  for sub = 1:length(subjects)
    [tf_RCR_gam,v_t_RCR_gam,v_f_RCR_gam,channels_RCR_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RCR_gam_zscore.lena']));
    [tf_RHSI_gam,v_t_RHSI_gam,v_f_RHSI_gam,channels_RHSI_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSI_gam_zscore.lena']));
    [tf_RHSC_gam,v_t_RHSC_gam,v_f_RHSC_gam,channels_RHSC_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSC_gam_zscore.lena']));

    tf_RCR_gam_all(sub,:,:,:) = permute(tf_RCR_gam,[3 2 1]);
    tf_RHSI_gam_all(sub,:,:,:) = permute(tf_RHSI_gam,[3 2 1]);
    tf_RHSC_gam_all(sub,:,:,:) = permute(tf_RHSC_gam,[3 2 1]);
%     tf_RCR_gam_all(sub,:,:) = squeeze(mean(tf_RCR_gam(:,:,ismember(elecNum,gruber_parietal)),3))';
%     tf_RHSI_gam_all(sub,:,:) = squeeze(mean(tf_RHSI_gam(:,:,ismember(elecNum,gruber_parietal)),3))';
%     tf_RHSC_gam_all(sub,:,:) = squeeze(mean(tf_RHSC_gam(:,:,ismember(elecNum,gruber_parietal)),3))';
    
    [tf_RCR_th,v_t_RCR_th,v_f_RCR_th,channels_RCR_th] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RCR_th_zscore.lena']));
    [tf_RHSI_th,v_t_RHSI_th,v_f_RHSI_th,channels_RHSI_th] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSI_th_zscore.lena']));
    [tf_RHSC_th,v_t_RHSC_th,v_f_RHSC_th,channels_RHSC_th] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSC_th_zscore.lena']));
    
    tf_RCR_th_all(sub,:,:,:) = permute(tf_RCR_th,[3 2 1]);
    tf_RHSI_th_all(sub,:,:,:) = permute(tf_RHSI_th,[3 2 1]);
    tf_RHSC_th_all(sub,:,:,:) = permute(tf_RHSC_th,[3 2 1]);
%     tf_RCR_th_all(sub,:,:) = squeeze(mean(tf_RCR_th(:,:,ismember(elecNum,gruber_frontal)),3))';
%     tf_RHSI_th_all(sub,:,:) = squeeze(mean(tf_RHSI_th(:,:,ismember(elecNum,gruber_frontal)),3))';
%     tf_RHSC_th_all(sub,:,:) = squeeze(mean(tf_RHSC_th(:,:,ismember(elecNum,gruber_frontal)),3))';
  end


  %% Plotting

  % could also try using plot_tf, which is included with fast_tf

  minmax_th = [0 5];
  minmax_gam = [0 1.2];

  % RCR
  figure(1)
  title_str = sprintf('RCR, Parietal (Gamma), GA');
  % %h = pcolor(v_t_RCR_gam,v_f_RCR_gam,squeeze(mean(tf_RCR_gam_all(:,:,ismember(elecNum,gruber_parietal)),3))');
  % h = pcolor(v_t_RCR_gam,v_f_RCR_gam,squeeze(mean(tf_RCR_gam_all,1)));
  % shading interp
  imagesc(v_t_RCR_gam,v_f_RCR_gam,squeeze(mean(mean(tf_RCR_gam_all(:,ismember(elecNum,gruber_parietal),:,:),2),1)),minmax_gam);
  %imagesc(v_t_RCR_gam,v_f_RCR_gam,squeeze(mean(tf_RCR_gam_all,1)),minmax_gam);
  axis xy
  colorbar
  title(title_str);
  figure(2)
  title_str = sprintf('RCR, Frontal (Theta), GA');
  % %h = pcolor(v_t_RCR_th,v_f_RCR_th,squeeze(mean(tf_RCR_th_all,1)));
  % shading interp
  imagesc(v_t_RCR_th,v_f_RCR_th,squeeze(mean(mean(tf_RCR_th_all(:,ismember(elecNum,gruber_frontal),:,:),2),1)),minmax_th);
  %imagesc(v_t_RCR_th,v_f_RCR_th,squeeze(mean(tf_RCR_th_all,1)),minmax_th);
  axis xy
  colorbar
  title(title_str);

  % RHSI
  title_str = sprintf('RHSI, Parietal (Gamma), GA');
  %[tf_RHSI_gam,v_t_RHSI_gam,v_f_RHSI_gam,channels_RHSI_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSI_gam_zscore.lena']));
  figure(3)
  % h = pcolor(v_t_RHSI_gam,v_f_RHSI_gam,squeeze(mean(tf_RHSI_gam_all,1)));
  % shading interp
  imagesc(v_t_RHSI_gam,v_f_RHSI_gam,squeeze(mean(mean(tf_RHSI_gam_all(:,ismember(elecNum,gruber_parietal),:,:),2),1)),minmax_gam);
  %imagesc(v_t_RHSI_gam,v_f_RHSI_gam,squeeze(mean(tf_RHSI_gam_all,1)),minmax_gam);
  axis xy
  colorbar
  title(title_str);
  figure(4)
  title_str = sprintf('RHSI, Frontal (Theta), GA');
  % h = pcolor(v_t_RHSI_th,v_f_RHSI_th,squeeze(mean(tf_RHSI_th_all,1)));
  % shading interp
  imagesc(v_t_RHSI_th,v_f_RHSI_th,squeeze(mean(mean(tf_RHSI_th_all(:,ismember(elecNum,gruber_frontal),:,:),2),1)),minmax_th);
  %imagesc(v_t_RHSI_th,v_f_RHSI_th,squeeze(mean(tf_RHSI_th_all,1)),minmax_th);
  axis xy
  colorbar
  title(title_str);

  % RHSC
  %[tf_RHSC_gam,v_t_RHSC_gam,v_f_RHSC_gam,channels_RHSC_gam] = read_lena(fullfile(saverooteeg,[subjects{sub},'_RHSC_gam_zscore.lena']));
  figure(5)
  title_str = sprintf('RHSC, Parietal (Gamma), GA');
  % h = pcolor(v_t_RHSC_gam,v_f_RHSC_gam,squeeze(mean(tf_RHSC_gam_all,1)));
  % shading interp
  imagesc(v_t_RHSC_gam,v_f_RHSC_gam,squeeze(mean(mean(tf_RHSC_gam_all(:,ismember(elecNum,gruber_parietal),:,:),2),1)),minmax_gam);
  %imagesc(v_t_RHSC_gam,v_f_RHSC_gam,squeeze(mean(tf_RHSC_gam_all,1)),minmax_gam);
  axis xy
  colorbar
  title(title_str);
  figure(6)
  title_str = sprintf('RHSC, Frontal (Theta), GA');
  % h = pcolor(v_t_RHSC_th,v_f_RHSC_th,squeeze(mean(tf_RHSC_th_all,1)));
  % shading interp
  imagesc(v_t_RHSC_th,v_f_RHSC_th,squeeze(mean(mean(tf_RHSC_th_all(:,ismember(elecNum,gruber_frontal),:,:),2),1)),minmax_th);
  %imagesc(v_t_RHSC_th,v_f_RHSC_th,squeeze(mean(tf_RHSC_th_all,1)),minmax_th);
  axis xy
  colorbar
  title(title_str);

  % gamma induced
  windowMS_gam = [210 + abs(offsetMS):330 + abs(offsetMS)];
  % convert milliseconds to samples
  windowSamp_gam = unique(fix(windowMS_gam * (resampledRate/1000)));
  % theta induced
  windowMS_th = [600 + abs(offsetMS):1200 + abs(offsetMS)];
  % convert milliseconds to samples
  windowSamp_th = unique(fix(windowMS_th * (resampledRate/1000)));

  freq_gam = [first_freq_gam:freq_step_gam:last_freq_gam];
  % 35--80
  freqWind_gam = [find(freq_gam == 35):find(freq_gam == 80)];
  freq_th = [first_freq_th:freq_step_th:last_freq_th];
  % 4--7.5
  freqWind_th = [find(freq_th == 4):find(freq_th == 7.5)];
  
  badchan = {};
  inducedGamma_par_allSub = [];
  inducedTheta_fr_allSub = [];
  for sub = 1:length(subjects)
    recHit_srcCor_sub = filterStruct(recHit_srcCor,'ismember(subject,varargin{1})',subjects{sub});
    % get the bad channels for this subject, if there are any
    [pathstr,eegfile] = fileparts(recHit_srcCor_sub(1).eegfile);
    badchan_file = fullfile(eegroot,subjects{sub},'session_0','eeg',[eegfile,'.bad_chan']);
    if exist(badchan_file,'file')
        badchan{sub} = textread(badchan_file,'%d');
    else
        badchan{sub} = [];
    end
    inducedGamma_par_allSub = [inducedGamma_par_allSub;...
                        mean(mean(mean(tf_RCR_gam_all(sub,ismember(elecNum,gruber_parietal) & ~ismember(elecNum,badchan{sub}),freqWind_gam,windowSamp_gam),4),3),2),...
                        mean(mean(mean(tf_RHSI_gam_all(sub,ismember(elecNum,gruber_parietal) & ~ismember(elecNum,badchan{sub}),freqWind_gam,windowSamp_gam),4),3),2),...
                        mean(mean(mean(tf_RHSC_gam_all(sub,ismember(elecNum,gruber_parietal) & ~ismember(elecNum,badchan{sub}),freqWind_gam,windowSamp_gam),4),3),2)];
    
    inducedTheta_fr_allSub = [inducedTheta_fr_allSub;...
                        mean(mean(mean(tf_RCR_th_all(sub,ismember(elecNum,gruber_frontal) & ~ismember(elecNum,badchan{sub}),freqWind_th,windowSamp_th),4),3),2),...
                        mean(mean(mean(tf_RHSI_th_all(sub,ismember(elecNum,gruber_frontal) & ~ismember(elecNum,badchan{sub}),freqWind_th,windowSamp_th),4),3),2),...
                        mean(mean(mean(tf_RHSC_th_all(sub,ismember(elecNum,gruber_frontal) & ~ismember(elecNum,badchan{sub}),freqWind_th,windowSamp_th),4),3),2)];
  end

  inducedGamma_par_err = std(inducedGamma_par_allSub,0,1)/sqrt(size(inducedGamma_par_allSub,1));
  inducedTheta_fr_err = std(inducedTheta_fr_allSub,0,1)/sqrt(size(inducedTheta_fr_allSub,1));
  
  figure(7);bar(inducedGamma_par_allSub);
  title(['Induced Gamma (',num2str(freq_gam(freqWind_gam(1))),':',num2str(freq_gam(freqWind_gam(end))),'Hz): Parietal, ',num2str(windowMS_gam(1) + offsetMS),':',num2str(windowMS_gam(end) + offsetMS),'ms']);
  %set(gca,'XTickLabel',{'CR','H-IS','H-CS'});
  ylabel('z-Transformed Power');
  figure(8);bar(inducedTheta_fr_allSub);
  title(['Induced Theta (',num2str(freq_th(freqWind_th(1))),':',num2str(freq_th(freqWind_th(end))),'Hz): Frontal, ',num2str(windowMS_th(1) + offsetMS),':',num2str(windowMS_th(end) + offsetMS),'ms']);
  %set(gca,'XTickLabel',{'CR','H-IS','H-CS'});
  ylabel('z-Transformed Power');

  [h,p] = ttest(inducedGamma_par_allSub(:,1),inducedGamma_par_allSub(:,2))
  [h,p] = ttest(inducedGamma_par_allSub(:,1),inducedGamma_par_allSub(:,3))
  [h,p] = ttest(inducedGamma_par_allSub(:,2),inducedGamma_par_allSub(:,3))
  [h,p] = ttest(inducedTheta_fr_allSub(:,1),inducedTheta_fr_allSub(:,2))
  [h,p] = ttest(inducedTheta_fr_allSub(:,1),inducedTheta_fr_allSub(:,3))
  [h,p] = ttest(inducedTheta_fr_allSub(:,2),inducedTheta_fr_allSub(:,3))

  inducedGamma_par = mean(inducedGamma_par_allSub,1);
  inducedTheta_fr = mean(inducedTheta_fr_allSub,1);

  figure(9);bar(inducedGamma_par);
  title(['Induced Gamma (',num2str(freq_gam(freqWind_gam(1))),':',num2str(freq_gam(freqWind_gam(end))),'Hz): Parietal, ',num2str(windowMS_gam(1) + offsetMS),':',num2str(windowMS_gam(end) + offsetMS),'ms']);
  set(gca,'XTickLabel',{'CR','H-IS','H-CS'});
  ylabel('z-Transformed Power');
  figure(10);bar(inducedTheta_fr);
  title(['Induced Theta (',num2str(freq_th(freqWind_th(1))),':',num2str(freq_th(freqWind_th(end))),'Hz): Frontal, ',num2str(windowMS_th(1) + offsetMS),':',num2str(windowMS_th(end) + offsetMS),'ms']);
  set(gca,'XTickLabel',{'CR','H-IS','H-CS'});
  ylabel('z-Transformed Power');
  
%   [h,p] = ttest(mean(mean(tf_RCR_gam_all(:,freqWind_gam,windowSamp_gam),3),2),mean(mean(tf_RHSI_gam_all(:,freqWind_gam,windowSamp_gam),3),2))
%   [h,p] = ttest(mean(mean(tf_RCR_gam_all(:,freqWind_gam,windowSamp_gam),3),2),mean(mean(tf_RHSC_gam_all(:,freqWind_gam,windowSamp_gam),3),2))
%   [h,p] = ttest(mean(mean(tf_RHSI_gam_all(:,freqWind_gam,windowSamp_gam),3),2),mean(mean(tf_RHSC_gam_all(:,freqWind_gam,windowSamp_gam),3),2))

%   [h,p] = ttest(mean(mean(tf_RCR_th_all(:,freqWind_th,windowSamp_th),3),2),mean(mean(tf_RHSI_th_all(:,freqWind_th,windowSamp_th),3),2))
%   [h,p] = ttest(mean(mean(tf_RCR_th_all(:,freqWind_th,windowSamp_th),3),2),mean(mean(tf_RHSC_th_all(:,freqWind_th,windowSamp_th),3),2))
%   [h,p] = ttest(mean(mean(tf_RHSI_th_all(:,freqWind_th,windowSamp_th),3),2),mean(mean(tf_RHSC_th_all(:,freqWind_th,windowSamp_th),3),2))
  
end % plotit
