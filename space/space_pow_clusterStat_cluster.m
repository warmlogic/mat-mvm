function space_pow_clusterStat_cluster(loadDir,loadFile,cfg_ft,cfg_ana)

%% load the data

load(fullfile(loadDir,loadFile));

[dirs] = mm_checkDirs(dirs);

%% decide who to kick out based on trial counts

% Subjects with bad behavior
exper.badBehSub = {{}};
% exper.badBehSub = {{'SPACE001','SPACE008','SPACE011','SPACE017','SPACE019','SPACE039'}};
% exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE039'}};
exper.badBehSub = {{'SPACE001','SPACE008','SPACE017','SPACE019','SPACE030'}};

% SPACE001 - low trial counts
% SPACE008 - did not do task correctly

% SPACE011 - noisy TF data

evToCheck = { ...
  { ...
  { ...
  'img_onePres' ...
  'img_RgH_rc_spac_p2','img_RgH_rc_mass_p2','img_RgH_fo_spac_p2','img_RgH_fo_mass_p2' ...
  %'img_RgH_rc_spac_p1','img_RgH_rc_mass_p1','img_RgH_fo_spac_p1','img_RgH_fo_mass_p1' ...
  %'img_RgM_spac_p1,'img_RgM_mass_p1'','img_RgM_spac_p2','img_RgM_mass_p2' ...
  } ...
  { ...
  'word_onePres' ...
  'word_RgH_rc_spac_p2','word_RgH_rc_mass_p2','word_RgH_fo_spac_p2','word_RgH_fo_mass_p2' ...
  %'word_RgH_rc_spac_p1','word_RgH_rc_mass_p1','word_RgH_fo_spac_p1','word_RgH_fo_mass_p1' ...
  %'word_RgM_spac_p1','word_RgM_mass_p1','word_RgM_spac_p2','word_RgM_mass_p2' ...
  } ...
  } ...
  };

% exclude subjects with low event counts
[exper,ana] = mm_threshSubs_multiSes(exper,ana,10,[],'vert',evToCheck);

%% run it

for lat = 1:size(cfg_ana.latencies,1)
  cfg_ft.latency = cfg_ana.latencies(lat,:);
  for fr = 1:size(cfg_ana.frequencies,1)
    cfg_ft.frequency = cfg_ana.frequencies(fr,:);
    
    [stat_clus] = mm_ft_clusterstatTFR(cfg_ft,cfg_ana,exper,ana,dirs,data_pow);
  end
end

