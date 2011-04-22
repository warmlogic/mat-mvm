function plotHead(z,z_range,skipelec,colors)
%PLOTHEAD - Plot scalp results on a head with color ranges.
%
%
%

%z_range = [-8.0 -10;
%           8.0 10];

% setup spline file

% from Tim
%splinefile = fullfile('~/eeg','Hydrocel_GSN_128.spl');
%setupfile = fullfile('~/eeg','Hydrocel_GSN_128_1.0_short.sfp');

% from EGI's FTP on April 3, 2009
splinefile = fullfile('~/eeg','GSN_HydroCel_129.spl');
%setupfile = fullfile('~/eeg','GSN_HydroCel_129_short.sfp');

%headplot_mod('setup',readlocs(setupfile,'filetype','sfp'),splinefile);

if ~exist('z','var') | isempty(z)
  % just make numbers only plots

  % left
  figure(1)
  clf
  reset(gcf)
  headplot_mod(zeros(129,1),splinefile,'view',[280 35],'labels',1);
  colormap(gray)
  
  % Right
  figure(2)
  clf 
  reset(gcf)
  headplot_mod(zeros(129,1),splinefile,'view',[80 35],'labels',1);  
  colormap(gray)
  
  return
end



if ~exist('skipelec','var')
  skipelec = [];
end

if ~exist('colors','var') | isempty(colors)
  colors.pos_col = [1 0 0];
  colors.pos_col_end = [1 .8 .8];
  colors.neg_col = [0 0 1];
  colors.neg_col_end = [.8 .8 1];
  colors.filler_col = [1 1 1];
end


% set up chans to remove so we don't spread into neck and face
%bad_chans = [127 22 17 14 126 8 1 125 121 120 114 108 100 95 89 82 74 69 63 56 49 44 39 128 33 26];
%chans_to_remove = [127 126 125 120 17 114 108 100 95 89 82 74 69 63 56 49 44 128];
chans_to_remove = [];
bad_chans = chans_to_remove;

z(bad_chans) = 0;
z(skipelec) = 0;

% make a colormap
cmap_len = 128;

if size(z_range,1) == 2
  % is two ranges
  % make z fall withing ranges
  z(z<z_range(1,2)) = z_range(1,2);
  z(z>z_range(2,2)) = z_range(2,2);

  % pos map
  col_len = ceil(abs(diff(z_range(1,:))/z_range(1,2)) * cmap_len/2)+1;
  pos_map = makecolormap(colors.pos_col,colors.pos_col_end,col_len);
  
  % neg map
  col_len = ceil(abs(diff(z_range(2,:))/z_range(2,2)) * cmap_len/2)+1;
  neg_map = makecolormap(colors.neg_col_end,colors.neg_col,col_len);
  
  % filler and total map
  filler_map = colors.filler_col(ones(cmap_len - (size(pos_map,1)+size(neg_map,1)),1),:);
  tot_map = [pos_map; filler_map; neg_map];
  
  % set the maplimit
  maplimit = [z_range(1,2) z_range(2,2)];
else
  % is one range
  
  % make z fall withing ranges
  z(z<z_range(1)) = z_range(1);
  z(z>z_range(2)) = z_range(2);

  shiftfix = .5;
  z(z<0) = z(z<0)-shiftfix;
  z(z>0) = z(z>0)+shiftfix;
  z_range = [z_range(1)-shiftfix -shiftfix; shiftfix z_range(2)+shiftfix];
%   % pos map
%   col_len = (cmap_len/2)-4;
%   pos_map = makecolormap(pos_col_end,pos_col,col_len);
  
%   % neg map
%   col_len = (cmap_len/2)-4;
%   neg_map = makecolormap(neg_col,neg_col_end,col_len);
  
%   % filler and total map
%   filler_map = filler_col(ones(cmap_len - (size(pos_map,1)+size(neg_map,1)),1),:);
%   tot_map = [neg_map; filler_map; pos_map];

  % pos map
  col_len = ceil(abs(diff(z_range(2,:))/z_range(2,2)) * cmap_len/2)+1;
  pos_map = makecolormap(colors.pos_col_end,colors.pos_col,col_len);
  
  % neg map
  col_len = ceil(abs(diff(z_range(1,:))/z_range(1,1)) * cmap_len/2)+1;
  neg_map = makecolormap(colors.neg_col,colors.neg_col_end,col_len);
  
  % filler and total map
  filler_map = colors.filler_col(ones(cmap_len - (size(pos_map,1)+size(neg_map,1)),1),:);
  tot_map = [neg_map; filler_map; pos_map];
  
  % set the maplimit
  %maplimit = z_range;
  maplimit = [z_range(1,1) z_range(2,2)];
end

%keyboard

%x = 1000;
%tot_map = floor(tot_map.*x)./x;
%for i = 1:length(tot_map(:))
%  tot_map(i) = str2num(sprintf('%0.5g',tot_map(i)));
%end
%tot_map(tot_map < .5) = tot_map(tot_map < .5) + .00001;
%tot_map(tot_map > .5) = tot_map(tot_map > .5) - .00001;

%tot_map = colormap(jet);
%tot_map(60:69,:) = 1;
  

% left
figure(1)
clf
reset(gcf)
%headplot_mod(z,splinefile,'maplimits',[z_range(1,2) z_range(2,2)],'colormap',tot_map,'view',[280 35],'elec_noplot',bad_chans);
headplot_mod(z,splinefile,'maplimits',maplimit,'view',[280 35],'elec_noplot',bad_chans,'colormap',tot_map);
%set(gcf, 'Colormap', tot_map);

% Right
figure(2)
clf 
reset(gcf)
%headplot_mod(z,splinefile,'maplimits',[z_range(1,2) z_range(2,2)],'colormap',tot_map,'view',[80 35],'elec_noplot',bad_chans);
headplot_mod(z,splinefile,'maplimits',maplimit,'view',[80 35],'elec_noplot',bad_chans,'colormap',tot_map);
%set(gcf, 'Colormap', tot_map);


