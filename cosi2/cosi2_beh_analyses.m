
%% initialize randomness

rng('shuffle','twister');

%% data

expName = 'cosi2';

exper.(expName).subjects = {...
  'COSI2008';
  'COSI2009';
  'COSI2010';
  'COSI2012';
  'COSI2013';
  'COSI2015';
  'COSI2016';
  'COSI2017';
  'COSI2018';
  'COSI2019';
  'COSI2020';
  'COSI2021';
  'COSI2022';
  'COSI2023';
  'COSI2024';
  'COSI2025';
  'COSI2026';
  'COSI2027';
  'COSI2028';
  'COSI2029';
  'COSI2030';
  'COSI2032';
  'COSI2033';
  'COSI2034';
  'COSI2035';
  'COSI2036';
  'COSI2037';
  'COSI2038';
  'COSI2039';
  'COSI2040';
  'COSI2042';
  'COSI2043';
  'COSI2044';
  'COSI2045';
  };


% exper.(expName).color.goodSub = {
%     'COSI2012';
%     'COSI2013';
%     'COSI2017';
%     'COSI2018';
%     'COSI2024';
%     'COSI2028';
%     'COSI2032';
%     'COSI2034';
%     'COSI2036';
%     'COSI2037';
%     'COSI2042';
%     'COSI2044';
%     'COSI2045';
%   };
% 
% exper.(expName).side.goodSub = {
%     'COSI2010';
%     'COSI2015';
%     'COSI2019';
%     'COSI2021';
%     'COSI2022';
%     'COSI2023';
%     'COSI2026';
%     'COSI2027';
%     'COSI2030';
%     'COSI2033';
%     'COSI2035';
%     'COSI2040';
%     'COSI2043';
%   };

% exclude for a number of reasons
exper.(expName).badSub = {'COSI2008','COSI2009','COSI2020','COSI2025','COSI2038','COSI2011','COSI2014','COSI2031','COSI2041','COSI2016','COSI2029','COSI2039'};

% % 8, 9, 20, 25: no F responses in one color/side SC/SI bin

% % 11, 14, 31, 41: no session_1

% % 38: potentially bad session_1 (puker)

% % 16, 29 have fewer than 15 trials for Side-SI

% % 39 has fewer than 15 trials for Color-SI

% % exclude based on automatic equating of accuracy
% exper.(expName).color.badSub = exper.(expName).subjects(~ismember(exper.(expName).subjects,exper.(expName).color.goodSub));
% exper.(expName).side.badSub = exper.(expName).subjects(~ismember(exper.(expName).subjects,exper.(expName).side.goodSub));

% color
COSI2_IHR_color = [0.7525, 0.95, 0.83, 0.9025, 0.8575, 0.9075, 0.965, 0.92, 0.765, 0.86, 0.9211, 0.945, 0.79, 0.7825, 0.935, 0.9475, 0.7325, 0.9125, 0.715, 0.975, 0.7975, 0.905, 0.6975, 0.9425, 0.52, 0.7625, 0.7925, 0.9018, 0.9825, 0.9175, 0.6575, 0.8625, 0.8725, 0.755];
COSI2_IFAR_color = [0.13, 0.05, 0.44, 0.155, 0.045, 0.36, 0.055, 0.265, 0.4, 0.59, 0.1217, 0.52, 0.51, 0.16, 0.505, 0.115, 0.25, 0.175, 0.21, 0.12, 0.21, 0.18, 0.04, 0.605, 0.08, 0.12, 0.24, 0.445, 0.07, 0.565, 0.11, 0.83, 0.59, 0.175];

COSI2_SHR_color = [0.6419, 0.9973, 0.6588, 0.7542, 0.7371, 0.5165, 0.9119, 0.7348, 0.7134, 0.569, 0.8187, 0.8021, 0.5671, 0.5513, 0.6649, 0.849, 0.4161, 0.7554, 0.6148, 0.7835, 0.6335, 0.8453, 0.5429, 0.6791, 0.7879, 0.6625, 0.6289, 0.7486, 0.9439, 0.5191, 0.7609, 0.4556, 0.6069, 0.6333];
COSI2_SFAR_color = [0.3595, 0.0026, 0.4259, 0.3516, 0.2262, 0.2099, 0.0933, 0.3422, 0.302, 0.5235, 0.1508, 0.3874, 0.5066, 0.465, 0.4262, 0.1283, 0.3264, 0.4199, 0.3775, 0.2857, 0.4177, 0.2376, 0.295, 0.3474, 0.5046, 0.2828, 0.4873, 0.3543, 0.0711, 0.3424, 0.472, 0.392, 0.4148, 0.4013];

COSI2_IDP_color = [1.8088, 3.2897, 1.1051, 2.3112, 2.7646, 1.684, 3.4101, 2.0331, 0.9758, 0.8528, 2.5788, 1.548, 0.7814, 1.7751, 1.5016, 2.8214, 1.2949, 2.2909, 1.3745, 3.135, 1.6391, 2.2259, 2.2679, 1.3098, 1.4552, 1.8894, 1.5214, 1.43, 3.5841, 1.2248, 1.6322, 0.1375, 0.9107, 1.6249];
COSI2_SDP_color = [0.7234, 5.5809, 0.596, 1.0686, 1.386, 0.8479, 2.6736, 1.0338, 1.0819, 0.1147, 1.9433, 1.1353, 0.1524, 0.2168, 0.6119, 2.1662, 0.238, 0.8939, 0.604, 1.35, 0.549, 1.7306, 0.6466, 0.8577, 0.7876, 0.994, 0.3608, 1.044, 3.0561, 0.4539, 0.7793, 0.1625, 0.4866, 0.5906];

% side
COSI2_IHR_side = [0.7175, 0.9225, 0.8025, 0.85, 0.86, 0.885, 0.9775, 0.8175, 0.72, 0.835, 0.9525, 0.9025, 0.7275, 0.6775, 0.91, 0.9625, 0.6025, 0.91, 0.735, 0.945, 0.7675, 0.8975, 0.5425, 0.9375, 0.47, 0.78, 0.6375, 0.7964, 0.97, 0.8925, 0.75, 0.84, 0.9225, 0.7525];
COSI2_IFAR_side = [0.105, 0.055, 0.325, 0.16, 0.04, 0.39, 0.095, 0.4, 0.485, 0.58, 0.08, 0.55, 0.475, 0.215, 0.575, 0.105, 0.285, 0.18, 0.26, 0.215, 0.27, 0.23, 0.09, 0.57, 0.07, 0.09, 0.3, 0.4112, 0.055, 0.755, 0.11, 0.71, 0.55, 0.185];

COSI2_SHR_side = [0.8129, 0.963, 0.811, 0.9096, 0.9266, 0.8914, 0.9646, 0.6747, 0.5608, 0.5839, 0.9184, 0.8523, 0.6667, 0.7197, 0.8503, 0.9588, 0.5938, 0.7542, 0.7368, 0.9194, 0.7347, 0.8971, 0.844, 0.7895, 0.5876, 0.8491, 0.5691, 0.7436, 0.9378, 0.5778, 0.9252, 0.7063, 0.7033, 0.8367];
COSI2_SFAR_side = [0.1892, 0.0778, 0.2739, 0.1149, 0.1317, 0.4302, 0.0311, 0.1429, 0.35, 0.3815, 0.0378, 0.4162, 0.4255, 0.2662, 0.3164, 0.1623, 0.1858, 0.1297, 0.1479, 0.0677, 0.1875, 0.0815, 0.3148, 0.3189, 0.1978, 0.0588, 0.3636, 0.2675, 0.0513, 0.3277, 0.183, 0.5341, 0.3476, 0.3182];

COSI2_IDP_side = [1.829, 3.0203, 1.3043, 2.0309, 2.831, 1.4797, 3.3152, 1.1592, 0.6204, 0.7722, 3.0747, 1.1703, 0.668, 1.2499, 1.1516, 3.034, 0.8279, 2.2561, 1.2714, 2.3874, 1.3435, 2.0063, 1.4475, 1.3577, 1.4005, 2.1129, 0.8762, 1.0535, 3.479, 0.5496, 1.901, 0.4411, 1.2964, 1.5789];
COSI2_SDP_side = [1.7697, 3.2063, 1.4826, 2.5392, 2.5688, 1.4101, 3.6724, 1.5205, 0.5383, 0.5133, 3.1705, 1.2578, 0.6185, 1.2063, 1.5154, 2.7215, 1.1305, 1.8154, 1.6792, 2.8938, 1.5142, 2.6603, 1.4934, 1.2753, 1.0709, 2.5971, 0.5229, 1.2748, 3.1693, 0.6425, 2.3447, 0.4569, 0.9257, 1.4539];

% RK accuracy
COSI2_RS_color = [0.6404, 0.9988, 0.7619, 0.8486, 0.8644, 0.8431, 0.9891, 0.9195, 0.8197, 0.7692, 0.8874, 0.7993, 0.48, 0.6439, 0.73, 0.9452, 0.7356, 0.7937, 0.7143, 0.7915, 0.8841, 0.8262, 0.7163, 0.8125, 0.9041, 0.8478, 0.5906, 0.7409, 0.9403, 0.717, 0.7589, 0.6154, 0.6238, 0.7081];
COSI2_RO_color = [0.7895, 0.9988, 0.4267, 0.5532, 0.6777, 0.5172, 0.7925, 0.5529, 0.4915, 0.4948, 0.5333, 0.4681, 0.5498, 0.4765, 0.4786, 0.5698, 0.5333, 0.5357, 0.6, 0.9167, 0.5371, 0.6364, 0.6296, 0.6135, 0.493, 0.6667, 0.9988, 0.5122, 0.0013, 0.4675, 0.6327, 0.4921, 0.5949, 0.5119];
COSI2_F_color = [0.4667, 0.9988, 0.5545, 0.5143, 0.5333, 0.5, 0.6441, 0.5102, 0.5781, 0.449, 0.0013, 0.4231, 0.5333, 0.4375, 0.5263, 0.9988, 0.4198, 0.5278, 0.5769, 0.5352, 0.52, 0.7143, 0.4643, 0.5882, 0.4844, 0.5625, 0.5427, 0.6047, 0.75, 0.5556, 0.5392, 0.5079, 0.5, 0.5088];

COSI2_RS_side = [0.8009, 0.9961, 0.9185, 0.984, 0.9836, 0.9744, 0.9903, 0.975, 0.6818, 0.7596, 0.9528, 0.804, 0.6176, 0.8872, 0.8732, 0.9657, 0.8866, 0.918, 0.9244, 0.9721, 0.992, 0.9315, 0.9266, 0.8833, 0.9538, 0.9643, 0.6695, 0.759, 0.9642, 0.7877, 0.981, 0.6818, 0.707, 0.8564];
COSI2_RO_side = [0.84, 0.8276, 0.6383, 0.6429, 0.7727, 0.5556, 0.9556, 0.6842, 0.6316, 0.5496, 0.7368, 0.5135, 0.6765, 0.5732, 0.6263, 0.7947, 0.64, 0.6438, 0.8, 0.9988, 0.7184, 0.8125, 0.6667, 0.6271, 0.619, 0.9431, 0.7, 0.7143, 0.9988, 0.5152, 0.6667, 0.5254, 0.6447, 0.5781];
COSI2_F_side = [0.9988, 0.8235, 0.5333, 0.6765, 0.5294, 0.4921, 0.7895, 0.5914, 0.5, 0.5051, 0.5, 0.5135, 0.4528, 0.5714, 0.6154, 0.9988, 0.5319, 0.5319, 0.5983, 0.62, 0.5063, 0.8, 0.5606, 0.5584, 0.4833, 0.7143, 0.5354, 0.6111, 0.625, 0.5063, 0.5873, 0.4797, 0.5405, 0.5714];

% response bias (Macmillan & Creelman, p. 29)
COSI2_IRB_color = [0.222, 0, -0.4016, -0.1404, 0.3131, -0.4835, -0.1069, -0.3885, -0.2346, -0.6539, -0.1228, -0.8242, -0.4157, 0.1069, -0.7633, -0.2104, 0.027, -0.2109, 0.1192, -0.3925, -0.0132, -0.1976, 0.6167, -0.9212, 0.6775, 0.2303, -0.0544, -0.5767, -0.3163, -0.7761, 0.4104, -1.0229, -0.6829, 0.1221];
COSI2_SRB_color = [-0.0018, 0.0051, -0.1112, -0.1534, 0.0584, 0.3826, -0.0159, -0.1105, -0.0223, -0.1164, 0.0612, -0.2816, -0.0927, -0.0205, -0.12, 0.0511, 0.3309, -0.2448, 0.0101, -0.1091, -0.0668, -0.1512, 0.2157, -0.0364, -0.4053, 0.0777, -0.1486, -0.1482, -0.0601, 0.179, -0.3194, 0.1927, -0.028, -0.0454];

COSI2_IRB_side = [0.3391, 0.0881, -0.1984, -0.021, 0.3352, -0.4605, -0.347, -0.3263, -0.2726, -0.588, -0.1323, -0.7108, -0.2713, 0.1642, -0.7649, -0.2634, 0.1541, -0.2127, 0.0077, -0.4045, -0.0589, -0.2643, 0.617, -0.8552, 0.7755, 0.2843, 0.0863, -0.3022, -0.1413, -0.9651, 0.276, -0.7739, -0.7739, 0.107];
COSI2_SRB_side = [-0.004, -0.183, -0.1402, -0.0689, -0.1662, -0.5291, 0.0288, 0.3073, 0.1161, 0.0449, 0.1911, -0.4173, -0.1215, 0.0212, -0.2799, -0.3757, 0.3281, 0.22, 0.2059, 0.0462, 0.13, 0.0647, -0.2645, -0.1669, 0.314, 0.2662, 0.0873, -0.0171, 0.0479, 0.1251, -0.2684, -0.314, -0.071, -0.2542];

%% averages and SEM


COSI2_RS_color_avg = mean(COSI2_RS_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_RO_color_avg = mean(COSI2_RO_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_F_color_avg = mean(COSI2_F_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

COSI2_RS_color_sem = std(COSI2_RS_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_RO_color_sem = std(COSI2_RO_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_F_color_sem = std(COSI2_F_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

COSI2_RS_side_avg = mean(COSI2_RS_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_RO_side_avg = mean(COSI2_RO_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_F_side_avg = mean(COSI2_F_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

COSI2_RS_side_sem = std(COSI2_RS_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_RO_side_sem = std(COSI2_RO_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_F_side_sem = std(COSI2_F_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));


% DP
COSI2_IDP_color_avg = mean(COSI2_IDP_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SDP_color_avg = mean(COSI2_SDP_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IDP_color_sem = std(COSI2_IDP_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SDP_color_sem = std(COSI2_SDP_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

COSI2_IDP_side_avg = mean(COSI2_IDP_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SDP_side_avg = mean(COSI2_SDP_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IDP_side_sem = std(COSI2_IDP_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SDP_side_sem = std(COSI2_SDP_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));


% HR
COSI2_IHR_color_avg = mean(COSI2_IHR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SHR_color_avg = mean(COSI2_SHR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IHR_color_sem = std(COSI2_IHR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SHR_color_sem = std(COSI2_SHR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

COSI2_IHR_side_avg = mean(COSI2_IHR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SHR_side_avg = mean(COSI2_SHR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IHR_side_sem = std(COSI2_IHR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SHR_side_sem = std(COSI2_SHR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

% FAR
COSI2_IFAR_color_avg = mean(COSI2_IFAR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SFAR_color_avg = mean(COSI2_SFAR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IFAR_color_sem = std(COSI2_IFAR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SFAR_color_sem = std(COSI2_SFAR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

COSI2_IFAR_side_avg = mean(COSI2_IFAR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SFAR_side_avg = mean(COSI2_SFAR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IFAR_side_sem = std(COSI2_IFAR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SFAR_side_sem = std(COSI2_SFAR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

% response bias (Macmillan & Creelman, p. 29)
COSI2_IRB_color_avg = mean(COSI2_IRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SRB_color_avg = mean(COSI2_SRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IRB_color_sem = std(COSI2_IRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SRB_color_sem = std(COSI2_SRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

COSI2_IRB_side_avg = mean(COSI2_IRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SRB_side_avg = mean(COSI2_SRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_IRB_side_sem = std(COSI2_IRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
COSI2_SRB_side_sem = std(COSI2_SRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)))/sqrt(sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));



%% ttest

chanceVec = 0.5*ones(1,sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));
zerosVec = zeros(1,sum(~ismember(exper.(expName).subjects,exper.(expName).badSub)));

% within source condition
fprintf('\nWIR versus chance\n');
[h,p,ci,stats] = ttest(COSI2_RS_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),chanceVec,0.05,'both');
fprintf('COSI2_RS_color (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_RS_color_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_RO_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),chanceVec,0.05,'both');
fprintf('COSI2_RO_color (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_RO_color_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_F_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),chanceVec,0.05,'both');
fprintf('COSI2_F_color (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_F_color_avg,stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(COSI2_RS_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),chanceVec,0.05,'both');
fprintf('COSI2_RS_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_RS_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_RO_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),chanceVec,0.05,'both');
fprintf('COSI2_RO_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_RO_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_F_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),chanceVec,0.05,'both');
fprintf('COSI2_F_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_F_side_avg,stats.df,stats.tstat,p);

fprintf('\nResponse bias versus zeros\n');
[h,p,ci,stats] = ttest(COSI2_IRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),zerosVec,0.05,'both');
fprintf('COSI2_IRB_color (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_IRB_color_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_SRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),zerosVec,0.05,'both');
fprintf('COSI2_SRB_color (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_SRB_color_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_IRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),zerosVec,0.05,'both');
fprintf('COSI2_IRB_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_IRB_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_SRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),zerosVec,0.05,'both');
fprintf('COSI2_SRB_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_SRB_side_avg,stats.df,stats.tstat,p);

% between source conditions

% WIR
fprintf('\nBetween source conditions\n');
fprintf('WIR\n');
[h,p,ci,stats] = ttest(COSI2_RS_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_RS_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_RS_color (M=%.2f) vs COSI2_RS_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_RS_color_avg,COSI2_RS_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_RO_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_RO_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_RO_color (M=%.2f) vs COSI2_RO_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_RO_color_avg,COSI2_RO_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_F_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_F_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_F_color (M=%.2f) vs COSI2_F_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_F_color_avg,COSI2_F_side_avg,stats.df,stats.tstat,p);

% DP
fprintf('DPrime\n');
[h,p,ci,stats] = ttest(COSI2_IDP_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_IDP_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_IDP_color (M=%.2f) vs COSI2_IDP_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_IDP_color_avg,COSI2_IDP_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_SDP_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_SDP_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_SDP_color (M=%.2f) vs COSI2_SDP_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_SDP_color_avg,COSI2_SDP_side_avg,stats.df,stats.tstat,p);

% HR
fprintf('HR\n');
[h,p,ci,stats] = ttest(COSI2_IHR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_IHR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_IHR_color (M=%.2f) vs COSI2_IHR_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_IHR_color_avg,COSI2_IHR_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_SHR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_SHR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_SHR_color (M=%.2f) vs COSI2_SHR_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_SHR_color_avg,COSI2_SHR_side_avg,stats.df,stats.tstat,p);

% FAR
fprintf('FAR\n');
[h,p,ci,stats] = ttest(COSI2_IFAR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_IFAR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_IFAR_color (M=%.2f) vs COSI2_IFAR_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_IFAR_color_avg,COSI2_IFAR_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_SFAR_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_SFAR_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_SFAR_color (M=%.2f) vs COSI2_SFAR_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_SFAR_color_avg,COSI2_SFAR_side_avg,stats.df,stats.tstat,p);

% response bias (Macmillan & Creelman, p. 29)
fprintf('Response Bias\n');
[h,p,ci,stats] = ttest(COSI2_IRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_IRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_IRB_color (M=%.2f) vs COSI2_IRB_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_IRB_color_avg,COSI2_IRB_side_avg,stats.df,stats.tstat,p);
[h,p,ci,stats] = ttest(COSI2_SRB_color(~ismember(exper.(expName).subjects,exper.(expName).badSub)),COSI2_SRB_side(~ismember(exper.(expName).subjects,exper.(expName).badSub)),0.05,'both');
fprintf('COSI2_SRB_color (M=%.2f) vs COSI2_SRB_side (M=%.2f): t(%d)=%.4f, p=%.10f\n',COSI2_SRB_color_avg,COSI2_SRB_side_avg,stats.df,stats.tstat,p);


%% plot WIR

bw_groupnames = {'Rem. Source';'Rem. Other';'Familiar'};
bw_title = 'Proportion of Source Correct responses';
bw_legend = {'Color','Location'};
bw_colormap = 'gray';
bw_data = [COSI2_RS_color_avg, COSI2_RS_side_avg; COSI2_RO_color_avg, COSI2_RO_side_avg; COSI2_F_color_avg, COSI2_F_side_avg];
bw_errors = [COSI2_RS_color_sem, COSI2_RS_side_sem; COSI2_RO_color_sem, COSI2_RO_side_sem; COSI2_F_color_sem, COSI2_F_side_sem];
bw_xlabel = 'RK Response';
bw_ylabel = 'Proportion Correct';

figure
h = barweb(bw_data,bw_errors,[],bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend);
set(h.legend,'Location','NorthEast');
axis([0.5 3.5 0 1]);
publishfig(gcf,0,[],[],[]);
hold on
plot([0.5 3.5], [0.5 0.5],'r--','LineWidth',2); % horiz chance line
%print(gcf,'-dpng','~/Desktop/COSI2_color_side_RS_RO_F_accuracy');
print(gcf,'-depsc2','~/Desktop/COSI2_color_side_RS_RO_F_accuracy');


%% plot DP

bw_groupnames = {'Item';'Source'};
bw_title = '';
bw_legend = {'Color','Location'};
bw_colormap = 'gray';
bw_data = [COSI2_IDP_color_avg, COSI2_IDP_side_avg; COSI2_SDP_color_avg, COSI2_SDP_side_avg];
bw_errors = [COSI2_IDP_color_sem, COSI2_IDP_side_sem; COSI2_SDP_color_sem, COSI2_SDP_side_sem];
bw_xlabel = 'Source Condition';
bw_ylabel = 'd''';

figure
h = barweb(bw_data,bw_errors,[],bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend);
set(h.legend,'Location','NorthEast');
axis([0.5 2.5 0 2.5]);
publishfig(gcf,0,[],[],[]);
%print(gcf,'-dpng','~/Desktop/COSI2_color_side_item_source_DP');
print(gcf,'-depsc2','~/Desktop/COSI2_color_side_item_source_DP');

%% plot HR

bw_groupnames = {'Item';'Source'};
bw_title = '';
bw_legend = {'Color','Location'};
bw_colormap = 'gray';
bw_data = [COSI2_IHR_color_avg, COSI2_IHR_side_avg; COSI2_SHR_color_avg, COSI2_SHR_side_avg];
bw_errors = [COSI2_IHR_color_sem, COSI2_IHR_side_sem; COSI2_SHR_color_sem, COSI2_SHR_side_sem];
bw_xlabel = 'Source Condition';
bw_ylabel = 'Hit Rate';

figure
h = barweb(bw_data,bw_errors,[],bw_groupnames,bw_title,bw_xlabel,bw_ylabel,bw_colormap,[],bw_legend);
set(h.legend,'Location','NorthEast');
axis([0.5 2.5 0 1]);
publishfig(gcf,0,[],[],[]);
%print(gcf,'-dpng','~/Desktop/COSI2_color_side_item_source_HR');
print(gcf,'-depsc2','~/Desktop/COSI2_color_side_item_source_HR');

%% try to equate

nTries = 10000;
dp_alpha = 0.1;
f_alpha = 0.05;

%foundIt = false;

for i = 1:nTries
  fprintf('.');
  subNums = randperm(length(exper.(expName).subjects));
  subNums = subNums(~ismember(exper.(expName).subjects(subNums),exper.(expName).badSub));
  
  colorHalf = subNums(1:length(subNums)/2);
  sideHalf = subNums((length(subNums)/2) + 1:end);
  
  [h,p,ci,stats] = ttest2(COSI2_SDP_color(colorHalf),COSI2_SDP_side(sideHalf),dp_alpha,'both');
  if h == 0
    [h_cs,p_cs,ci_cs,stats_cs] = ttest2(COSI2_F_color(colorHalf),COSI2_F_side(sideHalf),f_alpha,'both');
    if p_cs < f_alpha
      fprintf('Found p > %.2f and p_cs < %.2f after %d tries.\n',dp_alpha,f_alpha,i);
      
      fprintf('ttest2: Color half source d''=%.2f vs Side half source d''=%.2f: t(%d)=%.2f, p=%.3f',mean(COSI2_SDP_color(colorHalf)),mean(COSI2_SDP_side(sideHalf)),stats.df,stats.tstat,p);
      
      fprintf('\nColor half:');
      sort(exper.(expName).subjects(colorHalf))
      fprintf('\nSide half:');
      sort(exper.(expName).subjects(sideHalf))
      
      [h,p,ci,stats] = ttest2(COSI2_F_color(colorHalf),COSI2_F_side(sideHalf),f_alpha,'both');
      fprintf('ttest2: Color half F=%.2f vs Side half F=%.2f: t(%d)=%.2f, p=%.3f\n',mean(COSI2_F_color(colorHalf)),mean(COSI2_F_side(sideHalf)),stats.df,stats.tstat,p);
      
      [h,p,ci,stats] = ttest(COSI2_F_color(colorHalf),repmat(0.5,size(colorHalf)),f_alpha,'both');
      fprintf('ttest: Color half F=%.2f vs chance (%.1f): t(%d)=%.2f, p=%.3f\n',mean(COSI2_F_color(colorHalf)),0.5,stats.df,stats.tstat,p);
      
      [h,p,ci,stats] = ttest(COSI2_F_side(sideHalf),repmat(0.5,size(sideHalf)),f_alpha,'both');
      fprintf('ttest: Side half F=%.2f vs chance (%.1f): t(%d)=%.2f, p=%.3f\n',mean(COSI2_F_side(sideHalf)),0.5,stats.df,stats.tstat,p);
      
      %foundIt = true;
      break
    end
  end
  if i == nTries
    fprintf('Hit the %d-try limit.\n',i);
  end
end

%% ERP voltage

% 300--500 ms

volt.(expName).CCR.FS = [-1.2884, 0, -4.5967, -3.6045, -4.9645, -4.4493, -4.7565, -1.6172, -1.6075, 0.8224, -9.1271, -2.0727, -4.7027, -5.0464, -0.9493, -2.602, -5.4221, -4.0472, -3.6459, -2.9304, -2.0879, 0.2395, -5.3175, -5.2533, -2.5133, -3.243, -2.006, -4.852, -9.8395, -5.9207, -2.2315, -0.3564, -5.8941, -5.2847];
volt.(expName).CSC.FS = [-0.4125, 0, -4.7028, -3.1506, -2.3343, -4.3484, -4.0236, -0.7816, -1.818, 1.0766, -8.4049, -1.9197, -5.002, -4.3702, -0.8544, -1.2249, -4.8364, -3.6872, -2.7809, -3.5655, -1.8578, 0.2407, -3.6284, -4.6321, -1.2411, -3.2208, -1.323, -2.8932, -8.4693, -4.6172, -1.6819, -0.8204, -5.67, -4.068];
volt.(expName).CSI.FS = [-2.1639, 0, -3.7523, -2.7392, -3.2629, -3.3555, -4.4823, -2.3271, -1.1857, 1.147, -7.994, -3.099, -4.9589, -3.5465, -0.6474, -0.6838, -5.669, -2.3199, -1.1612, -2.7583, -1.2883, 0.8482, -4.1347, -4.481, -0.7869, -3.0424, -1.2367, -5.5774, -6.8555, -5.2343, -0.7741, -1.067, -5.7053, -4.7028];
volt.(expName).SCR.FS = [-2.2348, 0, -4.1371, -3.8239, -3.6971, -5.281, -5.2243, -2.1983, -1.2162, 0.0901, -9.4498, -2.7957, -7.4339, -5.1025, -1.8331, -2.1296, -6.6773, -2.9329, -3.57, -4.2448, -1.6552, 0.1186, -4.7742, -5.0761, -2.9338, -3.4451, -0.8937, -5.4631, -9.8632, -6.891, -2.384, -0.6825, -5.5755, -5.5535];
volt.(expName).SSC.FS = [-0.9526, 0, -3.9944, -3.4687, -3.1287, -3.2492, -4.0786, -1.8437, -1.6603, 0.9484, -8.5376, -2.4587, -5.5855, -3.5385, -1.1717, -1.2255, -4.9392, -3.8986, -3.5868, -1.8096, -1.3977, 0.3479, -4.3321, -4.7083, -1.839, -3.3204, -0.6242, -4.2756, -8.553, -5.2318, -1.3143, -0.4296, -4.4697, -4.119];
volt.(expName).SSI.FS = [-2.9904, 0, -4.8704, -3.3873, -3.9378, -3.1685, -3.638, -2.644, -1.758, 0.398, -10.9458, -1.907, -6.3302, -4.6442, -2.9286, -1.6045, -7.2576, -3.0487, -1.9105, -3.817, -1.5237, -0.9332, -4.5767, -4.7417, 0.1426, -6.0763, -1.6397, -3.1922, -11.0631, -5.3712, -1.4436, -0.9802, -5.3382, -3.9048];

volt.(expName).CCR.LAS = [-2.3704, 0, -4.7612, -3.166, -4.3364, -2.741, -3.0887, -1.2462, -1.1104, 1.2818, -7.2874, -4.4556, -4.7677, -5.6545, -0.7217, -0.8159, -5.9479, -3.0676, -2.3662, -1.6586, -1.965, 0.3622, -4.7877, -3.837, -2.6546, -3.9564, -0.9657, -4.1773, -6.5221, -5.4589, -1.8349, -1.0701, -4.5129, -5.2169];
volt.(expName).CSC.LAS = [-1.6846, 0, -4.265, -2.8017, -2.464, -2.5411, -3.1584, -0.9732, -0.8768, 1.2416, -6.4733, -3.5314, -4.9153, -4.435, -1.1063, 0.0239, -5.2881, -2.7717, -1.77, -1.1585, -1.474, 0.2448, -3.6054, -3.7667, -1.174, -4.1261, -0.8872, -2.7491, -6.2609, -3.6693, -1.4551, -0.2456, -4.3046, -3.9339];
volt.(expName).CSI.LAS = [-3.2284, 0, -2.9407, -2.6857, -3.2663, -1.937, -4.3706, -1.721, -0.2578, 1.6395, -6.5646, -4.2303, -4.5757, -4.1004, -1.1293, 0.241, -5.5866, -1.9829, -0.8604, -0.7257, -1.1981, 1.1524, -4.0427, -3.6849, -1.0709, -3.2334, -0.8661, -4.9184, -5.6263, -4.8262, -0.535, -0.5819, -4.1688, -3.8508];
volt.(expName).SCR.LAS = [-3.3876, 0, -3.424, -3.4442, -3.4786, -3.8791, -3.6211, -1.8671, 0.5658, 1.1511, -8.3585, -4.8209, -5.7531, -4.8516, -1.354, -0.8104, -6.5587, -1.8735, -2.7888, -2.6381, -1.2827, -0.3651, -3.9824, -3.9677, -3.015, -4.2109, -0.6681, -4.8378, -6.8589, -7.3409, -2.1384, -0.6326, -4.9367, -4.9845];
volt.(expName).SSC.LAS = [-2.2436, 0, -3.2699, -3.0166, -2.6649, -2.0689, -2.8712, -1.4446, -0.3684, 1.5829, -7.1498, -3.6239, -5.127, -3.9927, -1.2039, -0.0233, -5.0274, -3.5325, -2.5183, -1.0693, -1.2153, 0.8244, -4.0668, -3.5659, -2.3057, -4.0695, -0.3017, -3.7981, -5.999, -5.2135, -1.0847, 0.3134, -3.1664, -4.3144];
volt.(expName).SSI.LAS = [-3.7123, 0, -4.2497, -2.3724, -3.6714, -2.4931, -3.1606, -2.7047, -1.6668, 0.9629, -8.3152, -3.429, -6.1974, -5.4625, -2.642, -0.3314, -7.3773, -3.6186, -1.6486, -2.7305, -1.352, -1.1519, -4.3005, -3.6512, 0.0559, -7.0909, -1.3403, -2.7021, -9.0301, -4.7278, -0.7173, -0.3261, -3.8796, -3.4446];

volt.(expName).CCR.RAS = [-1.1603, 0, -3.4999, -2.8888, -4.0868, -4.6336, -5.134, -1.424, -1.8284, -0.2208, -9.5013, -0.0312, -4.336, -3.6886, -1.1362, -2.9074, -4.7592, -3.2447, -3.9604, -3.2535, -1.7257, 0.1992, -4.6973, -5.5669, -2.3279, -1.9477, -2.0774, -4.4966, -8.6623, -5.948, -2.343, 0.779, -5.7858, -5.1287];
volt.(expName).CSC.RAS = [-0.685, 0, -3.7955, -2.2762, -1.8661, -4.9067, -3.4942, 0.0637, -2.331, -0.1462, -8.9603, -0.5322, -4.7842, -3.4335, -0.5954, -1.6989, -3.7059, -3.0143, -3.0423, -3.7739, -1.5151, 0.1734, -2.807, -4.4164, -1.4185, -1.8918, -1.5111, -3.1792, -7.3922, -4.4589, -1.9669, -0.5467, -5.3547, -4.4736];
volt.(expName).CSI.RAS = [-2.1203, 0, -3.1455, -2.1793, -2.268, -3.5396, -2.5765, -1.5997, -1.4318, 0.1186, -7.5694, -1.3584, -4.388, -2.7465, -0.1606, -0.9713, -4.7093, -1.9265, -1.3276, -2.8908, -1.002, 0.7766, -3.2568, -4.2206, -0.633, -2.3798, -1.211, -5.3137, -6.195, -5.4321, -1.2045, -1.0117, -5.9825, -5.0553];
volt.(expName).SCR.RAS = [-1.7463, 0, -3.4487, -2.9667, -3.1483, -4.9059, -5.2402, -1.4232, -2.2818, -1.7275, -9.5674, -1.0071, -6.9514, -4.1822, -1.7381, -2.4764, -5.7616, -2.4242, -3.3949, -3.9527, -1.3951, 0.2706, -4.6508, -4.9816, -2.0163, -2.704, -1.3403, -5.2054, -9.3307, -6.573, -2.4337, 0.1711, -5.0959, -5.3899];
volt.(expName).SSC.RAS = [-1.0986, 0, -3.4888, -2.6191, -2.5518, -3.586, -4.0424, -1.0986, -2.1411, -0.5603, -9.0557, -1.0202, -5.1746, -2.7676, -0.8979, -1.6175, -4.2464, -2.5399, -3.5213, -2.521, -1.0823, 0.2099, -3.9128, -5.117, -1.2039, -2.3916, -0.9073, -4.3373, -7.9547, -4.7458, -1.5089, -0.6275, -4.9842, -4.5086];
volt.(expName).SSI.RAS = [-3.0398, 0, -3.4858, -3.2789, -3.4427, -3.008, -3.1826, -0.8995, -1.9622, -1.1283, -11.8653, -0.7917, -6.0095, -3.8122, -2.0289, -1.6166, -5.7875, -1.2912, -2.1203, -3.8622, -1.7483, -0.0054, -3.2848, -4.557, -0.6255, -4.257, -1.7174, -3.0666, -9.5797, -5.2969, -2.116, -1.0505, -5.9518, -4.2224];

% 500--800 ms

volt.(expName).CCR.LPS = [3.2454, 0, 5.2366, 5.5052, -0.0053, 0.9116, 4.2207, -1.176, 1.8677, 2.8175, -0.9458, 2.4252, -0.4645, 0.9494, 0.6199, 3.1727, 1.363, 4.3548, 0.7483, 4.1609, 1.0522, -0.7597, 3.6512, -0.757, 0.7816, 5.241, 4.6183, 0.9924, 9.9473, 1.0189, 2.5051, -0.1141, 2.9704, 3.9843];
volt.(expName).CSC.LPS = [4.3031, 0, 4.9954, 7.9285, 1.9241, 3.0933, 5.7211, -1.3879, 4.629, 4.2621, 1.9173, 5.1054, 0.575, 1.7827, 1.964, 3.8226, 1.6188, 5.6388, 2.237, 4.8375, 2.2835, 0.8848, 5.7535, -0.0026, 1.4965, 6.8814, 5.8296, 2.4618, 11.6729, 0.4694, 4.4676, 0.0165, 4.1454, 4.7236];
volt.(expName).CSI.LPS = [3.9067, 0, 4.4189, 6.6989, 0.1825, 2.8606, 2.3126, -2.004, 4.2705, 4.9003, 0.9013, 4.991, 0.6473, 2.7823, 0.8809, 3.6395, 1.8338, 4.2848, 2.7927, 3.9566, 2.0556, 0.2453, 5.6698, -0.5212, 0.2568, 5.3051, 5.8786, 1.1744, 13.7209, 1.6611, 4.4086, 0.9569, 4.0286, 4.5095];
volt.(expName).SCR.LPS = [3.4385, 0, 3.7707, 5.9768, -0.2595, 0.9479, 4.2411, -1.2419, 1.9033, 3.4283, -0.1312, 3.7108, -1.4177, 2.0477, -0.1593, 3.1432, 0.1917, 5.0211, 0.7121, 3.3494, 0.9334, -0.8825, 3.4164, -0.0681, 1.2047, 5.1369, 4.7138, 2.7153, 7.6948, -0.4795, 3.079, 0.3711, 2.6515, 4.2437];
volt.(expName).SSC.LPS = [3.422, 0, 6.556, 7.3797, 1.4836, 2.8426, 5.1429, -2.2452, 4.4615, 3.6947, 2.0886, 5.1251, 2.1672, 2.2368, 3.026, 5.2274, 2.4721, 6.1562, 2.1751, 3.488, 2.1332, 1.0206, 5.8796, -0.3226, 0.8158, 6.9593, 5.4919, 3.8618, 12.1235, 0.4414, 3.8851, 0.7156, 5.1252, 4.3754];
volt.(expName).SSI.LPS = [1.9853, 0, 5.8126, 7.0341, -1.1511, 1.3361, 6.4331, -1.7755, 2.8993, 3.0148, 0.2629, 3.6428, 0.4031, 0.6158, 1.4944, 3.9347, 0.6437, 4.5707, 1.491, 1.0137, 0.9121, -1.5886, 2.5714, -0.3446, 3.4975, 5.8301, 4.2175, 2.0409, 9.8337, 0.2895, 4.1621, 0.2051, 3.6962, 4.3072];

volt.(expName).CCR.RPS = [2.9502, 0, 5.3268, 5.4556, 1.1675, 1.0688, 0.3387, -1.5958, 0.1042, 2.675, -1.1688, 2.3774, -0.0741, 1.5087, 1.1749, 1.7503, 1.0584, 3.1437, 0.2136, 3.2278, 0.7706, -0.229, 4.4832, -0.8728, 0.4042, 2.8157, 3.4157, 0.4338, 8.8154, 3.1994, 1.341, 0.5965, 3.6773, 3.8603];
volt.(expName).CSC.RPS = [3.5573, 0, 4.1308, 6.9712, 4.2638, 1.8031, 1.4305, -1.6955, 1.8633, 3.6165, 0.4958, 3.2605, 0.8949, 1.7414, 3.0916, 2.8617, 1.3562, 3.9001, 2.5535, 3.548, 1.1216, 1.4779, 5.0875, 0.2962, 0.293, 4.9011, 4.2729, 1.2936, 9.2613, 3.0756, 2.2157, -0.3791, 4.0964, 4.156];
volt.(expName).CSI.RPS = [3.5878, 0, 3.5326, 6.5799, 3.5781, 1.3451, 0.0621, -2.1706, 1.5114, 4.7678, 0.5776, 3.6852, 1.7101, 3.1778, 2.2761, 1.8124, 1.0027, 3.0667, 2.9318, 2.7448, 1.3093, 1.2115, 6.1119, -0.4605, -1.9536, 3.763, 4.301, 0.0153, 11.043, 3.5026, 3.2004, -0.1021, 3.8109, 4.8271];
volt.(expName).SCR.RPS = [2.9992, 0, 3.2336, 6.4121, 1.9757, 0.1498, 1.1759, -1.4428, 0.3401, 3.0655, 0.2048, 3.7725, -0.2219, 1.7433, 1.1652, 1.8891, -0.3466, 2.5914, 0.2816, 3.042, 0.1578, 0.236, 3.627, 0.4351, 0.0715, 2.5819, 3.4316, 1.2884, 5.2236, 1.9573, 0.9754, -0.808, 3.9427, 4.5991];
volt.(expName).SSC.RPS = [3.3974, 0, 5.9698, 6.9839, 3.4782, 1.0775, 1.4722, -2.4847, 0.8171, 2.8095, 1.3575, 3.2246, 1.882, 2.5234, 3.7938, 3.0727, 2.1271, 5.4789, 1.4379, 3.5162, 1.282, 1.1881, 4.9821, 0.3201, 0.0476, 3.6448, 3.6772, 2.5452, 8.9924, 2.8356, 2.2072, -0.2349, 4.7905, 4.4903];
volt.(expName).SSI.RPS = [1.363, 0, 5.1858, 5.6884, 1.9584, 1.2621, 2.7422, -1.2095, 1.8516, 2.3296, -2.19, 2.7046, 0.7871, 1.0299, 1.88, 2.4482, -0.3276, 5.639, 0.7248, 2.4306, 0.6561, 1.1662, 3.1022, 1.8071, 1.1917, 3.3202, 3.0473, 1.5377, 8.1767, 1.7525, 2.267, -1.0603, 2.8775, 4.3575];


%% Set up the rest of the ANOVA

% IMPORTANT: This requires an equal number of subjects in each experiment

% IV1 (between-subjects factor): Experiment
exper.name = {'cosi2'};
% IV2 (within-subjects RM factor): ROI
exper.roi = {'LAS','RAS'};
% IV3 (within-subjects RM factor): Accuracy Condition
exper.cond = {'CCR','CSC','CSI'};
%exper.cond = {'SCR','SSC','SSI'};

% anovamat - data matrix (Size of matrix must be n-by-5;dependent variable=column 1;
%            independent variable 1=column 2;independent variable 2 (within-subjects)=column 3;
%            independent variable 3 (within-subjects)=column 4; subject=column 5 [Be careful
%            how the code for the subjects are entered.]) 

anovamat = [];

for e = 1:length(exper.name)
  goodSubInd = 0;
  for s = 1:length(exper.(exper.name{e}).subjects)
    if ~ismember(exper.(exper.name{e}).subjects(s),exper.(exper.name{e}).badSub)
      goodSubInd = goodSubInd + 1;
      for r = 1:length(exper.roi)
        for c = 1:length(exper.cond)
          anovamat = [anovamat; volt.(exper.name{e}).(exper.cond{c}).(exper.roi{r})(s) e r c goodSubInd];
        end
      end
    end
  end
end

alpha = 0.05;
%showtable = 1;
%calcGGHF = 0;

%[P] = RMAOV32_mod(anovamat,alpha,showtable,calcGGHF);
[P] = RMAOV32(anovamat,alpha);


%% dependent-samples t-test (Within an experiment)

alpha = 0.05;

rois = {'LAS','RAS'};
expNames = {'cosi2'};
condNames = {'color','side'};
%condNames = {};
condPairsAll = {{{'CSC','CSI'},{'CSC','CCR'},{'CSI','CCR'}},{{'SSC','SSI'},{'SSC','SCR'},{'SSI','SCR'}}};
%condPairs = {{'SSC','CSI'},{'SSC','SCR'},{'SSI','SCR'}};

for e = 1:length(expNames)
  expName = expNames{e};
  for r = 1:length(rois)
    roi = rois{r};
    for cn = 1:length(condNames)
      condPairs = condPairsAll{cn};
      for c = 1:length(condPairs)
        conds = condPairs{c};
        [h,p,ci,stats] = ttest(...
          volt.(expName).(conds{1}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).(condNames{cn}).badSub)),...
          volt.(expName).(conds{2}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).(condNames{cn}).badSub)),...
          alpha,'both');
        fprintf('ttest: %s: %s: %s: %s (M=%.2f) vs %s (M=%.2f): t(%d)=%.4f, p=%.10f',roi,expName,(condNames{cn}),...
          conds{1},mean(volt.(expName).(conds{1}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).(condNames{cn}).badSub))),...
          conds{2},mean(volt.(expName).(conds{2}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).(condNames{cn}).badSub))),...
          stats.df,stats.tstat,p);
        if p < alpha
          fprintf(' *');
        end
        fprintf('\n');
      end % c
    end % cn
  end
end

%% do an nchoosek count of F comparisons




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOCO/SOSI method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose based on source d' top (SOCO) or bottom (SOSI) nSubsubset subjects
nSubSubset = 15;

dp_alpha = 0.1;
f_alpha = 0.05;

% excluded based on trial counts, etc.
exper.(expName).badBehSub = {'COSI2008','COSI2009','COSI2020','COSI2025', 'COSI2011', 'COSI2014', 'COSI2031', 'COSI2041', 'COSI2016', 'COSI2029', 'COSI2039', 'COSI2038'};

% worst to best
[~,SOCO_sort] = sort(SOCO_SDP);
[~,SOSI_sort] = sort(SOSI_SDP);

% change SOCO to best to worst
SOCO_sort = fliplr(SOCO_sort);

% order the subjects from best to worst
exper.soco.subSort = exper.soco.subjects(SOCO_sort);
% order the subjects from worst to best
exper.sosi.subSort = exper.sosi.subjects(SOSI_sort);

% exclude the bad behavior ones
exper.soco.subSort = exper.soco.subSort(~ismember(exper.soco.subSort,exper.soco.badBehSub));
exper.sosi.subSort = exper.sosi.subSort(~ismember(exper.sosi.subSort,exper.sosi.badBehSub));

% best to worst
exper.soco.subSubset = exper.soco.subSort(1:nSubSubset);
% worst to best
exper.sosi.subSubset = exper.sosi.subSort(1:nSubSubset);

% figure out which subjects to exclude from the subset
exper.soco.badSub = exper.soco.subjects(~ismember(exper.soco.subjects,exper.soco.subSubset));
exper.sosi.badSub = exper.sosi.subjects(~ismember(exper.sosi.subjects,exper.sosi.subSubset));

% get just the data we want
SOCO_SDP_SUBSET = SOCO_SDP(ismember(exper.soco.subjects,exper.soco.subSubset));
SOSI_SDP_SUBSET = SOSI_SDP(ismember(exper.sosi.subjects,exper.sosi.subSubset));
SOCO_F_WIR_SUBSET = SOCO_F_WIR(ismember(exper.soco.subjects,exper.soco.subSubset));
SOSI_F_WIR_SUBSET = SOSI_F_WIR(ismember(exper.sosi.subjects,exper.sosi.subSubset));

%% start at 1 fewer than the number of subjects in the subset

%nSub = nSubSubset - 1;
nSub = 15;
minSub = 10;
nIter = (nSub - minSub + 1);

results = struct;
results.nSubSubset = nSubSubset;

counter = 0;

while nSub >= minSub
  counter = counter + 1;
  results.nSub(counter) = nSub;
  
  % calculate the unique combinations of subjects
  allCombos = single(nchoosek(1:nSubSubset,nSub));
  
  fprintf('Testing %d combinations for nSub=%d (out of the %d best SOCO & worst SOSI subjects) (%d of %d)...\n',size(allCombos,1),nSub,nSubSubset,counter,nIter);
  
  results.nCombos(counter) = size(allCombos,1);
  
  % initialize
%   results.DPh{counter} = false(size(allCombos,1),1);
%   results.DPp{counter} = nan(size(allCombos,1),1);
%   results.DPt{counter} = nan(size(allCombos,1),1);
%   results.Fh{counter} = false(size(allCombos,1),1);
%   results.Fp{counter} = nan(size(allCombos,1),1);
%   results.Ft{counter} = nan(size(allCombos,1),1);
  
  results.DPh{counter} = false(size(allCombos,1),size(allCombos,1));
  results.DPp{counter} = nan(size(allCombos,1),size(allCombos,1));
  results.DPt{counter} = nan(size(allCombos,1),size(allCombos,1));
  results.Fh{counter} = false(size(allCombos,1),size(allCombos,1));
  results.Fp{counter} = nan(size(allCombos,1),size(allCombos,1));
  results.Ft{counter} = nan(size(allCombos,1),size(allCombos,1));
  
  for i = 1:size(allCombos,1)
  for j = 1:size(allCombos,1)
    
    if mod(i,1000) == 0
      fprintf('.');
    end
    
    % test the d' difference
    [h,p,ci,stats] = ttest2(SOCO_SDP_SUBSET(allCombos(i,:)),SOSI_SDP_SUBSET(allCombos(j,:)),dp_alpha,'both');
    % store the results
    results.DPh{counter}(i) = h;
    results.DPp{counter}(i) = p;
    results.DPt{counter}(i) = stats.tstat;
    
    % if d' was equal
    if h == 0
      
      % test the F accuracy difference
      [h,p,ci,stats] = ttest2(SOCO_F_WIR_SUBSET(allCombos(i,:)),SOSI_F_WIR_SUBSET(allCombos(j,:)),f_alpha,'both');
      
      % store the results
      results.Fh{counter}(i) = h;
      results.Fp{counter}(i) = p;
      results.Ft{counter}(i) = stats.tstat;
      
    end
  end
  end
  fprintf('\nDone with nSub=%d.\n',nSub);
  nSub = nSub - 1;
end

save('~/Desktop/soco_sosi_permute_results.m','results');

%%

for i = 1:nIter
  fprintf('%d top (color) and bottom (side) subjects:\n',results.nSub(i));
  FtestedInd = ~isnan(results.Fp{i});
  %FlessInd = results.Fp{i}(~isnan(results.Fp{i})) < f_alpha;
  FlessInd = results.Fh{i};
  fprintf('\tColor vs Side source d'': %d of %d t-tests had p>%.2f (%.2f%%)\n',sum(FtestedInd),results.nCombos(i),dp_alpha,(sum(FtestedInd)/results.nCombos(i)*100));
  fprintf('\tColor vs Side F accuracy: %d of %d the null d'' tests had p<%.2f (%.2f%%)\n',sum(FlessInd),sum(FtestedInd),f_alpha,(sum(FlessInd)/sum(FtestedInd)*100));
end
