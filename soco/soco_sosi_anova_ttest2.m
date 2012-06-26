% soco sosi anova 2

exper.sosi.subjects = {
  'SOSI001';
  'SOSI002';
  'SOSI003';
  'SOSI004';
  'SOSI005';
  'SOSI006';
  'SOSI007';
  'SOSI008';
  'SOSI009';
  'SOSI010';
  'SOSI011';
  'SOSI012';
  'SOSI013';
  'SOSI014';
  'SOSI015';
  'SOSI016';
  'SOSI017';
  'SOSI018';
  'SOSI019';
  'SOSI020';
  'SOSI021';
  'SOSI022';
  'SOSI023';
  'SOSI024';
  'SOSI025';
  'SOSI026';
  'SOSI027';
  'SOSI028';
  'SOSI029';
  'SOSI030';
  };

exper.soco.subjects = {
  'SOCO001';
  'SOCO002';
  'SOCO003';
  'SOCO004';
  'SOCO005';
  'SOCO006';
  'SOCO007';
  'SOCO008';
  'SOCO009';
  'SOCO010';
  'SOCO011';
  'SOCO012';
  'SOCO013';
  'SOCO014';
  'SOCO015';
  'SOCO016';
  'SOCO017';
  'SOCO018';
  'SOCO019';
  'SOCO020';
  'SOCO021';
  'SOCO022';
  'SOCO023';
  'SOCO024';
  'SOCO025';
  'SOCO026';
  'SOCO027';
  'SOCO028';
  'SOCO029';
  'SOCO030';
  };

%% SC/SI/CR voltages

% SOSI

% 300--500 ms
expName = 'sosi';

% volt.(expName).SC.FC = [];
% volt.(expName).SI.FC = [];
% volt.(expName).CR.FC = [];

volt.(expName).SC.LAS = [-2.6655 -5.3282 2.1336 -0.3436 -2.9893 -0.3194 -3.9798 0.4425 -0.5791 -2.8256 -3.8347 -3.8586 -2.2305 2.0262 -2.2197 -2.7662 -1.6507 2.0416 -4.5649 -1.0875 -0.6422 0.5849 1.5534 -4.9000 -0.6134 -5.4975 -6.1148 0.3261 -1.5831 -6.1013];
volt.(expName).SI.LAS = [-2.9597 -6.2460 3.8133 -0.9311 -4.6546 -0.1542 -4.5073 -1.7020 -4.5547 -4.0719 -1.5209 -5.4462 -3.2213 1.0312 -2.8913 -2.4850 -2.2418 0.7065 -4.8392 -4.0836 -2.8619 0.3222 1.9493 -4.3831 -1.3142 -6.0216 -7.5182 -0.5924 -2.8436 -6.8067];
volt.(expName).CR.LAS = [-3.1463 -4.8294 1.6097 -0.6662 -4.1419 -0.7456 -5.4915 -0.9406 0.2197 -4.3825 -4.7072 -3.6813 -2.6704 0.8253 -2.6492 -4.0503 -2.2382 0.6886 -3.6315 -4.2025 -3.6350 0.5786 -0.6658 -4.8902 -1.0352 -6.7800 -6.9270 -1.4271 -2.1143 -7.1925];

volt.(expName).SC.RAS = [0.3377 -2.5180 0.2310 -4.1359 -3.9904 -0.7922 -3.3995 -0.8920 1.4175 -4.8029 -1.8326 -4.6273 -2.0758 0.2796 -3.2238 -3.4643 -0.6256 3.9810 -6.2437 -1.6409 -0.3349 1.1296 2.0991 -2.5555 -2.6606 -3.5950 -7.6283 -0.7605 -2.9484 -6.6137];
volt.(expName).SI.RAS = [-0.0958 -3.1581 2.7844 -4.4328 -5.6543 -0.2700 -4.1383 -1.8957 -0.0755 -5.0850 0.4226 -5.8351 -2.8625 -0.7956 -4.6455 -3.8094 -0.7139 2.1801 -5.5619 -3.3836 -1.6135 0.6990 2.3016 -2.3992 -2.9074 -4.0964 -8.0409 -2.0038 -3.8746 -7.0142];
volt.(expName).CR.RAS = [-0.7809 -2.3393 -0.0921 -4.6353 -7.0797 -1.7618 -5.2233 -1.9346 0.1590 -6.2723 -2.0139 -4.5546 -2.5343 -1.1607 -3.7629 -4.2927 -1.9780 2.5086 -5.5591 -2.8521 -1.5926 0.3806 0.2042 -3.4856 -2.4147 -5.6158 -8.0076 -2.4421 -3.2616 -9.0359];

volt.(expName).SC.LASRAS = [-1.1639 -3.9231 1.1823 -2.2397 -3.4899 -0.5558 -3.6897 -0.2247 0.4192 -3.8142 -2.8336 -4.2429 -2.1531 1.1529 -2.7217 -3.1152 -1.1381 3.0113 -5.4043 -1.3642 -0.4885 0.8572 1.8262 -3.7278 -1.6370 -4.5462 -6.8715 -0.2172 -2.2657 -6.3575];
volt.(expName).SI.LASRAS = [-1.5278 -4.7020 3.2988 -2.6819 -5.1545 -0.2121 -4.3228 -1.7989 -2.3151 -4.5785 -0.5492 -5.6406 -3.0419 0.1178 -3.7684 -3.1472 -1.4778 1.4433 -5.2006 -3.7336 -2.2377 0.5106 2.1255 -3.3912 -2.1108 -5.0590 -7.7795 -1.2981 -3.3591 -6.9104];
volt.(expName).CR.LASRAS = [-1.9636 -3.5844 0.7588 -2.6508 -5.6108 -1.2537 -5.3574 -1.4376 0.1893 -5.3274 -3.3605 -4.1180 -2.6023 -0.1677 -3.2061 -4.1715 -2.1081 1.5986 -4.5953 -3.5273 -2.6138 0.4796 -0.2308 -4.1879 -1.7249 -6.1979 -7.4673 -1.9346 -2.6880 -8.1142];

% SOCO - collapsed across 2 and 6 colors

% 300--500 ms
expName = 'soco';

% volt.(expName).SC.FC = [];
% volt.(expName).SI.FC = [];
% volt.(expName).CR.FC = [];

volt.(expName).SC.LAS = [-2.7565 -3.6308 -5.9637 -3.3498 -1.2388 -4.5147 -3.8627 -2.8872 0.5251 -5.4992 -5.0357 -2.2605 -2.8184 -2.1762 1.0264 -3.1328 -6.0156 1.5492 1.0797 -3.0264 -1.5892 -3.3534 -1.1989 -1.4894 -3.3634 -3.4462 -3.4250 -2.6073 -3.5232 -7.2087];
volt.(expName).SI.LAS = [-1.9033 -4.7366 -5.8680 -3.1973 -2.5325 -3.3691 -4.2032 -3.0596 -0.7879 -5.6840 -5.5624 -2.5038 -2.1400 -2.6527 1.7467 -4.1106 -7.1018 1.8665 2.4174 -2.5501 -1.5922 -3.1052 -2.1297 -0.8837 -2.6686 -2.1771 -3.1450 -2.0352 -4.0539 -5.3037];
volt.(expName).CR.LAS = [-3.4932 -2.5028 -5.6636 -4.1748 -2.5788 -6.0703 -4.6090 -3.7020 -0.5305 -6.4419 -5.5472 -3.8731 -4.3403 -4.5491 -0.9059 -3.5152 -6.4087 1.1601 0.0664 -4.3241 -1.8569 -4.0746 -2.1387 -2.2385 -4.7074 -4.3602 -3.6576 -1.4023 -3.8580 -7.8712];

volt.(expName).SC.RAS = [-2.9855 -3.2061 -5.6483 -3.4658 -1.3036 -3.3059 -3.6594 -2.3578 -0.5969 -5.3966 -5.4973 -4.1053 -3.1737 -3.8028 1.9464 -1.6682 -5.6184 1.8356 -1.2644 -3.1800 -2.3861 -1.0301 -0.5682 -1.4015 -3.4621 -3.6398 -3.2717 -2.8758 -3.2425 -3.9889];
volt.(expName).SI.RAS = [-2.6287 -4.2339 -5.4389 -3.1846 -2.7872 -3.2653 -0.8183 -3.0445 -1.9271 -4.7211 -5.8858 -4.1216 -3.2702 -4.3057 1.5639 -2.9869 -5.6080 1.5375 -1.3503 -2.4112 -1.9287 -1.5895 -0.8568 -0.4276 -3.4844 -3.6513 -1.5690 -1.7975 -4.3514 -3.4067];
volt.(expName).CR.RAS = [-2.9785 -3.0468 -5.7704 -4.5362 -3.2676 -4.7638 -4.9850 -3.2916 -2.9014 -6.2291 -5.9487 -4.2107 -4.7915 -5.1943 -0.2468 -2.5021 -5.0909 1.2033 -1.5196 -4.9821 -2.6463 -1.8886 -0.3001 -1.5136 -4.2296 -3.7829 -4.0574 -1.4579 -3.7634 -5.3834];

volt.(expName).SC.LASRAS = [-2.8710 -3.4184 -5.8060 -3.4078 -1.2712 -3.9103 -3.7611 -2.6225 -0.0359 -5.4479 -5.2665 -3.1829 -2.9961 -2.9895 1.4864 -2.4005 -5.8170 1.6924 -0.0924 -3.1032 -1.9876 -2.1918 -0.8836 -1.4455 -3.4128 -3.5430 -3.3483 -2.7416 -3.3829 -5.5988];
volt.(expName).SI.LASRAS = [-2.2660 -4.4852 -5.6535 -3.1910 -2.6598 -3.3172 -2.5108 -3.0520 -1.3575 -5.2025 -5.7241 -3.3127 -2.7051 -3.4792 1.6553 -3.5487 -6.3549 1.7020 0.5335 -2.4806 -1.7604 -2.3473 -1.4932 -0.6556 -3.0765 -2.9142 -2.3570 -1.9163 -4.2026 -4.3552];
volt.(expName).CR.LASRAS = [-3.2359 -2.7748 -5.7170 -4.3555 -2.9232 -5.4171 -4.7970 -3.4968 -1.7160 -6.3355 -5.7480 -4.0419 -4.5659 -4.8717 -0.5763 -3.0087 -5.7498 1.1817 -0.7266 -4.6531 -2.2516 -2.9816 -1.2194 -1.8761 -4.4685 -4.0716 -3.8575 -1.4301 -3.8107 -6.6273];

%% familiar FSC/FSI voltages

% SOSI

% 300--500 ms
expName = 'sosi';

volt.(expName).FSC.FC = [-4.0480 -5.0186 1.8212 -2.7003 -2.9184 -2.5009 2.5018 -2.5446 -0.4273 -1.7806 NaN -7.2002 -4.1178 2.7274 -1.9137 -4.5764 -2.9637 -0.0295 -6.7435 -0.8514 -1.3947 1.1006 -1.9086 -3.1241 -1.4716 -3.0964 -8.2993 -3.3915 -3.6727 NaN];
volt.(expName).FSI.FC = [2.4367 -5.5014 4.3549 -3.4854 -0.9532 -2.7307 NaN -3.4793 -1.5675 -2.8826 NaN -6.1934 -3.6526 -0.9732 -4.8026 -3.3881 -1.3236 1.7334 -7.2388 -3.0935 -2.3035 0.1551 1.6114 1.6231 -2.4092 -3.5421 -8.0983 -2.9892 -6.7195 NaN];
volt.(expName).CR.FC = [-0.8646 -4.0959 2.5585 -3.0653 -6.1732 -2.6924 -6.4574 -2.4376 0.3988 -3.9768 -2.9483 -4.7155 -3.6402 0.3930 -3.2524 -4.4515 -3.1205 1.4211 -5.8296 -2.6456 -2.8390 0.6694 -0.4649 -5.0050 -2.2152 -5.8483 -8.2590 -3.3991 -2.7241 -8.4500];

% SOCO - collapsed across 2 and 6 colors

% 300--500 ms
expName = 'soco';

volt.(expName).FSC.FC = [-4.5105 -2.9554 -7.6117 -2.6742 -3.3733 -0.7459 -6.0970 -2.6611 -1.0510 -5.5221 -7.8663 -5.5707 -0.9317 -5.1036 0.4520 -4.8213 -3.4921 3.9350 0.7436 -6.1692 -0.5915 -3.4525 -1.5437 -4.8453 -6.7256 -3.4791 -8.2976 -3.4591 -6.2893 -4.7023];
volt.(expName).FSI.FC = [-2.9845 -2.7535 -7.0479 -4.4131 0.7676 -2.0707 0.5528 0.9121 -2.2745 -3.8917 -8.2680 -7.0659 -1.1652 -3.7401 -0.0782 -3.2067 -3.9667 -1.3659 1.2499 -1.3317 -0.8806 -2.7858 -1.5819 0.8934 -9.6882 -3.3385 -3.5237 -1.1667 -6.4639 -2.3024];
volt.(expName).CR.FC = [-3.7688 -3.4448 -6.1004 -5.4201 -2.9426 -6.3042 -4.5275 -3.9689 -2.1893 -7.4815 -6.3442 -4.8952 -3.5511 -4.9506 -0.4946 -3.6649 -3.8392 1.3088 -1.5914 -5.1497 -1.6142 -2.6199 -1.2392 -2.9047 -4.3003 -3.4864 -5.3930 -1.1361 -4.9784 -6.3902];

%% independent-samples t-test for the main SC/SI FN400 ttest2s

% test the experiment X accuracy condition voltage difference interaction at one ROI

% who's bad
exper.sosi.badSub = {'SOSI001','SOSI007','SOSI011','SOSI030'};
exper.soco.badSub = {'SOCO018','SOCO026'};

roi = 'LAS';
conds = {'SC','SI'};
[h,p,ci,stats] = ttest2((volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (%s-%s) (M=%.2f) vs SOCO (%s-%s) (M=%.2f): t(%d)=%.4f, p=%.10f\n',...
  roi,...
  conds{1},conds{2},...
  mean(volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  conds{1},conds{2},...
  mean(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),...
  stats.df,stats.tstat,p);

roi = 'RAS';
conds = {'SC','SI'};
[h,p,ci,stats] = ttest2((volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (%s-%s) (M=%.2f) vs SOCO (%s-%s) (M=%.2f): t(%d)=%.4f, p=%.10f\n',...
  roi,...
  conds{1},conds{2},...
  mean(volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  conds{1},conds{2},...
  mean(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),...
  stats.df,stats.tstat,p);

% collapse across hemisphere
roi = 'LASRAS';
conds = {'SC','SI'};
[h,p,ci,stats] = ttest2(...
  (volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  (volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (%s-%s) (M=%.2f) vs SOCO (%s-%s) (M=%.2f): t(%d)=%.4f, p=%.10f\n',...
  roi,conds{1},conds{2},...
  mean(volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  conds{1},conds{2},...
  mean(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),...
  stats.df,stats.tstat,p);

roi = 'LASRAS';
conds = {'SC','CR'};
[h,p,ci,stats] = ttest2(...
  (volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  (volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (%s-%s) (M=%.2f) vs SOCO (%s-%s) (M=%.2f): t(%d)=%.4f, p=%.10f\n',...
  roi,conds{1},conds{2},...
  mean(volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  conds{1},conds{2},...
  mean(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),...
  stats.df,stats.tstat,p);

roi = 'LASRAS';
conds = {'SI','CR'};
[h,p,ci,stats] = ttest2(...
  (volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  (volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (%s-%s) (M=%.2f) vs SOCO (%s-%s) (M=%.2f): t(%d)=%.4f, p=%.10f\n',...
  roi,conds{1},conds{2},...
  mean(volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  conds{1},conds{2},...
  mean(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),...
  stats.df,stats.tstat,p);

%% independent-samples t-test for the familiar FSC/FSI ttest2s

% who's bad

exper.sosi.badSub = {
  'SOSI001'
  'SOSI003'
  'SOSI005'
  'SOSI007'
  'SOSI009'
  'SOSI011'
  'SOSI012'
  'SOSI017'
  'SOSI018'
  'SOSI022'
  'SOSI023'
  'SOSI024'
  'SOSI030'
  };

exper.soco.badSub = {
  'SOCO001'
  'SOCO002'
  'SOCO003'
  'SOCO004'
  'SOCO005'
  'SOCO006'
  'SOCO007'
  'SOCO008'
  'SOCO009'
  'SOCO010'
  'SOCO014'
  'SOCO016'
  'SOCO018'
  'SOCO019'
  'SOCO020'
  'SOCO021'
  'SOCO022'
  'SOCO024'
  'SOCO025'
  'SOCO026'
  'SOCO027'
  };

roi = 'FC';
conds = {'FSC','FSI'};
[h,p,ci,stats] = ttest2((volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (%s-%s) (M=%.2f) vs SOCO (%s-%s) (M=%.2f): t(%d)=%.4f, p=%.10f\n',...
  roi,...
  conds{1},conds{2},...
  mean(volt.sosi.(conds{1}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.(conds{2}).(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),...
  conds{1},conds{2},...
  mean(volt.soco.(conds{1}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.(conds{2}).(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),...
  stats.df,stats.tstat,p);
