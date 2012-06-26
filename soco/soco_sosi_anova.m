%% subject info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use a subject number style inclusion setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOSI019 and SOSI020 are out of order, but it's OK because of how we index
% in the excel spreadsheet. if we're getting voltages from anywhere else
% then they need to be in order.

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
  'SOSI020'; % out of order
  'SOSI019'; % out of order
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


% % good subjects found in automatic equating
%
% exper.soco.goodSub = {
%   'SOCO003';
%   'SOCO004';
%   'SOCO007';
%   'SOCO009';
%   'SOCO010';
%   'SOCO016';
%   'SOCO019';
%   'SOCO020';
%   'SOCO030';
%   };
% 
% exper.sosi.goodSub = {
%   'SOSI002';
%   'SOSI003';
%   'SOSI008';
%   'SOSI016';
%   'SOSI017';
%   'SOSI018';
%   'SOSI024';
%   'SOSI025';
%   'SOSI029';
%   };
%
% % exclude based on automatic equating of accuracy
% exper.soco.badSub = exper.soco.subjects(~ismember(exper.soco.subjects,exper.soco.goodSub));
% exper.sosi.badSub = exper.sosi.subjects(~ismember(exper.sosi.subjects,exper.sosi.goodSub));


% excluded based on trial counts, etc.
exper.sosi.badSub = {'SOSI001','SOSI007','SOSI011','SOSI030'};
exper.soco.badSub = {'SOCO018','SOCO026'};
% % just trying something out: equating numbers
% exper.soco.badSub = {'SOCO001','SOCO002','SOCO018','SOCO026'};

% % excluded to equate source d' (20 lowest side vs 20 highest color)
% exper.sosi.badSub = {'SOSI001','SOSI007','SOSI011','SOSI030','SOSI008','SOSI023','SOSI028','SOSI027','SOSI018','SOSI009'};
% exper.soco.badSub = {'SOCO018','SOCO026','SOCO022','SOCO001','SOCO006','SOCO029','SOCO024','SOCO012','SOCO015','SOCO030'};

% % excluded to equate source d' (15 lowest side vs 15 highest color)
% exper.sosi.badSub = {'SOSI001','SOSI007','SOSI011','SOSI030','SOSI008','SOSI023','SOSI028','SOSI027','SOSI018','SOSI009','SOSI004','SOSI015','SOSI013','SOSI025','SOSI029'};
% exper.soco.badSub = {'SOCO018','SOCO026','SOCO022','SOCO001','SOCO006','SOCO029','SOCO024','SOCO012','SOCO015','SOCO030','SOCO020','SOCO028','SOCO014','SOCO027','SOCO017'};



% exper.sosi.badSub = [1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
% exper.soco.badSub = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];
% % even it up
% exper.soco.badSub = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];

% %% These voltages are for Net Station + EP + NS preprocessing
% 
% % SOSI
% 
% % 300--500 ms
% 
% volt.sosi.RCR.FS = [-2.0967 -3.8552  0.2780 -3.1306 -3.9418 -1.4288 -5.4773 -1.8390  0.1859 -6.1592 -4.2967 -6.8965 -3.6632 -0.3500 -3.6431 -5.2454 -2.6920  2.1320 -4.7472 -5.2297 -2.8625  0.5506 -0.5511 -4.7820 -1.9545 -6.8076 -8.5398 -0.8121 -2.9559 -7.7213];
% volt.sosi.RHSC.FS = [-1.8586 -3.9610  1.5776 -2.0099 -3.5015 -0.5671 -4.9660 -0.7308  0.9708 -4.0975 -2.5098 -5.7959 -1.6103  1.6809 -2.5177 -3.5936 -1.1317  3.5929 -2.5763 -5.2777 -0.2098  0.4877  1.6143 -3.3757 -2.3006 -5.1245 -7.8329  0.3209 -2.7091 -6.4133];
% volt.sosi.RHSI.FS = [-2.1722 -4.0760  2.3633 -2.4579 -4.8204  0.4342 -4.8582 -2.0847 -1.6530 -4.0514 -3.4606 -5.7041 -3.0563  0.4130 -4.1497 -3.8431 -2.1402  1.5856 -4.8478 -6.2604 -2.5481  0.2427  2.5469 -2.8635 -2.3399 -5.3022 -8.9823 -0.3249 -4.0458 -7.5102];
% 
% volt.sosi.RCR.LAS = [-3.6343 -4.9876  0.5313 -1.3635 -1.6718 -0.8996 -4.5494 -0.9361 -0.4074 -4.1059 -5.2698 -6.5684 -3.1130  0.5545 -3.0541 -4.4020 -1.9895  0.6449 -3.9992 -3.5721 -3.4103  0.6435 -0.9608 -4.9836 -0.7034 -6.9694 -7.0104 -0.7275 -2.5144 -7.1506];
% volt.sosi.RHSC.LAS = [-3.6429 -5.2683  1.6674 -0.1701 -1.8575 -0.2446 -4.5573  0.0842 -0.5038 -2.8244 -3.6307 -4.8766 -1.3853  2.0416 -1.9544 -2.9401 -1.3059  1.8283 -1.2919 -3.8417 -0.2401 -0.4465  1.2029 -3.7441 -0.8174 -5.8618 -6.2737  0.3427 -1.7567 -5.9320];
% volt.sosi.RHSI.LAS = [-3.9067 -5.5388  3.1181 -0.3870 -2.4034  0.2772 -4.3981 -1.6188 -3.2428 -3.6879 -4.0634 -4.9036 -2.9848  0.8786 -2.9823 -2.6717 -2.3664  0.0713 -4.6378 -5.0076 -2.7900 -0.2780  1.8273 -2.8185 -1.2569 -5.7553 -7.1238 -0.2628 -3.0962 -7.9065];
% 
% volt.sosi.RCR.RAS = [-1.5316 -2.7553 -1.0681 -4.6711 -5.1588 -1.8174 -4.5874 -2.2206  0.3705 -7.3105 -2.3667 -5.9091 -2.9707 -1.0661 -3.8570 -5.1215 -1.6633  2.6286 -3.2765 -5.4804 -1.5595 -0.0625 -0.2402 -3.8614 -2.3104 -5.5373 -7.8201 -1.9532 -3.4330 -8.1835];
% volt.sosi.RHSC.RAS = [-1.6477 -2.4646 -0.0756 -4.3286 -3.6775 -0.9126 -3.9043 -1.2330  1.3822 -5.7476 -1.7558 -5.8615 -1.3068  0.4732 -2.9708 -2.8528 -0.1661  3.9909 -2.3085 -5.4341 -0.2759  0.6845  1.7312 -1.9848 -2.9603 -3.9287 -7.9835 -0.7554 -3.3818 -6.5891];
% volt.sosi.RHSI.RAS = [-0.7376 -2.5449  0.9046 -4.4421 -4.8673  0.3592 -3.7991 -1.7546 -0.2275 -4.6733 -2.7591 -5.5627 -3.0013 -0.6326 -4.6172 -4.4857 -0.8014  2.3457 -3.2002 -5.6092 -2.0604  0.4191  2.4630 -1.3916 -2.4760 -4.1783 -7.9587 -1.0994 -4.2752 -7.1604];
% 
% % 500--800 ms
% 
% volt.sosi.RCR.LPS = [5.8344  1.1240  7.5554  8.6566  8.5816  2.0974  2.5554  2.2148  0.7782  9.5262  2.0923  3.1186  3.5499  3.0801  2.2200  1.5636  0.8872  0.7725  5.7253  1.8833 -0.2981  3.2952  4.2548  3.6496  3.4887  1.6152  4.9915  1.0820  1.4559  2.0023];
% volt.sosi.RHSC.LPS = [9.0865  3.2846  7.7397  9.5748  7.0898  1.5327  2.1764  4.4227  3.7799 13.8807  4.0215  3.2167  7.5117  7.1600  3.1100  3.3258  1.7892  1.3158  7.6005  4.2802  1.0361  4.9406  4.7082  7.5131  3.6111  5.3059  6.6378  0.1270  2.3289  3.1684];
% volt.sosi.RHSI.LPS = [7.9848  1.1251  6.5587  8.4441  6.6771  1.2723  1.8165  3.0797  2.3192  9.1071  5.1915  0.3536  5.5413  3.7872  1.7552  2.3160  1.1248 -0.3985  3.9858  2.0816 -0.8110  4.0997  3.9381  2.4545  4.2815  4.7160  4.4488  0.5204  2.2861  1.0885];
% 
% volt.sosi.RCR.RPS = [6.1935  1.7732  5.9138  8.0010  3.3768  0.4649  0.2004  0.9880  3.1322  5.5834  5.5237  4.8885  1.7173  3.7109  2.1782  1.1966  1.2943  1.7888  6.9754  1.1239  2.2144  3.4816  3.2205  3.7692  1.7887  5.1551  4.7342 -0.6787  3.3575  1.9823];
% volt.sosi.RHSC.RPS = [9.3837  3.4628  5.6411  7.9116  3.4154  0.9035  1.5146  2.1896  6.7320  8.3790  7.8571  5.2242  3.7516  6.3837  2.3712  2.9660  2.1784  2.6127  7.1372  2.8614  3.5548  5.4465  4.1379  6.2661  0.9509  8.3370  5.3559 -1.6366  3.6053  3.7904];
% volt.sosi.RHSI.RPS = [8.8992  3.0665  4.3277  5.1487  3.7766 -0.0867  0.5299  1.6148  3.3566  5.2393  6.5508  1.8102  2.2773  3.1184  0.8914  1.9934  2.4873  0.8622  5.1467  0.9051  1.6737  5.7387  3.7846  4.8813  1.1686  6.8099  4.5022 -1.8476  4.5669  0.7449];
% 
% % SOCO - collapsed across 2 and 6 colors
% 
% % 300--500 ms
% 
% volt.soco.RCR.FS = [-4.4263 -4.2673 -6.6186 -4.5971 -3.8143 -4.5957 -4.3548 -3.7414 -0.8475 -6.5189 -6.5741 -3.4129 -4.8162 -5.0403  0.2152 -4.7303 -7.2774  1.2787 -1.5048 -5.3335 -2.2161 -2.7852 -1.1357 -2.0568 -4.6755 -4.0821 -3.5454 -1.8357 -5.2976 -7.4624];
% volt.soco.RHSC.FS = [-2.9110 -2.9149 -5.5456 -3.2689 -0.9778 -3.9752 -2.7713 -2.7332  0.2843 -5.6769 -5.6243 -2.6081 -3.4443 -3.5205  1.7614 -3.0127 -6.4316  2.0104 -0.3915 -3.6063 -2.2112 -3.2948 -1.0587 -0.6655 -3.2918 -3.3427 -2.7125 -2.8635 -3.9659 -6.4667];
% volt.soco.RHSI.FS = [-3.2068 -5.3709 -7.0747 -2.9279 -3.9713 -1.6485 -3.2922 -3.0829 -0.2152 -5.1956 -6.0994 -3.9005 -3.6132 -3.3564  3.0049 -4.6471 -6.6406  2.1379 -0.8339 -2.6294 -1.6958 -2.0932 -2.1832 -0.2426 -3.3735 -2.8939 -2.6085 -1.2968 -5.3664 -5.7779];
% 
% volt.soco.RCR.LAS = [-4.4359 -4.1655 -6.2993 -3.8230 -3.0300 -5.6394 -4.3226 -3.7157 -0.1519 -6.4638 -5.8869 -3.2843 -4.1194 -4.2291 -0.5452 -4.5968 -6.9352  1.1979 -0.1700 -4.6333 -1.8058 -3.6170 -1.4559 -2.6691 -4.6998 -4.0381 -3.6790 -1.6964 -4.4041 -7.4604];
% volt.soco.RHSC.LAS = [-2.1980 -3.0140 -5.2978 -3.2271 -0.8926 -4.1549 -3.0246 -2.9911 -0.0280 -5.3447 -5.2147 -1.5640 -2.7873 -2.6128  0.7637 -3.1485 -5.7939  1.4899  1.0202 -3.1459 -1.9198 -4.0407 -1.1269 -0.8285 -3.0559 -2.4839 -3.5533 -2.3288 -3.8195 -7.0774];
% volt.soco.RHSI.LAS = [-2.2941 -4.5857 -6.6682 -3.0684 -3.8804 -3.0375 -3.7908 -3.0343 -0.0654 -5.4578 -5.2609 -2.3068 -2.4181 -2.3566  2.0692 -3.5805 -6.6544  2.2695  1.4761 -2.6270 -1.5020 -2.7171 -2.8684 -0.1772 -2.9073 -2.4623 -4.1479 -1.2172 -4.6994 -5.7134];
% 
% volt.soco.RCR.RAS = [-3.2683 -3.1634 -6.1817 -4.5750 -3.4853 -4.3062 -4.9868 -3.3102 -2.6011 -5.9281 -6.1761 -3.6294 -4.7115 -4.5934 -0.2647 -2.9737 -5.6323  1.1940 -1.7951 -4.9945 -2.6143 -1.4058 -0.3990 -1.3465 -4.3735 -3.3166 -3.3796 -1.5520 -3.8310 -5.4894];
% volt.soco.RHSC.RAS = [-2.9455 -2.5359 -5.1129 -3.1449 -1.0485 -3.6446 -3.3463 -2.3662 -1.0588 -5.4484 -5.6652 -3.4596 -3.4131 -4.2676  1.7710 -1.4503 -5.4304  2.1303 -1.3179 -3.2143 -2.6935 -1.9219 -0.6594 -0.7153 -3.0623 -3.3107 -2.7316 -2.9642 -3.5525 -4.2572];
% volt.soco.RHSI.RAS = [-2.6473 -4.8319 -6.4007 -2.6204 -3.5545 -2.2085 -3.1485 -2.6378 -1.1201 -4.4697 -6.0656 -4.3144 -4.0242 -3.9972  1.8370 -3.4952 -5.2298  1.7744 -2.3388 -2.0484 -1.9279 -1.1872 -1.6088 -0.0845 -3.4262 -2.7280 -2.2538 -1.5297 -4.3416 -3.9163];
% 
% % 500--800 ms
% 
% volt.soco.RCR.LPS = [2.3292  3.1496  2.2735  0.6134  6.5524 -1.8759  3.0379  2.8124  1.0913 -0.6665  1.8014 -1.9390  4.3506  1.8640  4.4874  3.8394  4.2173  1.9340  3.9473  1.2643  5.2025  3.3855  1.8822 -0.7907  0.6523  2.7305  4.0701  2.5680  7.2843  3.4942];
% volt.soco.RHSC.LPS = [3.5345  2.5549  5.4269  1.1981  5.4972  1.6571  3.6243  3.3002 -0.0261  1.3720  3.5107 -0.6485  5.7275  3.6323  4.5930  4.9151  5.6490  2.4607  4.8145  2.5266  7.4115  2.9644  1.9704 -1.1152  1.1712  2.0772  1.5710  3.1455  7.7597  5.4456];
% volt.soco.RHSI.LPS = [3.5835  0.6268  3.1033 -0.4569  6.6642  0.9385  3.3154  1.9844 -0.5433  0.6589  1.9879 -0.6051  4.9629  2.3140  4.1435  4.7042  3.6278  2.5208  6.3592  2.3497  5.6473  2.9868  0.7567 -1.5841  0.9219  2.6947  1.9878  2.4869  8.2765  5.5157];
% 
% volt.soco.RCR.RPS = [0.5161  2.2974  0.0630  2.3052  5.3540  0.1067  4.2324  2.1232  2.4620  0.7514  1.3584 -0.9257  4.3063  2.4857  5.0936  3.3943  5.6716  1.5421  3.9601  1.2974  4.5005  4.3178  1.3971 -0.3253  2.4071  2.4970  4.2882  0.7146  2.3048  5.7924];
% volt.soco.RHSC.RPS = [0.5476  2.1180  2.6979  2.0726  4.9440  2.0899  4.8951  2.4230  3.0598  1.8572  2.1275 -1.1786  5.2131  2.6355  5.3861  5.4183  5.2102  3.2207  3.3577  1.8260  4.5445  4.1381  0.3213 -1.3986  3.0958  1.9461  3.0467  1.5816  1.7353  6.4500];
% volt.soco.RHSI.RPS = [0.8942  0.7133  0.6086  0.5652  4.6588  0.1798  3.4899  1.0253  1.5261  1.7710  0.3986 -1.4599  4.1189  0.8971  5.6317  3.7248  3.6879  2.4155  4.6601  1.9218  3.1553  3.4092  1.1035 -1.1680  2.1753  2.4036  3.4592  0.7305  3.5277  5.7347];
% 
% exper.badSub.sosi = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
% exper.badSub.soco = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];
% 
% % excluding SOSI007 and the last one from SOCO
% %exper.badSub.sosi = [0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
% %exper.badSub.soco = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 1];

%% These voltages are for EP-only preprocessing

% SOSI

% 300--500 ms

volt.sosi.RCR.FS = [-1.9172 -3.6134  1.1656 -2.6714 -7.0622 -1.1755 -6.2031 -1.7065  0.3363 -5.7751 -3.6106 -4.7340 -3.1307 -0.0866 -3.4406 -4.6947 -3.1345  2.1766 -4.5867 -5.4899 -2.8110  0.5739 -0.2673 -4.9143 -2.1985 -6.6682 -8.5591 -1.3930 -2.6110 -8.1833];
volt.sosi.RHSC.FS = [-0.8210 -3.9413  1.9309 -1.9211 -4.6313 -0.5552 -4.3391 -0.3853  1.0811 -3.6108 -2.7197 -4.9705 -2.4766  1.5901 -2.7753 -3.5386 -1.6044  3.7495 -2.0187 -6.1196 -0.3907  1.3128  1.9814 -4.2286 -2.0634 -4.8216 -7.5430  0.3538 -2.1756 -6.6609];
volt.sosi.RHSI.FS = [-1.3303 -4.6636  3.5354 -2.7475 -6.8748 -0.0654 -5.1770 -2.3431 -2.5520 -4.6431 -0.1733 -6.2553 -3.1005  0.2735 -3.9740 -3.4187 -2.1023  1.8624 -4.9144 -6.0207 -2.3075  0.7339  2.4351 -4.3312 -2.8802 -5.3028 -9.1620 -1.1075 -3.7090 -6.5819];

volt.sosi.RCR.LAS = [-3.1463 -4.8294  1.6097 -0.6662 -4.1419 -0.7456 -5.4915 -0.9406  0.2197 -4.3825 -4.7072 -3.6813 -2.6704  0.8253 -2.6492 -4.0503 -2.2382  0.6886 -4.2025 -3.6315 -3.6350  0.5786 -0.6658 -4.8902 -1.0352 -6.7800 -6.9270 -1.4271 -2.1143 -7.1925];
volt.sosi.RHSC.LAS = [-2.6655 -5.3282  2.1336 -0.3436 -2.9893 -0.3194 -3.9798  0.4425 -0.5791 -2.8256 -3.8347 -3.8586 -2.2305  2.0262 -2.2197 -2.7662 -1.6507  2.0416 -1.0875 -4.5649 -0.6422  0.5849  1.5534 -4.9000 -0.6134 -5.4975 -6.1148  0.3261 -1.5831 -6.1013];
volt.sosi.RHSI.LAS = [-2.9597 -6.2460  3.8133 -0.9311 -4.6546 -0.1542 -4.5073 -1.7020 -4.5547 -4.0719 -1.5209 -5.4462 -3.2213  1.0312 -2.8913 -2.4850 -2.2418  0.7065 -4.0836 -4.8392 -2.8619  0.3222  1.9493 -4.3831 -1.3142 -6.0216 -7.5182 -0.5924 -2.8436 -6.8067];

volt.sosi.RCR.RAS = [-0.7809 -2.3393 -0.0921 -4.6353 -7.0797 -1.7618 -5.2233 -1.9346  0.1590 -6.2723 -2.0139 -4.5546 -2.5343 -1.1607 -3.7629 -4.2927 -1.9780  2.5086 -2.8521 -5.5591 -1.5926  0.3806  0.2042 -3.4856 -2.4147 -5.6158 -8.0076 -2.4421 -3.2616 -9.0359];
volt.sosi.RHSC.RAS = [0.3377 -2.5180  0.2310 -4.1359 -3.9904 -0.7922 -3.3995 -0.8920  1.4175 -4.8029 -1.8326 -4.6273 -2.0758  0.2796 -3.2238 -3.4643 -0.6256  3.9810 -1.6409 -6.2437 -0.3349  1.1296  2.0991 -2.5555 -2.6606 -3.5950 -7.6283 -0.7605 -2.9484 -6.6137];
volt.sosi.RHSI.RAS = [-0.0958 -3.1581  2.7844 -4.4328 -5.6543 -0.2700 -4.1383 -1.8957 -0.0755 -5.0850  0.4226 -5.8351 -2.8625 -0.7956 -4.6455 -3.8094 -0.7139  2.1801 -3.3836 -5.5619 -1.6135  0.6990  2.3016 -2.3992 -2.9074 -4.0964 -8.0409 -2.0038 -3.8746 -7.0142];

% 500--800 ms

volt.sosi.RCR.LPS = [8.0339  1.6708  8.4304  8.7470  6.6991  2.0551  0.6166  2.8274  0.8179  9.4023  4.0488  3.6365  2.7293  3.0192  2.3687  1.8727  0.5140  0.8704  5.2175  2.5338 -0.0281  4.5451  4.0820  3.4593  3.4349  2.3289  4.4006  0.8257  1.7610  2.7585];
volt.sosi.RHSC.LPS = [9.0101  3.3680  8.1597  9.2373  5.9488  1.1885  2.6046  4.6178  3.6128 14.0967  4.0872  3.3684  6.7940  7.2088  3.1416  2.8266  1.6972  1.4475  7.5061  4.4319  1.0907  4.2695  4.7877  6.4367  3.6027  5.6785  6.8569 -0.2538  2.5972  3.8283];
volt.sosi.RHSI.LPS = [9.5862  1.7304  8.9801  8.5303  5.1543  1.1443  1.6817  2.9670  1.6302  8.5536  5.0330 -0.0054  3.7300  3.5351  1.5802  2.5283  1.0624  0.0939  4.8635  2.0944 -0.2918  3.5437  4.9531  1.4768  3.9292  5.0245  4.5797  0.0056  3.1976  2.5198];

volt.sosi.RCR.RPS = [7.8762  2.4607  6.1702  7.9303  2.2398  0.9817  0.0915  1.5484  3.0683  5.8761  7.6780  5.9739  0.8782  3.2918  1.5491  1.7494  0.7541  1.7296  6.2320  1.8531  2.8408  4.1840  3.4381  3.1396  2.2652  5.3833  4.4840 -1.0768  3.4694  2.1439];
volt.sosi.RHSC.RPS = [8.9773  3.4611  5.3617  7.5223  3.6210  0.3892  1.7419  2.4078  5.8791  9.6743  7.9206  4.6559  3.3945  6.2583  2.8010  2.4762  2.0791  2.6199  7.1277  3.1889  3.5042  4.7100  4.3878  6.3996  0.7347  8.2772  5.7980 -1.7670  4.3726  3.8734];
volt.sosi.RHSI.RPS = [10.7434  3.7136  7.4469  6.7333  1.6565 -0.1880  0.1962  2.0835  4.1152  5.4711  9.2001  1.9817  1.8150  3.2624  0.6188  2.3932  2.7812  1.0423  5.4733  1.3645  1.7063  5.0340  3.5125  3.6335  0.9418  7.8968  4.2344 -2.0905  4.8719  0.3723];

% SOCO - collapsed across 2 and 6 colors

% 300--500 ms

volt.soco.RCR.FS = [-3.9205 -3.1249 -6.0301 -4.7829 -3.3941 -5.3200 -4.8946 -3.7566 -1.1472 -6.6141 -6.2055 -3.9671 -4.9439 -5.5673 -0.1029 -3.7480 -6.6291  1.2282 -1.2566 -5.1571 -2.2468 -3.4239 -1.2527 -2.0062 -4.9494 -4.5880 -4.0856 -1.6536 -4.8376 -7.4799];
volt.soco.RHSC.FS = [-3.3205 -3.4120 -6.1208 -3.5268 -1.2989 -3.9480 -3.5266 -2.7005  0.8722 -5.6924 -5.4331 -3.3403 -3.3843 -3.0702  1.9759 -3.0397 -6.5875  1.6634 -0.2522 -3.5166 -1.7420 -2.2895 -0.9272 -1.3739 -3.7332 -3.9090 -3.1647 -3.1196 -3.7959 -6.6669];
volt.soco.RHSI.FS = [-2.8285 -5.0599 -6.0249 -3.3651 -2.7374 -2.5646 -2.6123 -3.2705 -0.7082 -5.4592 -6.0465 -4.0106 -3.0487 -3.6378  2.8581 -4.5646 -7.1191  1.7527  0.5536 -2.7155 -1.7629 -2.5810 -1.4520 -0.6134 -3.2613 -3.4571 -1.9316 -1.9017 -4.7256 -5.0517];

volt.soco.RCR.LAS = [-3.4932 -2.5028 -5.6636 -4.1748 -2.5788 -6.0703 -4.6090 -3.7020 -0.5305 -6.4419 -5.5472 -3.8731 -4.3403 -4.5491 -0.9059 -3.5152 -6.4087  1.1601  0.0664 -4.3241 -1.8569 -4.0746 -2.1387 -2.2385 -4.7074 -4.3602 -3.6576 -1.4023 -3.8580 -7.8712];
volt.soco.RHSC.LAS = [-2.7565 -3.6308 -5.9637 -3.3498 -1.2388 -4.5147 -3.8627 -2.8872  0.5251 -5.4992 -5.0357 -2.2605 -2.8184 -2.1762  1.0264 -3.1328 -6.0156  1.5492  1.0797 -3.0264 -1.5892 -3.3534 -1.1989 -1.4894 -3.3634 -3.4462 -3.4250 -2.6073 -3.5232 -7.2087];
volt.soco.RHSI.LAS = [-1.9033 -4.7366 -5.8680 -3.1973 -2.5325 -3.3691 -4.2032 -3.0596 -0.7879 -5.6840 -5.5624 -2.5038 -2.1400 -2.6527  1.7467 -4.1106 -7.1018  1.8665  2.4174 -2.5501 -1.5922 -3.1052 -2.1297 -0.8837 -2.6686 -2.1771 -3.1450 -2.0352 -4.0539 -5.3037];

volt.soco.RCR.RAS = [-2.9785 -3.0468 -5.7704 -4.5362 -3.2676 -4.7638 -4.9850 -3.2916 -2.9014 -6.2291 -5.9487 -4.2107 -4.7915 -5.1943 -0.2468 -2.5021 -5.0909  1.2033 -1.5196 -4.9821 -2.6463 -1.8886 -0.3001 -1.5136 -4.2296 -3.7829 -4.0574 -1.4579 -3.7634 -5.3834];
volt.soco.RHSC.RAS = [-2.9855 -3.2061 -5.6483 -3.4658 -1.3036 -3.3059 -3.6594 -2.3578 -0.5969 -5.3966 -5.4973 -4.1053 -3.1737 -3.8028  1.9464 -1.6682 -5.6184  1.8356 -1.2644 -3.1800 -2.3861 -1.0301 -0.5682 -1.4015 -3.4621 -3.6398 -3.2717 -2.8758 -3.2425 -3.9889];
volt.soco.RHSI.RAS = [-2.6287 -4.2339 -5.4389 -3.1846 -2.7872 -3.2653 -0.8183 -3.0445 -1.9271 -4.7211 -5.8858 -4.1216 -3.2702 -4.3057  1.5639 -2.9869 -5.6080  1.5375 -1.3503 -2.4112 -1.9287 -1.5895 -0.8568 -0.4276 -3.4844 -3.6513 -1.5690 -1.7975 -4.3514 -3.4067];

% 500--800 ms

volt.soco.RCR.LPS = [2.8924  3.8101  1.7287  0.3071  6.5817 -2.0555  2.9560  3.4058  0.7136 -0.6473  2.0259 -0.6280  5.1139  1.8156  4.6818  4.0246  4.5595  2.0213  3.9732  1.7989  4.9311  3.3360  2.2893 -0.6935  0.9412  2.4642  3.7828  2.9876  7.4062  4.1694];
volt.soco.RHSC.LPS = [3.7843  2.1577  5.3155  1.2655  7.1510  2.3315  3.1379  3.5085 -0.1639  1.3020  3.5647 -0.2225  6.1060  3.3134  4.5410  4.9297  6.0817  2.2880  4.4297  2.5439  7.6099  3.2244  2.5253 -0.4668  2.2161  1.4638  1.7016  2.6983  7.8951  6.0626];
volt.soco.RHSI.LPS = [3.4285  1.0817  4.4072  0.2154  5.4047  0.9614  4.1145  2.4727 -1.3779  0.5809  2.7547 -2.0327  4.7529  2.1999  3.7410  4.9065  4.8680  2.6822  6.7414  3.1787  5.4317  2.9002  1.6427 -1.4947  0.2848  2.2848  2.6063  2.6434  7.3232  5.6512];

volt.soco.RCR.RPS = [0.6786  2.1016 -0.8064  2.5133  5.2293  0.0164  4.1365  2.7280  2.4664  0.4018  1.3715  0.0626  4.6633  1.8379  5.5632  4.0205  5.4317  1.4834  3.8311  1.8693  4.7029  4.0216  1.9371 -0.6026  3.0211  2.1661  3.9080  1.7482  2.2918  6.2597];
volt.soco.RHSC.RPS = [1.4764  2.0859  2.7921  2.1787  5.5486  2.6662  4.1562  2.6387  2.3127  2.0167  2.4047 -0.9494  5.8477  2.3772  5.3358  5.3154  6.0888  2.8274  3.0142  1.8944  4.8745  4.5741  0.7194 -1.2650  4.0831  1.8921  3.1111  1.3888  2.5963  7.6789];
volt.soco.RHSI.RPS = [1.2703  1.0392  1.7856  0.0202  3.7864  0.0714  7.0961  1.3875  0.5329  2.0841  1.4192 -2.3599  4.3566  1.0700  5.6117  4.7704  5.3504  2.9853  5.2909  2.6676  2.9246  2.9459  1.3643 -1.4242  1.5422  0.6922  3.8726  1.6917  2.4800  6.0162];

%exper.badSub.sosi = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
%exper.badSub.soco = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];




%% Set up the rest of the ANOVA

% IMPORTANT: This requires an equal number of subjects in each experiment

% IV1 (between-subjects factor): Experiment
exper.name = {'sosi','soco'};
% IV2 (within-subjects RM factor): ROI
exper.roi = {'LAS','RAS'};
% IV3 (within-subjects RM factor): Accuracy Condition
exper.cond = {'RCR','RHSC','RHSI'};
%exper.cond = {'RHSC','RHSI'};

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

% I need to use R or SPSS to run:
%
% At a single ROI, do a 2-way ANOVA with one between-subjects factor
% (experiment) and one within-subjects factor (Accuracy Condition)



%% write it out so we can run unequal numbers of subjects

%exper.badSub.sosi = [1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
%exper.badSub.soco = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];

%nSub = sum(~exper.badSub.soco) + sum(~exper.badSub.sosi);

% IV1 (between-subjects factor): Experiment
exper.name = {'sosi','soco'};
% IV2 (within-subjects RM factor): ROI
exper.roi = {'LAS','RAS'};
% IV3 (within-subjects RM factor): Accuracy Condition
exper.cond = {'RCR','RHSC','RHSI'};

anova_write = [];

for e = 1:length(exper.name)
  goodSubInd = 0;
  for s = 1:length(exper.(exper.name{e}).subjects)
    if ~ismember(exper.(exper.name{e}).subjects(s),exper.(exper.name{e}).badSub)
      goodSubInd = goodSubInd + 1;
      for r = 1:length(exper.roi)
        for c = 1:length(exper.cond)
          %anova_write = [anova_write; goodSubInd sprintf('exp%d',e) sprintf('roi%d',r) sprintf('cnd%d',c) volt.(exper.name{e}).(exper.cond{c}).(exper.roi{r})(s)];
          anova_write = [anova_write; goodSubInd e r c volt.(exper.name{e}).(exper.cond{c}).(exper.roi{r})(s)];
        end
      end
    end
  end
end

outfile = fullfile(getenv('HOME'),'Desktop','soco_sosi_volt.txt');
fid = fopen(outfile,'w+');

% write the header
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n','sub','exp','roi','cond','volt');
% write the data
for i = 1:size(anova_write,1)
  fprintf(fid,'%d\t%d\t%d\t%d\t%.4f\n',anova_write(i,1),anova_write(i,2),anova_write(i,3),anova_write(i,4),anova_write(i,5));
end
% close it
fclose(fid);


%% write it out for SPSS

% exper.badSub.sosi = [1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
% exper.badSub.soco = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];
% % exper.badSub.soco = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];

%nSub = sum(~exper.badSub.soco) + sum(~exper.badSub.sosi);

% IV1 (between-subjects factor): Experiment
exper.name = {'sosi','soco'};
% IV2 (within-subjects RM factor): ROI
exper.roi = {'LAS','RAS'};
% IV3 (within-subjects RM factor): Accuracy Condition
exper.cond = {'RCR','RHSC','RHSI'};

anova_write = [];

for e = 1:length(exper.name)
  goodSubInd = 0;
  for s = 1:length(exper.badSub.(exper.name{e}))
    if ~ismember(exper.(exper.name{e}).subjects(s),exper.(exper.name{e}).badSub)
      goodSubInd = goodSubInd + 1;
      %anova_write = [anova_write; goodSubInd sprintf('exp%d',e) sprintf('roi%d',r) sprintf('cnd%d',c) volt.(exper.name{e}).(exper.cond{c}).(exper.roi{r})(s)];
      anova_write = [anova_write; goodSubInd, e,...
        volt.(exper.name{e}).(exper.cond{1}).(exper.roi{1})(s),...
        volt.(exper.name{e}).(exper.cond{2}).(exper.roi{1})(s),...
        volt.(exper.name{e}).(exper.cond{3}).(exper.roi{1})(s),...
        volt.(exper.name{e}).(exper.cond{1}).(exper.roi{2})(s),...
        volt.(exper.name{e}).(exper.cond{2}).(exper.roi{2})(s),...
        volt.(exper.name{e}).(exper.cond{3}).(exper.roi{2})(s)];
    end
  end
end



outfile = fullfile(getenv('HOME'),'Desktop','soco_sosi_volt_spss.txt');
fid = fopen(outfile,'w+');
% write the header
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','sub','exp','LAS_CR','LAS_HSC','LAS_HSI','RAS_CR','RAS_HSC','RAS_HSI');
% write the data
for i = 1:size(anova_write,1)
  fprintf(fid,'s%d\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',...
    anova_write(i,1),...
    exper.name{anova_write(i,2)},...
    anova_write(i,3),anova_write(i,4),anova_write(i,5),anova_write(i,6),anova_write(i,7),anova_write(i,8));
end
% close it
fclose(fid);



%% independent-samples t-test

% test the experiment X accuracy condition voltage difference interaction at one ROI

roi = 'LAS';
[h,p,ci,stats] = ttest2((volt.sosi.RHSC.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.RHSI.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),(volt.soco.RHSC.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.RHSI.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (RHSC-RHSI) (M=%.2f) vs SOCO (RHSC-RHSI) (M=%.2f): t(%d)=%.4f, p=%.10f\n',roi,mean(volt.sosi.RHSC.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.RHSI.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(volt.soco.RHSC.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.RHSI.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

roi = 'RAS';
[h,p,ci,stats] = ttest2((volt.sosi.RHSC.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.RHSI.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),(volt.soco.RHSC.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.RHSI.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('%s: SOSI (RHSC-RHSI) (M=%.2f) vs SOCO (RHSC-RHSI) (M=%.2f): t(%d)=%.4f, p=%.10f\n',roi,mean(volt.sosi.RHSC.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub)) - volt.sosi.RHSI.(roi)(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(volt.soco.RHSC.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub)) - volt.soco.RHSI.(roi)(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);


volt_sosi_RHSC = mean([volt.sosi.RHSC.('LAS')(~ismember(exper.sosi.subjects,exper.sosi.badSub));volt.sosi.RHSC.('RAS')(~ismember(exper.sosi.subjects,exper.sosi.badSub))],1);
volt_sosi_RHSI = mean([volt.sosi.RHSI.('LAS')(~ismember(exper.sosi.subjects,exper.sosi.badSub));volt.sosi.RHSI.('RAS')(~ismember(exper.sosi.subjects,exper.sosi.badSub))],1);
volt_sosi_RCR = mean([volt.sosi.RCR.('LAS')(~ismember(exper.sosi.subjects,exper.sosi.badSub));volt.sosi.RCR.('RAS')(~ismember(exper.sosi.subjects,exper.sosi.badSub))],1);
volt_soco_RHSC = mean([volt.soco.RHSC.('LAS')(~ismember(exper.soco.subjects,exper.soco.badSub));volt.soco.RHSC.('RAS')(~ismember(exper.soco.subjects,exper.soco.badSub))],1);
volt_soco_RHSI = mean([volt.soco.RHSI.('LAS')(~ismember(exper.soco.subjects,exper.soco.badSub));volt.soco.RHSI.('RAS')(~ismember(exper.soco.subjects,exper.soco.badSub))],1);
volt_soco_RCR = mean([volt.soco.RCR.('LAS')(~ismember(exper.soco.subjects,exper.soco.badSub));volt.soco.RCR.('RAS')(~ismember(exper.soco.subjects,exper.soco.badSub))],1);

[h,p,ci,stats] = ttest2((volt_sosi_RHSC - volt_sosi_RHSI),(volt_soco_RHSC - volt_soco_RHSI),0.05,'both');
fprintf('LAS+RAS: SOSI (RHSC-RHSI) (M=%.2f) vs SOCO (RHSC-RHSI) (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(volt_sosi_RHSC - volt_sosi_RHSI),mean(volt_soco_RHSC - volt_soco_RHSI),stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest2((volt_sosi_RHSC - volt_sosi_RCR),(volt_soco_RHSC - volt_soco_RCR),0.05,'both');
fprintf('LAS+RAS: SOSI (RHSC-RCR) (M=%.2f) vs SOCO (RHSC-RCR) (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(volt_sosi_RHSC - volt_sosi_RCR),mean(volt_soco_RHSC - volt_soco_RCR),stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest2((volt_sosi_RHSI - volt_sosi_RCR),(volt_soco_RHSI - volt_soco_RCR),0.05,'both');
fprintf('LAS+RAS: SOSI (RHSI-RCR) (M=%.2f) vs SOCO (RHSI-RCR) (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(volt_sosi_RHSI - volt_sosi_RCR),mean(volt_soco_RHSI - volt_soco_RCR),stats.df,stats.tstat,p);

%% dependent-samples t-test (Within an experiment)

alpha = 0.05;

rois = {'LAS','RAS'};
expNames = {'soco','sosi'};
condPairs = {{'RHSC','RHSI'},{'RHSC','RCR'},{'RHSI','RCR'}};

for e = 1:length(expNames)
  expName = expNames{e};
  for r = 1:length(rois)
    roi = rois{r};
    for c = 1:length(condPairs)
      conds = condPairs{c};
      [h,p,ci,stats] = ttest(...
        volt.(expName).(conds{1}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub)),...
        volt.(expName).(conds{2}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub)),...
        alpha,'both');
      fprintf('ttest: %s: %s: %s (M=%.2f) vs %s (M=%.2f): t(%d)=%.4f, p=%.10f',roi,expName,...
        conds{1},mean(volt.(expName).(conds{1}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub))),...
        conds{2},mean(volt.(expName).(conds{2}).(roi)(~ismember(exper.(expName).subjects,exper.(expName).badSub))),...
        stats.df,stats.tstat,p);
      if p < alpha
        fprintf(' *');
      end
      fprintf('\n');
    end
  end
end

%% behavioral accuracy test

%exper.badSub.sosi = [1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
%exper.badSub.soco = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0];

% SOCO separate colors
SOCO_C2_RS_WIR = [0.6304 0.7093 0.9323 0.9412 0.8333 0.4957 0.8769 0.8299 0.9322 0.7662 0.7875 0.7143 0.7674 0.5345 0.9167 0.8661 0.9975 0.6 0.9263 0.6263 0.6739 0.5455 0.9737 0.8182 0.7164 0.9975 0.7812 0.7333 0.5484 0.7857];
SOCO_C2_RO_WIR = [0.3269 0.566 0.6176 0.9 0.4386 0.6923 0.5269 0.5455 0.4466 0.5455 0.5 0.5312 0.475 0.45 0.4615 0.5694 0.5849 0.5192 0.4615 0.5696 0.5882 0.5412 0.5957 0.5385 0.5122 0.6 0.6019 0.625 0.4565 0.4714];
SOCO_C2_F_WIR = [0.5 0.7 0.5 0.5 0.7273 0.9975 0.5135 0.7143 0.5 0.0025 0.7 0.5806 0.561 0.6471 0.5 0.6667 0.6353 0.5 0.6 0.5 0.5217 0.6 0.5312 0.3571 0.5417 0.4091 0.25 0.5773 0.4894 0.5217];

SOCO_C6_RS_WIR = [0.6119 0.8 0.9048 0.8333 0.6714 0.5091 0.9059 0.8168 0.9143 0.7415 0.9519 0.7576 0.7286 0.8904 0.8148 0.9391 0.963 0.5882 0.9024 0.6306 0.6496 0.3793 0.9556 0.88 0.8056 0.6 0.76 0.7742 0.5974 0.6364];
SOCO_C6_RO_WIR = [0.4795 0.4082 0.5345 0.7561 0.5 0.4773 0.6479 0.6 0.598 0.5238 0.62 0.5283 0.5263 0.6216 0.6143 0.5745 0.5417 0.551 0.4528 0.5588 0.6207 0.5263 0.6386 0.5738 0.4643 0.5676 0.4321 0.6 0.4872 0.5];
SOCO_C6_F_WIR = [0.375 0.5833 0.3571 0.5769 0.6364 0.4286 0.6111 0.75 0.4615 0.3333 0.6562 0.5147 0.4091 0.5 0.3684 0.6957 0.5068 0.625 0.2105 0.0025 0.4545 0.4545 0.4828 0.4583 0.8 0.551 0.7143 0.5242 0.5581 0.525];

% SOCO collapsed across colors
SOCO_RS_WIR = mean([SOCO_C2_RS_WIR;SOCO_C6_RS_WIR],1);
SOCO_RO_WIR = mean([SOCO_C2_RO_WIR;SOCO_C6_RO_WIR],1);
SOCO_F_WIR = mean([SOCO_C2_F_WIR;SOCO_C6_F_WIR],1);

% SOSI
SOSI_RS_WIR = [0.9441 0.8333 0.9444 0.9302 0.7884 0.8621 0.5323 0.9349 0.9891 0.9878 0.9121 0.8398 0.9121 0.9789 0.8958 0.809 0.764 0.971 0.9298 0.817 0.8846 0.7308 0.9412 0.641 0.988 0.9679 0.9942 0.9767 0.8951 0.8339];
SOSI_RO_WIR = [0.7333 0.8636 0.9987 0.6667 0.7143 0.5354 0.4393 0.45 0.7391 0.8088 0.575 0.7273 0.6739 0.6792 0.625 0.6441 0.8681 0.6087 0.625 0.6053 0.6304 0.4898 0.6533 0.726 0.6825 0.5942 0.899 0.8395 0.625 0.5417];
SOSI_F_WIR = [0.7727 0.5741 0.5778 0.6094 0.6389 0.58 0.4286 0.4444 0.6842 0.5 0.0013 0.5135 0.5303 0.541 0.5 0.5067 0.6279 0.7037 0.5118 0.4667 0.5795 0.5882 0.4074 0.7368 0.5238 0.4468 0.5977 0.6 0.625 0.0013];

% between experiments
[h,p,ci,stats] = ttest2(SOSI_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub)),SOCO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub)),0.05,'both');
fprintf('ttest2: F WIR: SOSI (M=%.2f) vs SOCO (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(SOSI_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(SOCO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

% within experiment
[h,p,ci,stats] = ttest(SOCO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub)),0.5*ones(1,sum(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('ttest: SOCO F_WIR (M=%.2f) vs chance: t(%d)=%.4f, p=%.10f\n',mean(SOCO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(SOSI_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub)),0.5*ones(1,sum(~ismember(exper.sosi.subjects,exper.sosi.badSub))),0.05,'both');
fprintf('ttest: SOSI F_WIR (M=%.2f) vs chance: t(%d)=%.4f, p=%.10f\n',mean(SOSI_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub))),stats.df,stats.tstat,p);


%% collapse across RO and F
SOCO_RO_F_WIR = mean([SOCO_RO_WIR;SOCO_F_WIR],1);
SOSI_RO_F_WIR = mean([SOSI_RO_WIR;SOSI_F_WIR],1);

% within experiment
[h,p,ci,stats] = ttest(SOCO_RO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub)),0.5*ones(1,sum(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('SOCO RO+F WIR (M=%.2f) vs chance: t(%d)=%.4f, p=%.10f\n',mean(SOCO_RO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(SOSI_RO_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub)),0.5*ones(1,sum(~ismember(exper.sosi.subjects,exper.sosi.badSub))),0.05,'both');
fprintf('SOSI RO+F WIR (M=%.2f) vs chance: t(%d)=%.4f, p=%.10f\n',mean(SOSI_RO_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub))),stats.df,stats.tstat,p);

% between experiments
[h,p,ci,stats] = ttest2(SOSI_RO_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub)),SOCO_RO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub)),0.05,'both');
fprintf('RO+F WIR: SOSI (M=%.2f) vs SOCO (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(SOSI_RO_F_WIR(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(SOCO_RO_F_WIR(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

%% more data

% d'
SOCO_ITEM_DP = [0.9877, 1.8521, 1.6612, 2.3903, 1.5652, 1.2715, 3.0669, 2.0165, 1.5582, 2.7002, 2.434, 0.5012, 1.0083, 1.4696, 2.5912, 1.9567, 1.1204, 1.4478, 1.7541, 0.9364, 1.4128, 0.4553, 1.7471, 1.7161, 1.3115, 0.9281, 2.5018, 0.6834, 1.2993, 2.1327];
SOCO_SOURCE_DP = [0.0979, 0.7479, 1.6706, 1.8041, 0.6517, 0.1094, 1.0921, 1.5188, 0.8406, 1.2167, 1.3084, 0.4363, 0.6433, 0.6117, 0.4553, 1.5573, 0.6422, 0.2452, 1.1511, 0.5393, 0.675, 0.042, 0.9744, 0.3648, 0.7799, 0.3452, 0.6405, 0.5407, 0.1554, 0.4839];

% 19 and 20 out of order (use this)
SOSI_ITEM_DP = [2.4055, 1.0515, 1.9568, 1.3392, 0.875, 1.2222, 0.185, 1.9591, 1.9172, 1.1959, 2.4869, 0.9108, 1.0227, 1.705, 1.5148, 0.9685, 0.8106, 2.0084, 1.3761, 1.1812, 1.1057, 1.1521, 2.297, 0.4484, 1.4344, 1.5941, 1.398, 1.7473, 1.5984, 1.2503];
SOSI_SOURCE_DP = [2.6106, 1.2027, 1.515, 1.6033, 1.3849, 1.0842, -0.0091, 1.8163, 2.6303, 1.1873, 2.3238, 1.3773, 1.6582, 1.5506, 1.6074, 0.8615, 1.5816, 2.6045, 1.1, 1.1413, 1.4382, 0.8814, 2.1249, 1.0363, 1.6679, 1.5655, 2.2779, 2.1689, 1.6759, 1.6488];

% % 19 and 20 in correct order
% SOSI_ITEM_DP = [2.4055, 1.0515, 1.9568, 1.3392, 0.875, 1.2222, 0.185, 1.9591, 1.9172, 1.1959, 2.4869, 0.9108, 1.0227, 1.705, 1.5148, 0.9685, 0.8106, 2.0084, 1.1812, 1.3761, 1.1057, 1.1521, 2.297, 0.4484, 1.4344, 1.5941, 1.398, 1.7473, 1.5984, 1.2503];
% SOSI_SOURCE_DP = [2.6106, 1.2027, 1.515, 1.6033, 1.3849, 1.0842, -0.0091, 1.8163, 2.6303, 1.1873, 2.3238, 1.3773, 1.6582, 1.5506, 1.6074, 0.8615, 1.5816, 2.6045, 1.1413, 1.1, 1.4382, 0.8814, 2.1249, 1.0363, 1.6679, 1.5655, 2.2779, 2.1689, 1.6759, 1.6488];

% response bias
SOCO_ITEM_C = [-0.5915, -0.2024, -1.0839, -0.4808, -0.0479, -0.1114, 0.1114, -0.3802, -0.5772, 0.3453, -0.2225, -0.609, -0.3023, -0.1068, -0.2592, -0.6199, -0.1064, 1.1569, 0.1382, -0.8133, -0.0161, 0.1712, 0.0229, -0.1677, -0.4024, 0.7109, 0.4445, -0.6735, -0.1396, -0.0719];
SOCO_SOURCE_C = [-0.3196, -0.245, -0.0832, -0.0179, -0.0403, -0.0626, -0.231, -0.2677, -0.1783, 0.128, -0.1683, -0.1723, -0.1488, -0.1191, -0.2937, 0.1208, 0.0541, 0.0466, -0.1624, 0.3533, -0.0809, -0.2229, -0.1294, 0.0087, -0.3034, -0.426, 0.2808, -0.0599, 0.0777, 0.1782];

SOSI_ITEM_C = [-0.6091, 0.2464, 0.9816, 0.2458, -0.4186, -0.0282, -0.7053, -0.1235, 0.3229, 0.1207, 0.1288, -0.1499, -0.1395, -0.4855, -0.6821, 0.3396, 0.2536, -0.1446, 0.2362, -0.4142, -0.2605, -0.1084, -0.8561, 0.6688, -0.5079, 0.3293, -0.5409, -0.0176, 0.4273, -0.7256];
SOSI_SOURCE_C = [0.1518, -0.1634, -0.3268, -0.2707, -0.1447, -0.1485, -0.0598, 0.1957, -0.3433, 0.0536, -0.218, 0.1723, -0.2209, -0.2229, -0.0029, -0.1235, 0.1766, 0.1119, 0.2809, -0.1177, -0.118, -0.2956, -0.228, 0.1352, -0.0861, -0.0071, 0.0919, 0.0322, 0.1459, 0.045];


%% everyone

fprintf('\n');

% item d'

% between experiments
[h,p,ci,stats] = ttest2(SOSI_ITEM_DP(~ismember(exper.sosi.subjects,exper.sosi.badSub)),SOCO_ITEM_DP(~ismember(exper.soco.subjects,exper.soco.badSub)),0.05,'both');
fprintf('ttest2: Item d'': SOSI (M=%.2f) vs SOCO (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(SOSI_ITEM_DP(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(SOCO_ITEM_DP(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

% source d'

% between experiments
[h,p,ci,stats] = ttest2(SOSI_SOURCE_DP(~ismember(exper.sosi.subjects,exper.sosi.badSub)),SOCO_SOURCE_DP(~ismember(exper.soco.subjects,exper.soco.badSub)),0.05,'both');
fprintf('ttest2: Source d'': SOSI (M=%.2f) vs SOCO (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(SOSI_SOURCE_DP(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(SOCO_SOURCE_DP(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

fprintf('\n');

% Item response bias (c)

% between experiments
[h,p,ci,stats] = ttest2(SOSI_ITEM_C(~ismember(exper.sosi.subjects,exper.sosi.badSub)),SOCO_ITEM_C(~ismember(exper.soco.subjects,exper.soco.badSub)),0.05,'both');
fprintf('ttest2: Item response bias (c): SOSI (M=%.2f) vs SOCO (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(SOSI_ITEM_C(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(SOCO_ITEM_C(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

% within experiment
[h,p,ci,stats] = ttest(SOCO_ITEM_C(~ismember(exper.soco.subjects,exper.soco.badSub)),zeros(1,sum(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('ttest: SOCO Item response bias (c) (M=%.2f) vs zero: t(%d)=%.4f, p=%.10f\n',mean(SOCO_ITEM_C(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(SOSI_ITEM_C(~ismember(exper.sosi.subjects,exper.sosi.badSub)),zeros(1,sum(~ismember(exper.sosi.subjects,exper.sosi.badSub))),0.05,'both');
fprintf('ttest: SOSI Item response bias (c) (M=%.2f) vs zero: t(%d)=%.4f, p=%.10f\n',mean(SOSI_ITEM_C(~ismember(exper.sosi.subjects,exper.sosi.badSub))),stats.df,stats.tstat,p);

% Source response bias (c)

% between experiments
[h,p,ci,stats] = ttest2(SOSI_SOURCE_C(~ismember(exper.sosi.subjects,exper.sosi.badSub)),SOCO_SOURCE_C(~ismember(exper.soco.subjects,exper.soco.badSub)),0.05,'both');
fprintf('ttest2: Source response bias (c): SOSI (M=%.2f) vs SOCO (M=%.2f): t(%d)=%.4f, p=%.10f\n',mean(SOSI_SOURCE_C(~ismember(exper.sosi.subjects,exper.sosi.badSub))),mean(SOCO_SOURCE_C(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

% within experiment
[h,p,ci,stats] = ttest(SOCO_SOURCE_C(~ismember(exper.soco.subjects,exper.soco.badSub)),zeros(1,sum(~ismember(exper.soco.subjects,exper.soco.badSub))),0.05,'both');
fprintf('ttest: SOCO Source response bias (c) (M=%.2f) vs zero: t(%d)=%.4f, p=%.10f\n',mean(SOCO_SOURCE_C(~ismember(exper.soco.subjects,exper.soco.badSub))),stats.df,stats.tstat,p);

[h,p,ci,stats] = ttest(SOSI_SOURCE_C(~ismember(exper.sosi.subjects,exper.sosi.badSub)),zeros(1,sum(~ismember(exper.sosi.subjects,exper.sosi.badSub))),0.05,'both');
fprintf('ttest: SOSI Source response bias (c) (M=%.2f) vs zero: t(%d)=%.4f, p=%.10f\n',mean(SOSI_SOURCE_C(~ismember(exper.sosi.subjects,exper.sosi.badSub))),stats.df,stats.tstat,p);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mix and match subjects from each experiment to equate source d' but
% differ on F accuracy


%% initialize randomness

rng('shuffle','twister');

%% try to equate

nTries = 10000;
p_thresh = 0.1;

nSub = 20;

foundIt = false;

while nSub > 5
  fprintf('nSub = %d\n',nSub);
  
  for i = 1:nTries
    fprintf('.');
    
    SOCO_subNums = randperm(length(exper.soco.subjects));
    SOCO_subNums = SOCO_subNums(~ismember(exper.soco.subjects(SOCO_subNums),exper.soco.badSub));
    SOCO_subNums = SOCO_subNums(1:nSub);
    SOSI_subNums = randperm(length(exper.sosi.subjects));
    SOSI_subNums = SOSI_subNums(~ismember(exper.sosi.subjects(SOSI_subNums),exper.sosi.badSub));
    SOSI_subNums = SOSI_subNums(1:nSub);
    
    [h,p,ci,stats] = ttest2(SOCO_SOURCE_DP(SOCO_subNums),SOSI_SOURCE_DP(SOSI_subNums),0.05,'both');
    if p > p_thresh
      [h_cs,p_cs,ci_cs,stats_cs] = ttest2(SOCO_F_WIR(SOCO_subNums),SOSI_F_WIR(SOSI_subNums),0.05,'both');
      if p_cs < .05
        fprintf('Found p > %.2f and p_cs < .05 after %d tries (nSub=%d).\n',p_thresh,i,nSub);
        
        fprintf('ttest2: SOCO source d''=%.2f vs SOSI source d''=%.2f: t(%d)=%.2f, p=%.3f',mean(SOCO_SOURCE_DP(SOCO_subNums)),mean(SOSI_SOURCE_DP(SOSI_subNums)),stats.df,stats.tstat,p);
        
        fprintf('\nSOCO:');
        sort(exper.soco.subjects(SOCO_subNums))
        fprintf('\nSOSI:');
        sort(exper.sosi.subjects(SOSI_subNums))
        
        [h,p,ci,stats] = ttest2(SOCO_F_WIR(SOCO_subNums),SOSI_F_WIR(SOSI_subNums),0.05,'both');
        fprintf('ttest2: SOCO F=%.2f vs SOSI F=%.2f: t(%d)=%.2f, p=%.3f\n',mean(SOCO_F_WIR(SOCO_subNums)),mean(SOSI_F_WIR(SOSI_subNums)),stats.df,stats.tstat,p);
        
        [h,p,ci,stats] = ttest(SOCO_F_WIR(SOCO_subNums),repmat(0.5,size(SOCO_subNums)),0.05,'both');
        fprintf('ttest: SOCO F=%.2f vs chance (%.1f): t(%d)=%.2f, p=%.3f\n',mean(SOCO_F_WIR(SOCO_subNums)),0.5,stats.df,stats.tstat,p);
        
        [h,p,ci,stats] = ttest(SOSI_F_WIR(SOSI_subNums),repmat(0.5,size(SOSI_subNums)),0.05,'both');
        fprintf('ttest: SOSI F=%.2f vs chance (%.1f): t(%d)=%.2f, p=%.3f\n',mean(SOSI_F_WIR(SOSI_subNums)),0.5,stats.df,stats.tstat,p);
        
        foundIt = true;
        break
      end
    end
    if i == nTries
      fprintf('Hit the %d-try limit.\n',i);
    end
  end
  if foundIt == true
    break
  else
    nSub = nSub - 1;
  end
end

%% do an nchoosek count of F comparisons

% choose based on source d' top (SOCO) or bottom (SOSI) nSubsubset subjects
nSubSubset = 15;

dp_alpha = 0.1;
f_alpha = 0.05;

% excluded based on trial counts, etc.
exper.sosi.badBehSub = {'SOSI001','SOSI007','SOSI011','SOSI030'};
exper.soco.badBehSub = {'SOCO018','SOCO026'};

% worst to best
[~,SOCO_sort] = sort(SOCO_SOURCE_DP);
[~,SOSI_sort] = sort(SOSI_SOURCE_DP);

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
SOCO_SOURCE_DP_SUBSET = SOCO_SOURCE_DP(ismember(exper.soco.subjects,exper.soco.subSubset));
SOSI_SOURCE_DP_SUBSET = SOSI_SOURCE_DP(ismember(exper.sosi.subjects,exper.sosi.subSubset));
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
    [h,p,ci,stats] = ttest2(SOCO_SOURCE_DP_SUBSET(allCombos(i,:)),SOSI_SOURCE_DP_SUBSET(allCombos(j,:)),dp_alpha,'both');
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
  fprintf('\tColor vs Side source d'': %d of %d t-tests had p>%.2f (%.2f%%)\n',sum(FtestedInd(:)),results.nCombos(i)^2,dp_alpha,((sum(FtestedInd(:))/results.nCombos(i)^2)*100));
  fprintf('\tColor vs Side F accuracy: %d of %d the null d'' tests had p<%.2f (%.2f%%)\n',sum(FlessInd(:)),sum(FtestedInd(:)),f_alpha,(sum(FlessInd(:))/sum(FtestedInd(:))*100));
end
