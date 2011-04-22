
% reaction times by source/correct
S_RH_SrcCor_rt = [592.2281 252.2021 433.837 305.386 392.5204 376.2593 288.1765 950.2547 289.7979 370.033 339.1471 319.5385 414.4032 406.7476 594.6076 420.2333 322 536.89 524.6029 439.3878 519.8571 376.2692 346.8911 481.9706];
S_RH_SrcInc_rt = [509.48 338.2973 415.25 347.75 397.96 475.55 300.0256 1174.8519 335.4444 440.3333 515.3421 603.1892 531.913 384.9 572.1379 389.5833 336.7273 542.8214 649 545.25 436.65 445.72 802 536.3571];
L_RH_SrcCor_rt = [648.3393 314.3647 443.05 324.6442 388.4255 439.1176 285.2778 1047.6897 298.5538 432.5857 432.25 405.6508 478.2 428.9623 457.3333 348.8718 338.0897 567.4528 622.7895 468.9552 442.3284 350.9012 341.25 499.6111];
L_RH_SrcInc_rt = [551.1628 369.5294 477.963 319.8 429.2069 570.6176 321.7857 1821.5789 318.0488 457.9091 414.1702 452.7778 368.0952 486.614 580.4211 454.4 390.2273 691.4531 539.4348 474.1562 460 442.5405 317.5135 512.3929];

nSub = length(S_RH_SrcCor_rt);

rt_src = {S_RH_SrcCor_rt S_RH_SrcInc_rt;...
  L_RH_SrcCor_rt L_RH_SrcInc_rt};

cfg.task = {'Size','Life'};
cfg.accur = {'Correct','Incorrect'};
cfg.alpha = .05;
cfg.showtable = 1;
cfg.printTable_tex = 1;
if length(cfg.task) > 2 || length(cfg.accur) > 2
  cfg.calcGGHF = 1;
else
  cfg.calcGGHF = 0;
end

fprintf('====================== RMAOV2_mod =========================\n');
fprintf('IV1: Task (%d;%s), IV2: Source Accuracy RT (%d;%s)\n',length(cfg.task),sprintf(repmat(' %s',1,length(cfg.task)),cfg.task{:}),length(cfg.accur),sprintf(repmat(' %s',1,length(cfg.accur)),cfg.accur{:}));

% col 1 = dependent variable (rate)
% col 2 = independent variable 1 (study task: size or life)
% col 3 = independent variable 2 (accuracy: correct or incorrect)
% col 4 = subject number

anovamat = [];
for t = 1:length(cfg.task)
  for a = 1:length(cfg.accur)
    for sub = 1:nSub
      anovamat = [anovamat; rt_src{t,a}(sub) t a sub];
    end
  end
end
p_vals = RMAOV2_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

%%%%%%%%%%%%%%%%%%%%

% reaction times overall
S_RH_rt = mean(cat(1,S_RH_SrcCor_rt,S_RH_SrcInc_rt),1);
L_RH_rt = mean(cat(1,L_RH_SrcCor_rt,L_RH_SrcInc_rt),1);
CR_rt = [530.2529 325.119 417.5547 360.6131 409.6667 358.3971 400.8058 1011.5932 289.8333 381.2555 335.5556 417.0741 512.6972 516.7246 579.9606 382.6978 316.7832 737.2333 469.3071 557.6466 434.1471 357.1318 395.1544 514.4603];

nSub = length(S_RH_rt);

rt_overall = {S_RH_rt,L_RH_rt,CR_rt};

cfg.task = {'Size','Life','CR'};
cfg.alpha = .05;
cfg.showtable = 1;
cfg.printTable_tex = 1;
if length(cfg.task) > 2
  cfg.calcGGHF = 1;
else
  cfg.calcGGHF = 0;
end

fprintf('====================== RMAOV1_mod =========================\n');
fprintf('IV1: Task (%d;%s)\n',length(cfg.task),sprintf(repmat(' %s',1,length(cfg.task)),cfg.task{:}));

% col 1 = dependent variable (rate)
% col 2 = independent variable 1 (study task: size or life)
% col 3 = subject number

anovamat = [];
for t = 1:length(cfg.task)
  for sub = 1:nSub
    anovamat = [anovamat; rt_overall{t}(sub) t sub];
  end
end
p_vals = RMAOV1_mod(anovamat,cfg.alpha,cfg.showtable,cfg.calcGGHF,cfg.printTable_tex);

