function [P,X] = RMAOV2_mod(X,alpha,showtable,calcGGHF,printTable_tex,colsubs)
% RMAOV2 Repeated Measures Two-way Analysis of Variance Test.
%   ANOVA with two within-subject variables is used to analyze the relationship 
%   between two independent variables and the dependent variable. This procedure
%   requires fewer participants but can give rise to research confounds such as 
%   order or practice effects and participant bias.
%   For a more detailed explanation on repeated measures designs, we suggest you 
%   read the help section of the RMAOV1 file which you can find on this Matlab File
%   Exchange site.
%   
%   Syntax: [P,X] = RMAOV2_mod(X,alpha,showtable,calcGGHF,printTable_tex,colsubs)
%      
%     Inputs:
%          X - data matrix (Size of matrix must be n-by-4;dependent variable=column 1;
%              independent variable 1=column 2;independent variable 2=column 3;
%              subject=column 4). 
%      alpha - significance level (default = 0.05).
%  showtable - 0 or 1 to show the table of results.
%    Outputs:
%            - Complete Analysis of Variance Table.
%            - Strength of the relationships.
%
%    Example: From the example of Schuyler W. Huck other's on-line resources Chapter 16 (within-subjects ANOVA)
%             of the Readings Statistics and Research book [http://www.readingstats.com/third/index.html]. 
%             Considering the following hypothetical experiment. A total of four rats are run in a maze 
%             three trials (T) a day for two days (D). The number of wrong turns on each trial are shown
%             on the below table. Use a significance level = 0.05.
%
%                                                      D1                         D2
%                                          ---------------------------------------------------
%                               Subject       T1       T2       T3       T1       T2       T3
%                             ----------------------------------------------------------------
%                                  1         10        9        7        8        7        5
%                                  2          9        8        5        5        5        3           
%                                  3          9        7        7        7        5        5
%                                  4          8        6        4        3        2        1
%                             ----------------------------------------------------------------
%                                       
%     Data matrix must be:
%     X=[10 1 1 1;9 1 1 2;9 1 1 3;8 1 1 4;9 1 2 1;8 1 2 2;7 1 2 3;6 1 2 4;7 1 3 1;5 1 3 2;7 1 3 3;4 1 3 4;
%     8 2 1 1;5 2 1 2;7 2 1 3;3 2 1 4;7 2 2 1;5 2 2 2;5 2 2 3;2 2 2 4;5 2 3 1;3 2 3 2;5 2 3 3;1 2 3 4];
%
%     Calling on Matlab the function: 
%             RMAOV2(X)
%
%       Answer is:
%
%    The number of IV1 levels are: 2
%
%    The number of IV2 levels are: 3
%
%    The number of subjects are: 4
%
%    Repeated Measures Two-Way Analysis of Variance Table.
%    ---------------------------------------------------------------------------
%    SOV                  SS          df           MS             F        P
%    ---------------------------------------------------------------------------
%    Subjects           43.458         3         14.486[       24.716   0.0000]
%    IV1                45.375         1         45.375        33.000   0.0105
%    Error(IV1)          4.125         3          1.375
%    IV2                30.333         2         15.167        24.818   0.0013
%    Error(IV2)          3.667         6          0.611
%    IV1xIV2             1.000         2          0.500         3.000   0.1250
%    Error(IV1xIV2)      1.000         6          0.167
%    [Error              8.792        15          0.586]
%    Total             128.958        23
%    ---------------------------------------------------------------------------
%    If the P result are smaller than 0.05
%    the corresponding Ho's tested result statistically significant. Otherwise, are not significative.
%    [Generally speaking, no Mean Square is computed for the variable "subjects" since it is assumed
%    that subjects differ from one another thus making a significance test of "subjects" superfluous.
%    However, for all the interested people we are given it anyway].
%  
%    The percentage of the variability in the DV associated with the IV1 is 91.67
%    (After the effects of individual differences have been removed).
%
%    The percentage of the variability in the DV associated with the IV2 is 89.22
%    (After the effects of individual differences have been removed).
%
%    Created by A. Trujillo-Ortiz, R. Hernandez-Walls and R.A. Trujillo-Perez
%               Facultad de Ciencias Marinas
%               Universidad Autonoma de Baja California
%               Apdo. Postal 453
%               Ensenada, Baja California
%               Mexico.
%               atrujo@uabc.mx
%
%    Copyright.July 25, 2004.
%
%    To cite this file, this would be an appropriate format:
%    Trujillo-Ortiz, A., R. Hernandez-Walls and R.A. Trujillo-Perez. (2004). RMAOV2:Two-way repeated
%      measures ANOVA. A MATLAB file. [WWW document]. URL http://www.mathworks.com/matlabcentral/
%      fileexchange/loadFile.do?objectId=5578
%
%    References:
%    Huck, S. W. (2000), Reading Statistics and Research. 3rd. ed. 
%             New-York:Allyn&Bacon/Longman Pub. Chapter 16.

% downloaded from:
% http://www.mathworks.com/matlabcentral/fileexchange/5578
%
% covered under the Simplified BSD License
%
% Copyright (c) 2009, Antonio Trujillo-Ortiz
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

if nargin < 6
  colsubs = 0;
  if nargin < 5
    printTable_tex = 1;
    if nargin < 4
      if max(X(:,2)) > 2 || max(X(:,3)) > 2
        calcGGHF = 1;
      else
        calcGGHF = 0;
      end
      if nargin < 3
        showtable = 1;
        if nargin < 2
          alpha = 0.05; %(default)
        end
      end
    end
  end
end

% parameters for printing the LaTeX table
printGG_tex = 1;
printHF_tex = 1;

if (alpha <= 0 || alpha >= 1)
   fprintf('Warning: significance level must be between 0 and 1\n');
   return;
end;

if nargin < 1, 
   error('Requires at least one input argument.');
end;

% number of IV1 levels
a = max(X(:,2));
% number of IV2 levels
b = max(X(:,3));
% number of subjects
s = max(X(:,4));

% check on IV levels for calculating GG/HF correction
if a <= 2 && b <= 2 && calcGGHF == 1
  fprintf('Turning off GG/HF correction because the number of IV1 levels is %d and IV2 levels is %d.\n\n',a,b);
  calcGGHF = 0;
end

if showtable
  fprintf('The number of IV1 levels are:%2i\n\n', a);
  fprintf('The number of IV2 levels are:%2i\n\n', b);
  fprintf('The number of subjects are:%2i\n\n', s);
end

% collapse across extra IVs that are implicitly present by having multiple
% DV values for each subject within the two DV here (i.e., average across
% each subject's DV values for those IVs).  This will only happen when
% you've set up the data matrix for an ANOVA that has more IVs (i.e.,
% columns) than you want to test here.
if colsubs
  Xcol = [];
  for iv1 = 1:a
    for iv2 = 1:b
      for sub = 1:s
        Xcol = [Xcol; mean(X(X(:,2) == iv1 & X(:,3) == iv2 & X(:,4) == sub,1)), iv1, iv2, sub];
      end
    end
  end
  X = Xcol;
end


indice = X(:,2);
for i = 1:a
    Xe = find(indice==i);
    eval(['A' num2str(i) '=X(Xe,1);']);
end;

indice = X(:,3);
for j = 1:b
    Xe = find(indice==j);
    eval(['B' num2str(j) '=X(Xe,1);']);
end;

indice = X(:,4);
for k = 1:s
    Xe = find(indice==k);
    eval(['S' num2str(k) '=X(Xe,1);']);
end;

C = (sum(X(:,1)))^2/length(X(:,1));  %correction term
SSTO = sum(X(:,1).^2)-C;  %total sum of squares
dfTO = length(X(:,1))-1;  %total degrees of freedom
   
%procedure related to the IV1 (independent variable 1).
A = [];
for i = 1:a
    eval(['x =((sum(A' num2str(i) ').^2)/length(A' num2str(i) '));']);
    A = [A,x];
end;
SSA = sum(A)-C;  %sum of squares for the IV1
dfA = a-1;  %degrees of freedom for the IV1
MSA = SSA/dfA;  %mean square for the IV1

%procedure related to the IV2 (independent variable 2).
B = [];
for j = 1:b
    eval(['x =((sum(B' num2str(j) ').^2)/length(B' num2str(j) '));']);
    B =[B,x];
end;
SSB = sum(B)-C;  %sum of squares for the IV2
dfB = b-1;  %degrees of freedom for the IV2
MSB = SSB/dfB;  %mean square for the IV2

%procedure related to the within-subjects.
S = [];
for k = 1:s
    eval(['x =((sum(S' num2str(k) ').^2)/length(S' num2str(k) '));']);
    S = [S,x];
end;
SSS = sum(S)-C;  %sum of squares for the within-subjects
dfS = k-1;  %degrees of freedom for the within-subjects
MSS = SSS/dfS;  %mean square for the within-subjects

%procedure related to the IV1-error.
for i = 1:a
    for k = 1:s
        Xe = find((X(:,2)==i) & (X(:,4)==k));
        eval(['IV1S' num2str(i) num2str(k) '=X(Xe,1);']);
    end;
end;
EIV1 = [];
for i = 1:a
    for k = 1:s
        eval(['x =((sum(IV1S' num2str(i) num2str(k) ').^2)/length(IV1S' num2str(i) num2str(k) '));']);
        EIV1 = [EIV1,x];
    end;
end;
SSEA = sum(EIV1)-sum(A)-sum(S)+C;  %sum of squares of the IV1-error
dfEA = dfA*dfS;  %degrees of freedom of the IV1-error
MSEA = SSEA/dfEA;  %mean square for the IV1-error


%procedure related to the IV2-error.
for j = 1:b
    for k = 1:s
        Xe = find((X(:,3)==j) & (X(:,4)==k));
        eval(['IV2S' num2str(j) num2str(k) '=X(Xe,1);']);
    end;
end;
EIV2 = [];
for j = 1:b
    for k = 1:s
        eval(['x =((sum(IV2S' num2str(j) num2str(k) ').^2)/length(IV2S' num2str(j) num2str(k) '));']);
        EIV2 = [EIV2,x];
    end;
end;
SSEB = sum(EIV2)-sum(B)-sum(S)+C;  %sum of squares of the IV2-error
dfEB = dfB*dfS;  %degrees of freedom of the IV2-error
MSEB = SSEB/dfEB;  %mean square for the IV2-error

%procedure related to the IV1 and IV2.
for i = 1:a
    for j = 1:b
        Xe = find((X(:,2)==i) & (X(:,3)==j));
        eval(['AB' num2str(i) num2str(j) '=X(Xe,1);']);
    end;
end;
AB = [];
for i = 1:a
    for j = 1:b
        eval(['x =((sum(AB' num2str(i) num2str(j) ').^2)/length(AB' num2str(i) num2str(j) '));']);
        AB = [AB,x];
    end;
end;
SSAB = sum(AB)-sum(A)-sum(B)+C;  %sum of squares of the IV1xIV2
dfAB = dfA*dfB;  %degrees of freedom of the IV1xIV2
MSAB = SSAB/dfAB;  %mean square for the IV1xIV2

%procedure related to the IV1xIV2-error.
SSEAB = SSTO-SSA-SSEA-SSB-SSEB-SSAB-SSS;  %sum of squares of the IV1xIV2-error
dfEAB = dfTO-dfA-dfEA-dfB-dfEB-dfAB-dfS;  %degrees of freedom of the IV1xIV2-error
MSEAB = SSEAB/dfEAB;  %mean square for the IV1xIV2-error

%procedure related to the within-subject error.
SSSE = SSEA+SSEB+SSEAB;
dfSE = dfEA+dfEB+dfEAB;
MSSE = SSSE/dfSE;

%F-statistics calculation
F1 = MSA/MSEA;
F2 = MSB/MSEB;
F3 = MSAB/MSEAB;
F4 = MSS/MSSE;

%degrees of freedom re-definition
v1 = dfA;
v2 = dfEA;
v3 = dfB;
v4 = dfEB;
v5 = dfAB;
v6 = dfEAB;
v7 = dfS;
v8 = dfSE;
v9 = dfTO;

%Probability associated to the F-statistics.
P1 = 1 - fcdf(F1,v1,v2);    
P2 = 1 - fcdf(F2,v3,v4);   
P3 = 1 - fcdf(F3,v5,v6);
P4 = 1 - fcdf(F4,v7,v8);

% eta^2
eta21 = SSA/SSTO;
eta22 = SSB/SSTO;

% partial eta^2
etap21 = SSA/(SSA+SSEA)*100;
etap22 = SSB/(SSB+SSEB)*100;

if P1 >= alpha;
    ds1 ='NS';
else
    ds1 =' S';
end;
if  P2 >= alpha;
    ds2 ='NS';
else
    ds2 =' S';
end;
if  P3 >= alpha;
    ds3 ='NS';
else
    ds3 =' S';
end;
if  P4 >= alpha;
    ds4 ='NS';
else
    ds4 =' S';
end;

% collect all the info into matrices for easy printing
P = [P1 P2 P3];
F = [F1 F2 F3];
SS = [SSA SSB SSAB];
SSE = [SSEA SSEB SSEAB];
MS = [MSA MSB MSAB];
MSE = [MSEA MSEB MSEAB];
% ds = decisions
ds = {ds1 ds2 ds3};
% v = degrees of freedom
v_num = [v1 v3 v5];
v_den = [v2 v4 v6];
% set up the strings
fact_str = {'IV1       ',...
  'IV2       ',...
  'IV1xIV2   '};
if calcGGHF
  fact_GG_str = {'IV1 GG    ',...
    'IV2 GG    ',...
    'IV1xIV2 GG'};
  fact_HF_str = {'IV1 HF    ',...
    'IV2 HF    ',...
    'IV1xIV2 HF'};
end

if calcGGHF
  % calculate Mauchly's values; convert Ps, etc. to GG and HF
  [EpsHF EpsList EpsGG Mau] = GenCalcHFEps(X(:,1),[],X(:,[2 3]),X(:,4));
  
  % set the levels for this particular ANOVA
  levels = {a b [a b]};
  
  % initialize
  P_GG = nan(1,length(EpsGG));
  df_num_GG = nan(1,length(EpsGG));
  df_denom_GG = nan(1,length(EpsGG));
  P_HF = nan(1,length(EpsHF));
  df_num_HF = nan(1,length(EpsHF));
  df_denom_HF = nan(1,length(EpsHF));
  
  for i = 1:length(EpsGG)
    if length(levels{i}) > 1
      k = prod(levels{i} - 1);
    else
      k = levels{i} - 1;
    end
    %adjusted numerator degrees of freedom
    df_num_GG(i) = EpsGG(i)*k;
    %adjusted denominator degrees of freedom
    df_denom_GG(i) = EpsGG(i)*k*(s-1);
    % calculate the GG ps
    P_GG(i) = 1 - fcdf(F(i),df_num_GG(i),df_denom_GG(i));
  end
  
  for i = 1:length(EpsHF)
    if length(levels{i}) > 1
      k = prod(levels{i} - 1);
    else
      k = levels{i} - 1;
    end
    %adjusted numerator degrees of freedom
    df_num_HF(i) = EpsHF(i)*k;
    %adjusted denominator degrees of freedom
    df_denom_HF(i) = EpsHF(i)*k*(s-1);
    % calculate the HF ps
    P_HF(i) = 1 - fcdf(F(i),df_num_HF(i),df_denom_HF(i));
  end
  
  % Print Mauchly's test; factors and interactions with less than 3 levels
  % are not printed
  if showtable
    disp('Mauchly''s Test of Sphericity')
    fprintf('----------------------------------------------\n');
    fprintf('SOV          W       Chi-sq     df        P\n')
    fprintf('----------------------------------------------\n');
    for ind = 1:length(P)
      if sum(levels{ind} > 2) > 0
        fprintf('%s %3.4f%10.4f%7.i%12.4f\n',fact_str{ind},Mau.w(ind),Mau.chiSq(ind),Mau.df(ind),Mau.p(ind));
      end
    end
    fprintf('--------------------------------------------------\n');
    fprintf('If P < %3.2f, assumption of sphericity has been violated.\n\n',alpha);
  end
end

if showtable
  disp('Repeated Measures Two-Way Analysis of Variance Table')
  if calcGGHF
    disp('with Greenhouse-Geisser (G and G, 1959) and Huynh-Feldt (Huynh, 1978) correction values')
  end
  fprintf('------------------------------------------------------------------------------------------\n');
  fprintf('SOV                       SS         df            MS            F       P      Conclusion\n');
  fprintf('------------------------------------------------------------------------------------------\n');
  fprintf('Subjects          %11.3f%10i%15.3f[%13.3f%9.4f]%8s\n\n',SSS,v7,MSS,F4,P4,ds4);
  
  for ind = 1:length(P)
    % IV: SS,v,MS,F,P,ds
    fprintf('%s        %11.3f%10i%15.3f%14.3f%9.4f%9s\n',fact_str{ind},SS(ind),v_num(ind),MS(ind),F(ind),P(ind),ds{ind});
    % Error(IV): SSE,v,MSE
    fprintf('Error: %s %11.3f%10i%15.3f\n',fact_str{ind},SSE(ind),v_den(ind),MSE(ind));
    
    if calcGGHF == 1 && sum(levels{ind} > 2) > 0
      % GG
      fprintf('%s        %11.3f%10.4f%15.3f%14.3f%9.4f\n',fact_GG_str{ind},SS(ind),df_num_GG(ind),MS(ind)/EpsGG(ind),F(ind),P_GG(ind));
      fprintf('Error: %s %11.3f%10.4f%15.3f\n',fact_GG_str{ind},SSE(ind),df_denom_GG(ind),MSE(ind)/EpsGG(ind));
      % HF
      fprintf('%s        %11.3f%10.4f%15.3f%14.3f%9.4f\n',fact_HF_str{ind},SS(ind),df_num_HF(ind),MS(ind)/EpsHF(ind),F(ind),P_HF(ind));
      fprintf('Error: %s %11.3f%10.4f%15.3f\n\n',fact_HF_str{ind},SSE(ind),df_denom_HF(ind),MSE(ind)/EpsHF(ind));
    else
      fprintf('\n');
    end
  end
  
  fprintf('[Error            %11.3f%10i%15.3f]\n\n',SSSE,v8,MSSE);
  fprintf('Total             %11.3f%10i\n\n',SSTO,v9);
  fprintf('------------------------------------------------------------------------------------------\n');

  fprintf('If the P results are smaller than %3.2f\n',alpha);
  fprintf('the corresponding Ho''s tested result statistically significant. Otherwise, are not significative.\n');
  fprintf('[Generally speaking, no Mean Square is computed for the variable "subjects" since it is assumed\n');
  fprintf('that subjects differ from one another thus making a significance test of "subjects" superfluous.\n');
  fprintf('However, for all the interested people we are given it anyway].\n\n');
  
  fprintf('eta^2 = %0.5f\n', eta21);
  fprintf('The percentage of the variability in the DV associated with IV1 (partial eta^2) is %3.2f\n', etap21);
  fprintf('(after the effects of individual differences have been removed).\n\n');
  
  fprintf('eta^2 = %0.5f\n', eta22);
  fprintf('The percentage of the variability in the DV associated with IV2 (partial eta^2) is %3.2f\n', etap22);
  fprintf('(after the effects of individual differences have been removed).\n\n');
  
  %% print a table for use in a LaTeX doc
  if printTable_tex
    fprintf('%% ANOVA summary table for use in a LaTeX document\n');
    fprintf('\\begin{table}\n');
    fprintf('  \\centering\n');
    fprintf('  \\begin{tabular}{lllll}\n');
    fprintf('    \\hline\n');
    fprintf('    Effect & d.f. & $F$ & M.S.E. & $p$ \\\\\n');
    fprintf('    \\hline\n');
    for ind = 1:length(P)
      % print out the normal ANOVA stats
      p_str = [];
      if P(ind) > alpha
        p_str = strrep(sprintf('%.3f',P(ind)),'0.','.');
      elseif P(ind) < 0.0000000001
        p_str = sprintf('<.0000000001');
      elseif P(ind) < 0.000000001
        p_str = sprintf('<.000000001');
      elseif P(ind) < 0.00000001
        p_str = sprintf('<.00000001');
      elseif P(ind) < 0.0000001
        p_str = sprintf('<.0000001');
      elseif P(ind) < 0.000001
        p_str = sprintf('<.000001');
      elseif P(ind) < 0.00001
        p_str = sprintf('<.00001');
      elseif P(ind) < 0.0001
        p_str = sprintf('<.0001');
      elseif P(ind) < 0.001
        p_str = sprintf('<.001');
      elseif P(ind) < 0.01
        p_str = sprintf('<.01');
      elseif P(ind) < 0.05
        p_str = sprintf('<.05');
      else
        p_str = sprintf('<%.2d',alpha);
      end
      fprintf('    %s & $%d,%d$ & $%.3f$ & $%.3f$ & $%s$ \\\\ %% No correction\n',strrep(fact_str{ind},' ',''),v_num(ind),v_den(ind),F(ind),MSE(ind),p_str);
      
      % if we want GG or HF correction, print those stats too
      if calcGGHF == 1 && sum(levels{ind} > 2) > 0
        if printGG_tex
          p_str = [];
          if P_GG(ind) > alpha
            p_str = strrep(sprintf('%.3f',P_GG(ind)),'0.','.');
          elseif P_GG(ind) < 0.0000000001
            p_str = sprintf('<.0000000001');
          elseif P_GG(ind) < 0.000000001
            p_str = sprintf('<.000000001');
          elseif P_GG(ind) < 0.00000001
            p_str = sprintf('<.00000001');
          elseif P_GG(ind) < 0.0000001
            p_str = sprintf('<.0000001');
          elseif P_GG(ind) < 0.000001
            p_str = sprintf('<.000001');
          elseif P_GG(ind) < 0.00001
            p_str = sprintf('<.00001');
          elseif P_GG(ind) < 0.0001
            p_str = sprintf('<.0001');
          elseif P_GG(ind) < 0.001
            p_str = sprintf('<.001');
          elseif P_GG(ind) < 0.01
            p_str = sprintf('<.01');
          elseif P_GG(ind) < 0.05
            p_str = sprintf('<.05');
          else
            p_str = sprintf('<%.2d',alpha);
          end
          fprintf('    %s & $%.2f,%.2f$ & $%.3f$ & $%.3f$ & $%s$ \\\\ %% GG\n',strrep(fact_GG_str{ind},' ',''),df_num_GG(ind),df_denom_GG(ind),F(ind),MSE(ind)/EpsGG(ind),p_str);
        end
        if printHF_tex
          p_str = [];
          if P_HF(ind) > alpha
            p_str = strrep(sprintf('%.3f',P_HF(ind)),'0.','.');
          elseif P_HF(ind) < 0.0000000001
            p_str = sprintf('<.0000000001');
          elseif P_HF(ind) < 0.000000001
            p_str = sprintf('<.000000001');
          elseif P_HF(ind) < 0.00000001
            p_str = sprintf('<.00000001');
          elseif P_HF(ind) < 0.0000001
            p_str = sprintf('<.0000001');
          elseif P_HF(ind) < 0.000001
            p_str = sprintf('<.000001');
          elseif P_HF(ind) < 0.00001
            p_str = sprintf('<.00001');
          elseif P_HF(ind) < 0.0001
            p_str = sprintf('<.0001');
          elseif P_HF(ind) < 0.001
            p_str = sprintf('<.001');
          elseif P_HF(ind) < 0.01
            p_str = sprintf('<.01');
          elseif P_HF(ind) < 0.05
            p_str = sprintf('<.05');
          else
            p_str = sprintf('<%.2d',alpha);
          end
          fprintf('    %s & $%.2f,%.2f$ & $%.3f$ & $%.3f$ & $%s$ \\\\ %% HF\n',strrep(fact_HF_str{ind},' ',''),df_num_HF(ind),df_denom_HF(ind),F(ind),MSE(ind)/EpsHF(ind),p_str);
        end
      end
      
    end
    fprintf('    \\hline\n');
    fprintf('  \\end{tabular}\n');
    fprintf('  \\caption{%d~$\\times$~%d Repeated Measures ANOVA results}\n',a,b);
    fprintf('  \\label{tab:_%dx%d}\n',a,b);
    fprintf('  %%\\ref{tab:_%dx%d}\n',a,b);
    fprintf('\\end{table}\n\n');
  end

end

return;