function [P,X] = RMAOV1_mod(X,alpha,showtable,calcGGHF,printTable_tex,colsubs)
% RMAOV1 Repeated Measures Single-Factor Analysis of Variance Test.
%   One-way repeated measures ANOVA is used to analyze the relationship 
%   between the independent variable and dependent variable when:(1) 
%   the dependent variable is quantitative in nature and is measured on
%   a level that at least approximates interval characteristics, (2) the
%   independent variable is within-subjects in nature, and (3) the
%   independent variable has three or more levels. It is an extension of
%   the correlated-groups t test, where the main advantage is controlling
%   the disturbance variables or individual differences that could influence
%   the dependent variable. 
%   In contrast to the independent groups ANOVA (between-subjects), the 
%   repeated measures procedure is generally more powerful (ie. leads to 
%   a greater likelihood of rejecting a false null hypothesis). Statistically,
%   this is the case because unwanted variability (error) due to individual 
%   differences among participants can be quantified and removed from the 
%   denominator of the F-ratio. This procedure requires fewer participants,
%   but it can give rise to research confounds such as order or practice effects
%   and participant bias.
%   As with the independent groups ANOVA, the total variability of a set of scores
%   is subdivided into sources that provide estimates of population variability. 
%   However with the repeated measures ANOVA, the process of subdividing the 
%   variability is a two step process outlined below.
%
%   Total variability broken down into two components:
%     -Between subjects. Variability in scores due to individual differences 
%      among participants. 
%     -Within subjects. 
%   The within subjects variability is subdivided into the following components: 
%     -Treatment. Variance among the treatment means (same as MS between in the 
%      independent groups ANOVA). 
%     -Residual. Leftover or unwanted error variability (can be thought of an 
%      inconsistencies in the effects of the treatments across the participants). 
%
%   The F-ratio is calculated by dividing the treatment variance (Mean Square 
%   treatment) by the residual variance (Mean Square residual). The numerator 
%   variance is the same as the MS between groups in the independent groups ANOVA
%   while the denominator variance (residual) is the leftover or unwanted variability
%   (error) but without variability due to individual differences among participants
%   which has been separated from it.
%
%   At first glance, type of experiment resambles a randomized complete block design.
%   However, the block design considering groups or sets of subjects (blocks), whereas 
%   the within-subjects design each of the n subjects are exposed to different experimental
%   conditions or trials.
%
%   Syntax: [P,X] = RMAOV1_mod(X,alpha,showtable,calcGGHF,printTable_tex,colsubs)
%      
%     Inputs:
%          X - data matrix (Size of matrix must be n-by-3;dependent variable=column 1,
%              independent variable=column 2;subject=column 3). 
%      alpha - significance level (default = 0.05).
%  showtable - 0 or 1 to show the table of results.
%    Outputs:
%            - Complete Analysis of Variance Table.
%            - Strength of the relationship.
%
%    Example: From the example given by Dr. Matthias Winkel* (http://www.stats.ox.ac.uk/~winkel/phs.html) 
%             on the relaxation therapy against migrane. Nine subjects participated in a relaxation therapy
%             with several weeks baseline frequency/duration recording (w1 and w2) and several weeks
%             therapy (w3 to w5). Its is of interest to test if there exist differences on the relaxation
%             therapy and within subjects with a significance level = 0.05.
%
%                                                           Weeks
%                             ------------------------------------------------------
%                              Subject        1       2       3       4       5
%                             ------------------------------------------------------
%                                  1         21      22       8       6       6
%                                  2         20      19      10       4       4           
%                                  3         17      15       5       4       5
%                                  4         25      30      13      12      17
%                                  5         30      27      13       8       6
%                                  6         19      27       8       7       4
%                                  7         26      16       5       2       5        
%                                  8         17      18       8       1       5       
%                                  9         26      24      14       8       9
%                             ------------------------------------------------------
%                                       
%     *Note: Due to a typing error, on the data table given by Dr. Winkel the value of subject 6 on week 4 must be 7, not 6.
%
%     Data matrix must be:
%     X=[21 1 1;20 1 2;17 1 3;25 1 4;30 1 5;19 1 6;26 1 7;17 1 8;26 1 9;
%     22 2 1;19 2 2;15 2 3;30 2 4;27 2 5;27 2 6;16 2 7;18 2 8;24 2 9;
%     8 3 1;10 3 2;5 3 3;13 3 4;13 3 5;8 3 6;5 3 7;8 3 8;14 3 9;
%     6 4 1;4 4 2;4 4 3;12 4 4;8 4 5;7 4 6;2 4 7;1 4 8;8 4 9;
%     6 5 1;4 5 2;5 5 3;17 5 4;6 5 5;4 5 6;5 5 7;5 5 8;9 5 9];
%
%     Calling on Matlab the function: 
%             RMAOV1(X)
%
%       Answer is:
%
%    The number of IV levels are: 5
%
%    The number of subjects are: 9
%
%    Analysis of Variance Table.
%    --------------------------------------------------------------
%    SOV            SS          df         MS         F        P
%    --------------------------------------------------------------
%    Subjects    486.711         8      60.839[     8.450   0.0000]
%    IV         2449.200         4     612.300     85.042   0.0000
%    Error       230.400        32       7.200
%    Total      3166.311        44
%    --------------------------------------------------------------
%    If the P result is smaller than 0.05
%    the Ho tested results statistically significant. Otherwise, it is not significative.
%    [Generally speaking, no Mean Square is computed for the variable "subjects" since it is assumed
%    that subjects differ from one another thus making a significance test of "subjects" superfluous.
%    However, for all the interested people we are given it anyway].
%
%    The percentage of the variability in the DV associated with the IV is 91.40
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
%    Copyright.July 21, 2004.
%
%    To cite this file, this would be an appropriate format:
%    Trujillo-Ortiz, A., R. Hernandez-Walls and R.A. Trujillo-Perez. (2004). RMAOV1:One-way repeated
%      measures ANOVA. A MATLAB file. [WWW document]. URL http://www.mathworks.com/matlabcentral/
%      fileexchange/loadFile.do?objectId=5576
%
%    References:
%    Huck, S. W. (2000), Reading Statistics and Research. 3rd. ed. 
%             New-York:Allyn&Bacon/Longman Pub. Chapter 16.
%    Winkel, M. http://www.stats.ox.ac.uk/~winkel/phs.html
%    Zar, J. H. (1999), Biostatistical Analysis. 4th. ed.   
%           New-Jersey:Upper Saddle River. p. 255-259.

% downloaded from:
% http://www.mathworks.com/matlabcentral/fileexchange/5576
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
      if max(X(:,2)) > 2
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

% number of IV levels
a = max(X(:,2));
% number of subjects
s = max(X(:,3));

% check on IV levels for calculating GG/HF correction
if a <= 2 && calcGGHF == 1
  fprintf('Turning off GG/HF correction because the number of IV levels is %d.\n\n',a);
  calcGGHF = 0;
end

if showtable
  fprintf('The number of IV levels are:%2i\n\n', a);
  fprintf('The number of subjects are:%2i\n\n', s);
end

% collapse across extra IVs that are implicitly present by having multiple
% DV values for each subject within the one DV here (i.e., average across
% each subject's DV values for those IVs).  This will only happen when
% you've set up the data matrix for an ANOVA that has more IVs (i.e.,
% columns) than you want to test here.
if colsubs
  Xcol = [];
  for iv = 1:a
    for sub = 1:s
      Xcol = [Xcol; mean(X(X(:,2) == iv & X(:,3) == sub,1)), iv, sub];
    end
  end
  X = Xcol;
end

%Analysis of Variance Procedure.
m=[];n=[];nn=[];A=[];
indice = X(:,2);
for i = 1:a
  % get the indices for this IV
  Xe = find(indice==i);
  % get the DV values for this IV
  eval(['X' num2str(i) '=X(Xe,1);']);
  % average across this DV
  eval(['m' num2str(i) '=mean(X' num2str(i) ');'])
  % get the length of this DV
  eval(['n' num2str(i) '=length(X' num2str(i) ') ;'])
  % square the length of this DV
  eval(['nn' num2str(i) '=(length(X' num2str(i) ').^2);'])
  % store these in new variables
  eval(['xm = m' num2str(i) ';'])
  eval(['xn = n' num2str(i) ';'])
  eval(['xnn = nn' num2str(i) ';'])
  % ?get the sum of squares for this DV?
  eval(['x =(sum(X' num2str(i) ').^2)/(n' num2str(i) ');']);
  m=[m;xm];n=[n;xn];nn=[nn,xnn];A=[A,x];
end;


S=[];
indice=X(:,3);
for j=1:s
  % get the indices for this subject
  Xe=find(indice==j);
  % get this subject's DV values across all levels of the IV
  eval(['S' num2str(j) '=X(Xe,1);']);
  % ?get the sum of squares for this subject?
  eval(['x =((sum(S' num2str(j) ').^2)/length(S' num2str(j) '));']);
  S=[S,x]; 
end;

C = (sum(X(:,1)))^2/length(X(:,1)); %correction term
SST = sum(X(:,1).^2)-C; %total sum of squares
dfT = length(X(:,1))-1; %total degrees of freedom

SSA = sum(A)-C; %IV sum of squares
v1 = a-1; %IV degrees of freedom
SSS = sum(S)-C; %within-subjects sum of squares
v2 = s-1; %within-subjects degrees of freedom
SSEA = SST-SSA-SSS; %error sum of squares
v3 = v1*v2; %error degrees of freedom
MSA = SSA/v1; %IV mean squares
MSS = SSS/v2; %within-subjects mean squares
MSEA = SSEA/v3; %error mean squares
F1 = MSA/MSEA; %IV F-statistic
F2 = MSS/MSEA; %within-subjects F-statistic

%Probability associated to the F-statistics.
P1 = 1 - fcdf(F1,v1,v3);    
P2 = 1 - fcdf(F2,v2,v3);   

eta2 = SSA/SST;

etap2 = SSA/(SSA+SSEA)*100;

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

P = P1;
F = F1;
SS = SSA;
SSE = SSEA;
MS = MSA;
MSE = MSEA;
% ds = decisions
ds = {ds1};
% v = degrees of freedom
v_num = v1;
v_den = v2;
% set up the strings
fact_str = {'IV1   '};
if calcGGHF
  fact_GG_str = {'IV1 GG'};
  fact_HF_str = {'IV1 HF'};
end

if calcGGHF
  % calculate Mauchly's values; convert Ps, etc. to GG and HF
  [EpsHF EpsList EpsGG Mau] = GenCalcHFEps(X(:,1),[],X(:,2),X(:,3));
  
  % set the levels for this particular ANOVA
  levels = {a};
  
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
    fprintf('-----------------------------------------------\n');
    fprintf('SOV           W       Chi-sq     df        P\n');
    fprintf('-----------------------------------------------\n');
    for ind = 1:length(P)
      if sum(levels{ind} > 2) > 0
        fprintf('%s %11.4f%10.4f%7.i%12.4f\n',fact_str{ind},Mau.w(ind),Mau.chiSq(ind),Mau.df(ind),Mau.p(ind));
      end
    end
    fprintf('-----------------------------------------------\n');
    fprintf('If P < %3.2f, assumption of sphericity has been violated.\n\n',alpha);
  end
end

if showtable
    disp('Repeated Measures One-Way Analysis of Variance Table')
  if calcGGHF
    disp('with Greenhouse-Geisser (G and G, 1959) and Huynh-Feldt (Huynh, 1978) correction values')
  end
  fprintf('---------------------------------------------------------------------------------------\n');
  fprintf('SOV                  SS          df           MS             F        P      Conclusion\n');
  fprintf('---------------------------------------------------------------------------------------\n');
  fprintf('Subjects      %11.3f%10i%15.3f[%13.3f%9.4f]%8s\n\n',SSS,v2,MSS,F2,P2,ds2);
  
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
  
  %fprintf('IV            %11.3f%10i%15.3f%14.3f%9.4f\n\n',SSA,v1,MSA,F1,P1);
  %fprintf('Error         %11.3f%10i%15.3f\n\n',SSEA,v3,MSEA);
  
  fprintf('Total         %11.3f%10i\n\n',SST,dfT);
  fprintf('---------------------------------------------------------------------------------------\n');

  fprintf('If the P result is smaller than% 3.2f\n', alpha );
  disp('the Ho tested results statistically significant. Otherwise, it is not significative.');
  disp('[Generally speaking, no Mean Square is computed for the variable "subjects" since it is assumed');
  disp('that subjects differ from one another thus making a significance test of "subjects" superfluous.');
  disp('However, for all the interested people we are given it anyway].');
  disp('  ');
  fprintf('eta^2 is% 0.5f\n', eta2);
  fprintf('The percentage of the variability in the DV associated with the IV (partial eta^2) is% 3.2f\n', etap2);
  disp('(After the effects of individual differences have been removed).');
  
  %% print a LaTeX table
  if printTable_tex
    fprintf('%% ANOVA summary table for use in a LaTeX document\n');
    fprintf('\\begin{table}\n');
    fprintf('  \\centering\n');
    fprintf('  \\begin{tabular}{lllll}\n');
    fprintf('    \\hline\n');
    fprintf('    Effect & d.f. & $F$ & M.S.E. & $p$ \\\\\n');
    fprintf('    \\hline\n');
    for ind = 1:length(P)
      % print out the normal stats
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
    fprintf('  \\caption{One-Way Repeated Measures ANOVA results}\n');
    fprintf('  \\label{tab:}\n');
    fprintf('  %%\\ref{tab:}\n');
    fprintf('\\end{table}\n\n');
  end
  
end

return;