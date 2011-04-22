function [P] = RMAOV33_mod(X,alpha,showtable,calcGGHF,printTable_tex)
% RMAOV33 Three-way Analysis of Variance With Repeated Measures on Three Factors Test.
%   This is a three-factor analysis of variance design in which all factors are within-
%   subjects variables. In repeated measures designs, the same participants are used
%   in all conditions. This is like an extreme matching. This allows for reduction of 
%   error variance due to subject factors. Fewer participants can be used in an repeated
%   measures design. Repeated measures designs make it easier to see an effect of the
%   independent variable on the dependent variable (if there is such an effect).
%   Due that there is no way to obtain an independent estimate of error component, for
%   we have only one score per cell, and therefore no within-cell variance. However,
%   each of the interactions with subjects can be shown to serve as a denominator for 
%   an F ratio. So, each effect to be tested has its own error term. Thus every effect
%   is tested by the interaction of that effect with the Subject effect.
%   The SS components are divided up for this design in a way that is best illustrated
%   in a SS tree, as shown:
%                                         
%                                 /    
%                                | SSS  
%                                |           /        / 
%                                |          |        | SSA 
%                                |          |   A   < 
%                                |          |        | [SSEA]
%                                |          |         \         
%                                |          |    
%                                |          |         / 
%                                |          |        | SSB 
%                                |          |   B   < 
%                                |          |        | [SSEB]
%                                |          |         \         
%                                |          |    
%                                |          |         / 
%                                |          |        | SSC 
%                                |          |   C   < 
%                                |          |        | [SSEC]
%                                |          |         \         
%                                |          |    
%                        SSTO   <           |         /
%                                |          |        | SSAxB
%                                | SSW-S   <   AB   <  
%                                |          |        | [SSEAB]
%                                |          |         \
%                                |          | 
%                                |          |         / 
%                                |          |        | SSAxC
%                                |          |  AC   < 
%                                |          |        | [SSEAC]
%                                |          |         \
%                                |          |
%                                |          |         / 
%                                |          |        | SSBxC 
%                                |          |  BC   <
%                                |          |        | [SSEBC]
%                                |          |
%                                |          |         / 
%                                |          |        | SSAxBxC 
%                                |          | ABC   <
%                                |          |        | [SSEABC]
%                                 \          \        \         
%                    
%   Syntax: function [RMAOV33] = RMAOV33(X,alpha) 
%      
%     Inputs:
%          X - data matrix (Size of matrix must be n-by-5;dependent variable=column 1;
%              independent variable 1 (within subjects)=column 2;independent variable 2
%              (within subjects)=column 3; independent variable 3 (within subjects)
%              =column 4; subject=column 5). 
%      alpha - significance level (default = 0.05).
%    Outputs:
%            - Complete Analysis of Variance Table.
%
%    Example: From the example of Howell (2002, p. 510-512). There he examine driver behavior as a 
%             function of two times of day, three types of course, and three models of cars. There
%             were three drivers, each of whom drove each model on each course at each time of day.
%             The dependent variable is the number of steering errors as shown on the below table.
%             Use a significance level = 0.05.
%
%                                  T1                                           T2                     
%             -----------------------------------------------------------------------------------------
%                   C1             C2             C3             C1             C2             C3      
%             -----------------------------------------------------------------------------------------
%    Subject   M1   M2   M3   M1   M2   M3   M1   M2   M3   M1   M2   M3   M1   M2   M3   M1   M2   M3 
%    --------------------------------------------------------------------------------------------------
%       1      10    8    6    9    7    5    7    6    3    5    4    3    4    3    3    2    2    1  
%       2       9    8    5   10    6    4    4    5    2    4    3    3    4    2    2    2    3    2         
%       3       8    7    4    7    4    3    3    4    2    4    1    2    3    3    2    1    0    1 
%    --------------------------------------------------------------------------------------------------
%                                       
%     Data matrix must be:
%     X=[10 1 1 1 1;9 1 1 1 2;8 1 1 1 3;8 1 1 2 1;8 1 1 2 2;7 1 1 2 3;6 1 1 3 1;5 1 1 3 2;4 1 1 3 3;9 1 2 1 1;
%     10 1 2 1 2;7 1 2 1 3;7 1 2 2 1;6 1 2 2 2;4 1 2 2 3;5 1 2 3 1;4 1 2 3 2;3 1 2 3 3;7 1 3 1 1;4 1 3 1 2;
%     3 1 3 1 3;6 1 3 2 1;5 1 3 2 2;4 1 3 2 3;3 1 3 3 1;2 1 3 3 2;2 1 3 3 3;5 2 1 1 1;4 2 1 1 2;4 2 1 1 3;
%     4 2 1 2 1;3 2 1 2 2;1 2 1 2 3;3 2 1 3 1;3 2 1 3 2;2 2 1 3 3;4 2 2 1 1;4 2 2 1 2;3 2 2 1 3;3 2 2 2 1;
%     2 2 2 2 2;3 2 2 2 3;3 2 2 3 1;2 2 2 3 2;2 2 2 3 3;2 2 3 1 1;2 2 3 1 2;1 2 3 1 3;2 2 3 2 1;3 2 3 2 2;
%     0 2 3 2 3;1 2 3 3 1;2 2 3 3 2;1 2 3 3 3];
%
%     Calling on Matlab the function: 
%             RMAOV33(X)
%
%     Answer is:
%
%    The number of IV1 levels are: 2
%    The number of IV2 levels are: 3
%    The number of IV3 levels are: 3
%    The number of subjects are:    3
%
%    Three-Way Analysis of Variance With Repeated Measures on Three Factors (Within-Subjects) Table.
%    ---------------------------------------------------------------------------------------------------
%    SOV                             SS          df           MS             F        P      Conclusion
%    ---------------------------------------------------------------------------------------------------
%    Between-Subjects              24.111         2
%
%    Within-Subjects              316.167        51
%    IV1                          140.167         1        140.167       120.143   0.0082        S
%    Error(IV1)                     2.333         2          1.167
%
%    IV2                           56.778         2         28.389      1022.000   0.0000        S
%    Error(IV2)                     0.111         4          0.028
%
%    IV3                           51.444         2         25.722        92.600   0.0004        S
%    Error(IV3)                     1.111         4          0.278
%
%    IV1xIV2                        5.444         2          2.722         2.085   0.2397       NS
%    Error(IV1xIV2)                 5.222         4          1.306
%
%    IV1xIV3                       16.778         2          8.389        37.750   0.0025        S
%    Error(IV1xIV3)                 0.889         4          0.222
%
%    IV2xIV3                        8.778         4          2.194         3.762   0.0524       NS
%    Error(IV2-IV3)                 4.667         8          0.583
%
%    IV1xIV2xIV3                    2.778         4          0.694         1.923   0.2000       NS
%    Error(IV1-IV2-IV3)             2.889         8          0.361
%    ---------------------------------------------------------------------------------------------------
%    Total                        323.500        53
%    ---------------------------------------------------------------------------------------------------
%    With a given significance level of: 0.05
%    The results are significant (S) or not significant (NS).
%  
%    Created by A. Trujillo-Ortiz, R. Hernandez-Walls and F.A. Trujillo-Perez
%               Facultad de Ciencias Marinas
%               Universidad Autonoma de Baja California
%               Apdo. Postal 453
%               Ensenada, Baja California
%               Mexico.
%               atrujo@uabc.mx
%
%    Copyright.January 10, 2006.
%
%    ---Special thanks are given to Georgina M. Blanc from the Vision Center Laboratory of the 
%       Salk Institute for Biological Studies, La Jolla, CA, for encouraging us to create
%       this m-file-- 
%
%    To cite this file, this would be an appropriate format:
%    Trujillo-Ortiz, A., R. Hernandez-Walls and F.A. Trujillo-Perez. (2006). RMAOV33: Three-way 
%      Analysis of Variance With Repeated Measures on Three Factors Test. A MATLAB file. [WWW document].
%      URL http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9638
%
%    References:
%    Howell, D. C. (2002), Statistical Methods for Psychology. 5th ed. 
%             Pacific Grove, CA:Duxbury Wadsworth Group.

% downloaded from:
% http://www.mathworks.com/matlabcentral/fileexchange/9638
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

if nargin < 5
  printTable_tex = 1;
  if nargin < 4
    if max(X(:,2)) > 2 || max(X(:,3)) > 2 || max(X(:,4)) > 2
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

% parameters for what to print in the LaTeX table
printGG_tex = 1;
printHF_tex = 1;

if (alpha <= 0 || alpha >= 1)
   fprintf('Warning: significance level must be between 0 and 1\n');
   return;
end;

if nargin < 1, 
   error('Requires at least one input argument.');
end;

a = max(X(:,2));
b = max(X(:,3));
c = max(X(:,4));
s = max(X(:,5));

% check on IV levels for calculating GG/HF correction
if a <= 2 && b <= 2 && c <= 2 && calcGGHF == 1
  fprintf('Turning off GG/HF correction because the number of IV1 levels is %d, IV2 levels is %d, and IV3 levels is %d.\n\n',a,b,c);
  calcGGHF = 0;
end

if showtable
  disp('   ');
  fprintf('The number of IV1 levels are:%2i\n', a);
  fprintf('The number of IV2 levels are:%2i\n', b);
  fprintf('The number of IV3 levels are:%2i\n', c);
  fprintf('The number of subjects are:  %2i\n\n', s);
end

CT = (sum(X(:,1)))^2/length(X(:,1));  %correction term
SSTO = sum(X(:,1).^2)-CT;  %total sum of squares
v16 = length(X(:,1))-1;  %total degrees of freedom
   
%procedure related to the subjects.
S = [];
indice = X(:,5);
for l = 1:s
    Xe = find(indice==l);
    eval(['S' num2str(l) '=X(Xe,1);']);
    eval(['x =((sum(S' num2str(l) ').^2)/length(S' num2str(l) '));']);
    S = [S,x];
end;

SSS = sum(S)-CT;
v15 = s-1;

%--Procedure Related to the Within-Subjects--
%procedure related to the IV1 (independent variable 1 [within-subjects]).
A = [];
indice = X(:,2);
for i = 1:a
    Xe = find(indice==i);
    eval(['A' num2str(i) '=X(Xe,1);']);
    eval(['x =((sum(A' num2str(i) ').^2)/length(A' num2str(i) '));']);
    A = [A,x];
end;
SSA = sum(A)-CT;  %sum of squares for the IV1
v1 = a-1;  %degrees of freedom for the IV1
MSA = SSA/v1;  %mean square for the IV1

%procedure related to the IV1-error.
EIV1 = [];
for i = 1:a
    for l = 1:s
        Xe = find((X(:,2)==i) & (X(:,5)==l));
        eval(['IV1S' num2str(i) num2str(l) '=X(Xe,1);']);
        eval(['x =((sum(IV1S' num2str(i) num2str(l) ').^2)/length(IV1S' num2str(i) num2str(l) '));']);
        EIV1 = [EIV1,x];
    end;
end;
SSEA = sum(EIV1)-sum(A)-sum(S)+CT;  %sum of squares of the IV1-error
v2 = v1*v15;  %degrees of freedom of the IV1-error
MSEA = SSEA/v2;  %mean square for the IV1-error

%F-statistics calculation.
F1 = MSA/MSEA;

%Probability associated to the F-statistics.
P1 = 1 - fcdf(F1,v1,v2);    

%procedure related to the IV2 (independent variable 2 [within-subjects]).
B = [];
indice = X(:,3);
for j = 1:b
    Xe = find(indice==j);
    eval(['B' num2str(j) '=X(Xe,1);']);
    eval(['x =((sum(B' num2str(j) ').^2)/length(B' num2str(j) '));']);
    B =[B,x];
end;
SSB = sum(B)-CT;  %sum of squares for the IV2
v3 = b-1;  %degrees of freedom for the IV2
MSB = SSB/v3;  %mean square for the IV2

%procedure related to the IV2-error.
EIV2 = [];
for j = 1:b
    for l = 1:s
        Xe = find((X(:,3)==j) & (X(:,5)==l));
        eval(['IV2S' num2str(j) num2str(l) '=X(Xe,1);']);
        eval(['x =((sum(IV2S' num2str(j) num2str(l) ').^2)/length(IV2S' num2str(j) num2str(l) '));']);
        EIV2 = [EIV2,x];
    end;
end;
SSEB = sum(EIV2)-sum(B)-sum(S)+CT;  %sum of squares of the IV2-error
v4 = v3*v15;  %degrees of freedom of the IV2-error
MSEB = SSEB/v4;  %mean square for the IV2-error

%F-statistics calculation.
F2 = MSB/MSEB;

%Probability associated to the F-statistics.
P2 = 1 - fcdf(F2,v3,v4);    

%procedure related to the IV3 (independent variable 3 [within-subject]).
C = [];
indice = X(:,4);
for k = 1:c
    Xe = find(indice==k);
    eval(['C' num2str(k) '=X(Xe,1);']);
    eval(['x =((sum(C' num2str(k) ').^2)/length(C' num2str(k) '));']);
    C =[C,x];
end;
SSC = sum(C)-CT;  %sum of squares for the IV3
v5 = c-1;  %degrees of freedom for the IV3
MSC = SSC/v5;  %mean square for the IV3

%procedure related to the IV3-error.
EIV3 = [];
for k = 1:c
    for l = 1:s
        Xe = find((X(:,4)==k) & (X(:,5)==l));
        eval(['IV3S' num2str(k) num2str(l) '=X(Xe,1);']);
        eval(['x =((sum(IV3S' num2str(k) num2str(l) ').^2)/length(IV3S' num2str(k) num2str(l) '));']);
        EIV3 = [EIV3,x];
    end;
end;
SSEC = sum(EIV3)-sum(C)-sum(S)+CT;  %sum of squares of the IV3-error
v6 = v5*v15;  %degrees of freedom of the IV3-error
MSEC = SSEC/v6;  %mean square for the IV3-error

%F-statistics calculation.
F3 = MSC/MSEC;

%Probability associated to the F-statistics.
P3 = 1 - fcdf(F3,v5,v6);    

%procedure related to the IV1 and IV2 (within- and within- subject).
AB = [];
for i = 1:a
    for j = 1:b
        Xe = find((X(:,2)==i) & (X(:,3)==j));
        eval(['AB' num2str(i) num2str(j) '=X(Xe,1);']);
        eval(['x =((sum(AB' num2str(i) num2str(j) ').^2)/length(AB' num2str(i) num2str(j) '));']);
        AB = [AB,x];
    end;
end;
SSAB = sum(AB)-sum(A)-sum(B)+CT;  %sum of squares of the IV1xIV2
v7 = v1*v3;  %degrees of freedom of the IV1xIV2
MSAB = SSAB/v7;  %mean square for the IV1xIV2

%procedure related to the IV1-IV2-error.
EIV12 = [];
for i = 1:a
    for j = 1:b
        for l = 1:s
            Xe = find((X(:,2)==i) & (X(:,3)==j) & (X(:,5)==l));
            eval(['IV12S' num2str(i) num2str(j) num2str(l) '=X(Xe,1);']);
            eval(['x =((sum(IV12S' num2str(i) num2str(j) num2str(l) ').^2)/length(IV12S' num2str(i) num2str(j) num2str(l) '));']);
            EIV12 = [EIV12,x];
        end;
    end;
end;
SSEAB = sum(EIV12)-sum(AB)-sum(EIV2)+sum(B)-sum(EIV1)+sum(A)+sum(S)-CT;
v8= v2*v3;  %degrees of freedom of the IV1-IV2-error
MSEAB = SSEAB/v8;  %mean square for the IV1-IV2-error

%F-statistics calculation
F4 = MSAB/MSEAB;

%Probability associated to the F-statistics.
P4 = 1 - fcdf(F4,v7,v8);

%procedure related to the IV1 and IV3 (between- and within- subject).
AC = [];
for i = 1:a
    for k = 1:c
        Xe = find((X(:,2)==i) & (X(:,4)==k));
        eval(['AC' num2str(i) num2str(k) '=X(Xe,1);']);
        eval(['x =((sum(AC' num2str(i) num2str(k) ').^2)/length(AC' num2str(i) num2str(k) '));']);
        AC = [AC,x];
    end;
end;
SSAC = sum(AC)-sum(A)-sum(C)+CT;  %sum of squares of the IV1xIV3
v9 = v1*v5;  %degrees of freedom of the IV1xIV3
MSAC = SSAC/v9;  %mean square for the IV1xIV3

%procedure related to the IV1-IV3-error.
EIV13 = [];
for i = 1:a
    for k = 1:c
        for l = 1:s
            Xe = find((X(:,2)==i) & (X(:,4)==k) & (X(:,5)==l));
            eval(['IV13S' num2str(i) num2str(k) num2str(l) '=X(Xe,1);']);
            eval(['x =((sum(IV13S' num2str(i) num2str(k) num2str(l) ').^2)/length(IV13S' num2str(i) num2str(k) num2str(l) '));']);
            EIV13 = [EIV13,x];
        end;
    end;
end;
SSEAC = sum(EIV13)-sum(AC)-sum(EIV3)+sum(C)-sum(EIV1)+sum(A)+sum(S)-CT;
v10 = v2*v5;  %degrees of freedom of the IV1-IV3-error
MSEAC = SSEAC/v10;  %mean square for the IV1-IV3-error

%F-statistics calculation
F5 = MSAC/MSEAC;

%Probability associated to the F-statistics.
P5 = 1 - fcdf(F5,v9,v10);

%procedure related to the IV2 and IV3 (within- and within- subject).
BC = [];
for j = 1:b
    for k = 1:c
        Xe = find((X(:,3)==j) & (X(:,4)==k));
        eval(['BC' num2str(j) num2str(k) '=X(Xe,1);']);
        eval(['x =((sum(BC' num2str(j) num2str(k) ').^2)/length(BC' num2str(j) num2str(k) '));']);
        BC = [BC,x];
    end;
end;
SSBC = sum(BC)-sum(B)-sum(C)+CT;  %sum of squares of the IV2xIV3
v11 = v3*v5;  %degrees of freedom of the IV2xIV3
MSBC = SSBC/v11;  %mean square for the IV2xIV3

%procedure related to the IV2-IV3-error.
EIV23 = [];
for j = 1:b
    for k = 1:c
        for l = 1:s
            Xe = find((X(:,3)==j) & (X(:,4)==k) & (X(:,5)==l));
            eval(['IV23S' num2str(j) num2str(k) num2str(l) '=X(Xe,1);']);
            eval(['x =((sum(IV23S' num2str(j) num2str(k) num2str(l) ').^2)/length(IV23S' num2str(j) num2str(k) num2str(l) '));']);
            EIV23 = [EIV23,x];
        end;
    end;
end;
SSEBC = sum(EIV23)-sum(BC)-sum(EIV3)+sum(C)-sum(EIV2)+sum(B)+sum(S)-CT;
v12 = v4*v5;  %degrees of freedom of the IV2-IV3-error
MSEBC = SSEBC/v12;  %mean square for the IV2-IV3-error

%F-statistics calculation
F6 = MSBC/MSEBC;

%Probability associated to the F-statistics.
P6 = 1 - fcdf(F6,v11,v12);

%procedure related to the IV1, IV2 and IV3 (within, within- and within- subject).
ABC = [];
for i = 1:a
    for j = 1:b
        for k = 1:c
            Xe = find((X(:,2)==i) & (X(:,3)==j) & (X(:,4)==k));
            eval(['AB' num2str(i) num2str(j) num2str(k) '=X(Xe,1);']);
            eval(['x =((sum(AB' num2str(i) num2str(j) num2str(k) ').^2)/length(AB' num2str(i) num2str(j) num2str(k) '));']);
            ABC = [ABC,x];
        end;
    end;
end;
SSABC = sum(ABC)+sum(A)+sum(B)+sum(C)-sum(AB)-sum(AC)-sum(BC)-CT;  %sum of squares of the IV1xIV2xIV3
v13 = v1*v3*v5;  %degrees of freedom of the IV1xIV2xIV3
MSABC = SSABC/v13;  %mean square for the IV1xIV2xIV3

%procedure related to the IV1-IV2-IV3-error.
EIV123 = [];
for i = 1:a
    for j = 1:b
        for k = 1:c
            for l = 1:s
                Xe = find((X(:,2)==i) &(X(:,3)==j) & (X(:,4)==k) & (X(:,5)==l));
                eval(['IV123S' num2str(i) num2str(j) num2str(k) num2str(l) '=X(Xe,1);']);
                eval(['x =((sum(IV123S' num2str(i) num2str(j) num2str(k) num2str(l) ').^2)/length(IV123S' num2str(i) num2str(j) num2str(k) num2str(l) '));']);
                EIV123 = [EIV123,x];
            end;
        end;
    end;
end;
SSEABC = sum(EIV123)-sum(ABC)-sum(EIV23)+sum(BC)-sum(EIV13)+sum(AC)+sum(EIV3)-sum(C)-sum(EIV12)+sum(AB)+sum(EIV2)-sum(B)+sum(EIV1)-sum(A)-sum(S)+CT;  %sum of squares of the IV1-IV2-IV3-error
v14 = v2*v3*v5;  %degrees of freedom of the IV1-IV2-IV3-error
MSEABC = SSEABC/v14;  %mean square for the IV1-IV2-IV3-error

%F-statistics calculation
F7 = MSABC/MSEABC;

%Probability associated to the F-statistics.
P7 = 1 - fcdf(F7,v13,v14);

SSWS = SSA+SSEA+SSB+SSEB+SSC+SSEC+SSAB+SSEAB+SSAC+SSAC+SSEAC+SSBC+SSEBC+SSABC+SSEABC;
vWS = v1+v2+v3+v4+v5+v6+v7+v8+v9+v10+v11+v12+v13+v14;

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
if P5 >= alpha;
    ds5 ='NS';
else
    ds5 =' S';
end;
if  P6 >= alpha;
    ds6 ='NS';
else
    ds6 =' S';
end;
if  P7 >= alpha;
    ds7 ='NS';
else
    ds7 =' S';
end;

% eta^2 values were not in the original file
eta21 = SSA/(SSA+SSEA)*100;
eta22 = SSB/(SSB+SSEB)*100;
eta23 = SSC/(SSC+SSEC)*100;

% collect all the info into matrices for easy printing
P = [P1 P2 P3 P4 P5 P6 P7];
F = [F1 F2 F3 F4 F5 F6 F7];
SS = [SSA SSB SSC SSAB SSAC SSBC SSABC];
SSE = [SSEA SSEB SSEC SSEAB SSEAC SSEBC SSEABC];
MS = [MSA MSB MSC MSAB MSAC MSBC MSABC];
MSE = [MSEA MSEB MSEC MSEAB MSEAC MSEBC MSEABC];
% ds = decisions
ds = {ds1 ds2 ds3 ds4 ds5 ds6 ds7};
% v = degrees of freedom
v_num = [v1 v3 v5 v7 v9 v11 v13];
v_den = [v2 v4 v6 v8 v10 v12 v14];
% set up the strings
fact_str = {'IV1           ',...
  'IV2           ',...
  'IV3           ',...
  'IV1xIV2       ',...
  'IV1xIV3       ',...
  'IV2xIV3       ',...
  'IV1xIV2xIV3   '};
if calcGGHF
  fact_GG_str = {'IV1 GG        ',...
    'IV2 GG        ',...
    'IV3 GG        ',...
    'IV1xIV2 GG    ',...
    'IV1xIV3 GG    ',...
    'IV2xIV3 GG    ',...
    'IV1xIV2xIV3 GG'};
  fact_HF_str = {'IV1 HF        ',...
    'IV2 HF        ',...
    'IV3 HF        ',...
    'IV1xIV2 HF    ',...
    'IV1xIV3 HF    ',...
    'IV2xIV3 HF    ',...
    'IV1xIV2xIV3 HF'};
end

if calcGGHF
  % calculate Mauchly's values; convert Ps, etc. to GG and HF
  [EpsHF EpsList EpsGG Mau] = GenCalcHFEps(X(:,1),[],X(:,[2 3 4]),X(:,5));
  
  % set the levels for this particular ANOVA
  levels = {a b c [a b] [a c] [b c] [a b c]};
  
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
    fprintf('--------------------------------------------------\n');
    fprintf('SOV              W       Chi-sq     df        P\n')
    fprintf('--------------------------------------------------\n');
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
  disp('Three-Way Analysis of Variance With Repeated Measures on Three Factors (Within-Subjects) Table')
  if calcGGHF
    disp('with Greenhouse-Geisser (G and G, 1959) and Huynh-Feldt (Huynh, 1978) correction values')
  end
  fprintf('----------------------------------------------------------------------------------------------\n');
  fprintf('SOV                           SS         df            MS            F       P      Conclusion\n');
  fprintf('----------------------------------------------------------------------------------------------\n');
  fprintf('Between-Subjects      %11.3f%10i\n\n',SSS,v15);
  fprintf('Within-Subjects       %11.3f%10i\n\n',SSWS,vWS);
  
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
  
  fprintf('----------------------------------------------------------------------------------------------\n');
  fprintf('Total                 %11.3f%10i\n',SSTO,v16);
  fprintf('----------------------------------------------------------------------------------------------\n');
  fprintf('With a given significance level of: %.2f\n', alpha);
  fprintf('The results are significant (S) or not significant (NS).\n\n');
  
  fprintf('The percentage of the variability in the DV associated with IV1 (partial eta^2) is %3.2f\n', eta21);
  fprintf('(after the effects of individual differences have been removed).\n\n');
  
  fprintf('The percentage of the variability in the DV associated with IV2 (partial eta^2) is %3.2f\n', eta22);
  fprintf('(after the effects of individual differences have been removed).\n\n');
  
  fprintf('The percentage of the variability in the DV associated with IV3 (partial eta^2) is %3.2f\n', eta23);
  fprintf('(after the effects of individual differences have been removed).\n\n');
  
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
    fprintf('  \\caption{%d~$\\times$~%d~$\\times$~%d Repeated Measures ANOVA results}\n',a,b,c);
    fprintf('  \\label{tab:_%dx%dx%d}\n',a,b,c);
    fprintf('  %%\\ref{tab:_%dx%dx%d}\n',a,b,c);
    fprintf('\\end{table}\n\n');
  end
  
end

return;