% Copy function from Hooftman et al. 2022. Only partially used for this
% work,

function [Outputs] = Accuracy_statistics_World(testArray,Parameters)
% clean data set
testArray(isinf(testArray)==1) = NaN;
orginal_order = 1:1:length(testArray(:,1));
orginal_order = orginal_order';
a1=find((isnan(testArray(:,1))==1));
b=find((isnan(testArray(:,2))==1));
all = [a1;b];
a(:,1) = unique(all);
testArray(a,:) = [];  %#ok<*FNDSB>
orginal_order(a,:) = [];
missing_order = a;
Outputs.datapoints = size(testArray,1);
%% Make log if applicable, but not for Ensembles
testVar = testArray;
%% Correlation stats: Spearman: Overall Accuracy
isnnx = length(find(isnan(testVar(:,1)))==1);
isnny = length(find(isnan(testVar(:,2)))==1);
if (exist('testVar','var') == 1) && (isempty(testVar)~= 1) &&  (isnnx ~= length(testVar(:,1))) &&  (isnny ~= length(testVar(:,2)))
%     MinX = min(testVar(:,1)); 
%     MinY = min(testVar(:,2)); 
    if  Parameters.make_log == 1
        testVar = log10(testVar+1);
    end
    [Outputs.RHO,Outputs.PVAL] = corr(testVar(:,1),testVar(:,2),'type','Spearman');
    Outputs.PVAL = single(Outputs.PVAL);
    Outputs.RHO = (round(Outputs.RHO.*Parameters.Precision(2)))./Parameters.Precision(2);
    %% Inverse deviance against a 1:1 line
    clear x_range y_range
    % Mean double normalised deviation: !Winzorisation!
    % Two sides and on range!!
    if Parameters.ensemble == 0
       [x_range,~] = WinsorFunction(testVar(:,1),Parameters.WinsorValue,Parameters); 
       [y_range,~] = WinsorFunction(testVar(:,2),Parameters.WinsorValue,Parameters); 
    elseif Parameters.ensemble == 1
        x_range= testVar(:,1);
        y_range = testVar(:,2);
        x_range(x_range<0) = 0;
        y_range(y_range<0) = 0;
        x_range(x_range>1) = 1;
        y_range(y_range>1) = 1;
    end
    Outputs.deviation_point= abs(y_range-x_range); %% Accuracy per point
    Outputs.mean_double_deviation = 1- ((sum(abs(y_range-x_range)))/Outputs.datapoints); %%Accuracy overall
    Outputs.LeastSquares = -sum((abs(y_range-x_range)).^2);
    Outputs.std_double_deviation = nanstd(1-(abs(y_range-x_range))); % Output for sensitivities
    Outputs.mean_double_deviation = (round( Outputs.mean_double_deviation.*Parameters.Precision(2)))./Parameters.Precision(2);
    Outputs.std_double_deviation  = (round(Outputs.std_double_deviation.*Parameters.Precision(2)))./Parameters.Precision(2);
    % %% recreate x and y ranges in orginal order, for individual model point
    % storage to create ensembles
Outputs.xes(orginal_order,1) = x_range;
Outputs.xes(missing_order,1) = NaN;
Outputs.yes(orginal_order,1) = y_range;
Outputs.yes(missing_order,1) = NaN;
Outputs.deviation_point(orginal_order,1) = Outputs.deviation_point;
Outputs.deviation_point(missing_order,1) = NaN;
else 
   Outputs.RHO = -9999;
end
end

function  [OutVar,prct] = WinsorFunction(InVar,percLow,Parameters)
InVar = reshape(InVar,((size(InVar,1)).*(size(InVar,2))),1);
InVar(InVar<0) = NaN;
prct =  prctile(InVar,percLow);
if prct < 0
    display('zero Values present; correct this first')
    cccc
end
InVar_norm = InVar - prct;
clear InVar
InVar_norm( InVar_norm<0) = 0;
if isfield(Parameters,'SmallerMargin') ~= 0
   if Parameters.SmallerMargin == 1
        % Remove all zeros from the test set
        InVar_org = InVar_norm;
        testtmp = find(InVar_norm == 0);
        upboud = ceil(((100-percLow)/100)*(length(testtmp)));
        if upboud >= 1
        testlist = testtmp(1:upboud);
        else
        testlist = testtmp;
        end
        InVar_norm(testlist) = [];
        percLow = percLow./2.5; %Elevate percHigh to 99%
    end
end
prct(2) = prctile(InVar_norm,(100-percLow));
if exist('InVar_org') == 1; %#ok<EXIST>
    InVar_norm = InVar_org;
end
OutVar = (InVar_norm./prct(2));
clear InVar_norm InVar_org testlist testtmp upboud
OutVar(OutVar>1) = 1;
end