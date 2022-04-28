% The function creating the weights for the weighted ensembles.
% This is a shortend version of Hooftman et al. 2022
% It is build as main function with multiple nested subfunctions below


function  [EnsemblePoints,Weighting] = Make_Ensembles_World(PointsIN,Weighting,Parameters)
%model combine module, based on Combination
warning off
%Reset parameters
Parameters.make_log = 0;
Parameters.ensemble = 1;
EnsemblePoints = [];
%HIMaxentStore = [];
%% All preparations
% Set data_set
mod_tmp = [];
for i=1:1: Parameters.data_set_max
    model_Values = PointsIN(:,i);
    mod_tmp = [mod_tmp,model_Values]; %#ok<*AGROW>
end
model_values_store = mod_tmp;
clear i model_Values val_tmp* mod_tmp*
%% Ensemble stats
for Ensemble = 1:Parameters.Nr_ensembles
    [ValOut,Weighting] = Actual_Ensembles(Weighting,Parameters,model_values_store,Ensemble);
     if strcmp('Yes',Parameters.CollectWinsor) == 1
         EnsemblePoints(:,Ensemble) = ValOut; 
     end
end
end
%%
function [ValOut,Weighting] = Actual_Ensembles(Weighting,Parameters,model_values_store,Ensemble)
name = Parameters.Ensemble_Names(Ensemble);
str = sprintf('Ensemble running = %s ',char(name));
disp(str)
clear InvModels new_weights model_values Weights NewWeights
if Ensemble == 1
    %PCA weights
    InvModels = model_values_store;
    [InvModels, ~] = ModelClean(InvModels, []);
    [coefs,~,~,~] = princomp(InvModels);
    PCA_weights = coefs(:,1);
    InvModels = model_values_store; % Restore full dataset
    [ValOut,~,new_weights] = Weightin_algo(InvModels,PCA_weights,Parameters);
    Weighting.PCA(:,1) = new_weights;
    clear PCA_weights
elseif Ensemble == 2 % Corrcoef
    %disp('      Correlation coefficient among ouputs')
    InvModels = model_values_store;
    [InvModels, ~] = ModelClean(InvModels, []);
    NewWeights = mean(corrcoef(InvModels));
    Weights = ShapeWeights(NewWeights,Parameters);
    InvModels = model_values_store; % Restore full dataset
    [ValOut,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
    Weighting.CorCoef(:,1) = Weights;
    clear Weights NewWeights
elseif Ensemble == 3 %Regression to the median
    InvModels = model_values_store;
    Predicted(:,1) = (nanmedian(InvModels,2));
    [InvModels, Predicted] = ModelClean(InvModels, Predicted);
    %disp('      Regressing fit to median')
    [Median_Models]  = RegresAlgo(Parameters,InvModels,Predicted);
    InvModels = model_values_store; % Restore full dataset
    [ValOut,~,new_weights] = Weightin_algo(InvModels,Median_Models,Parameters);
    Weighting.RegresstoMedian(:,1) = new_weights; %Store to compare
    clear  Median_Models Weights NewWeights
elseif Ensemble == 4 %Leave one out regression
    InvModels = model_values_store;
    %Step 1: Estimate all one to all fits, restrict weights and calcuate deviance
    step = 1;
    nrM = size(model_values_store,2);
    for i = 1:nrM
        NewWeights(:,i) = LeaveOutRegress(model_values_store,i,nrM,Parameters);
    end   
    for i = 1:1:nrM
        take_array = RemvOne(nrM,i);
        Weights(i,1) = nanmean(NewWeights(i,take_array));
    end
    Weights = ShapeWeights(Weights,Parameters);
    InvModels = model_values_store; % Restore full dataset
    [ValOut,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
    Weighting.LeaveOneOut(:,1) = new_weights;
     clear Weights NewWeights
elseif Ensemble == 5 % GridSize
    %disp('      Gridsize differences among maps')
    NewWeights = reshape((1./(log10(Parameters.GridSizes+1))),length(Parameters.GridSizes),1);
    Weights = ShapeWeights(NewWeights,Parameters);
    InvModels = model_values_store; % Restore full dataset
    [ValOut,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
    Weighting.GridSize(:,1) = Weights;
    clear Weights NewWeights
elseif Ensemble == 6 % Mean Ensemble
    ValOut = (nanmean(model_values_store,2)); %The recalculation method
elseif Ensemble == 7 % Median Ensemble
    ValOut = (nanmedian(model_values_store,2)); %The recalculation method
elseif Ensemble == 8 % Variation among models
    ValOut = (nanstd(model_values_store,1,2))./sqrt(Parameters.NumberofModels.NumberofModels); %The recalculation method
end
end
%% Subfunctions
%%
function WeightOut = LeaveOutRegress(model_values_store,i,nrM,Parameters)
InvModels = model_values_store;
take_array = RemvOne(nrM,i);
Predicted(:,1) = InvModels(:,i);
InvModels =  InvModels(:,take_array);
[InvModels, Predicted] = ModelClean(InvModels, Predicted);
%         str = sprintf('      Leave one out Regression = %s ',char(Parameters.SetNames(i)));
%         disp(str)
[weightings]  = RegresAlgo(Parameters,InvModels,Predicted);
[~,~,new_weights] = Weightin_algo(InvModels,weightings,Parameters);
WeightOut(take_array,1) = new_weights;
if i == nrM
    WeightOut(nrM,1) = 0;
end
end

%%
function [ValOut,StdOut,WdS] = Weightin_algo(InvM,WdS,Parameters)
for k = 1:1:length(InvM)
    WdS = ShapeWeights(WdS,Parameters);
    nomi = 0;
    denomi = 0;
    for j = 1:1:length(WdS)
        if isnan(InvM(k,j))~= 1
            nomi = nomi + (InvM(k,j).* WdS(j));
            denomi = denomi +  WdS(j);
        end
    end
    % denomi keeps needed since sometimes values drop out as NaN for indiviudal data points
    % so sum of weights < 1 so correction is needed
    ValOut(k,1) = nomi/denomi;
    clear nomi
    nomi = 0; %denomi is same as above
    for j = 1:1:length(WdS)
        if isnan(InvM(k,j))~= 1
            nomi = nomi + (((InvM(k,j)- ValOut(k))^2).*WdS(j));
        end
    end
    StdOut(k,1) =  nomi/denomi;
    clear nomi denomi
end
end

%%
function  InV = ShapeWeights(InV,Parameters)
InV = reshape(InV,(length(InV)),1);
InV= max(InV,0);
InV = InV./sum(InV);
InV= min(max(InV,0),1);
InV = (round(InV.*Parameters.Precision(1)))./Parameters.Precision(1);
end


function [Inv,InvP] = ModelClean(Inv,InvP)
num = [];
for f = size(Inv,1):-1:1
    Inv(f,(isnan(Inv(f,:))==1)) = nanmean(Inv(f,:));
    testArray = find(isnan(Inv(f,:))~=1);
    if isempty(testArray==1)
        num = [num;f];
    end
end
if isempty(InvP)~=1
    num = unique([num;(find(isnan(InvP)==1))]);
    InvP(num) = [];
end
Inv(num,:) =[];
end

function take_array = RemvOne(nrM,i)
take_array = 1:nrM;
take_array(take_array==i) = [];
end

%% The regression functions
function [WghtOut] = RegresAlgo(Parameters,InParY,InParX)
[modelFunw,~, prior,~,~] = CreateModelFun(InParY);
Options = statset('FunValCheck','off','Display','off','MaxIter',200,...
    'TolFun',1.0000e-4, 'TolX',1.0e-4, 'DerivStep', 6.0555e-06,'Robust',...
    'off', 'WgtFun', 'bisquare');
[WghtOut,~,~,~,~]  = nlinfit(InParY,InParX,modelFunw,prior,'Options',Options);
WghtOut = ShapeWeights(WghtOut,Parameters);
end

function [modelFunw,series, prior,Options,Group] = CreateModelFun(InvModels)
series = [1:size(InvModels,2)];
prior = (ones((size(InvModels,2)),1)/(size(InvModels,2)))';
if (size(InvModels,2)) == 11
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9))+(b(10).*InvModels(:,10))+(b(11).*InvModels(:,11));
elseif (size(InvModels,2)) == 10
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9))+(b(10).*InvModels(:,10));
elseif (size(InvModels,2)) == 9
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9));
elseif (size(InvModels,2)) == 8
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8));
elseif (size(InvModels,2)) == 7
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7));
elseif (size(InvModels,2)) == 6
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6));
elseif (size(InvModels,2)) == 5
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5));
elseif (size(InvModels,2)) == 4
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) ;
elseif (size(InvModels,2)) == 3
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3));
elseif (size(InvModels,2)) == 12
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9))+ (b(10).*InvModels(:,10)) + (b(11).*InvModels(:,11)) + ...
        (b(12).*InvModels(:,12));
elseif (size(InvModels,2)) == 13
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9))+ (b(10).*InvModels(:,10)) + (b(11).*InvModels(:,11)) + ...
        (b(12).*InvModels(:,12)) + (b(13).*InvModels(:,13));
elseif (size(InvModels,2)) == 14
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9))+ (b(10).*InvModels(:,10)) + (b(11).*InvModels(:,11)) + ...
        (b(12).*InvModels(:,12)) + (b(13).*InvModels(:,13))+(b(14).*InvModels(:,14));
end
modelFunw = @(b,InvModels)modelFun(b,InvModels);
Options = statset('FunValCheck','on','Display','off','MaxIter',100,...
    'TolFun',1.0000e-2, 'TolX',1.0e-2, 'Jacobian','off', 'DerivStep', 6.0555e-03, 'OutputFcn',' ');
Group = grp2idx(InvModels(:,1));
end