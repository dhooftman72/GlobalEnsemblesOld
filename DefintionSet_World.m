function     [Parameters, DatasetUse,PointsToTest] = DefintionSet_World(Parameters,RunOverall)
for withinRun = 1:Parameters.JoinedRuns
    PointsToTest(:,withinRun) = sort(randperm(Parameters.TotalPoints,Parameters.NrPoints))'; %#ok<AGROW>
end
DatasetUse = [];
cd(Parameters.InputFolder)
for LoopIn = 1:length(Parameters.filenames)
    [DatasetUse,Parameters] = CollectData(Parameters,LoopIn,PointsToTest,DatasetUse,RunOverall);
end
% save('all')
% ccc
cd(Parameters.currentFolder)
end

function [DatasetUse,Parameters] = CollectData(Parameters,LoopIn,PointsToTest,DatasetUse,RunOverall)
% load data
load(char(Parameters.filenames(LoopIn)));
vars = whos('-file', char(Parameters.filenames(LoopIn)));
Data =  eval(vars(1).name);
for withinRun = 1: Parameters.JoinedRuns
    DatasetUse(:,LoopIn,withinRun) = Data(PointsToTest(:,withinRun));
    % replace missing values with mean value
    if strcmp('Yes',Parameters.CollectWinsor) ~= 1
    List = find(isnan(DatasetUse(:,LoopIn))==1);
    DatasetUse(List,LoopIn,withinRun)= nanmean(DatasetUse(:,LoopIn,withinRun)); %#ok<FNDSB>
    end
end
Parameters.SetNames(LoopIn) = {char(vars(1).name)};
if RunOverall == 1
   [~, Parameters.FullPercentiles(:,LoopIn)] = WinsorFunction(Data,Parameters.WinsorValue,Parameters);
    if strcmp('Yes',Parameters.CollectWinsor) == 1
        TmpX = abs(randn(length(Data),1));
        testArray = [TmpX,Data]; % X, Y
        cd (Parameters.currentFolder)
        [Outputs] = Accuracy_statistics_World(testArray,Parameters);
        if  exist('allWinsorPoints.mat') ~= 0
            load('allWinsorPoints.mat')
        end
        WinsorPoints(:,LoopIn) = Outputs.yes;
        save('allWinsorPoints','WinsorPoints','-v7.3');
        cd(Parameters.InputFolder)
    end
    clear Datatmp
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