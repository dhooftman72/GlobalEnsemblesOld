% Main module for generating weighted ensembles from gridded Global model
% outputs following Hooftman et al. 2022. This module calculates the
% weights that will be used for input numbers in arcpy code files.
% Mean and median ensembles are calculated within arcgis

%Inputs are defined in the parameter files (example below added)
% Output is the weights to be excported to arcpy codes (see there)
% Via Parameters.CollectWinsor = 'Yes'), the individual outputs can be
% saved, whioch is done for accumulated water

function Main_TheRuns_WorldPublish
warning off
clc
close all hidden
if matlabpool('size') ~= 0
    matlabpool close
end
warning off
if exist('Output_Dir','dir') ~= 7
    mkdir('Output_Dir')
end

% define parameters, amount of runs and model collect model dara
Parameters = [];
Parameters = DefineParameters(Parameters);
% Double loop in tranches, selecting bootstrapped random grid cells.
% Afterwards those are run with parfor within tranches
for RunOverall= 1:Parameters.runOverall
    display(['Running = ',int2str(RunOverall),'th overall run'])
    display ('Collecting data')
    %Collect data and parameters
    cd(Parameters.currentFolder)
    % Collect bootstrapped grid cell sets per tranch
    % See seperate DefinitionSet_World codes
    [Parameters, DatasetUse,~] = DefintionSet_World(Parameters,RunOverall);
    if (ceil(Parameters.JoinedRuns./2)) <=10
        eval(['matlabpool open Full ',int2str(ceil(Parameters.JoinedRuns./2))])
    else
        matlabpool open Full
    end
    display('Running paralel runs')
    parfor withinRun = 1:Parameters.JoinedRuns
        run = ((RunOverall-1).*(Parameters.JoinedRuns))+ withinRun; %#ok<PFBNS>
        Actual_Runs(Parameters,run, DatasetUse(:,:,withinRun)); % Codes Below
    end
    if matlabpool('size') ~= 0
        matlabpool close
    end
    clear DatasetUse PointsTested
end
% Combining all runs via JoinFunc below
display ('Combining all Runs')
display ('  ')
[FinalWeights, AllRunWeights] = JoinFunc(Parameters); %#ok<NASGU
if strcmp('Yes',Parameters.CollectWinsor) == 1
    % This for water supply, so the values per catchment can be transported to arcgis
    % In that case the number of runs in the parameter file is set as 1,
    % so it takes all points
    load('allEnsemblesPoints')
    [WaterAccuEnsembles,WaterAccuModels] = WaterAccumulatedCodes(FinalWeights,WinsorPoints); %#ok<NASGU,ASGLU>
    save(Parameters.output_file,'WaterAccuEnsembles','WaterAccuModels','FinalWeights','AllRunWeights','Parameters');
    delete('allEnsemblePoints.mat')
else
    %for grid outputs
    save(Parameters.output_file,'FinalWeights','AllRunWeights','Parameters');
end
str = sprintf('Ready with genereating Ensembles %s ',Parameters.output_file);
disp(str)
end

%% Subfunction Actual Runs 
% Asking for seperate functions: 
%  Accuracy_statistics_World
%  Make_Ensembles_World
function Actual_Runs(Parameters,run,DatasetIn)
Parameters.data_set_max =  size(DatasetIn,2);
TmpX = abs(randn(Parameters.NrPoints,1));
for data_set = 1:1:Parameters.data_set_max
    if data_set == 1 % Initiate outputs
        Points = [];
        Weighting= [];
    end
    clear testArray Outputs % make sure they are not lingering around
    testArray = [TmpX,DatasetIn(:,data_set)]; % X, Y
    % Calculate accuracy statistics and normalised model points per boots
    % via seperate Accuracy_statistics_World codes. The accuracy statistics
    % are not removed from the code, although not used here.
    [Outputs] = Accuracy_statistics_World(testArray,Parameters);
    Points.Models(:,data_set) = Outputs.yes;
    Points.Validation(:,data_set) = Outputs.xes;
    clear testArray Outputs Sizes
    %%  Ensembles run on results from individual model runs after 
    % all models have been finished. Seperate code file: Make_Ensembles_World
    if data_set ==  Parameters.data_set_max
        [EnsemblePoints,Weighting] = Make_Ensembles_World(Points.Models,Weighting,Parameters); %#ok<ASGLU>
        % Save the data away in a file
        cd('Output_Dir')
        Output_file = [Parameters.output_file,'_',int2str(run)];
        save(Output_file,'Weighting','Parameters');
        cd ..
        if strcmp('Yes',Parameters.CollectWinsor) == 1 % For water supply
            WinsorPoints = Points.Models;  %#ok<NASGU>
            save('allEnsemblesPoints','EnsemblePoints','WinsorPoints','-v7.3');
        end
    end
    clear  testArray
end
clear DatasetUse TmpX Results Points Weighting
end

%% Subfunction combine runs for weights
function [FinalWeights, Weights] = JoinFunc(Parameters)
cd('Output_Dir')
% Loop through all output files
for run = 1:(Parameters.runOverall.*Parameters.JoinedRuns)
    Output_file = [Parameters.output_file,'_',int2str(run),'.mat'];
    load(Output_file,'Weighting','UsedPoints');
    for i = 1:5
    Weights.(genvarname(char(Parameters.Ensemble_Names(i))))(:,run) =  Weighting.(genvarname(char(Parameters.Ensemble_Names(i))))(:,1);
    end
end
cd ..

for i = 1:5
FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))) = dataset(Parameters.SetNames','Varnames',char({'Names'}));
FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))).Mean(:,1) = nanmean(Weights.(genvarname(char(Parameters.Ensemble_Names(i)))),2);
FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))).Median(:,1) = nanmedian(Weights.(genvarname(char(Parameters.Ensemble_Names(i)))),2);
FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))).STD(:,1) = nanstd(Weights.(genvarname(char(Parameters.Ensemble_Names(i)))),1,2);
FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))).CV(:,1) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))).STD./...
                                                                         FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))).Mean;
FinalWeights.(genvarname(char(Parameters.Ensemble_Names(i)))).Range(:,1) = prctile(Weights.(genvarname(char(Parameters.Ensemble_Names(i)))),97.5,2) -...
                                                                           prctile(Weights.(genvarname(char(Parameters.Ensemble_Names(i)))),2.5,2);
end
end
%% Subfunction Parameter define
function Parameters = DefineParameters(Parameters)
Parameters.currentFolder = pwd;
display ('Setting Parameters for:')
Parameters = ParametersCarbon; % Below as example
%cd('ParameterFiles')
%ParametersWaterAccuGRDC
%ParametersWaterAccuHydro
%ParametersGrazing
%ParametersFuelWood
%ParametersRecreation
end


%% Collate Water Ensembles
function [WaterAccuEnsembles,WaterAccuModels] = WaterAccumulatedCodes(FinalWeights,WinsorPoints, Parameters)
WaterAccuModels = dataset(WinsorPoints(:,1),'Varnames',Parameters.filenames(1));
for i = 2:8
    WaterAccuModels.(genvarname(char(Parameters.filenames(i)))) = WinsorPoints(:,i);
end
Mean_Ensemble(:,1) = nanmean(WinsorPoints,2);
Mean_Ensemble = Mean_Ensemble./max(Mean_Ensemble); % It is normalised, so don't redo the 2.5% procedure
Median_Ensemble(:,1) = nanmedian(WinsorPoints,2);
Median_Ensemble = Median_Ensemble./max(Median_Ensemble);  % It is normalised, so don't redo the 2.5% procedure
WaterAccuEnsembles = dataset(Mean_Ensemble,'Varnames','MeanEnsemble');
WaterAccuEnsembles.MedianEnsemble = Median_Ensemble;

for x = 1:1:5
    WinsorWeights (1:length(WinsorPoints),1) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(1);
    WinsorWeights (1:length(WinsorPoints),2) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(2);
    WinsorWeights (1:length(WinsorPoints),3) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(3);
    WinsorWeights (1:length(WinsorPoints),4) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(4);
    WinsorWeights (1:length(WinsorPoints),5) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(5);
    WinsorWeights (1:length(WinsorPoints),6) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(6);
    WinsorWeights (1:length(WinsorPoints),7) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(7);
    WinsorWeights (1:length(WinsorPoints),8) = FinalWeights.(genvarname(char(Parameters.Ensemble_Names(x)))).Mean(8);
    list = find(isnan(WinsorPoints));
    WinsorWeights(list) = 0; %#ok<*AGROW>
    WinsorUse = WinsorPoints;
    WinsorUse(list) = 0;
    for i = 1:length(WinsorWeights)
        tot = sum(WinsorWeights(i,:));
        WinsorWeight(1,:) = WinsorWeights(i,:)./tot;
        EnsembleT(i,1) = (WinsorWeight(1).*WinsorUse(i,1)) +...
            (WinsorWeight(2).*WinsorUse(i,2)) +...
            (WinsorWeight(3).*WinsorUse(i,3)) +...
            (WinsorWeight(4).*WinsorUse(i,4)) +...
            (WinsorWeight(5).*WinsorUse(i,5)) +...
            (WinsorWeight(6).*WinsorUse(i,6)) +...
            (WinsorWeight(7).*WinsorUse(i,7)) +...
            (WinsorWeight(8).*WinsorUse(i,8));
    end
    WaterAccuEnsembles.(genvarname(char(Parameters.Ensemble_Names(x)))) = EnsembleT./max(EnsembleT);
    clear EnsembleT
end
end

%% Example parameter file
% Carbon
function Parameters = ParametersCarbon
Parameters.output_file= 'Carbon_tons';
Parameters.InputFolder = 'CarbonData';
Parameters.CollectWinsor = 'No'; %'Yes' for water supply
display(Parameters.output_file)
display ('  ')
Parameters.ServiceName = 'Carbon_tons';
% The Model output grid files
Parameters.filenames = ...
    {'1_AriesCarbonOnly.mat';...
    '2_InVestPlusOnly.mat';...
    '3_CNCarbonOnly.mat';...
    '4_TEEBCarbonOnly.mat';...
    '5_LPJ2015Only.mat';...
    '6_NoonConIntOnly.mat';...
    '7_KindermannOnly.mat';...
    '8_BaredoworldOnly.mat';...
    '9_HarrisForestOnly.mat';...
    '10_JRCForestOnly.mat';...
    '11_AvitableCarbonOnly.mat';...
    '12_GEOCarbonOnly.mat';...
    '13_ORLDAACOnly.mat';...
    '14_ESACIIOnly.mat'};
% Relative original grid sizes rounded at one decimal
Parameters.GridSizes = [9.4, 9.4, 9.4, 3.1, 562.6, 3, 9.4, 9.4, 2.8, 0.9, 9.4, 11.3, 3.1, 1];
Parameters.make_log = 0;
NewFolder = [currentFolder,'\',Parameters.InputFolder];
cd(NewFolder)
load(char(Parameters.filenames(1)));
vars = whos('-file', char(Parameters.filenames(1)));
Data =  eval(vars(1).name);
Parameters.TotalPoints = length(Data);
clear Data vars
cd(currentFolder)
%%
Parameters.Ensemble_Names = {'PCA';'CorCoef';'RegresstoMedian';'LeaveOneOut';'GridSize'};
% Parameters.testRun = 1;
Parameters.runMax = 250;
Parameters.JoinedRuns = 15;
Parameters.Nr_ensembles = length(Parameters.Ensemble_Names);
Parameters.Precision = 0.000001; % rounding factors
Parameters.Precision(2) = Parameters.Precision;%.*10;
Parameters.Precision = (1./Parameters.Precision);
Parameters.ensemble = 0;
Parameters.NrPoints = 1000000;
Parameters.runOverall = ceil(Parameters.runMax./Parameters.JoinedRuns);
Parameters.SmallerMargin = 0; % 1 for recreation to remove more zeros
Parameters.WinsorValue = 2.5;
end