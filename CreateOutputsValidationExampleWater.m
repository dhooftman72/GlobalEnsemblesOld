function CreateOutputsValidationExampleWater
warning off
load('D:\Data\Dropbox\Lactuca\Projects\Global Ensembles\Calculations\WaterData\WaterDataAccumulatedGRDC\WaterGRDCSheds.mat')
Dataset = WaterGRDCSheds;
Validation = Dataset.NormalisedValidation;
Transfer.NrModels = 8;
Transfer.OtherEnsemblesTestd = 4;
Transfer.AccWater = 'Yes'; % winsorise the data
% if 'minus' it is one sided as well as two sided
Transfer.TestedDirection = 'Minus';

% Validation data, models and Water ensembles are normalised; 
% The cross-ensembles are NOT pre-normalised!
% Everything is area corrected through taking the mean (so is per 0.001
% degree squared) or maximum divided by area (WW and Aquestat)

tester = get(Dataset);
Outfile ='WaterValidation&StatsCore';
Transfer.EnsemblesTested = 2;
Transfer.Bonferroni = 5;
Transfer.DoRegres = 1; % Doing the regressions
Transfer.BaseRun = 1; % Storing the baserun;
% Mean: ObsList = [tester.VarNames(13:20)';'MeanAmongModels';tester.VarNames(21)';tester.VarNames(22)';tester.VarNames(29)';tester.VarNames(31:33)'];
ObsList = [tester.VarNames(13:20)';'MedianAmongModels';tester.VarNames(21)';tester.VarNames(22)';tester.VarNames(34:37)'];
[Transfer] = MainbodyCodeValidation(Outfile,Dataset,Validation,Transfer,ObsList); 

display('Starting second part')
Outfile ='WaterValidationOnlyFull';
Transfer.EnsemblesTested = 7;
Transfer.Bonferroni = 11;
Transfer.DoRegres = 0; %No stats
Transfer.BaseRun = 0; %Repearting the baseruns
ObsList = [tester.VarNames(13:20)';'MedianAmongModels';tester.VarNames(21:27)';tester.VarNames(34:37)'];
[~]=MainbodyCodeValidation(Outfile,Dataset,Validation,Transfer,ObsList); 
end
