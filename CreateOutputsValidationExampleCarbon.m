function CreateOutputsValidationExampleCarbon 
warning off
load('D:\Data\Dropbox\Lactuca\Projects\Global Ensembles\Calculations\CarbonData\PlotsCarbon.mat');
Dataset = PlotsCarbon;
Validation = Dataset.NormalisedValidation;
Transfer.NrModels = 14;
Transfer.OtherEnsemblesTestd = 4;
Transfer.AccWater = 'No';
% if 'minus' it is one sided as well as two sided
Transfer.TestedDirection = 'Minus';

% All data are prenormalised
% Everything is area corrected through taking the mean (so is per 0.001
% degree squared)

tester = get(Dataset);

Outfile ='CarbonValidation&StatsCore';
Transfer.EnsemblesTested = 2;
Transfer.Bonferroni = 5;
Transfer.DoRegres = 1; % Doing the regressions
Transfer.BaseRun = 1; % Storing the baserun;
% Mean: ObsList = [tester.VarNames(15:28)';'MeanAmongModels';tester.VarNames(29)';tester.VarNames(30)';tester.VarNames(37)';tester.VarNames(39:41)'];
ObsList = [tester.VarNames(15:28)';'MedianAmongModels';tester.VarNames(29)';tester.VarNames(30)';tester.VarNames(42:45)'];
[Transfer] = MainbodyCodeValidation(Outfile,Dataset,Validation,Transfer,ObsList); 

display('Starting second part')
Outfile ='CarbonValidationOnlyFull';
Transfer.EnsemblesTested = 7;
Transfer.Bonferroni = 11;
Transfer.DoRegres = 0; %No stats
Transfer.BaseRun = 0; %Repearting the baseruns
ObsList = [tester.VarNames(15:28)';'MedianAmongModels';tester.VarNames(29:35)';tester.VarNames(42:45)'];
[~] = MainbodyCodeValidation(Outfile,Dataset,Validation,Transfer,ObsList); 
end


