% Main body calculation for 1) comparative accuracy calculation and 2) statistics
% against wealth metrics. 1) run as bootstrap runs on 10% of the data
% (with min N = 100) and 2) as convergence loop until Sun of Squares is
% stable

function [Transfer] = MainbodyCodeValidationPublication(Outfile,Dataset,Validation,Transfer,ObsList)
if matlabpool('size') ~= 0 
       matlabpool close
end
Transfer.TotalObs = length(ObsList);
Transfer.maxBoots = 1000; % the number of accuracy bootstraps
% Regression convergence parameters
Transfer.Thresh = 0.0005; %(0.05% stats difference threshold)
Transfer.DFtarget = 178; % target number of total df's per bootstrap = Number of Recreation DFDs
Transfer.CloseTarget = 25; % how many times it needs to be equal for convergence
Transfer.MaxCounter = 10000; % fail safe against everlasting loops
Transfer.IncludeAuto = 'Yes'; %Include Spatial autocorellation (Yes/No)

OutputsList = dataset(ObsList,'ObsNames',ObsList,'Varnames','Names');
Transfer.UnValidValidation =  find(isnan(Validation)==1);
Transfer.ObsList = ObsList;
for x = 1:Transfer.NrModels
    display(char(ObsList(x)));
    Transfer.count = x;
    InVar = double(Dataset.(genvarname(char(ObsList(x)))));
    if  strcmpi(Transfer.AccWater,'Yes') == 1
        [OutVar] = WinsorFunction(InVar,2.5);
    else
        OutVar = InVar ;
    end
    testArray = [Validation,OutVar];
    [OutputsList,~,Transfer]  = CreateOutPutList(OutputsList,Dataset,testArray, Transfer);   %#ok<*ASGLU>
end
All = double(OutputsList(1:Transfer.NrModels,2:5));
OutputsList.Rho((Transfer.NrModels+1),1) = nanmedian(All(:,1));
OutputsList.PVal((Transfer.NrModels+1),1) =  NaN;
OutputsList.Rho_STD((Transfer.NrModels+1),1) =  NaN;
OutputsList.InversedDeviation((Transfer.NrModels+1),1) =nanmedian(All(:,4));
OutputsList.InversedDeviationSTD((Transfer.NrModels+1),1) =  NaN;
OutputsList.ValidDataPoints((Transfer.NrModels+1),1) =  NaN;
Transfer.IndividualRuns.Deviation.(genvarname(char(Transfer.ObsList(Transfer.NrModels+1)))) =...
    reshape((nanmean(double(Transfer.IndividualRuns.Deviation),2)),[],1);
Transfer.IndividualRuns.RHO.(genvarname(char(Transfer.ObsList(Transfer.NrModels+1)))) =...
    reshape((nanmean(double(Transfer.IndividualRuns.RHO),2)),[],1);

clear x InVar OutVar testArray
save(Outfile, 'OutputsList')

Transfer.first = 0;
Transfer.ToCrossVerification = 100;
for x = 1:Transfer.EnsemblesTested
    Transfer.count = x+(Transfer.NrModels+1);
    display(char(ObsList(Transfer.count)));
    InVar = double(Dataset.(genvarname(char(ObsList(Transfer.count)))));
    if  strcmpi(Transfer.AccWater,'Yes') == 1
        [OutVar] = WinsorFunction(InVar,2.5);
    else
        OutVar = InVar ;
    end
    testArray = [Validation,OutVar];
    [OutputsList,Regression,Transfer]  = CreateOutPutList(OutputsList,Dataset,testArray, Transfer);  %#ok<*ASGLU>
    if Transfer.count == (Transfer.NrModels+2) || Transfer.count == (Transfer.NrModels+3)
        Statistics.(genvarname(char(ObsList(Transfer.count)))) = Regression; %#ok<*STRNU>
    end
    save(Outfile, 'OutputsList','Statistics') %,'
end
clear x InVar OutVar testArray

Transfer.ToCrossVerification = Transfer.count;
for x = 1:Transfer.OtherEnsemblesTestd
    Transfer.count = x+(Transfer.ToCrossVerification);
    display(char(ObsList(Transfer.count)));
    InVar = double(Dataset.(genvarname(char(ObsList(Transfer.count))))); % Input is already area corrected
    [OutVar] = WinsorFunction(InVar,2.5);
    testArray = [Validation,OutVar];
    [OutputsList,~,Transfer]  = CreateOutPutList(OutputsList,Dataset,testArray, Transfer);   %#ok<*ASGLU>
    save(Outfile, 'OutputsList','Statistics')
end
% Note is is all pairwise, because of repeated bootstrap numbers!!
IndividualAccuracyRuns =  Transfer.IndividualRuns;
NonModelsLength = Transfer.TotalObs-Transfer.NrModels;
emptylist = zeros(NonModelsLength,1);
emptylist(:,:) = NaN;

RankingError = Transfer.RankingError;
DifferenceSignificanceMatrix.DeviationBonferoni = dataset(emptylist,'ObsNames',Transfer.ObsList((Transfer.NrModels+1):Transfer.TotalObs),'VarNames',{char(Transfer.ObsList(Transfer.NrModels+1))});
DifferenceSignificanceMatrix.RhoBonferoni = DifferenceSignificanceMatrix.DeviationBonferoni;
DifferenceSignificanceMatrix.DeviationAbs = DifferenceSignificanceMatrix.DeviationBonferoni;
DifferenceSignificanceMatrix.RhoAbs = DifferenceSignificanceMatrix.DeviationBonferoni;
DifferenceSignificanceMatrix.DeviationHoghberg = DifferenceSignificanceMatrix.DeviationBonferoni;
DifferenceSignificanceMatrix.RhoHoghberg = DifferenceSignificanceMatrix.DeviationBonferoni;

ImprovementMatrix.Deviation.Mean = dataset(emptylist,'ObsNames',Transfer.ObsList((Transfer.NrModels+1):Transfer.TotalObs),'VarNames',{char(Transfer.ObsList(Transfer.NrModels+1))});
ImprovementMatrix.Deviation.Std = ImprovementMatrix.Deviation.Mean;
ImprovementMatrix.Rho.Mean = ImprovementMatrix.Deviation.Mean;
ImprovementMatrix.Rho.Std = ImprovementMatrix.Deviation.Mean;
for x = 1: (NonModelsLength)
    for y =  1: (NonModelsLength)
        if x ~= y
            TestVar1 = double(IndividualAccuracyRuns.Deviation.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))));
            TestVar2 = double(IndividualAccuracyRuns.Deviation.(genvarname(char(Transfer.ObsList(y+Transfer.NrModels)))));
            [p, PAbs] = ttester (TestVar1,TestVar2,Transfer.Bonferroni);
            Difference.Deviation(y,1) = p;
            Difference.DeviationAbs(y,1) = PAbs;
            %save('all')
            Improvement.Deviation(y,1) = nanmean(rdivide(TestVar2,TestVar1)-1);
            ImprovementSTD.Deviation(y,1) = nanstd(rdivide(TestVar2,TestVar1)-1);
            TestVar1 = double(IndividualAccuracyRuns.RHO.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))));
            TestVar2 = double(IndividualAccuracyRuns.RHO.(genvarname(char(Transfer.ObsList(y+Transfer.NrModels)))));
            [p, PAbs] = ttester (TestVar1,TestVar2,Transfer.Bonferroni);
            Difference.Rho(y,1) = p;
            Difference.RhoAbs(y,1) = PAbs;
            Improvement.Rho(y,1) =  nanmean(rdivide(TestVar2,TestVar1)-1);
            ImprovementSTD.Rho(y,1) = nanstd(rdivide(TestVar2,TestVar1)-1);
        else
            Difference.Deviation(y,1) = NaN;
            Difference.DeviationAbs(y,1) = NaN;
            Improvement.Deviation(y,1) = NaN;
            Difference.Rho(y,1) = NaN;
            Difference.RhoAbs(y,1) = NaN;
            Improvement.Rho(y,1) = NaN;
        end
    end
    
    DifferenceSignificanceMatrix.DeviationBonferoni.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = Difference.Deviation;
    DifferenceSignificanceMatrix.DeviationAbs.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = Difference.DeviationAbs;
    StarSignificant = HoghbergFunc(Difference.DeviationAbs);
    DifferenceSignificanceMatrix.DeviationHoghberg.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = StarSignificant;
    
    DifferenceSignificanceMatrix.RhoBonferoni.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = Difference.Rho;
    DifferenceSignificanceMatrix.RhoAbs.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = Difference.RhoAbs;
    StarSignificant = HoghbergFunc(Difference.RhoAbs);
    DifferenceSignificanceMatrix.RhoHoghberg.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = StarSignificant;
    
    ImprovementMatrix.Deviation.Mean.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = Improvement.Deviation;
    ImprovementMatrix.Deviation.Std.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = ImprovementSTD.Deviation;
    ImprovementMatrix.Rho.Mean.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = Improvement.Rho;
    ImprovementMatrix.Rho.Std.(genvarname(char(Transfer.ObsList(x+Transfer.NrModels)))) = ImprovementSTD.Rho;
end
if matlabpool('size') ~= 0
    matlabpool close
end
save(Outfile, 'OutputsList', 'DifferenceSignificanceMatrix','IndividualAccuracyRuns','ImprovementMatrix', 'Statistics','RankingError','Transfer')
display(OutputsList)
end

function  [p, PAbs] = ttester (TestVar1,TestVar2,Bonferroni)
[h,p] = ttest2(TestVar1,TestVar2); % ttest would provide a paired sample. The draws are paired already!
PAbs = p;
p = p * Bonferroni;
if p < 0.000001
    p = 0;
end
if p >0.05
    p = -9;
end
end

function StarSignificant = HoghbergFunc(VarIn)
Nr = find((isnan(VarIn))==1);
[TrueOuts(:,1)]=Hoghberg(VarIn,0.05);
[TrueOuts(:,2)]=Hoghberg(VarIn,0.01);
[TrueOuts(:,3)]=Hoghberg(VarIn,0.001);
StarSignificant = sum(TrueOuts,2);
StarSignificant(Nr) = NaN;
end

function [OutputsList,Regression,Transfer]  = CreateOutPutList(OutputsList,Variables,testArray, Transfer) %#ok<*STOUT>
[Outputs,Transfer] = Accuracy_statistics_Validation(testArray,Transfer);
OutputsList.Rho(Transfer.count,1) = Outputs.RHO;
OutputsList.Rho_STD(Transfer.count,1) = Outputs.RHOSTD;
OutputsList.PVal(Transfer.count,1) =  Outputs.PVAL;
if OutputsList.PVal(Transfer.count,1)  < 0.000001
    OutputsList.PVal(Transfer.count,1)  = 0;
end
OutputsList.InversedDeviation(Transfer.count,1) = Outputs.mean_double_deviation;
OutputsList.InversedDeviationSTD(Transfer.count,1) = Outputs.mean_double_deviationSTD;
OutputsList.ValidDataPoints(Transfer.count,1) = Outputs.Datapoints;

if  Transfer.count == 1
    Transfer.IndividualRuns.Deviation = dataset([1;1],'VarNames',(char(Transfer.ObsList(Transfer.count))));
    Transfer.IndividualRuns.RHO = dataset([1;1],'VarNames',(char(Transfer.ObsList(Transfer.count))));
end
Transfer.IndividualRuns.Deviation.(genvarname(char(Transfer.ObsList(Transfer.count)))) =  reshape(Outputs.AllMeanDeviation,[],1);
Transfer.IndividualRuns.RHO.(genvarname(char(Transfer.ObsList(Transfer.count)))) =  reshape(Outputs.AllRho,[],1);
if Transfer.count == (Transfer.NrModels+2) ||  Transfer.count == (Transfer.NrModels+3)
    Transfer.RankingError.Rho.(genvarname(char(Transfer.ObsList(Transfer.count)))) =  Outputs.Rho_point;
    Transfer.RankingError.Deviation.(genvarname(char(Transfer.ObsList(Transfer.count)))) =  Outputs.deviation_point;
end
Regression = NaN;
if  (Transfer.count == (Transfer.NrModels+2) || Transfer.count == (Transfer.NrModels+3)) && Transfer.DoRegres == 1 
    %(Transfer.count > (Transfer.NrModels+1) && Transfer.count <= Transfer.ToCrossVerification)
    Size = length(Outputs.deviation_point);
    btsfact = 5; %5 degrees maximum interaction, else autocorrelation coefficient = 0
    if Transfer.first == 0;
        [wij] = Morans(Variables.XCoor,Variables.YCoor,5,(Size.*btsfact));
        save('wij.mat','wij');
        Transfer.first = 1;
    end
    if isempty(strmatch('PlotTypesDummy',get(Variables,'VarNames')))== 1
        Variables.PlotTypesDummy(:,1) = 1;
    end
    
    for j = 1:2
        warning off
        if j == 1
            VarIn = Outputs.deviation_point;
        elseif j == 2
            VarIn = Outputs.Rho_point;
        elseif j == 3
            VarIn = Variables.SEMModels;
        end
        Type = j;
        AutoCorrelationCoef = AutoCorFunc(VarIn,Size);
        [AnovaTable] = InteractionModel(VarIn,Variables,AutoCorrelationCoef,Transfer,Type);
        if j == 1
            Regression.DeviationBespoke = AnovaTable;
        elseif j == 2 %&& (Transfer.count == (Transfer.NrModels+2))
            Regression.RhoPointBespoke = AnovaTable;
        elseif j == 3 %&& (Transfer.count == (Transfer.NrModels+2))
            Regression.SEMBespoke = AnovaTable;
        end
    end
end
end

function AutoCorrelationCoef = AutoCorFunc(InVariation,Size)
wij = load('wij.mat');
parfor i = 1:Size
    Autot = 0;
    for j = 1:Size
        if i ~= j
            if isnan(InVariation(j))~= 1
                Autot = Autot + (wij.wij(i,j).*InVariation(j));  %#ok<*NODEF,*PFIIN>
            end
        end
    end
    if Autot > 0
        AutoCorrelationCoef(i,1) = Autot./nansum(wij.wij(i,:)); %#ok<*PFBNS,*AGROW>
    else
        AutoCorrelationCoef(i,1) = 0;
    end
end
end


function [Outputs,Transfer] = Accuracy_statistics_Validation(testArray,Transfer)
% clean data set to determine true N and where top put NaN;
Precision = 1./(0.00001);
testArray(isinf(testArray)==1) = NaN;
testArrayOrg = testArray;
Outputs.DatapointsOrg = size(testArray,1);
[orginal_order,NaNPoints,testArray]  = CleanOutNaN(testArray);
Outputs.Datapoints = size(testArray,1);
BootTarget = max(100,(ceil(Outputs.Datapoints./10)));

%% Make log if applicable, but not for Ensembles
testVar = testArray;
%% Correlation stats: Spearman: Overall Accuracy
isnnx = length(find(isnan(testVar(:,1)))==1);
isnny = length(find(isnan(testVar(:,2)))==1);
if (exist('testVar','var') == 1) && (isempty(testVar)~= 1) &&  (isnnx ~= length(testVar(:,1))) &&  (isnny ~= length(testVar(:,2)))
    for boots = 1:(Transfer.maxBoots+1)
        if boots == 1 % first create the full list
            testSet = testVar;
            %Datapoints = Outputs.Datapoints;
            OutPoints = CorrFunc(testSet, Precision, boots);
            % %% recreate deviation and rho in orginal order
            Outputs.deviation_point(orginal_order,1) = OutPoints.deviation_point;
            Outputs.deviation_point(NaNPoints,1) = NaN;
            Outputs.Rho_point(orginal_order,1) = OutPoints.Rho_point;
            Outputs.Rho_point(NaNPoints,1) = NaN;
            testVar = testArrayOrg; % so this include the NaN list
        else
            if Transfer.count == 1 && Transfer.BaseRun == 1
                tmper = randperm(Outputs.DatapointsOrg);
                tmperList = tmper(1:BootTarget);
                Transfer.AllbootsList((boots-1),:) = tmperList;
            else % so note it is repeating so pairwise comparisons; 
                % as well as repeated for both short and long run by
                % storing AllBootsList
                tmperList = Transfer.AllbootsList((boots-1),:);
                if boots == 2
                display('Repeated Bootstrap')
                end
            end
            testSet = testVar(tmperList,:);
            OutPoints = CorrFunc(testSet, Precision, boots);
            Out.mean_double_deviation(boots-1) = OutPoints.mean_double_deviation;
            Out.RHO(boots-1) = OutPoints.RHO;
            Out.PVAL(boots-1) = single(OutPoints.PVAL);
        end
        
    end
    Outputs.mean_double_deviation = nanmean(Out.mean_double_deviation);
    Outputs.AllMeanDeviation = Out.mean_double_deviation;
    Outputs.RHO = nanmean(Out.RHO);
    Outputs.AllRho = Out.RHO;
    Outputs.PVAL = nanmean(Out.PVAL);
    Outputs.mean_double_deviationSTD = nanstd(Out.mean_double_deviation);
    Outputs.RHOSTD = nanstd(Out.RHO);
else
    Outsputs.RHO = -9999;
end
    function  Outs = CorrFunc(testSet, Precision, Boots)
        if Boots  == 1
            RankVali = tiedrank(testSet(:,1));
            RankY = tiedrank(testSet(:,2));
            Dsingle(:,1) = (abs((RankVali-RankY)))/(length(testSet(:,1))./3);
            % 1) the expected difference under a random is 1/3 times length
            % (I have tested that with 1,000,000 random draws)
            % 2) If the ranking difference increases (so less accuracy) the value
            % of Dsingle increases, so ranking error is good (low) to
            % bad(high) and normalised (below)
            Outs.Rho_point(:,1) = Dsingle(:,1)./max(Dsingle(:,1));
            %% Inverse deviance against a 1:1 line
            clear x_range y_range
            x_range = testSet(:,1);
            y_range = testSet(:,2);
            Outs.deviation_point= (abs(y_range-x_range)); %% Accuracy per point
        else
            c1=find((isnan(testSet(:,1))==1));
            d=find((isnan(testSet(:,2))==1));
            alL = [c1;d];
            A(:,1) = unique(alL);
            testSet(A,:) = [];  %#ok<*FNDSB>
            Datapoint = size(testSet,1);
            [Outs.RHO,Outs.PVAL] = corr(testSet(:,1),testSet(:,2),'type','Spearman');
            Outs.PVAL = single(Outs.PVAL);
            Outs.RHO = (round(Outs.RHO.*Precision))./Precision;
            % Inverse deviance against a 1:1 line
            clear x_range y_range
            x_range = testSet(:,1);
            y_range = testSet(:,2);
            Outs.deviation_point= abs(y_range-x_range); %% Accuracy per point
            Outs.mean_double_deviation = 1- (nansum(Outs.deviation_point)/Datapoint); %%Accuracy overall
            Outs.mean_double_deviation = (round( Outs.mean_double_deviation.*Precision))./Precision;
        end
    end
end

function [Anova_Table] = InteractionModel(VarIn,Variables,AutoCorrelationCoef,Transfer,Type)
clc
% Interaction convergence model with Type 1 splits
ValidList = 1:1:length(VarIn(:,1));
ValidList = ValidList';
ValidEnsem =  find(isnan(VarIn)==1);
all = [Transfer.UnValidValidation;ValidEnsem];
UnValid(:,1) = unique(all);
ValidList(UnValid,:) = [];

%Clean out non valid Validation points
VarInValid = VarIn(ValidList);
ValidLength = length(VarInValid);
VariablesValid = Variables(ValidList,:);
AutoValid = AutoCorrelationCoef(ValidList);
if ValidLength <= Transfer.DFtarget
    Transfer.DFtarget = ValidLength; % just in case there are too many NaN;
end
% Anova Table = dataset
[Anova_Table,Coeffs,Rsquared] = runRegressMany(VarInValid,VariablesValid,AutoValid,Transfer,ValidLength,Type);

%F-values
for para = 1:12
    F_Value(para,1) = (Anova_Table.MeanSq(para))./(Anova_Table.MeanSq(13));
    P_Value(para,1) =  1- fcdf(F_Value(para),cell2mat(Anova_Table.DFs(para)),cell2mat(Anova_Table.DFs(13)));
end
Anova_Table.F = F_Value;
Anova_Table.Pvalue = P_Value;
Anova_Table.Direction(13,1) = 0;
Anova_Table.RSquared(13,1) = 0;
Anova_Table.NrStarsSigni(:,1) = NaN;
clc
if strcmpi(Transfer.TestedDirection,'Minus') == 1
   for para = 3:11
        fprintf('Recalulating F values one-sided for parameter nr: %i for %s \n', (para-2), (char(Transfer.ObsList(Transfer.count))))
        OldF = Anova_Table.F(para,1);
       if Coeffs(para-2) > 0
           Anova_Table.Pvalue(para,1) = 1-P_Value(para,1);
           reach = 0;
           Iter = 0;
           while reach == 0
               TestedF = F_Value(para,1) - Iter;
               Tester =  1- fcdf(TestedF,cell2mat(Anova_Table.DFs(para)),cell2mat(Anova_Table.DFs(12)));
               if Tester >=  Anova_Table.Pvalue(para,1) ||  TestedF <= 0
                   reach = 1;
                   if TestedF <= 0
                       TestedF = 0;
                   end
                   Anova_Table.F(para,1) = TestedF;
               else
                   Iter = Iter + 0.0001;
               end
           end  
       end
       fprintf('Old F value %1.4f with P %1.4f becomes new F value %1.4f with P %1.4f  \n',  OldF,P_Value(para,1), Anova_Table.F(para,1), Anova_Table.Pvalue(para,1))
   end
end
%Hoghberg correction for significane
StarSignificant = HoghbergFunc(Anova_Table.Pvalue(3:12));
Anova_Table.NrStarsSigni(3:12,1) = StarSignificant;
for NrVariable = 1:7
    %if Anova_Table.Pvalue((NrVariable+2),1) <= 0.1
        Anova_Table.Direction((NrVariable+2),1) = Coeffs(NrVariable);
        Anova_Table.RSquared((NrVariable+2),1) = Rsquared(NrVariable);
    %end       
end
Anova_Table.Direction(10,1) = Coeffs(8); 
Anova_Table.Direction(11,1) = Coeffs(9); 
Anova_Table.Direction(11,1) = Coeffs(10); 

% Two sided Statistics
Anova_Table.TwoSideF = F_Value;
Anova_Table.TwoSideP = P_Value;
Anova_Table.TwoSideNrStarsSigni(:,1) = NaN;
StarSignificant = HoghbergFunc(Anova_Table.TwoSideP(3:12));
Anova_Table.TwoSideNrStarsSigni(3:12,1) = StarSignificant;
%Anova_Table.Direction(Anova_Table.Direction == 0) = NaN;
end

%
function   [Anova_Table,Coeffs,RSquared] = runRegressMany(VarInValid,VariablesValid,AutoTransIn,Transfer,ValidLength,Type)
for NrVariable = 1:8
    clc
    %Iniate parameters
    ConVer = 0;
    Counter = 1;
    CloseCount = 0;
    SumSquares = zeros(3,1);
    if NrVariable == 8
    SumSquares = zeros(5,1); 
    end
    SSError = 0;
    DFError = 0;
    RsquaredSum = 0;
    CoefList =  zeros(1,1);
    DifAverage = 1;
    while ConVer == 0
        % pick random selection
        tmper = randperm(ValidLength);
        tmperList = sort(tmper(1:Transfer.DFtarget));
        VarInTrans = VarInValid(tmperList);
        VariablesTrans = VariablesValid(tmperList,:);
        AutoTrans = AutoTransIn(tmperList);
        %Data & arcsin corrections 
        Variable(:,1)  = log10(VariablesTrans.IncomePerHead+1);
        Variable(:,2) = asin(sqrt(double(VariablesTrans.GiniCoef)));
        Variable(:,3) = asin(sqrt(double(VariablesTrans.UrbanPropPeople)));
        Variable(:,4) = double(VariablesTrans.HDI);
        Variable(:,5) =  asin(sqrt((VariablesTrans.PeopleinR_D)./100000));      
        Variable(:,6) = asin(sqrt((double(VariablesTrans.PropR_D)./100)));  
        Variable(:,7) = asin(sqrt((double(VariablesTrans.GiniEquality)./100)));
        NameList = [{'PlotType (Dummy)'};{'AutoCorrelation'};{'IncomePerhead'};{'LandUseGini'};{'UrbanProp'};{'HumanDevelopIndex'};{'PeopleinR&D'};{'ProportionR&D'};...
            {'IncomeEquality'};{'IncomePerhead (x Equality)'};{'IncomeEquality (x Income)'};{'Interaction'};{'Error'}];
        ATable = dataset({'Dummy'},'Varnames',char('Source'));
        
        if strcmpi(Transfer.IncludeAuto,'Yes') == 1
            [~,outstmp,statstmp] = anovan(VarInTrans,{VariablesTrans.PlotTypesDummy, AutoTrans},'sstype',1,...
                'model',[1 0 ; 0 1],'continuous', [2],'random', [1], 'display', 'off','varnames', {'PlotType (Dummy)','AutoCorrelation'});
            ATable.Source(1:2,1) = outstmp(2:3,1);
            ATable.SumSq(1:2,1) = outstmp(2:3,2);
            ATable.DFs(1:2,1) = outstmp(2:3,3);
            SecondTierVar = statstmp.resid';
        else
            ATable.Source(1:2,1) =  {'PlotType (Dummy)';'AutoCorrelation'};
            ATable.SumSq(1:2,1) = {0};
            ATable.DFs(1:2,1) = {0};
            SecondTierVar = VarInTrans;
        end
        for i = NrVariable;
            if NrVariable <= 7
            [~,outs,stats] = anovan(SecondTierVar,{Variable(:,i)},'sstype',1,...
                'model',[1],'continuous', [1], 'display', 'off',...
                'varnames', {'Variable'});
            ATable.Source((3:4),1) =  {NameList(i+2);'Error'};
            ATable.SumSq((3:4),1) = outs(2:3,2);
            ATable.DFs((3:4),1) = outs((2:3),3);
            coeffs = stats.coeffs(2);
            else
             [~,outs,stats] = anovan(SecondTierVar,{Variable(:,1),Variable(:,7)},'sstype',3,...
                'model',[1 0 ; 0 1; 1 1],'continuous', [1 2], 'display', 'off','varnames', {'Income','GiniCoef'});
            ATable.Source((3:6),1) =  {'Income','GiniCoef','Interaction','Error'};
            ATable.SumSq((3:6),1) = outs(2:5,2);
            ATable.DFs((3:6),1) = outs((2:5),3);
            coeffs = stats.coeffs(2:4);
            end
        end
        %save('all')
        Points_Corrected =  reshape(SecondTierVar,[],1) - reshape(stats.resid,[],1);    
        Rsquared = RsquaredFunc(Points_Corrected,reshape(SecondTierVar,[],1));
        clear Points_Corrected
        % test for proportional differences compared to all previous runs
        if NrVariable <= 7
            SumSquaresTest =  SumSquares./(Counter-1);
  
            Difs = zeros(3,1);
            SumSquaresNew =  (SumSquares + cell2mat(ATable.SumSq(1:3)))./(Counter);
            for j = 1:3
                if isinf(SumSquaresTest(j))~= 1
                    Dif(j) = abs(SumSquaresTest(j)-SumSquaresNew(j));
                    if Dif(j) > 0
                        Dif(j) = Dif(j)./SumSquaresTest(j);
                    end
                    if Dif(j) < Transfer.Thresh && isinf(Dif(j))~= 1
                        Difs(j) = 1;
                    end
                end
                DifAverage = max(Dif);
                clear Dif
            end
            if sum(Difs) ~= 3|| ((Transfer.DFtarget*Counter) < ValidLength*10) || Counter < 10
                SumSquares =  SumSquares + cell2mat(ATable.SumSq(1:3));          
                RsquaredSum =  RsquaredSum + cell2mat(Rsquared);
                CoefList = CoefList + coeffs;
                SSError = SSError + cell2mat(ATable.SumSq(4));
                DFError = DFError + cell2mat(ATable.DFs(4));
                Counter = Counter + 1;
                CloseCount = 0;
                if (Counter/100 == floor(Counter/100)) || Counter == 10
                    fprintf('Convergence # run within %s type %i Variable %i: %i (Max. convergence level left: %2.2f percent) \n', (char(Transfer.ObsList(Transfer.count))), Type, NrVariable, Counter,(DifAverage.*100))
                end
            elseif  sum(Difs) == 3
                CloseCount = CloseCount + 1;
                if CloseCount == Transfer.CloseTarget
                    ConVer = 1;
                end
            else
                display('Fatal error interaction Model')
                break
            end
            if Counter >= Transfer.MaxCounter; % Fail safe counter
                ConVer = 1;
            end
        else 
            SumSquaresTest =  SumSquares./(Counter-1);
            Difs = zeros(5,1);
            SumSquaresNew =  (SumSquares + cell2mat(ATable.SumSq(1:5)))./(Counter);
            for j = 1:5
                if isinf(SumSquaresTest(j))~= 1
                    Dif(j) = abs(SumSquaresTest(j)-SumSquaresNew(j));
                    if Dif(j) > 0
                        Dif(j) = Dif(j)./SumSquaresTest(j);
                    end
                    if Dif(j) < Transfer.Thresh && isinf(Dif(j))~= 1
                        Difs(j) = 1;
                    end
                end
                DifAverage = max(Dif);
            end
            if sum(Difs) ~= 5|| ((Transfer.DFtarget*Counter) < ValidLength*10) || Counter < 10
                SumSquares =  SumSquares + cell2mat(ATable.SumSq(1:5));
                RsquaredSum =  RsquaredSum + cell2mat(Rsquared);
                CoefList = CoefList + coeffs;
                SSError = SSError + cell2mat(ATable.SumSq(5));
                DFError = DFError + cell2mat(ATable.DFs(5));
                Counter = Counter + 1;
                CloseCount = 0;
                if (Counter/100 == floor(Counter/100)) || Counter == 10
                    fprintf('Convergence # run within %s type %i Variable %i: %i (Max. convergence level left: %2.2f percent) \n', (char(Transfer.ObsList(Transfer.count))), Type, NrVariable, Counter,(DifAverage.*100))
                end
            elseif  sum(Difs) == 5
                CloseCount = CloseCount + 1;
                if CloseCount == Transfer.CloseTarget
                    ConVer = 1;
                end
            else
                display('Fatal error interaction Model')
                break
            end
            if Counter >= Transfer.MaxCounter; % Fail safe counter
                ConVer = 1;
            end
        end
    end
    % Create AnovaTable
    if NrVariable == 1
        Anova_Table = dataset(NameList,'Varnames',char('Source'));
        Anova_Table.SumSq(1:2,1) = SumSquaresTest(1:2,1);
        Anova_Table.DFs(1:2,1)= ATable.DFs(1:2,1);
    end
    if NrVariable <= 7
    Anova_Table.SumSq(NrVariable+2,1) = SumSquaresTest(3);
    Anova_Table.DFs(NrVariable+2,1)= ATable.DFs(3);
    Coeffs(NrVariable) = CoefList./(Counter-1);
    RSquared(NrVariable) = RsquaredSum./(Counter-1);
    if NrVariable == 7
        Anova_Table.SumSq(13) = SSError./(Counter-1);
        Anova_Table.DFs(13) = mat2cell(ceil(DFError./(Counter-1)));
    end
    elseif NrVariable == 8
    Anova_Table.SumSq(NrVariable+2,1) = SumSquaresTest(3);
    Anova_Table.SumSq(NrVariable+3,1) = SumSquaresTest(4);
    Anova_Table.SumSq(NrVariable+4,1) = SumSquaresTest(5);
    Anova_Table.DFs(NrVariable+2,1)= ATable.DFs(3);
    Anova_Table.DFs(NrVariable+3,1)= ATable.DFs(4);
    Anova_Table.DFs(NrVariable+4,1)= ATable.DFs(5);
    Anova_Table.MeanSq(1:13,1) = Anova_Table.SumSq(1:13,1)./ cell2mat(Anova_Table.DFs(1:13,1));
    Coeffs(NrVariable) = CoefList(1)./(Counter-1);
    Coeffs(NrVariable+1) = CoefList(2)./(Counter-1);
    Coeffs(NrVariable+2) = CoefList(3)./(Counter-1);
    Coeffs(NrVariable+1) = CoefList(2)./(Counter-1);
    RSquared(NrVariable) = RsquaredSum./(Counter-1);
    RSquared(NrVariable+1) = RsquaredSum./(Counter-1);
    RSquared(NrVariable+2) = RsquaredSum./(Counter-1);
    end
end
end


function [wij] = Morans(X,Y,MaxD,bts) %#ok<INUSD>
if length(X) ~= length(Y)
    display('Fatal Error, unequal grid size')
    return
end
Size = length(X);
if matlabpool('size') == 0 && Size > 1000
    matlabpool open 8
end
parfor i = 1:Size
    for j = 1:Size
        if i == j
            wij(i,j) = 0;  %#ok<*AGROW>
        else
            Dist =  log10((sqrt ((((X(i)-X(j))^2) + (Y(i)-Y(j))^2)))+1);
            wijtmp = (log10(MaxD)-Dist)./log10(MaxD);
            wijtmp(wijtmp<0) = 0;
            wij(i,j) = wijtmp;
        end
    end
end
if matlabpool('size') ~= 0 
       matlabpool close
end
end

function [OutVarOrg,OutVargMis,ArrayOut] = CleanOutNaN(ArrayIn)
OutVarOrg = 1:1:length(ArrayIn(:,1));
OutVarOrg = OutVarOrg';
a1=find((isnan(ArrayIn(:,1))==1));
b=find((isnan(ArrayIn(:,2))==1));
all = [a1;b];
a(:,1) = unique(all);
ArrayOut = ArrayIn;
ArrayOut(a,:) = [];  
OutVarOrg(a,:) = [];
OutVargMis = a;
end

function  [OutVar] = WinsorFunction(InVar,percLow)
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
prct(2) = prctile(InVar_norm,(100-percLow));
if exist('InVar_org') == 1; %#ok<EXIST>
    InVar_norm = InVar_org;
end
OutVar = (InVar_norm./prct(2));
clear InVar_norm InVar_org testlist testtmp upboud
OutVar(OutVar>1) = 1;
end

function Rsquared = RsquaredFunc(ExpecIn,ObservedIn)
% ExpecIn = ExpecIn - (min(ExpecIn));
% ObservedIn =ObservedIn - (min(ObservedIn));
List = find((isnan(ExpecIn)== 1));
List2 = find((isnan(ObservedIn)== 1));
List = [(reshape(List,[],1));(reshape(List2,[],1))];
ExpecIn(List) = [];
ObservedIn(List) = [];
meanY = mean(ObservedIn);
Size = length(ExpecIn);
for t = 1:1:Size
    ssres(t) = ((ExpecIn(t)-ObservedIn(t)).^2);
    sstot(t) =  ((ObservedIn(t)-meanY).^2);
end
Rsquared= {1- ((sum(ssres)) /(sum(sstot)))};
end

%% To this work tailored and shortend version of 
% https://nl.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh
function [TrueOuts]=Hoghberg(pvals,alpha)
if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1)),
        error('Some p-values are greater than 1.');
    end
end
if nargin<2,
    alpha=.05;
end
s=size(pvals);
if (length(s)>2) || s(1)>1,
    p_sorted=sort(reshape(pvals,[],1));
else
    %p-values are already a row vector
    p_sorted=sort(pvals);
end
m=length(p_sorted); %number of tests
    %BH procedure for independence or positive dependence
    % see also https://www.statisticshowto.com/benjamini-hochberg-procedure/
    baseL=(1:m)*alpha/m;
rej=p_sorted<=baseL;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id),
    TrueOuts=pvals*0;
else
    crit_p=p_sorted(max_id);
    TrueOuts=pvals<=crit_p;
end
TrueOuts = reshape(TrueOuts,[],1);
end
