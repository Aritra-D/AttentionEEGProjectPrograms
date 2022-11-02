function displayPopulationLDA_FlickerStimuli(regFlag,transformType,numFolds,useEqualStimRepsFlag)

close all; % closes any open Figures

% Set FolderSourceString automatically; add files & folders to MATLAB path
fileNameWithPath = mfilename('fullpath');
[filePath,~,~] = fileparts(fileNameWithPath);
cd(filePath); cd ..
addpath(genpath(pwd));
folderSourceString = fullfile(pwd,'analyzedData');

if ~exist('gridType','var');            gridType = 'EEG';                                       end
if ~exist('protocolType','var');        protocolType = 'SRC-Long';                              end
if ~exist('badTrialStr','var');         badTrialStr = 'v10';                                    end

% regFlag = 2;
% transformType = 3;
% numFolds = 2;
% useEqualStimRepsFlag = 1;

tapers = [1 1];
fileName = fullfile(folderSourceString,['powerData_',protocolType,'_',gridType,'_allSubjects_N26_tapers_' num2str(tapers(2)) '_badTrial_' badTrialStr '.mat']);
if exist(fileName,'file')
    load(fileName,'powerData_AllElecs','powerData_GroupWise');
else
    error('Datafile Not Available! Check Data Path!')
end

[~, data_ElecGroupWise] = getDataForLDA(1:26,powerData_AllElecs,powerData_GroupWise);

dPrimeVals = getProjection_v2(data_ElecGroupWise,regFlag,transformType,numFolds,useEqualStimRepsFlag);
nanFlag = 'omitnan'; % NaN Values are processed using this flag during averaging

% Display options
hFig = figure;
set(hFig,'units','normalized','outerPosition',[0 0 1 1]);
hPlot = getPlotHandles(2,3,[0.1 0.1, 0.8 0.8],0.1,0.14,0); linkaxes(hPlot)
fontSize = 14;
tickLength = get(hPlot(1,1),'TickLength');

colors = {'k','r','c','m'};
neuralMeasureLabels = {'alpha','gamma','SSVEP','all'};
elecGroups = {'Parieto-Occipital','Frontal','Centro-Parietal','Fronto-Central','Temporal','all'};

for i= 1:6
    clear data mBars eBars
    data = squeeze(dPrimeVals(:,i,:));
    for iBar = 1:size(dPrimeVals,3)
        mBars(iBar) = mean(data(:,iBar),1,nanFlag); %#ok<*SAGROW>
        eBars(iBar) = std(data(:,iBar),[],1,nanFlag)./sqrt(size(data,1));
        if i<=3
            subplot(hPlot(1,i)); hold(hPlot(1,i),'on');
            barPlot = bar(iBar,mBars(iBar));
            barPlot.FaceColor = colors{iBar};
        else
            subplot(hPlot(2,i-3)); hold(hPlot(2,i-3),'on');
            barPlot = bar(iBar,mBars(iBar));
            barPlot.FaceColor = colors{iBar};
        end
    end
    if i<=3
        errorbar(hPlot(1,i),1:length(mBars),mBars,eBars,'.','color','k');
        set(hPlot(1,i),'xTick',1:4,'xTickLabel',neuralMeasureLabels,'XTickLabelRotation',30,'yTick',-0.3:0.1:0.3,'yTickLabel',-0.3:0.1:0.3,'fontSize',fontSize,'TickDir','out','TickLength',2*tickLength)
        title(hPlot(1,i),elecGroups{i})
    else
        errorbar(hPlot(2,i-3),1:length(mBars),mBars,eBars,'.','color','k');
        set(hPlot(2,i-3),'xTick',1:4,'xTickLabel',neuralMeasureLabels,'XTickLabelRotation',30,'yTick',-0.3:0.1:0.3,'yTickLabel',-0.3:0.1:0.3,'fontSize',fontSize,'TickDir','out','TickLength',2*tickLength)
        title(hPlot(2,i-3),elecGroups{i})
    end
end

xlim(hPlot(1,1),[0 length(neuralMeasureLabels)+1]);
ylim(hPlot(1,1),[-0.3 0.3]);
ylabel(hPlot(1,1),{'Neural d'' (Att-Ign)'})
ylabel(hPlot(2,1),{'Neural d'' (Att-Ign)'})
end


% Accessory Functions
function [data_allElecs, data_ElecGroupWise] = getDataForLDA(subjectIdx,powerData_AllElecs,powerData_GroupWise) % Process power data for LDA
neuralMeasures = {'alpha','gamma','SSVEP'};
elecGroups = {'PO','F','CP','FC','T'};
TFs =[1 2]; % 1- Hits 2 - Miss
attLoc = [1 2];  % 1- Right 2 - Left

for iSub = 1:length(subjectIdx)
    powerDataTG_allElecs = powerData_AllElecs{subjectIdx(iSub)}.powerValsTG;
    powerDataTG_GroupWise = powerData_GroupWise{subjectIdx(iSub)}.powerValsTG;
    for iElecGroup = 1: length(elecGroups)
        elecs = getElectrodeGroupInformation(elecGroups{iElecGroup});
        for iNeuralMeasure = 1:length(neuralMeasures)
            for iCount = 1:length(attLoc)*length(TFs)
                switch iCount
                    % here elec_side refers to the row of the cell array of powerData_Groupwise matrix - 1st row is for all left elecs and 2nd row for right elecs
                    case 1; att_Condition = 3; ign_Condition = 6;  elec_side = 2; elecList = elecs{elec_side}; SSVEPPos = 3; % att_condition - HL 12Hz ign_condition - HR 16Hz
                    case 2; att_Condition = 4; ign_Condition = 5;  elec_side = 1; elecList = elecs{elec_side}; SSVEPPos = 3; % att_condition - HR 12Hz ign_condition - HL 16Hz
                    case 3; att_Condition = 5; ign_Condition = 4;  elec_side = 2; elecList = elecs{elec_side}; SSVEPPos = 4; % att_condition - HL 16Hz ign_condition - HR 12Hz
                    case 4; att_Condition = 6; ign_Condition = 3;  elec_side = 1; elecList = elecs{elec_side}; SSVEPPos = 4; % att_condition - HR 16Hz ign_condition - HL 12Hz
                end
                
                % Processing the data for alpha and gamma power
                if strcmp(neuralMeasures{iNeuralMeasure},'alpha')||strcmp(neuralMeasures{iNeuralMeasure},'gamma')
                    % Processing relevant electrode data from all 64 EEG electrodes
                    hits_AttData_allElecs{iSub,iElecGroup}{iCount,iNeuralMeasure} = powerDataTG_allElecs{att_Condition}(elecList,:,iNeuralMeasure)'; %#ok<*AGROW>
                    hits_IgnData_allElecs{iSub,iElecGroup}{iCount,iNeuralMeasure} = powerDataTG_allElecs{ign_Condition}(elecList,:,iNeuralMeasure)';
                    % Processing Group-wise electrode (averaged across electrodes in an electrode-group) data
                    hits_AttData_GroupWise{iSub}{iCount,iNeuralMeasure}(iElecGroup,:) = powerDataTG_GroupWise{elec_side,att_Condition}(iElecGroup,:,iNeuralMeasure)';
                    hits_IgnData_GroupWise{iSub}{iCount,iNeuralMeasure}(iElecGroup,:) = powerDataTG_GroupWise{elec_side,ign_Condition}(iElecGroup,:,iNeuralMeasure)';
                    
                    % Processing the data for SSVEP power since SSVEP freq should be selected based on the att/ign location for the set of contralateral electrodes
                elseif strcmp(neuralMeasures{iNeuralMeasure},'SSVEP')
                    hits_AttData_allElecs{iSub,iElecGroup}{iCount,iNeuralMeasure} = powerDataTG_allElecs{att_Condition}(elecList,:,SSVEPPos)';
                    hits_IgnData_allElecs{iSub,iElecGroup}{iCount,iNeuralMeasure} = powerDataTG_allElecs{ign_Condition}(elecList,:,SSVEPPos)';
                    hits_AttData_GroupWise{iSub}{iCount,iNeuralMeasure}(iElecGroup,:) = powerDataTG_GroupWise{elec_side,att_Condition}(iElecGroup,:,SSVEPPos)';
                    hits_IgnData_GroupWise{iSub}{iCount,iNeuralMeasure}(iElecGroup,:) = powerDataTG_GroupWise{elec_side,ign_Condition}(iElecGroup,:,SSVEPPos)';
                end
            end
        end
    end
end

data_allElecs.attData = hits_AttData_allElecs;
data_allElecs.ignData = hits_IgnData_allElecs;
data_ElecGroupWise.attData = hits_AttData_GroupWise;
data_ElecGroupWise.ignData = hits_IgnData_GroupWise;
end

% function dPrimeVals = getProjection(data_allElecs,regFlag,transformType,numFolds,useEqualStimRepsFlag)
% TFs =[1 2]; % 1- Hits 2 - Miss
% attLoc = [1 2];  % 1- Right 2 - Left
% neuralMeasures = {'alpha','gamma','SSVEP'};
% dPrimeVals = zeros(size(data_allElecs.attData,1),size(data_allElecs.attData,2),length(neuralMeasures)+1);
%
% % Computing d' values for all Subjects x all elec Groups x all Neural
% % measures where LDA was performed on n-electrode dataset for all the
% % trials
% for iSub = 1:size(data_allElecs.attData,1)
%     disp(['Processing data for Sub: ' num2str(iSub)])
%     for iElecGroup = 1:size(data_allElecs.attData,2)
%         for iNeuralMeasure = 1:length(neuralMeasures)
%             dPrimeTMP = zeros(1,length(TFs)*length(attLoc));
%             for iCond = 1:length(TFs)*length(attLoc)
%                 data1 = data_allElecs.attData{iSub,iElecGroup}{iCond,iNeuralMeasure}';
%                 data2 = data_allElecs.ignData{iSub,iElecGroup}{iCond,iNeuralMeasure}';
%                 data1 = data1(all(~isnan(data1),2),:); % removing Bad Electrode data
%                 data2 = data2(all(~isnan(data2),2),:); % removing Bad Electrode data
%                 N1 = size(data1,2);
%                 N2 = size(data2,2);
%
%                 % if att/ign data for a subject is empty because all
%                 % electrodes in a group are bad electrodes
%                 if isempty(data1)||isempty(data2)
%                     dPrimeTMP(iCond) = NaN;
%                 else
%                     [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag);
%                     [~,p1{iCond},p2{iCond},allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,regFlag);
%                     dPrimeTMP(iCond) = getDPrime(p1{iCond}(allIndices1),p2{iCond}(allIndices2));
%                 end
%             end
%             dPrimeVals(iSub,iElecGroup,iNeuralMeasure) = mean(dPrimeTMP,2,'omitnan');
%         end
%     end
% end
%
% % computing d' values for all Neural measures combined (n-electrode data
% % has been pooled and LDA was performed on this pooled dataset for all the
% % trials
%
% for iSub = 1:size(data_allElecs.attData,1)
%     disp(['Processing Combined data for Sub: ' num2str(iSub)])
%     for iElecGroup = 1:size(data_allElecs.attData,2)
%         for iCond = 1:length(TFs)*length(attLoc)
%             clear attDataTMP ignDataTMP
%             for iNeuralMeasure = 1:length(neuralMeasures)
%                 attDataTMP(iNeuralMeasure,:,:) = data_allElecs.attData{iSub,iElecGroup}{iCond,iNeuralMeasure}'; %#ok<*AGROW>
%                 ignDataTMP(iNeuralMeasure,:,:) = data_allElecs.ignData{iSub,iElecGroup}{iCond,iNeuralMeasure}';
%             end
%
%             % concatenating n-electrode dataset
%             data1 = reshape(attDataTMP,[length(neuralMeasures)*size(attDataTMP,2) size(attDataTMP,3)]);
%             data2 = reshape(ignDataTMP,[length(neuralMeasures)*size(attDataTMP,2) size(ignDataTMP,3)]);
%             data1 = data1(all(~isnan(data1),2),:); % removing Bad Electrode data
%             data2 = data2(all(~isnan(data2),2),:); % removing Bad Electrode data
%             N1 = size(data1,2);
%             N2 = size(data2,2);
%
%             % if att/ign data for a subject is empty because all
%             % electrodes in a group are bad electrodes
%             if isempty(data1)||isempty(data2)
%                 dPrimeTMP(iCond) = NaN;
%             else
%                 [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag);
%                 [~,p1,p2,allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,regFlag);
%                 dPrimeTMP(iCond) = getDPrime(p1(allIndices1),p2(allIndices2));
%             end
%         end
%         dPrimeVals(iSub,iElecGroup,iNeuralMeasure+1) = mean(dPrimeTMP,2,'omitnan');
%     end
% end
%
% end
function dPrimeVals = getProjection_v2(data_ElecGroupWise,regFlag,transformType,numFolds,useEqualStimRepsFlag)
TFs =[1 2]; % 1- Hits 2 - Miss
attLoc = [1 2];  % 1- Right 2 - Left
neuralMeasures = {'alpha','gamma','SSVEP'};
elecGroups = {'PO','F','CP','FC','T'};
dPrimeVals = zeros(size(data_ElecGroupWise.attData,2),length(elecGroups)+1,length(neuralMeasures)+1);

% Computing d' values for all Subjects x all elec Groups x all Neural
% measures
for iSub = 1:size(data_ElecGroupWise.attData,2)
    disp(['Processing data for Sub: ' num2str(iSub)])
    for iElecGroup = 1:length(elecGroups)
        for iNeuralMeasure = 1:length(neuralMeasures)
            dPrimeTMP = zeros(1,length(TFs)*length(attLoc));
            for iCond = 1:length(TFs)*length(attLoc)
                clear data1 data2
                data1 = data_ElecGroupWise.attData{iSub}{iCond,iNeuralMeasure}(iElecGroup,:);
                data2 = data_ElecGroupWise.ignData{iSub}{iCond,iNeuralMeasure}(iElecGroup,:);
                dPrimeTMP(iCond) = getDPrime(data1,data2);
            end
            dPrimeVals(iSub,iElecGroup,iNeuralMeasure) = mean(dPrimeTMP,2,'omitnan');
        end
    end
end

% computing d' values for all Neural measures pooled together
% and LDA was performed on this pooled dataset for all the
% trials and saving them in all Subjects x all elec Groups x
% allNeuralMeasures +1

for iSub = 1:size(data_ElecGroupWise.attData,2)
    disp(['Processing Combined data for Sub: ' num2str(iSub)])
    for iElecGroup = 1:length(elecGroups)
        for iCond = 1:length(TFs)*length(attLoc)
            clear attDataTMP ignDataTMP data1 data2
            for iNeuralMeasure = 1:length(neuralMeasures)
                attDataTMP(iNeuralMeasure,:) = data_ElecGroupWise.attData{iSub}{iCond,iNeuralMeasure}(iElecGroup,:); %#ok<*AGROW>
                ignDataTMP(iNeuralMeasure,:) = data_ElecGroupWise.ignData{iSub}{iCond,iNeuralMeasure}(iElecGroup,:);
            end
            
            data1 = attDataTMP;
            data2 = ignDataTMP;
            data1 = data1(all(~isnan(data1),2),:); % removing Bad Electrode data
            data2 = data2(all(~isnan(data2),2),:); % removing Bad Electrode data
            N1 = size(data1,2);
            N2 = size(data2,2);
            
            % if att/ign data for a subject is empty because all
            % electrodes in a group are bad electrodes
            if isempty(data1)||isempty(data2)
                dPrimeTMP(iCond) = NaN;
            else
                [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag);
                [~,p1,p2,allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,regFlag);
                dPrimeTMP(iCond) = getDPrime(p1(allIndices1),p2(allIndices2));
            end
        end
        dPrimeVals(iSub,iElecGroup,length(neuralMeasures)+1) = mean(dPrimeTMP,2,'omitnan');
    end
end

% computing d' values for all electrode groups combined and
% LDA was performed on this pooled dataset for all the trials
% and saving them in all Subjects x all elec Groups+1 x
% allNeuralMeasures

for iSub = 1:size(data_ElecGroupWise.attData,2)
    disp(['Processing Combined data for Sub: ' num2str(iSub)])
    for iNeuralMeasure = 1:length(neuralMeasures)
        clear attDataTMP ignDataTMP clear data1 data2
        for iCond = 1:length(TFs)*length(attLoc)
            attDataTMP = data_ElecGroupWise.attData{iSub}{iCond,iNeuralMeasure}; %#ok<*AGROW>
            ignDataTMP = data_ElecGroupWise.ignData{iSub}{iCond,iNeuralMeasure};
            
            data1 = attDataTMP;
            data2 = ignDataTMP;
            data1 = data1(all(~isnan(data1),2),:); % removing Bad Electrode data
            data2 = data2(all(~isnan(data2),2),:); % removing Bad Electrode data
            N1 = size(data1,2);
            N2 = size(data2,2);
            
            % if att/ign data for a subject is empty because all
            % electrodes in a group are bad electrodes
            if isempty(data1)||isempty(data2)
                dPrimeTMP(iCond) = NaN;
            else
                [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag);
                [~,p1,p2,allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,regFlag);
                dPrimeTMP(iCond) = getDPrime(p1(allIndices1),p2(allIndices2));
            end
        end
        dPrimeVals(iSub,length(elecGroups)+1,iNeuralMeasure) = mean(dPrimeTMP,2,'omitnan');
    end
end

% computing d' values for all electrode groups and all neural measures
% pooled together and LDA was performed on this pooled dataset for all
% the trials and saving them in all Subjects x all elec Groups+1 x
% allNeuralMeasures+1
for iSub = 1:size(data_ElecGroupWise.attData,2)
    disp(['Processing Combined data for Sub: ' num2str(iSub)])
    for iCond = 1:length(TFs)*length(attLoc)
            clear attDataTMP ignDataTMP clear data1 data2            
        for iNeuralMeasure = 1:length(neuralMeasures)
            attDataTMP(iNeuralMeasure,:,:) = data_ElecGroupWise.attData{iSub}{iCond,iNeuralMeasure}; %#ok<*AGROW>
            ignDataTMP(iNeuralMeasure,:,:) = data_ElecGroupWise.ignData{iSub}{iCond,iNeuralMeasure};
        end
        % concatenating neuralMeasures and electrodeGroups dataset
        data1 = reshape(attDataTMP,[length(neuralMeasures)*length(elecGroups) size(attDataTMP,3)]);
        data2 = reshape(ignDataTMP,[length(neuralMeasures)*length(elecGroups) size(ignDataTMP,3)]);
        data1 = data1(all(~isnan(data1),2),:); % removing Bad Electrode data
        data2 = data2(all(~isnan(data2),2),:); % removing Bad Electrode data
        N1 = size(data1,2);
        N2 = size(data2,2);
        
        % if att/ign data for a subject is empty because all
        % electrodes in a group are bad electrodes
        if isempty(data1)||isempty(data2)
            dPrimeTMP(iCond) = NaN;
        else
            [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag);
            [~,p1,p2,allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,regFlag);
            dPrimeTMP(iCond) = getDPrime(p1(allIndices1),p2(allIndices2));
        end
    end
    dPrimeVals(iSub,length(elecGroups)+1,length(neuralMeasures)+1) = mean(dPrimeTMP,2,'omitnan');
end

end

% Aceessory Functions used from displayPopulationLDA2 from
% MayoProjectPrograms repository by Supratim Ray
function [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag)

allIndices1 = randperm(N1);
allIndices2 = randperm(N2);

if useEqualStimRepsFlag
    N1=min(N1,N2);
    N2=N1;
end

allIndices1 = sort(allIndices1(1:N1));
allIndices2 = sort(allIndices2(1:N2));
allIndices = [allIndices1 allIndices2];

testingIndices1 = cell(1,numFolds);
testingIndices2 = cell(1,numFolds);

if numFolds==1 % No cross Validation
    testingIndices1{1} = allIndices1;
    testingIndices2{1} = allIndices2;
else
    Y = [zeros(N1,1) ; ones(N2,1)];
    cvp = cvpartition(Y,'KFold',numFolds); % Partition data
    
    for i=1:numFolds
        testingIDs = find(cvp.test(i)==1);
        testingIndices1{i} = allIndices(testingIDs(testingIDs<=N1));
        testingIndices2{i} = allIndices(testingIDs(testingIDs>N1));
    end
end
end
function [weightVector,projections1,projections2,fullSetIndices1,fullSetIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,regFlag)

numFolds = length(testingIndices1);
fullSetIndices1=[]; fullSetIndices2=[];
for i=1:numFolds
    fullSetIndices1 = cat(2,fullSetIndices1,testingIndices1{i});
    fullSetIndices2 = cat(2,fullSetIndices2,testingIndices2{i});
end
fullSetIndices1 = sort(fullSetIndices1);
fullSetIndices2 = sort(fullSetIndices2);

projections1 = zeros(size(data1,2),1);
projections2 = zeros(size(data2,2),1);
D = size(data1,1);

weightVectorTMP = zeros(D,numFolds);

for i=1:numFolds
    t1 = testingIndices1{i};
    t2 = testingIndices2{i};
    
    if numFolds==1 % No cross validation. Train and test on the same data
        train1 = t1;
        train2 = t2;
    else
        train1 = setdiff(fullSetIndices1,t1);
        train2 = setdiff(fullSetIndices2,t2);
    end
    d1 = data1(:,train1); d2 = data2(:,train2);
    m1 = size(d1,2); m2 = size(d2,2);
    
    if transformType==1 % LDA for uncorrelated case
        weightVectorTMP(:,i) = 1/size(d1,1);
        
    elseif transformType==2 % LDA for uncorrelated case
        meanDiff = mean(d1,2) - mean(d2,2);
        var1 = var(d1,[],2);
        var2 = var(d2,[],2);
        pooledVar = ((m1-1)*var1 + (m2-1)*var2)/(m1+m2-2); %Pooled Variance
        
        weightVectorTMP(:,i) = meanDiff./pooledVar;
        
    elseif transformType==3 % LDA
        
        label1 = repmat({'Att'},m1,1);
        label2 = repmat({'Ign'},m2,1);
        labelAll = cat(1,label1,label2);
        dataAll = cat(2,d1,d2);
        
        if regFlag==0 % No regularization
            Mdl = fitcdiscr(dataAll',labelAll);
            
        elseif regFlag==1 % Only optimize gamma
            Mdl = fitcdiscr(dataAll',labelAll);
            [err,gamma,~,~] = cvshrink(Mdl,'NumGamma',20);
            [~,minIndex] = min(err);
            if minIndex>1
                Mdl.Gamma = gamma(minIndex); % Changing gamma changes the model weights also
            end
            
        elseif regFlag==2 % Optimize gamma and delta
            Mdl = fitcdiscr(dataAll',labelAll);
            [err,gamma,delta,~] = cvshrink(Mdl,'NumGamma',20,'numDelta',20);
            minerr = min(err(:));
            [x,y] = find(err==minerr);
            if x(1)>1 || y(1)>1
                Mdl.Gamma = gamma(x(1)); % Take the smallest gamma
                Mdl.Delta = delta(x(1),y(1)); % Take the smallest delta
            end
            
        elseif regFlag==3 % Hyper optimize gamma and delta
            rng(1)
            myOpts.AcquisitionFunctionName = 'expected-improvement-plus';
            myOpts.ShowPlots = 0;
            myOpts.Verbose = 0;
            Mdl = fitcdiscr(X,Y,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',myOpts);
            
        end
        weightVectorTMP(:,i) = Mdl.Coeffs(1,2).Linear;
    end
    projections1(t1) = data1(:,t1)' *  weightVectorTMP(:,i);
    projections2(t2) = data2(:,t2)' *  weightVectorTMP(:,i);
end

weightVector = mean(weightVectorTMP,2);
end

function elecs = getElectrodeGroupInformation(elecGroup)
if strcmp(elecGroup,'PO')
    elecs{1} = [24 29 57 61];
    elecs{2} = [26 31 58 63];
elseif strcmp(elecGroup,'F')
    elecs{1} = [1 33 34 3 37 4];
    elecs{2} = [2 35 36 6 40 7];
elseif strcmp(elecGroup,'CP')
    elecs{1} = [18 19 52 23 56];
    elecs{2} = [20 21 27 54 59];
elseif strcmp(elecGroup,'FC')
    elecs{1} = [8 9 13 43 47 48];
    elecs{2} = [10 11 15 44 49 50];
elseif strcmp(elecGroup,'T')
    elecs{1} = [12 17 41 42 51];
    elecs{2} = [16 22 45 46 55];
end
end

function d = getDPrime(x1,x2)
n1 = length(x1);    n2 = length(x2);
stdVal = sqrt(((n1-1)*var(x1,'omitnan')+(n2-1)*var(x2,'omitnan'))/(n1+n2-2));
d = (mean(x1,'omitnan')- mean(x2,'omitnan'))/stdVal;
end