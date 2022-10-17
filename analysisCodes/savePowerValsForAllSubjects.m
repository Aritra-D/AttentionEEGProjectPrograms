% This Program saves the power values for Attention Dataset (26 human
% individuals) for all 64 electrodes segregated in  12 cell arrays,
% where the cells correspond to the following Stimulus/Attention 
% Conditions as well as saves the power data for different electrode 
% groups for both left and right hemispheric EEG electrodes
% The 12 cell arrays correspond to the following Stimulus/Attention Conditions: 
% {'H0V_0Hz','H1V_0Hz','H0V_12Hz','H1V_12Hz','H0V_16Hz','H1V_16Hz',
%  'M0V_0Hz','M1V_0Hz','M0V_12Hz','M1V_12Hz','M0V_16Hz','M1V_16Hz',}


function [powerData_AllElecs, powerData_GroupWise, badElecs, freqRanges_SubjectWise] = savePowerValsForAllSubjects(protocolType, gridType, refType, dataFolderSourceString, tapers, badTrialStr)
if ~exist('folderSourceString','var');  folderSourceString ='E:\data\human\SRCLong';          end
if ~exist('gridType','var');            gridType = 'EEG';               end
if ~exist('protocolType','var');        protocolType = 'SRC-Long';      end
if ~exist('badTrialStr','var');         badTrialStr = 'v10';            end
if ~exist('saveDataFlag','var');        saveDataFlag = 1;               end

capType = 'actiCap64';
folderSave = strtok(dataFolderSourceString,'\');
saveFileName = fullfile(folderSave, ['powerData_' protocolType '_' gridType '_allSubjects_N26' '_tapers_' num2str(tapers(2)) '_badTrial_' badTrialStr '.mat']);
if exist(saveFileName,'file')
    load(saveFileName)
else
    [powerData_AllElecs, powerData_GroupWise, badElecs, freqRanges_SubjectWise] = getDataForAllSubjects(protocolType, gridType, capType, refType, dataFolderSourceString, tapers, badTrialStr);
    save(saveFileName,'powerData_AllElecs','powerData_GroupWise','badElecs','freqRanges_SubjectWise')
end
end

%% Accessory Nested Functions

% get Data for allSubjects
function [powerData_AllElecs, powerData_GroupWise, badElecs, freqRanges_SubjectWise] = getDataForAllSubjects(protocolType, gridType, capType, refType, dataFolderSourceString, tapers, badTrialStr)

% Loading data/protocol Information
[subjectNames,expDates,protocolNames,dataFolderSourceString] = dataInformationSRCProtocols_HumanEEG(gridType,protocolType,dataFolderSourceString);

for iSub = 1:length(subjectNames)
    tic
    subject = subjectNames{iSub};
    expDate = expDates{iSub};
    protocolName = protocolNames{iSub};
    disp(['SUBJECT_ID: ' num2str(iSub) ', SUBJECTNAME: ' subjectNames{iSub} ', EXPDATE: ' num2str(expDates{iSub}) ', PROTOCOL: ' num2str(protocolNames{iSub})])
    
    if iSub<8 % First Set of Recording- Nov-Dec 2021
        freqRanges{1} = [8 12];    % alpha
        freqRanges{2} = [25 70];   % gamma
        freqRanges{3} = [23 23];   % SSVEP Left Stim; Flicker Freq moved by 0.5 Hz due one extra blank Frame
        freqRanges{4} = [31 31];   % SSVEP Right Stim; Flicker Freq moved by 0.5 Hz due one extra blank Frame
    else % Second Set of Recording- Jan-Mar 2022
        freqRanges{1} = [8 12];    % alpha
        freqRanges{2} = [25 70];   % gamma
        freqRanges{3} = [24 24];   % SSVEP Left Stim; Flicker Freq bug Fixed
        freqRanges{4} = [32 32];   % SSVEP Right Stim; Flicker Freq bug Fixed
    end
    
    numFreqs = length(freqRanges); %#ok<*NASGU>
    freqRanges_SubjectWise{iSub} = freqRanges; %#ok<*AGROW>
    
    [powerData_AllElecs{iSub}, powerData_GroupWise{iSub}, badElecs{iSub}]= getDataForIndividualSubject(subject, expDate, protocolName, dataFolderSourceString, gridType, capType, refType, tapers, badTrialStr, freqRanges);
    toc
end
end

% get Data for Individual Subjects
function [powerData_AllElecs, powerData_GroupWise, badElecs] = getDataForIndividualSubject(subject, expDate, protocolName, dataFolderSourceString, gridType, capType, refType, tapers, badTrialStr, freqRanges)

electrodes = getElectrodeList(capType, refType, 1);
folderName = fullfile(dataFolderSourceString,'data',subject,gridType,expDate,protocolName);
folderExtract= fullfile(folderName,'extractedData');
folderSegment= fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');

% Get Parameter Combinations for SRC-Long Protocols
[parameterCombinations,cValsUnique,tValsUnique,eValsUnique,aValsUnique,sValsUnique] = loadParameterCombinations(folderExtract);

% Segregating Trials to Static and Flickering trials
trialIDs_Static = parameterCombinations.stimOnset{1,1,3,3,1};
trialIDs_Flickering = unique([parameterCombinations.stimOnset{1,2,3,3,1} parameterCombinations.stimOnset{1,3,3,3,1}]);
totalTrialNums = length(trialIDs_Static)+ length(trialIDs_Flickering);

trialTFs = zeros(1,totalTrialNums);
trialTFs(1,trialIDs_Flickering) = 1; % 0: Static 1: Flickering

% load LFP Info
[analogChannelsStored, timeVals, ~, ~] = loadlfpInfo(folderLFP);

% Timing Related Information
Fs = round(1/(timeVals(2)- timeVals(1)));
blRange = [-1 0];
tgRange = [-1 0];
range = blRange;
rangePos = round(diff(range)*Fs);
blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
tgPos = find(timeVals>=tgRange(1),1)+ (1:rangePos);


% Set up params for MT
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 250];
params.trialave = 0;

photoDiodeChannels = [65 66];
% Get bad trials and Electrodes
badTrialFile = fullfile(folderSegment,['badTrials_' badTrialStr '.mat']);
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badElecs = []; badTrials=[];
else
    [badTrials,badElectrodes,badTrialsUnique] = loadBadTrials(badTrialFile);
    badElecsAll = unique([badElectrodes.badImpedanceElecs; badElectrodes.noisyElecs; badElectrodes.flatPSDElecs; badElectrodes.flatPSDElecs]);
    disp([num2str(length(union(badTrials,badTrialsUnique.badEyeTrials))) ' bad trials']);
    
    badElecs = intersect(setdiff(analogChannelsStored,photoDiodeChannels),badElecsAll);
    disp(['Unipolar, all bad elecs: ' num2str(length(badElecsAll))]);
end

contrast = length(cValsUnique);
tempFreq = length(tValsUnique);
eotCodes = length(eValsUnique);
attLoc = length(aValsUnique);
stimType = length(sValsUnique);

goodPos_all = setdiff(parameterCombinations.targetOnset{contrast,tempFreq+1,eotCodes+1,attLoc+1,stimType},union(badTrials,badTrialsUnique.badEyeTrials));

% Segmenting data according to timePos
hW1 = waitbar(0,'collecting data...');
for iElec = 1: length(electrodes)
    waitbar((iElec-1)/length(electrodes),hW1,['collecting data from electrode: ' num2str(iElec) ' of ' num2str(length(electrodes))]);
    clear x
    x = load(fullfile(folderLFP,['elec' num2str(electrodes{iElec}{1}) '.mat']));
    dataBL = x.analogData.stimOnset(:,blPos)';
    dataTG = x.analogData.targetOnset(:,tgPos)';
    
    % power spectral density estimation
    [tmpEBL,~] = mtspectrumc(dataBL,params);
    [tmpETG,freqVals] = mtspectrumc(dataTG,params);
    
    
    for iTrial = 1:size(dataTG,2)
        if trialTFs(iTrial) == 0
            tfVal1 = 0; tfVal2 =0;
        elseif trialTFs(iTrial) == 1
            tfVal1 = unique(freqRanges{3}/2);
            tfVal2 = unique(freqRanges{4}/2);
        end
        for iFreqRange=1:length(freqRanges)
            if iFreqRange == 3||iFreqRange == 4
                remove_NthHarmonicOnwards = 3;
            else
                remove_NthHarmonicOnwards = 2;
            end
            deltaF_LineNoise = 2; deltaF_tfHarmonics = 0;
            badFreqPos = getBadFreqPos(freqVals,deltaF_LineNoise,deltaF_tfHarmonics,remove_NthHarmonicOnwards,tfVal1,tfVal2);
            if find(iElec == badElecs)
                powerValsBL(iElec,iTrial,iFreqRange) = NaN;
                powerValsTG(iElec,iTrial,iFreqRange) = NaN;
            else
                powerValsBL(iElec,iTrial,iFreqRange) = getMeanEnergyForAnalysis(tmpEBL(:,iTrial),freqVals,freqRanges{iFreqRange},badFreqPos);
                powerValsTG(iElec,iTrial,iFreqRange) = getMeanEnergyForAnalysis(tmpETG(:,iTrial),freqVals,freqRanges{iFreqRange},badFreqPos);
            end
        end
    end
end
close(hW1);
nanFlag = 'omitnan';

% Grouping actiCAP-64 EEG electrodes into separate groups according to
% location

% Parieto-Occipital Group % electrode numbers 24 26 57 58 are originally included in centro-parietal electrode list from electrodePositionOnGrid func.
showOccipitalElecsUnipolarLeft = [24 29 57 61]; 
showOccipitalElecsUnipolarRight = [26 31 58 63];

% Frontal Group
showFrontalElecsUnipolarLeft = [1 33 34 3 37 4];
showFrontalElecsUnipolarRight = [2 35 36 6 40 7];

% Centro-parietal Group
showCentroParietalElecsUnipolarLeft = [18 19 52 23 56];
showCentroParietalElecsUnipolarRight = [20 21 27 54 59];

% Fronto-central Group
showFrontoCentralElecsUnipolarLeft = [8 9 13 43 47 48];
showFrontoCentralElecsUnipolarRight = [10 11 15 44 49 50];

% Temporal Group
showTemporalElecsUnipolarLeft = [12 17 41 42 51];
showTemporalElecsUnipolarRight = [16 22 45 46 55];

% Segregate the data in  12 cell array of , where the cells
% correspond to the following Stimulus/Attention Conditions
% {'H0V_0Hz','H1V_0Hz','H0V_12Hz','H1V_12Hz','H0V_16Hz','H1V_16Hz',
%  'M0V_0Hz','M1V_0Hz','M0V_12Hz','M1V_12Hz','M0V_16Hz','M1V_16Hz',}

for iCondition = 1:12
    clear c tf eotCode attLoc s
    switch iCondition
        case 1;  c = 1; tf = 1; eotCode = 1; attLoc = 2;  s=1;
        case 2;  c = 1; tf = 1; eotCode = 1; attLoc = 1;  s=1;
        case 3;  c = 1; tf = 2; eotCode = 1; attLoc = 2;  s=1;
        case 4;  c = 1; tf = 2; eotCode = 1; attLoc = 1;  s=1;
        case 5;  c = 1; tf = 3; eotCode = 1; attLoc = 2;  s=1;
        case 6;  c = 1; tf = 3; eotCode = 1; attLoc = 1;  s=1;
        case 7;  c = 1; tf = 1; eotCode = 2; attLoc = 2;  s=1;
        case 8;  c = 1; tf = 1; eotCode = 2; attLoc = 1;  s=1;
        case 9;  c = 1; tf = 2; eotCode = 2; attLoc = 2;  s=1;
        case 10; c = 1; tf = 2; eotCode = 2; attLoc = 1;  s=1;
        case 11; c = 1; tf = 3; eotCode = 2; attLoc = 2;  s=1;
        case 12; c = 1; tf = 3; eotCode = 2; attLoc = 1;  s=1;
    end
    goodPos{iCondition} = setdiff(parameterCombinations.targetOnset{c,tf,eotCode,attLoc,s},union(badTrials,badTrialsUnique.badEyeTrials));
    powerValsBL_ConditionWise{iCondition} = powerValsBL(:,goodPos{iCondition},:);
    powerValsTG_ConditionWise{iCondition} = powerValsTG(:,goodPos{iCondition},:);
end

powerData_AllElecs.powerValsBL = powerValsBL_ConditionWise;
powerData_AllElecs.powerValsTG = powerValsTG_ConditionWise;


% Segregate the data in 2 x 12 cell arrays where the first cell
% dimension represents the Left (Index 1) and Right (Index 2)
% Electrodes respectively and 12 cells correspond to the following
% Stimulus/Attention Conditions
% {'H0V_0Hz','H1V_0Hz','H0V_12Hz','H1V_12Hz','H0V_16Hz','H1V_16Hz',
%  'M0V_0Hz','M1V_0Hz','M0V_12Hz','M1V_12Hz','M0V_16Hz','M1V_16Hz',}

powerValsBL_ElecGroupAndConditionWise = cell(2,12);
powerValsTG_ElecGroupAndConditionWise = cell(2,12);
elecGroupLabels = {'Occipital','Frontal','CentroParietal','FrontoCentral','Temporal'};

for iElecSide = 1:2
    switch iElecSide
        case 1 % Left hemispheric EEG electrodes
            elecGroups{1} = showOccipitalElecsUnipolarLeft;
            elecGroups{2} = showFrontalElecsUnipolarLeft;
            elecGroups{3} = showCentroParietalElecsUnipolarLeft;
            elecGroups{4} = showFrontoCentralElecsUnipolarLeft;
            elecGroups{5} = showTemporalElecsUnipolarLeft;
        case 2 % Right hemispheric EEG electrodes
            elecGroups{1} = showOccipitalElecsUnipolarRight;
            elecGroups{2} = showFrontalElecsUnipolarRight;
            elecGroups{3} = showCentroParietalElecsUnipolarRight;
            elecGroups{4} = showFrontoCentralElecsUnipolarRight;
            elecGroups{5} = showTemporalElecsUnipolarRight;
    end
    
    for iCondition = 1:12
        for iElecGroup = 1:length(elecGroups)
            powerValsBL_ElecGroupAndConditionWise{iElecSide,iCondition}(iElecGroup,:,:) = squeeze(mean(powerValsBL_ConditionWise{iCondition}(elecGroups{iElecGroup},:,:),1,nanFlag));
            powerValsTG_ElecGroupAndConditionWise{iElecSide,iCondition}(iElecGroup,:,:) = squeeze(mean(powerValsTG_ConditionWise{iCondition}(elecGroups{iElecGroup},:,:),1,nanFlag));
        end
    end
end

powerData_GroupWise.powerValsBL = powerValsBL_ElecGroupAndConditionWise;
powerData_GroupWise.powerValsTG = powerValsTG_ElecGroupAndConditionWise;

end

% Load LFP Info
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored); %#ok<NODEF>
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

% Get parameter combinations
function [parameterCombinations,cValsUnique,tValsUnique,eValsUnique,...
    aValsUnique,sValsUnique] = ...
    loadParameterCombinations(folderExtract)
load(fullfile(folderExtract,'parameterCombinations.mat')); %#ok<*LOAD>
end

% Get Bad Trials
function [badTrials,badElecs,badTrialsUnique] = loadBadTrials(badTrialFile) %#ok<*STOUT>
load(badTrialFile);
end

% Get MeanEnergy for different frequency bands
function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange,badFreqPos)

posToAverage = setdiff(intersect(find(freq>=freqRange(1)),find(freq<=freqRange(2))),badFreqPos);
eValue   = sum(mEnergy(posToAverage));
end

function badFreqPos = getBadFreqPos(freqVals,deltaF_LineNoise,deltaF_TFHarmonics,remove_NthHarmonicOnwards,tf1,tf2)
% During this Project, line Noise was at
% 51 Hz for 1 Hz Freq Resolution and

if nargin<2
    deltaF_LineNoise = 1; deltaF_TFHarmonics = 0; tf1 = 0; tf2 = 0;
end

if tf1>0 && tf2>0 % Flickering Stimuli
    badFreqs = 51:51:max(freqVals);
    tfHarmonics1 = remove_NthHarmonicOnwards*tf1:tf1:max(freqVals); % remove nth SSVEP harmonic and beyond
    tfHarmonics2 = remove_NthHarmonicOnwards*tf2:tf2:max(freqVals); % remove nth SSVEP harmonic and beyond
    tfHarmonics = unique([tfHarmonics1 tfHarmonics2]);
elseif tf1==0 && tf2==0 % Static Stimuli
    badFreqs = 51:51:max(freqVals);
end

badFreqPos = [];
for i=1:length(badFreqs)
    badFreqPos = cat(2,badFreqPos,intersect(find(freqVals>=badFreqs(i)-deltaF_LineNoise),find(freqVals<=badFreqs(i)+deltaF_LineNoise)));
end

if exist('tfHarmonics','var')
    freqPosToRemove =  [];
    for i=1:length(badFreqs)
        freqPosToRemove = cat(2,freqPosToRemove,intersect(find(freqVals>=tfHarmonics(i)-deltaF_TFHarmonics),find(freqVals<=tfHarmonics(i)+deltaF_TFHarmonics)));
    end
    badFreqPos = unique([badFreqPos freqPosToRemove]);
end
end
