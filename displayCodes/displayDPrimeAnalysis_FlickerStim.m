% d-prime analysis results for flicker stimuli for
% three neural measures -  alpha, gamma, and SSVEP

clear; % clears the Workspace
close all; % closes any open Figures

% Set FolderSourceString automatically; add files & folders to MATLAB path
fileNameWithPath = mfilename('fullpath');
[filePath,name,ext] = fileparts(fileNameWithPath);
cd(filePath); cd ..
addpath(genpath(pwd));
folderSourceString = fullfile(pwd,'analyzedData');

if ~exist('gridType','var');            gridType = 'EEG';                                       end
if ~exist('protocolType','var');        protocolType = 'SRC-Long';                              end
if ~exist('badTrialStr','var');         badTrialStr = 'v10';                                    end

tapers = [1 1];
fileName = fullfile(folderSourceString,['powerData_',protocolType,'_',gridType,'_allSubjects_N26_tapers_' num2str(tapers(2)) '_badTrial_' badTrialStr '.mat']);
if exist(fileName,'file')
    load(fileName,'powerData_AllElecs','powerData_GroupWise');
else
    error('Datafile Not Available! Check Data Path!')
end

nanFlag = 'omitnan'; % NaN Values are processed using this flag during averaging

% Display options
hFig = figure;
set(hFig,'units','normalized','outerPosition',[0 0 1 1]);
hPlot = getPlotHandles(2,5,[0.2 0.1, 0.75 0.7],0.05,0.1,0);
fontSize = 14;
subjectIdx = 1:26; % Numbers 1-26 are subject IDs.
colors = {'k','r','c'}; % Barplot colors for alpha, gamma and SSVEP respectively
ElecGroup  = 1:5 ;% 1: Parieto-Occipital , 2: Frontal, 3: Centro-Parietal Group, 4: Fronto-Central Group, 5: Temporal Group
freqRangeLabels = {'alpha','gamma','SSVEP'}; % alpha, gamma, and SSVEP

% Compute dPrime Values for Att
for iElecGroup =1:length(ElecGroup)
    for iFreqRange = 1:length(freqRangeLabels)
        clear dPrime_Hit dPrime_Miss
        switch iFreqRange
            case 1;  neuralMeasure = 'alpha';
            case 2;  neuralMeasure = 'gamma';
            case 3;  neuralMeasure = 'SSVEP';
        end
        [dPrime_Hit,dPrime_Miss] = compute_dPrimeVals_flickerStim(powerData_GroupWise,subjectIdx,neuralMeasure,iElecGroup,iFreqRange);
        dPrime_Hit_AttMInusIgn{iElecGroup,iFreqRange} = dPrime_Hit; %#ok<*SAGROW>
        dPrime_Miss_AttMinusIgn{iElecGroup,iFreqRange} = dPrime_Miss;
    end
    
    % Plotting Bar Plots of d-Prime Values
    for iAttend = 1:2   %1 - contralateral attend ; 2- ipsilateral attend
        switch iAttend
            case 1;  dPrime = dPrime_Hit_AttMInusIgn;
            case 2;  dPrime = dPrime_Miss_AttMinusIgn;
        end
        for iBar=1:length(freqRangeLabels)
            mBar = mean(dPrime{iElecGroup,iBar},1,nanFlag); %mean across subjects
            errorBar = std(dPrime{iElecGroup,iBar},[],1,nanFlag)./sqrt(length(dPrime{iElecGroup,iBar}));
            mBars(iBar) =   mBar;
            eBars(iBar) = errorBar;
            subplot(hPlot(iAttend,iElecGroup));hold(hPlot(iAttend,iElecGroup),'on');
            barPlot = bar(iBar,mBar);
            barPlot.FaceColor = colors{iBar};
        end
        errorbar(hPlot(iAttend,iElecGroup),1:length(mBars),mBars,eBars,'.','color','k');
    end
end

%setting up legends and labels
elecGroup = [{'Parieto-Occipital'}, {'      Frontal'},{'Centro-Parietal'},{'Fronto-Central'},{'Temporal'}];
for iElec =1:length(elecGroup)
    subplot(hPlot(1,iElec))
    p=get(gca);
    annotation('textbox',[1*(p.Position(1,1)) 0.84 0.094 0.055],'EdgeColor','none','String',elecGroup{iElec},'fontSize',16,'fontWeight','bold');
end
subplot(hPlot(1,1))
annotation('textbox',[0.06 0.64 0.084 0.058],'EdgeColor','none','String','Hits','fontSize',20,'fontWeight','bold');
subplot(hPlot(2,1))
annotation('textbox',[0.06 0.198 0.084 0.058],'EdgeColor','none','String','Misses','fontSize',20,'fontWeight','bold');
Datalabels = {'alpha','gamma','SSVEP-TbT'};
tickLength = 2*get(hPlot(1,1),'TickLength');
for i=1:5
    if i>1
        set(hPlot(1,i),'xTick',1:3,'xTickLabel',[],'yLim',[-0.3 0.3],'yTick',-0.3:0.1:0.3,'TickDir','out','TickLength',tickLength,'fontSize',fontSize)
        set(hPlot(2,i),'xTick',1:3,'xTickLabel',[],'yLim',[-0.3 0.3],'yTick',-0.3:0.1:0.3,'TickDir','out','TickLength',tickLength,'fontSize',fontSize)
    else
        set(hPlot(1,i),'xTick',1:3,'xTickLabel',Datalabels,'XTickLabelRotation',30,'yLim',[-0.3 0.3],'yTick',-0.3:0.1:0.3,'TickDir','out','TickLength',tickLength,'fontSize',fontSize)
        set(hPlot(2,i),'xTick',1:3,'xTickLabel',Datalabels,'XTickLabelRotation',30,'yLim',[-0.3 0.3],'yTick',-0.3:0.1:0.3,'TickDir','out','TickLength',tickLength,'fontSize',fontSize)
    end
end
subplot(hPlot(2,1))
xlim=get(gca,'XLim');
ylim=get(gca,'YLim');
ht = text(xlim(1)-2,ylim(1)+0.3*ylim(2),'Neural d'' (Att-Ign)', 'fontSize',14);
set(ht,'Rotation',90)

% Accessory function to calculate d Prime Values for attendVsIgn condition
% for hits and miss
function [dPrime_Hits_AttMinusIgnore,dPrime_Miss_AttMinusIgnore] = compute_dPrimeVals_flickerStim(powerData,subjectIdx,neuralMeasure,elecGroup,freqRange)

TFs =[1 2]; % 1- Hits 2 - Miss
attLoc = [1 2];  % 1- Right 2 - Left

if strcmp(neuralMeasure,'alpha')||strcmp(neuralMeasure,'gamma')
    for iSub=1:length(subjectIdx)
        clear topoDataTG 
        topoDataTG = powerData{subjectIdx(iSub)}.powerValsTG;
        % Each condition refers to separate attended location and
        % flicker freuency of either 12 or 16 Hz frequency. For each condition we
        % calculate the d prime for attend - ignored from the hits condition
        for iCount=1:length(attLoc)*length(TFs)
            %        clear hitsData_Attend hitsData_Ignored
            switch iCount
                % here elec_side refers to the row of the cell array of powerData_Groupwise matrix - 1st row is for all left elecs and 2nd row for right elecs
                case 1; att_Condition = 3; ign_Condition = 6;  elec_side = 2;   % att_condition - HL 12Hz ign_condition - HR 16Hz
                case 2; att_Condition = 4; ign_Condition = 5;  elec_side = 1;    % att_condition - HR 12Hz ign_condition - HL 16Hz
                case 3; att_Condition = 5; ign_Condition = 4;  elec_side = 2;  % att_condition - HL 16Hz ign_condition - HR 12Hz
                case 4; att_Condition = 6; ign_Condition = 3;  elec_side = 1;  % att_condition - HR 16Hz ign_condition - HL 12Hz
            end
            hitsData_Attend = topoDataTG{elec_side,att_Condition}(elecGroup,:,freqRange)';
            hitsData_Ignored = topoDataTG{elec_side,ign_Condition}(elecGroup,:,freqRange)';
            dPrime(iCount,:) = getDPrime(hitsData_Attend,hitsData_Ignored); %#ok<*AGROW>
        end
        % averaging actoss the attend conditions for each subject
        dPrime_Hits_AttMinusIgnore(iSub,:) = mean(dPrime,1);
        
        % Each condition refers to separate attended location and
        % flicker freuency of either 12 or 16 Hz frequency. For each condition we
        % calculate the d prime for attend - ignored from the miss condition
        for iCount=1:length(attLoc)*length(TFs)
            clear hitsData_Attend hitsData_Ignored
            switch iCount
                % here elec_side refers to the row of the cell array of powerData_Groupwise matrix - 1st row is for all left elecs and 2nd row for right elecs
                case 1; att_Condition = 9; ign_Condition =12;  elec_side = 2;   % att_condition - ML 12Hz ign_condition - MR 16Hz
                case 2; att_Condition = 10; ign_Condition =11;  elec_side = 1;    % att_condition - MR 12Hz ign_condition - ML 16Hz
                case 3; att_Condition = 11; ign_Condition =10;  elec_side = 2;  % att_condition - ML 16Hz ign_condition - MR 12Hz
                case 4; att_Condition = 12; ign_Condition =9;  elec_side = 1;  % att_condition - MR 16Hz ign_condition - ML 12Hz
            end
            
            if length(size(topoDataTG{elec_side,att_Condition}))==2 || length(size(topoDataTG{elec_side,ign_Condition}))==2
                hitsData_Attend = topoDataTG{elec_side,att_Condition}(elecGroup,freqRange)';
                hitsData_Ignored = topoDataTG{elec_side,ign_Condition}(elecGroup,freqRange)';
            elseif length(size(topoDataTG{elec_side,att_Condition}))==3 && length(size(topoDataTG{elec_side,ign_Condition}))==3
                hitsData_Attend = topoDataTG{elec_side,att_Condition}(elecGroup,:,freqRange)';
                hitsData_Ignored = topoDataTG{elec_side,ign_Condition}(elecGroup,:,freqRange)';
            end
            dPrime2(iCount,:) = getDPrime(hitsData_Attend,hitsData_Ignored); %#ok<*AGROW>
        end
        dPrime_Miss_AttMinusIgnore(iSub,:) = mean(dPrime2,1);
    end
    
elseif strcmp(neuralMeasure,'SSVEP')
    
    for iSub=1:length(subjectIdx)
        clear topoDataTG topoDataBL
        topoDataTG = powerData{subjectIdx(iSub)}.powerValsTG;
        % Attend Contra - Each condition refers to separate attend side and
        % attending to either 12 or 16 Hz frequency. For each condition we
        % calculate the d prime for hits and miss from the contralateral
        % electrodes
        for iCount=1:length(attLoc)*length(TFs)
            clear hitsData_Attend hitsData_Ignored
            switch iCount
                % here elec_side refers to the row of the cell array of powerData_Groupwise matrix - 1st row is for all left elecs and 2nd row for right elecs
                case 1; SSVEPFreqPos = 3; att_Condition = 3; ign_Condition =6;  elec_side = 2;   % att_condition - HL 12Hz ign_condition - HR 16Hz
                case 2; SSVEPFreqPos = 3; att_Condition = 4; ign_Condition =5;  elec_side = 1;    % att_condition - HR 12Hz ign_condition - HL 16Hz
                case 3; SSVEPFreqPos = 4; att_Condition = 5; ign_Condition =4;  elec_side = 2;  % att_condition - HL 16Hz ign_condition - HR 12Hz
                case 4; SSVEPFreqPos = 4; att_Condition = 6; ign_Condition =3;  elec_side = 1;  % att_condition - HR 16Hz ign_condition - HL 12Hz
            end
            hitsData_Attend = topoDataTG{elec_side,att_Condition}(elecGroup,:,SSVEPFreqPos)';
            hitsData_Ignored = topoDataTG{elec_side,ign_Condition}(elecGroup,:,SSVEPFreqPos)';
            dPrime(iCount,:) = getDPrime(hitsData_Attend,hitsData_Ignored); %#ok<*AGROW>
        end
        
        dPrime_Hits_AttMinusIgnore(iSub,:) = mean(dPrime,1);
        
        % Attend Ipsi - Each condition refers to separate attend side and
        % attending to either 12 or 16 Hz frequency. For each condition we
        % calculate the d prime for hits and miss from the ipsilateral
        % electrodes
        for iCount=1:length(attLoc)*length(TFs)
            clear hitsData_Attend hitsData_Ignored
            switch iCount
                %here elec_side refers to the row of the cell array of powerData_Groupwise matrix - 1st row is for all left elecs and 2nd row for right elecs
                case 1; SSVEPFreqPos = 3;  att_Condition = 9; ign_Condition =12;  elec_side = 2; % att_condition - ML 12Hz ign_condition - MR 16Hz
                case 2; SSVEPFreqPos = 3; att_Condition = 10; ign_Condition =11;  elec_side = 1;    % att_condition - MR 12Hz ign_condition - ML 16Hz
                case 3; SSVEPFreqPos = 4; att_Condition = 11; ign_Condition =10;  elec_side = 2;  % att_condition - ML 16Hz ign_condition - MR 12Hz
                case 4; SSVEPFreqPos = 4; att_Condition = 12; ign_Condition =9;  elec_side = 1;  % att_condition - MR 16Hz ign_condition - ML 12Hz
            end
            if length(size(topoDataTG{elec_side,att_Condition}))==2 || length(size(topoDataTG{elec_side,ign_Condition}))==2 % if only one trial is available for either att or ign condition
                hitsData_Attend = topoDataTG{elec_side,att_Condition}(elecGroup,SSVEPFreqPos)';
                hitsData_Ignored = topoDataTG{elec_side,ign_Condition}(elecGroup,SSVEPFreqPos)';
            elseif length(size(topoDataTG{elec_side,att_Condition}))==3 && length(size(topoDataTG{elec_side,ign_Condition}))==3 % if more than one trial is available for both att and ign condition
                hitsData_Attend = topoDataTG{elec_side,att_Condition}(elecGroup,:,SSVEPFreqPos)';
                hitsData_Ignored = topoDataTG{elec_side,ign_Condition}(elecGroup,:,SSVEPFreqPos)';
            end
            dPrime2(iCount,:) = getDPrime(hitsData_Attend,hitsData_Ignored); %#ok<*AGROW>
        end
        dPrime_Miss_AttMinusIgnore(iSub,:) = mean(dPrime2,1);
    end
end
end

function d = getDPrime(x1,x2)
n1 = length(x1);    n2 = length(x2);
stdVal = sqrt(((n1-1)*var(x1,'omitnan')+(n2-1)*var(x2,'omitnan'))/(n1+n2-2));
d = (mean(x1,'omitnan')- mean(x2,'omitnan'))/stdVal;
end
