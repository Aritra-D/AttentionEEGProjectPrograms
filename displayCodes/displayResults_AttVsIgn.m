function displayResults_AttVsIgn(folderSourceString,subjectIdx,stimType,badTrialStr,colorMap,topoplot_style)
if ~exist('folderSourceString','var');  folderSourceString ='E:\data\human\SRCLong';            end
if ~exist('gridType','var');            gridType = 'EEG';                                       end
if ~exist('protocolType','var');        protocolType = 'SRC-Long';                              end
if ~exist('badTrialStr','var');         badTrialStr = 'v10';                                    end

% close all;

tapers = [1 1];
fileName = fullfile(folderSourceString,['powerData_',protocolType,'_',gridType,'_allSubjects_N26_tapers_' num2str(tapers(2)) '_badTrial_' badTrialStr '.mat']);
if exist(fileName,'file')
    load(fileName,'powerData_AllElecs');
else
    error('Datafile Not Available! Check Data Path!')
end
nanFlag = 'omitnan'; % NaN Values are processed using this flag during averaging

% Figure Handles
if strcmp(stimType,'static')
    hFig1 = figure;colormap(colorMap)
    set(hFig1,'units','normalized','outerPosition',[0 0 1 1]);
    hPlot1 = getPlotHandles(2,3,[0.05 0.5, 0.5 0.4],0.005,0.07,0);
    hPlot2 = getPlotHandles(1,2,[0.61 0.5, 0.35 0.3],0.02,0.07,0);
    tickPlotLength = get(hPlot1(1,1),'TickLength');
    freqRanges = [1 2];% alpha and gamma
    
elseif strcmp(stimType,'flicker')
    hFig1 = figure; colormap(colorMap)
    set(hFig1,'units','normalized','outerPosition',[0 0 1 1]);
    hPlot1 = getPlotHandles(3,3,[0.07 0.07, 0.5 0.85],0.02,0.04,0);
    hPlot2 = getPlotHandles(1,2,[0.63 0.5, 0.35 0.3],0.02,0.07,0);
    tickPlotLength = get(hPlot1(1,1),'TickLength');
    freqRanges = [1 2 3];  % alpha, gamma and SSVEP
end

% Electrodes to mark in topoplot
showOccipitalElecsUnipolarLeft = [24 29 57 61]; showOccipitalElecsUnipolarRight = [26 31 58 63];
showOccipitalElecsUnipolar = [showOccipitalElecsUnipolarLeft showOccipitalElecsUnipolarRight];
showFrontalElecsUnipolarLeft = [1 33 34 3 37 4]; showFrontalElecsUnipolarRight = [2 35 36 6 40 7];
showFrontalElecsUnipolar = [showFrontalElecsUnipolarLeft showFrontalElecsUnipolarRight];
showElecIDs =[showOccipitalElecsUnipolar showFrontalElecsUnipolar];

% Get the electrode list
capLayout = {'actiCap64'};
cL_Unipolar = load(fullfile('C:\Users\RayLabPC-Aritra\Dropbox\Lab Workbench','Programs\ProgramsMAP','Montages',...
    'Layouts',capLayout{1},[capLayout{1} '.mat']));
chanlocs = cL_Unipolar.chanlocs;
fontSize = 14;
colors = {'k','r','c'}; % Colors for Barplots
elecGroups = {'PO','F','CP','FC','T'}; % PO: Parieto-Occipital; F: Frontal; CP: Centro-Parietal, FC: Frontro-Central; T: Temporal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% static stimuli %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(stimType,'static') % Static gratings -- Generates Figure 5
    att_Condition=1; ign_Condition=2; % LeftVsRight Topoplot
    % Getting Topoplot data
    [attDataTG,ignDataTG,ignDataBL]= getAttendVsIgnored_TopoPlotPowerData_staticStim(powerData_AllElecs,freqRanges,att_Condition,ign_Condition,subjectIdx,nanFlag)  ;
    
    for iFreq= 1:length(freqRanges)
        switch iFreq
            case 1 
                cLims{1} = [-2 2]; % range in dB
                cLimsDiff{1} = [-0.5 2]; % range in dB
            case 2  
                cLims{2} = [-2 2]; % range in dB
                cLimsDiff{2} = [-1 2]; % range in dB
        end
        cLim = cLims{iFreq}; cLimDiff = cLimsDiff{iFreq};
        topoPlot_Attended = 10*(attDataTG{iFreq} - ignDataBL{iFreq});
        topoPlot_Ignored = 10*(ignDataTG{iFreq} - ignDataBL{iFreq});
        topoPlot_AttendMinusIgnored = topoPlot_Attended -  topoPlot_Ignored;
        
        subplot(hPlot1(iFreq,1)); cla; hold on;
        topoplot_murty(topoPlot_Attended,chanlocs,'electrodes','on','style',...
            topoplot_style,'drawaxis','off','nosedir','+X','emarkercolors',...
            topoPlot_Attended); caxis(cLim);
        cBar_Att = colorbar; cTicks = [cLim(1) 0 cLim(2)];
        set(cBar_Att,'Ticks',cTicks,'tickLength',3*tickPlotLength(1),'TickDir','out','fontSize',fontSize)
        topoplot_murty([],chanlocs,'electrodes','on','style','blank',...
            'drawaxis','off','nosedir','+X','plotchans',showElecIDs);
        
        subplot(hPlot1(iFreq,2)); cla; hold on;
        topoplot_murty(topoPlot_Ignored,chanlocs,'electrodes','on','style',...
            topoplot_style,'drawaxis','off','nosedir','+X','emarkercolors',...
            topoPlot_Ignored); caxis(cLim); cTicks = [cLim(1) 0 cLim(2)];
        cBar_Ign = colorbar; set(cBar_Ign,'Ticks',cTicks,'tickLength',3*tickPlotLength(1),'TickDir','out','fontSize',12);
        topoplot_murty([],chanlocs,'electrodes','on','style','blank',...
            'drawaxis','off','nosedir','+X','plotchans',showElecIDs);
        
        subplot(hPlot1(iFreq,3)); cla; hold on;
        topoplot_murty(topoPlot_AttendMinusIgnored,chanlocs,...
            'electrodes','on','style',topoplot_style,'drawaxis','off',...
            'nosedir','+X','emarkercolors',topoPlot_AttendMinusIgnored);
        caxis(cLimDiff); cTicks = [cLimDiff(1) 0 cLimDiff(2)];
        cBar_Diff = colorbar; set(cBar_Diff,'Ticks',cTicks,'tickLength',3*tickPlotLength(1),'TickDir','out','fontSize',12);
        topoplot_murty([],chanlocs,'electrodes','on','style','blank',...
            'drawaxis','off','nosedir','+X','plotchans',showElecIDs);
        if iFreq == 2
            cBar_Diff.Label.String ='\Delta Power (dB)'; cBar_Diff.Label.FontSize = 12;
        end
    end
    
    % Processing & plotting data for barplot in frontal and occipital electrodes
    for iElecGroup =1:2 % 1: Occipital , 2: Frontal
        clear elecsLeft elecRight
        [elecsLeft,elecsRight] = getElectrodeGroupInformation(elecGroups{iElecGroup});
        [attDataTG_barPlot,ignDataTG_barPlot,ignDataBL_barPlot] = getAttendVsIgnored_BarPlotData_StaticStimuli(powerData_AllElecs,elecsLeft,elecsRight,subjectIdx,nanFlag);

        for iBar=1:length(freqRanges)
            attData = (attDataTG_barPlot{iBar}- ignDataBL_barPlot{iBar});
            ignData = (ignDataTG_barPlot{iBar}- ignDataBL_barPlot{iBar});
            diffData = 10*(attData-ignData);
            mBar = mean(diffData,1,nanFlag); % mean across subjects
            errorBar = std(diffData,[],1,nanFlag)./sqrt(length(diffData));
            mBars(iBar) = mBar; %#ok<*AGROW>
            eBars(iBar) = errorBar;
            subplot(hPlot2(1,iElecGroup));hold(hPlot2(1,iElecGroup),'on');
            barPlot = bar(iBar,mBar);
            barPlot.FaceColor = colors{iBar};
            ylim(hPlot2(1,iElecGroup),[-1 1]);
        end
        errorbar(hPlot2(1,iElecGroup),1:length(mBars),mBars,eBars,'.','color','k');
    end
    
    % Plot Titles
    annotation('textbox',[0.08 0.97 0.1 0.0241],'EdgeColor','none','String','Attend Left','fontSize',14,'fontWeight','bold');
    annotation('textbox',[0.13+ 0.12 0.97 0.1 0.0241],'EdgeColor','none','String','Attend Right','fontSize',14,'fontWeight','bold');
    annotation('textbox',[0.29+ 0.1 0.97 0.3 0.0241],'EdgeColor','none','String','Attend Left - Attend Right','fontSize',14,'fontWeight','bold');
    
    % Neural Measure Labels
    annotation('textbox',[0.001 0.82 0.08 0.0241],'EdgeColor','none','HorizontalAlignment','center','String',{'Alpha' '(8-12 Hz)'},'fontSize',14);
    annotation('textbox',[0.001 0.82-0.23 0.08 0.0241],'EdgeColor','none','HorizontalAlignment','center','String',{'Gamma' '(25-70 Hz)'},'fontSize',14);
    
    % Drawing ellipses indicating Attended Locations (filled) and ignored locations (blank) 
    for iAttLoc=1:2
    plotStimDisks_Static(hPlot1(1,iAttLoc),iAttLoc)
    plotStimDisks_Static(hPlot1(2,iAttLoc),iAttLoc)
    end
    
    % Formatting axes properties of barplots
    Datalabels = {'alpha','gamma'};
    set(hPlot2(1,1),'xTick',1:2,'xTickLabel',Datalabels,'XTickLabelRotation',30,'yTick',-1:0.5:1,'fontSize',fontSize,'box','off','tickLength',3*tickPlotLength,'TickDir','out')
    ylabel(hPlot2(1,1),'Change in Power (dB)')
    set(hPlot2(1,2),'xTick',1:2,'xTickLabel',Datalabels,'XTickLabelRotation',30,'yTick',-1:0.5:1,'yTickLabel',[],'fontSize',fontSize,'box','off','tickLength',3*tickPlotLength,'TickDir','out')
    title(hPlot2(1,1),{'Parieto-Occipital' 'Electrodes'})
    title(hPlot2(1,2),{'Frontal' 'Electrodes'})

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% flicker stimuli %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(stimType,'flicker')
    for iFreq=1:length(freqRanges)
        switch iFreq
            case 1; neuralMeasure ='alpha';  cLimsRaw = [-2 2];  cLimsDiff = [-1 2];
            case 2; neuralMeasure = 'gamma'; cLimsRaw = [-2 2];  cLimsDiff = [-1 2];
            case 3; neuralMeasure = 'SSVEP'; cLimsRaw = [-1 2];  cLimsDiff = [0 1] ;
        end
        %getting data for topoplots for flicker stimuli
        [attDataTG_allCond,ignDataTG_allCond,ignDataBL_allCond]= getAttendVsIgnored_TopoPlotPowerData_FlickerStim(powerData_AllElecs,neuralMeasure,iFreq,subjectIdx,nanFlag);
        attDataTG{iFreq}= attDataTG_allCond;
        ignDataTG{iFreq}= ignDataTG_allCond;
        ignDataBL{iFreq}= ignDataBL_allCond;
        
        % plotting topoplots for flicker stimuli
        topoPlot_Attended =  10*(attDataTG{iFreq}-ignDataBL{iFreq});
        topoPlot_Ignored  =  10*(ignDataTG{iFreq}-ignDataBL{iFreq});
        topoPlot_AttendedMinusIgnored = topoPlot_Attended-topoPlot_Ignored;
        
        subplot(hPlot1(iFreq,1)); cla; hold on;
        topoplot_murty(topoPlot_Attended,chanlocs,'electrodes','on',...
            'style',topoplot_style,'drawaxis','off','nosedir','+X',...
            'emarkercolors',topoPlot_Attended);
        caxis(cLimsRaw);    cBar_Att = colorbar; cTicks = [cLimsRaw(1) 0 cLimsRaw(2)];
        set(cBar_Att,'Ticks',cTicks,'tickLength',3*tickPlotLength(1),'TickDir','out','fontSize',fontSize);
        
        topoplot_murty([],chanlocs,'electrodes','on',...
            'style','blank','drawaxis','off','nosedir','+X','plotchans',showElecIDs);
        
        subplot(hPlot1(iFreq,2)); cla; hold on;
        topoplot_murty(topoPlot_Ignored,chanlocs,'electrodes','on','style',topoplot_style,'drawaxis','off','nosedir','+X','emarkercolors',topoPlot_Ignored);
        caxis(cLimsRaw); cBar_Ign = colorbar; cTicks = [cLimsRaw(1) 0 cLimsRaw(2)];
        set(cBar_Ign,'Ticks',cTicks,'tickLength',3*tickPlotLength(1),'TickDir','out','fontSize',fontSize);
        
        topoplot_murty([],chanlocs,'electrodes','on','style','blank','drawaxis','off','nosedir','+X','plotchans',showElecIDs);
        
        subplot(hPlot1(iFreq,3)); cla; hold on;
        topoplot_murty(topoPlot_AttendedMinusIgnored,chanlocs,...
            'electrodes','on','style',topoplot_style,'drawaxis','off',...
            'nosedir','+X','emarkercolors',topoPlot_AttendedMinusIgnored);
        caxis(cLimsDiff);   cBar_Diff = colorbar; cTicks = [unique([cLimsDiff(1) 0]) cLimsDiff(2)];
        set(cBar_Diff,'Ticks',cTicks,'tickLength',3*tickPlotLength(1),'TickDir','out','fontSize',fontSize);
        topoplot_murty([],chanlocs,'electrodes','on',...
            'style','blank','drawaxis','off','nosedir','+X','plotchans',showElecIDs);
        if iFreq == 3
            cBar_Diff.Label.String ='\Delta Power (dB)'; cBar_Diff.Label.FontSize = 14;
        end
        
        % Getting data for bar plots for flicker stimuli for 
        % parieto-occipital and Frontal elec groups
        for iElecGroup=1:2
            clear elecsLeft elecRight
            [elecsLeft,elecsRight] = getElectrodeGroupInformation(elecGroups{iElecGroup});
            [attDataTG_allCond_barPlot,ignDataTG_allCond_barPlot,ignDataBL_allCond_barPlot]= getAttendVsIgnored_BarPlotData_FlickerStimuli(powerData_AllElecs,neuralMeasure,iFreq,elecsLeft,elecsRight,subjectIdx,nanFlag);
            attDataTG_barPlot{iElecGroup,iFreq}= attDataTG_allCond_barPlot;
            ignDataTG_barPlot{iElecGroup,iFreq}= ignDataTG_allCond_barPlot;
            ignDataBL_barPlot{iElecGroup,iFreq}= ignDataBL_allCond_barPlot;
        end
    end
    % Plotting bar plot data for 2 electrode groups
    for iElecGroup=1:2  % 1: Parieto-Occipital 2: Frontal
        for iBar=1:length(freqRanges)
            attData = (attDataTG_barPlot{iElecGroup,iBar}- ignDataBL_barPlot{iElecGroup,iBar});
            ignData = (ignDataTG_barPlot{iElecGroup,iBar}- ignDataBL_barPlot{iElecGroup,iBar});
            diffData = 10*(attData-ignData);
            mBar = mean(diffData,1,nanFlag);
            errorBar = std(diffData,[],1,nanFlag)./sqrt(length(diffData));
            mBars(iBar) = mBar; %#ok<*AGROW>
            eBars(iBar) = errorBar;
            subplot(hPlot2(1,iElecGroup));hold(hPlot2(1,iElecGroup),'on');
            barPlot = bar(iBar,mBar);
            barPlot.FaceColor = colors{iBar};
            ylim(hPlot2(1,iElecGroup),[-2 2]);
        end
        errorbar(hPlot2(1,iElecGroup),1:length(mBars),mBars,eBars,'.','color','k');
    end
    % Stim TF: Left 12 Hz; Right 16 Hz; Attended Left: 12 Hz; Ignored Left: 12 Hz;
    attLoc{1,1} = 1;
    attLoc{1,2} = 2;
    attLoc{2,1} = 1;
    attLoc{2,2} = 2;
    attLoc{3,1} = 1;
    attLoc{3,2} = 2;
    
    for iRow =1:3
        for iColumn = 1:2
            plotStimDisks_Flicker(hPlot1(iRow,iColumn),attLoc{iRow,iColumn})
        end
    end
    
    % Titles for topoplots
    annotation('textbox',[0.09 0.97 0.1 0.0241],'EdgeColor','none','String','Attend Left','fontSize',14,'fontWeight','bold');
    annotation('textbox',[0.15+ 0.11 0.97 0.1 0.0241],'EdgeColor','none','String','Attend Right','fontSize',14,'fontWeight','bold');
    annotation('textbox',[0.31+ 0.11 0.97 0.3 0.0241],'EdgeColor','none','String','Attend Left - Attend Right','fontSize',14,'fontWeight','bold');
   
    % Labeling Neural Measures for Topoplots
    stringLabels1 = {'alpha' '(8-12 Hz)'};
    stringLabels2 = {'gamma' '(25-70 Hz)'};
    stringLabels3 = {'SSVEP' 'Trial-by-trial'};
    annotation('textbox',[0.001 0.81 0.08 0.0241],'EdgeColor','none','HorizontalAlignment','center','String',stringLabels1,'fontSize',14);
    annotation('textbox',[0.001 0.8-0.29 0.08 0.0241],'EdgeColor','none','HorizontalAlignment','center','String',stringLabels2,'fontSize',14);
    annotation('textbox',[0.001 0.8-2*0.295 0.08 0.0241],'EdgeColor','none','HorizontalAlignment','center','String',stringLabels3,'fontSize',14);
    
    % Formatting axes properties of barplots
    Datalabels = {'alpha','gamma','SSVEP'};
    set(hPlot2(1,1),'xTick',1:3,'xTickLabel',Datalabels,'XTickLabelRotation',30,'yTick',-2:1:2,'fontSize',fontSize,'box','off','tickLength',3*tickPlotLength,'TickDir','out')
    ylabel(hPlot2(1,1),'Change in Power (dB)')
    set(hPlot2(1,2),'xTick',1:3,'xTickLabel',Datalabels,'XTickLabelRotation',30,'yTick',-2:1:2,'yTickLabel',[],'fontSize',fontSize,'box','off','tickLength',3*tickPlotLength,'TickDir','out')
    title(hPlot2(1,1),{'Parieto-Occipital' 'Electrodes'})
    title(hPlot2(1,2),{'Frontal' 'Electrodes'})
end
end

% Processes TopoPower Data for static stimuli
function [attDataTG,ignDataTG,ignDataBL]=getAttendVsIgnored_TopoPlotPowerData_staticStim(powerData,freqRanges,att_Condition,ign_Condition,subjectIdx,nanFlag)
% alpha and gamma
for iSub=1:length(subjectIdx)
    clear topoDataTG topoDataBL
    topoDataTG = powerData{subjectIdx(iSub)}.powerValsTG;
    topoDataBL = powerData{subjectIdx(iSub)}.powerValsBL;
    %mean across trials for attend and ignored condition
    attDataTGTMP(iSub,:,:,:,:) = squeeze(mean(topoDataTG{att_Condition},2,nanFlag));
    ignDataTGTMP(iSub,:,:,:,:) = squeeze(mean(topoDataTG{ign_Condition},2,nanFlag));
    ignDataBLTMP(iSub,:,:,:,:) = squeeze(mean(topoDataBL{ign_Condition},2,nanFlag));
    
end
attDataTG = cell(1,length(freqRanges));
ignDataTG = cell(1,length(freqRanges));
ignDataBL = cell(1,length(freqRanges));

% Saves topoplot data for alpha and gamma powerVals after 
% converting to individual subject data in log10 scale and 
% subsequently averaged across subjects
for iFreq = 1:length(freqRanges)
    attDataTG{iFreq} = squeeze(mean(log10(attDataTGTMP(:,:,iFreq)),1,nanFlag));
    ignDataTG{iFreq} = squeeze(mean(log10(ignDataTGTMP(:,:,iFreq)),1,nanFlag));
    ignDataBL{iFreq} = squeeze(mean(log10(ignDataBLTMP(:,:,iFreq)),1,nanFlag));
end
end

% Processes BarPlot Data for static stimuli
function [attDataTG,ignDataTG,ignDataBL] = getAttendVsIgnored_BarPlotData_StaticStimuli(powerData,elecsLeft,elecsRight,subjectIdx,nanFlag)
freqRanges =[1 2] ; % alpha and gamma
for iSub=1:length(subjectIdx)
    clear topoDataTG topoDataBL
    topoDataTG = powerData{subjectIdx(iSub)}.powerValsTG;
    topoDataBL = powerData{subjectIdx(iSub)}.powerValsBL;
    
    for iCount=1:2
        switch iCount
            case 1 ; att_Condition = 2;  ign_Condition = 1; elecs = elecsLeft;
            case 2 ; att_Condition = 1;  ign_Condition = 2; elecs = elecsRight;
        end
        % Mean across trials for attend and ignored condition
        attDataTGTMP(iSub,iCount,:,:) = squeeze(mean(topoDataTG{att_Condition}(elecs,:,:),2,nanFlag));
        ignDataTGTMP(iSub,iCount,:,:) = squeeze(mean(topoDataTG{ign_Condition}(elecs,:,:),2,nanFlag));
        ignDataBLTMP(iSub,iCount,:,:) = squeeze(mean(topoDataBL{ign_Condition}(elecs,:,:),2,nanFlag));
    end
end
% Saves power data for alpha and gamma powerVals after taking average across 
% both attend-Left and attend-Right conditions for individual subject
% and subsequently converting to log10 values
for iFreq= 1:length(freqRanges)
    attDataTG{iFreq} = squeeze(log10(mean(mean(attDataTGTMP(:,:,:,iFreq),2,nanFlag),3,nanFlag)));
    ignDataTG{iFreq} = squeeze(log10(mean(mean(ignDataTGTMP(:,:,:,iFreq),2,nanFlag),3,nanFlag)));
    ignDataBL{iFreq} = squeeze(log10(mean(mean(ignDataBLTMP(:,:,:,iFreq),2,nanFlag),3,nanFlag)));
end
end


% Processes TopoPower Data for flicker stimuli
function [attDataTG,ignDataTG,ignDataBL]=getAttendVsIgnored_TopoPlotPowerData_FlickerStim(powerData,neuralMeasure,freqRange,subjectIdx,nanFlag)
SSVEPFreq = [3 4]; % SSVEPFreq; 3- 24 Hz; 4- 32 Hz

% Averaged across 12 & 16 Hz flicker conditions (Attend-Left & Attend-Right)
% and log-averaged over all subjects for alpha and gamma power 
if strcmp(neuralMeasure,'alpha')||strcmp(neuralMeasure,'gamma')
    for iSub = 1:length(subjectIdx)
        clear topoDataTG topoDataBL
        topoDataTG = powerData{subjectIdx(iSub)}.powerValsTG;
        topoDataBL = powerData{subjectIdx(iSub)}.powerValsBL;
        
        for iCount = 1:2 % Attend-Left Conditions averaged over cued 12 Hz and 16 Hz Counterphase Flickering conditions
            switch iCount
                case 1; att_Condition = 3; ign_Condition =6; % 3-H0V_12Hz; 6-H1V_16Hz
                case 2; att_Condition = 5; ign_Condition =4; % 5-H0V_16Hz; 4-H1V_12Hz
            end
            % Averaged across trials
            attDataTMP_perSub = squeeze(mean(topoDataTG{att_Condition},2,nanFlag));
            ignDataTMP_perSub = squeeze(mean(topoDataTG{ign_Condition},2,nanFlag));
            ignDataBLTMP_perSub = squeeze(mean(topoDataBL{ign_Condition},2,nanFlag));
            % Saves the powerVals for particular frequency range
            attDataTMP(iSub,iCount,:) = attDataTMP_perSub(:,freqRange)';
            ignDataTMP(iSub,iCount,:) = ignDataTMP_perSub(:,freqRange)';
            ignDataBLTMP(iSub,iCount,:) = ignDataBLTMP_perSub(:,freqRange)';
        end
    end
    
    attDataTG = squeeze(mean(log10(mean(attDataTMP,2,nanFlag)),1,nanFlag))';
    ignDataTG = squeeze(mean(log10(mean(ignDataTMP,2,nanFlag)),1,nanFlag))';
    ignDataBL = squeeze(mean(log10(mean(ignDataBLTMP,2,nanFlag)),1,nanFlag))';
    
elseif strcmp(neuralMeasure,'SSVEP')
    % Averaged across all four stimulus and Attention conditions 
    % and log-averaged over all subjects with AttR conditions mirrored with z-line
    for iSub = 1:length(subjectIdx)
        clear topoDataTG topoDataBL
        topoDataTG = powerData{subjectIdx(iSub)}.powerValsTG;
        topoDataBL = powerData{subjectIdx(iSub)}.powerValsBL;
        for iCount = 1:4
            switch iCount
                case 1; att_Condition = 3; ign_Condition =6; SSVEPFreq = 3; % 24 Hz
                case 2; att_Condition = 4; ign_Condition =5; SSVEPFreq = 3; % 24 Hz
                case 3; att_Condition = 5; ign_Condition =4; SSVEPFreq = 4; % 32 Hz
                case 4; att_Condition = 6; ign_Condition =3; SSVEPFreq = 4; % 32 Hz
            end
            % Averaged across trials
            attDataTMP_perSub = squeeze(mean(topoDataTG{att_Condition},2,nanFlag));
            ignDataTMP_perSub = squeeze(mean(topoDataTG{ign_Condition},2,nanFlag));
            ignDataBLTMP_perSub = squeeze(mean(topoDataBL{ign_Condition},2,nanFlag));
            % Saves the powerVals for particular frequency range
            if iCount==2 || iCount==4
                attDataTMP(iSub,iCount,:) = mirrorTopoplotData(attDataTMP_perSub(:,SSVEPFreq)','');
                ignDataTMP(iSub,iCount,:) = mirrorTopoplotData(ignDataTMP_perSub(:,SSVEPFreq)','');
                ignDataBLTMP(iSub,iCount,:) = mirrorTopoplotData(ignDataBLTMP_perSub(:,SSVEPFreq)','');
            else
                attDataTMP(iSub,iCount,:) = attDataTMP_perSub(:,SSVEPFreq)';
                ignDataTMP(iSub,iCount,:) = ignDataTMP_perSub(:,SSVEPFreq)';
                ignDataBLTMP(iSub,iCount,:) = ignDataBLTMP_perSub(:,SSVEPFreq)';
            end
        end
    end
    % Saving topo power data for 64 elecs after taking average across 
    % different stimulus/attention conditions, converting to log and 
    % then taking average across subjects
    attDataTG = squeeze(mean(log10(mean(attDataTMP,2,nanFlag)),1,nanFlag))';
    ignDataTG = squeeze(mean(log10(mean(ignDataTMP,2,nanFlag)),1,nanFlag))';
    ignDataBL = squeeze(mean(log10(mean(ignDataBLTMP,2,nanFlag)),1,nanFlag))';
end
end

% Processes BarPlot Data for flicker stimuli
function [attDataTG,ignDataTG,ignDataBL] = getAttendVsIgnored_BarPlotData_FlickerStimuli(powerData,neuralMeasure,iFreq,elecsLeft,elecsRight,subjectIdx,nanFlag)
if strcmp(neuralMeasure,'alpha')||strcmp(neuralMeasure,'gamma')
    for iSub=1:length(subjectIdx)
        clear dataTG dataBL
        dataTG = powerData{subjectIdx(iSub)}.powerValsTG;
        dataBL = powerData{subjectIdx(iSub)}.powerValsBL;
        for iCondition=1:4
            switch iCondition
                case 1; att_Condition = 3; ign_Condition =6; elecs = elecsRight;
                case 2; att_Condition = 4; ign_Condition =5; elecs = elecsLeft;
                case 3; att_Condition = 5; ign_Condition =4; elecs = elecsRight;
                case 4; att_Condition = 6; ign_Condition =3; elecs = elecsLeft;
            end
            % Averaging across trials for each electrode group
            attDataTMP_perSub = squeeze(mean(dataTG{att_Condition}(elecs,:,:),2,nanFlag));
            ignDataTMP_perSub = squeeze(mean(dataTG{ign_Condition}(elecs,:,:),2,nanFlag));
            ignDataBLTMP_perSub = squeeze(mean(dataBL{ign_Condition}(elecs,:,:),2,nanFlag));
            % Aaves the powerVals for particular frequency range
            attDataTMP(iSub,:,iCondition) = attDataTMP_perSub(:,iFreq);
            ignDataTMP(iSub,:,iCondition) = ignDataTMP_perSub(:,iFreq);
            ignDataBLTMP(iSub,:,iCondition) = ignDataBLTMP_perSub(:,iFreq);
        end
    end
    
    % Saves the powerData for bar plot by taking average across the
    % conditions and subsequently log-averaging across subjects
    attDataTG = squeeze(log10(mean(mean(attDataTMP,3,nanFlag),2,nanFlag)));
    ignDataTG = squeeze(log10(mean(mean(ignDataTMP,3,nanFlag),2,nanFlag)));
    ignDataBL = squeeze(log10(mean(mean(ignDataBLTMP,3,nanFlag),2,nanFlag)));
    
elseif strcmp(neuralMeasure,'SSVEP')
    
    for iSub=1:length(subjectIdx)
        clear topoDataTG topoDataBL
        dataTG = powerData{subjectIdx(iSub)}.powerValsTG;
        dataBL = powerData{subjectIdx(iSub)}.powerValsBL;
        for iCondition=1:4
            switch iCondition
                case 1; att_Condition = 3; ign_Condition =6; elecs = elecsRight; SSVEPFreq = 3; %24 Hz
                case 2; att_Condition = 4; ign_Condition =5; elecs = elecsLeft; SSVEPFreq = 3; %24 Hz
                case 3; att_Condition = 5; ign_Condition =4; elecs = elecsRight; SSVEPFreq = 4; %32 Hz
                case 4; att_Condition = 6; ign_Condition =3; elecs = elecsLeft; SSVEPFreq = 4; %32 Hz
            end
            % Averaged across trials for each electrode group
            attDataTMP_perSub = squeeze(mean(dataTG{att_Condition}(elecs,:,:),2,nanFlag));
            ignDataTMP_perSub = squeeze(mean(dataTG{ign_Condition}(elecs,:,:),2,nanFlag));
            ignDataBLTMP_perSub = squeeze(mean(dataBL{ign_Condition}(elecs,:,:),2,nanFlag));
            
            % Saves the powerVals for particular frequency range
            attDataTMP(iSub,:,iCondition) = attDataTMP_perSub(:,SSVEPFreq);
            ignDataTMP(iSub,:,iCondition) = ignDataTMP_perSub(:,SSVEPFreq);
            ignDataBLTMP(iSub,:,iCondition) = ignDataBLTMP_perSub(:,SSVEPFreq);
        end
    end
    attDataTG = squeeze(log10(mean(mean(attDataTMP,3,nanFlag),2,nanFlag)));
    ignDataTG = squeeze(log10(mean(mean(ignDataTMP,3,nanFlag),2,nanFlag)));
    ignDataBL = squeeze(log10(mean(mean(ignDataBLTMP,3,nanFlag),2,nanFlag)));
end
end

function [elecsLeft,elecsRight] = getElectrodeGroupInformation(elecGroup)

if strcmp(elecGroup,'PO')
    elecsLeft = [24 29 57 61];
    elecsRight = [26 31 58 63];
elseif strcmp(elecGroup,'F')
    elecsLeft = [1 33 34 3 37 4];
    elecsRight = [2 35 36 6 40 7];
elseif strcmp(elecGroup,'CP')
    elecsLeft = [18 19 52 23 56];
    elecsRight = [20 21 27 54 59];
elseif strcmp(elecGroup,'FC')
    elecsLeft = [8 9 13 43 47 48];
    elecsRight = [10 11 15 44 49 50];
elseif strcmp(elecGroup,'T')
    elecsLeft = [12 17 41 42 51];
    elecsRight = [16 22 45 46 55];
end
end

function mirrored_topoData = mirrorTopoplotData(data,refType)
if ~exist(refType,'var')
    refType = 'unipolar';
end
if strcmp(refType,'unipolar')
    mirror_elecNums = [2 1	7	6	5	4	3	11	10	9	8	16	15	14 ...
        13	12	22	21	20	19	18	17	27	26	25	24	23	32	31	30	29 ...
        28	36	35	34	33	40	39	38	37	46	45	44	43	42	41	50	49 ...
        48	47	55	54	53	52	51	59	58	57	56	64	63	62	61	60];
elseif strcmp(refType,'bipolar')
    mirror_elecNums = [ 4	3	2	1	6	5	14	13	12	11	10	9	8	7 ...
        22	21	20	19	18	17	16	15	30	29	28	27	26	25	24	23	38	...
        37	36	35	34	33	32	31	46	45	44	43	42	41	40	39	54	53	...
        52	51	50	49	48	47	63	62	61	60	59	58	57	56	55	73	72	...
        71	70	69	68	67	66	65	64	82	81	80	79	78	77	76	75	74	...
        90	89	88	87	86	85	84	83	99	98	97	96	95	94	93	92	91	...
        103	102	101	100	110	109	108	107	106	105	104	112	111];
end
mirrored_topoData = data(:,mirror_elecNums,:,:,:);
end

function plotStimDisks_Static(hPlot,attLoc)
stimDiskDistanceFromMidline = 0.01;
% textStartPosGapFromMidline = 0.001;
ellipseYGap = 0.17;
ellipseWidth = 0.015;
ellipseHeight = 0.012;
% textYGap = 0.19;
% textWidth = 0.04;
% textHeight = 0.0241;


AttendPlotPos = get(hPlot,'Position');
AttendPlotMidline = AttendPlotPos(1)+ AttendPlotPos(3)/2;
elpsL = annotation('ellipse',[AttendPlotMidline-ellipseWidth- stimDiskDistanceFromMidline AttendPlotPos(2)+ellipseYGap ellipseWidth ellipseHeight],'units','normalized');
elpsR = annotation('ellipse',[AttendPlotMidline + stimDiskDistanceFromMidline AttendPlotPos(2)+ellipseYGap ellipseWidth ellipseHeight]);

if attLoc==1
    elpsL.FaceColor = 'k'; elpsR.FaceColor = 'none';
elseif attLoc==2
    elpsL.FaceColor = 'none'; elpsR.FaceColor = 'k';
end

% if ~isempty(ssvepFreqs)
%     annotation('textbox',[AttendPlotMidline-(textWidth+textStartPosGapFromMidline) AttendPlotPos(2)+textYGap textWidth textHeight],...
%         'EdgeColor','none','String',[num2str(ssvepFreqs(1)) ' Hz'],'fontSize',10,'EdgeColor','none','FitBoxToText','on',...
%         'HorizontalAlignment','center');
%
%     annotation('textbox',[AttendPlotMidline+textStartPosGapFromMidline AttendPlotPos(2)+textYGap textWidth textHeight],...
%         'EdgeColor','none','String',[num2str(ssvepFreqs(2)) ' Hz'],'fontSize',10,...
%         'EdgeColor','none','FitBoxToText','on','HorizontalAlignment','center');
% end
end

function plotStimDisks_Flicker(hPlot,attLoc)
stimDiskDistanceFromMidline = 0.01;
% textStartPosGapFromMidline = 0.001;
ellipseYGap = 0.25;
ellipseWidth = 0.02;
ellipseHeight = 0.016;
% textYGap = 0.19;
% textWidth = 0.04;
% textHeight = 0.0241;


AttendPlotPos = get(hPlot,'Position');
AttendPlotMidline = AttendPlotPos(1)+ AttendPlotPos(3)/2;
elpsL = annotation('ellipse',[AttendPlotMidline-ellipseWidth- stimDiskDistanceFromMidline AttendPlotPos(2)+ellipseYGap ellipseWidth ellipseHeight],'units','normalized');
elpsR = annotation('ellipse',[AttendPlotMidline + stimDiskDistanceFromMidline AttendPlotPos(2)+ellipseYGap ellipseWidth ellipseHeight]);

if attLoc==1
    elpsL.FaceColor = 'k'; elpsR.FaceColor = 'none';
elseif attLoc==2
    elpsL.FaceColor = 'none'; elpsR.FaceColor = 'k';
end

% if ~isempty(ssvepFreqs)
%     annotation('textbox',[AttendPlotMidline-(textWidth+textStartPosGapFromMidline) AttendPlotPos(2)+textYGap textWidth textHeight],...
%         'EdgeColor','none','String',[num2str(ssvepFreqs(1)) ' Hz'],'fontSize',10,'EdgeColor','none','FitBoxToText','on',...
%         'HorizontalAlignment','center');
%
%     annotation('textbox',[AttendPlotMidline+textStartPosGapFromMidline AttendPlotPos(2)+textYGap textWidth textHeight],...
%         'EdgeColor','none','String',[num2str(ssvepFreqs(2)) ' Hz'],'fontSize',10,...
%         'EdgeColor','none','FitBoxToText','on','HorizontalAlignment','center');
% end
end