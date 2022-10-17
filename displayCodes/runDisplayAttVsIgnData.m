% This Program is the executive code that displays the Attended versus
% Ignored topoplot and bar plot Results for attention EEG dataset 
% collected from human individuals

clear; % clears the Workspace
close all; % closes any open Figures

% Set FolderSourceString automatically; add files & folders to MATLAB path
fileNameWithPath = mfilename('fullpath');
[filePath,name,ext] = fileparts(fileNameWithPath);
cd(filePath); cd ..
addpath(genpath(pwd));
folderSourceString = fullfile(pwd,'analyzedData');

subjectIdx = 1:26; % Select more than one Subject. Numbers 1-26 are subject IDs.
badTrialStr = 'v10'; % bad Trial Version used
colorMap ='jet'; % Topoplot color scheme;
topoplot_style = 'both'; % Draws Contour Maps and Electrodes

stimType= 'static'; % Use either 'static' or 'flicker' to generate respective plots/figures
displayResults_AttVsIgn(folderSourceString,subjectIdx,stimType,badTrialStr,colorMap,topoplot_style)

stimType= 'flicker'; 
displayResults_AttVsIgn(folderSourceString,subjectIdx,stimType,badTrialStr,colorMap,topoplot_style)