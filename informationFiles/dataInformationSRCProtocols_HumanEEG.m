function [subjectNames,expDates,protocolNames,dataFolderSourceString] = dataInformationSRCProtocols_HumanEEG(dataFolderSourceString,gridType,protocolType)

if ~exist('dataFolderSourceString','var');  dataFolderSourceString='';  end
if ~exist('gridType','var');         gridType = 'EEG';                  end
if ~exist('protocolType','var');     protocolType = 'SRC-Long';         end

% FolderSourceString for extracted dataset
if isempty(dataFolderSourceString)
    dataFolderSourceString = 'N:\Projects\Aritra_AttentionEEGProject\data\SRCLong\';
end

[allSubjectNames,allExpDates,allProtocolNames,~,~,~] = eval(['allProtocolsAttention' gridType 'Project']);

% stimList = cell(1,4);
if strcmp(protocolType, 'SRC-Long_Psychophysics')
    protocolList{1} =  247;  % M
    protocolList{2} =  249;  % F
    protocolList{3} =  251;  % M
    protocolList{4} =  260;  % M
    protocolList{5} =  264;  % M
    protocolList{6} =  266;  % M
    protocolList{7} =  271;  % F
    protocolList{8} =  276;  % F
    protocolList{9} =  287;  % M
    protocolList{10} = 293;  % M
    protocolList{11} = 299;  % M
    protocolList{12} = 305;  % M
    protocolList{13} = 314;  % F
    protocolList{14} = 323;  % M
    protocolList{15} = 333;  % F
    protocolList{16} = 340;  % M
    protocolList{17} = 348;  % F
    protocolList{18} = 355;  % F
    protocolList{19} = 361;  % F
    protocolList{20} = 382;  % M
    protocolList{21} = 386;  % M
    protocolList{22} = 408;  % M
    protocolList{23} = 414;  % F
    protocolList{24} = 422;  % F
    protocolList{25} = 444;  % M
    protocolList{26} = 459;  % F

elseif strcmp(protocolType, 'SRC-Long')
    protocolList{1} = 248;  
    protocolList{2} = 250;  
    protocolList{3} = 252;  
    protocolList{4} = 263;  
    protocolList{5} = 265;  
    protocolList{6} = 270;  
    protocolList{7} = 273;  
    protocolList{8} = 277;  
    protocolList{9} = 288;  
    protocolList{10} = 294;  
    protocolList{11} = 300;  
    protocolList{12} = 306;  
    protocolList{13} = 317;  
    protocolList{14} = 325;  
    protocolList{15} = 334;  
    protocolList{16} = 342;  
    protocolList{17} = 350;  
    protocolList{18} = 356;  
    protocolList{19} = 362;  
    protocolList{20} = 383;  
    protocolList{21} = 403;  
    protocolList{22} = 409;  
    protocolList{23} = 417;  
    protocolList{24} = 425;  
    protocolList{25} = 446;  
    protocolList{26} = 460;  

end

numSubjects = length(protocolList);
subjectNames = cell(numSubjects,length(protocolList{1}));
expDates = cell(numSubjects,length(protocolList{1}));
protocolNames = cell(numSubjects,length(protocolList{1}));

for i=1:numSubjects
    subjectNames{i} = allSubjectNames{protocolList{i}};
    expDates{i} = allExpDates{protocolList{i}};
    protocolNames{i} = allProtocolNames{protocolList{i}};
end

