function [subjectNames,expDates,protocolNames,stimType,deviceName,capLayout] = allProtocolsAttentionEEGProject

%% Recorded in Rig 4: Stim type 2: Sampling Rate: 1000 Hz, timeStartFromBaseLineList(2) = -1.148; deltaTList(2) = 2.048;
clear index; index=1; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '221015'; protocolNames{index} = 'SRC_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=2; deviceName{index} = 'BP'; subjectNames{index} = 'MD'; expDates{index} = '160116'; protocolNames{index} = 'SRC_0001'; capLayout{index} = 'brainCap64'; stimType{index} = 2;

%% Recorded in RIG 2: Stim type 4: timeStartFromBaseLineList(4) = -0.8192;deltaTList(4)=1.6384; % in seconds- with fs as 2500samp/sec you get 4096 samples with this 
clear index; index=3; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '020416'; protocolNames{index} = 'SRC_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=4; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '020416'; protocolNames{index} = 'SRC_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=5; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '020416'; protocolNames{index} = 'SRC_0003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=6; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '020416'; protocolNames{index} = 'SRC_0004'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=7; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '100416'; protocolNames{index} = 'SRC_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=8; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '100416'; protocolNames{index} = 'SRC_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=9; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '230516'; protocolNames{index} = 'SRC_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=10; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '230516'; protocolNames{index} = 'SRC_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=11; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '230516'; protocolNames{index} = 'SRC_0003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=12; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '260516'; protocolNames{index} = 'GRF_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=13; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '260516'; protocolNames{index} = 'GRF_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=14; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '260516'; protocolNames{index} = 'GRF_0003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=15; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '270516'; protocolNames{index} = 'GRF_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=16; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '270516'; protocolNames{index} = 'GRF_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=17; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '270516'; protocolNames{index} = 'GRF_0003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=18; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '050616'; protocolNames{index} = 'GRF_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=19; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '050616'; protocolNames{index} = 'GRF_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=20; deviceName{index} = 'BP'; subjectNames{index} = 'AK'; expDates{index} = '080616'; protocolNames{index} = 'GRF_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=21; deviceName{index} = 'BP'; subjectNames{index} = 'AK'; expDates{index} = '080616'; protocolNames{index} = 'GRF_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=22; deviceName{index} = 'BP'; subjectNames{index} = 'DJ'; expDates{index} = '140616'; protocolNames{index} = 'GRF_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=23; deviceName{index} = 'BP'; subjectNames{index} = 'DJ'; expDates{index} = '140616'; protocolNames{index} = 'GRF_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=24; deviceName{index} = 'BP'; subjectNames{index} = 'DJ'; expDates{index} = '140616'; protocolNames{index} = 'GRF_0003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=25; deviceName{index} = 'BP'; subjectNames{index} = 'MD'; expDates{index} = '160616'; protocolNames{index} = 'GRF_0001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=26; deviceName{index} = 'BP'; subjectNames{index} = 'MD'; expDates{index} = '160616'; protocolNames{index} = 'GRF_0002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=27; deviceName{index} = 'BP'; subjectNames{index} = 'MD'; expDates{index} = '160616'; protocolNames{index} = 'GRF_0003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=28; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '140916'; protocolNames{index} = 'GRF_0003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

%% Preliminary Data Recorded in RIG 2 in 2021 after Covid-19 lockdown: Stim type 4; RiG 2 shifted back from CNS to PRL: timeStartFromBaseLineList(4) = -0.848;deltaTList(4)=2.048; % in seconds- with fs as 1000 samples/sec you get 2048 samples with this 

clear index; index=29; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '200121'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=30; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '200121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=31; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '220121'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=32; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '220121'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=33; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '280121'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=34; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '280121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=35; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '280121'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=36; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '280121'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=37; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '290121'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=38; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=39; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=40; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=41; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=42; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=43; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=44; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=45; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=46; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=47; deviceName{index} = 'BP'; subjectNames{index} = 'Aritra'; expDates{index} = '120221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=48; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '150221'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=49; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '150221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=50; deviceName{index} = 'BP'; subjectNames{index} = 'test'; expDates{index} = '150221'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=51; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=52; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=53; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=54; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=55; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=56; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=57; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=58; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=59; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=60; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '260221'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=61; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=62; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=63; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=64; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=65; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=66; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=67; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=68; deviceName{index} = 'BP'; subjectNames{index} = 'SVP'; expDates{index} = '030321'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=69; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open  
clear index; index=70; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=71; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=72; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=73; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 2, ele = 2, diameter = 2;
clear index; index=74; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 3, ele = 3, diameter = 3;
clear index; index=75; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 4, ele = 4, diameter = 4;
clear index; index=76; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 5, ele = 5, diameter = 5;
clear index; index=77; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 6, ele = 6, diameter = 6;
clear index; index=78; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 7, ele = 7, diameter = 7;
clear index; index=79; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 8, ele = 8, diameter = 8;
clear index; index=80; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 6.5, ele = 0, diameter = 6.5;
clear index; index=81; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % azi = 8.5, ele = 0, diameter = 8.5;
clear index; index=82; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % 
clear index; index=83; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '120421'; protocolNames{index} = 'GRF_015'; capLayout{index} = 'actiCap64'; stimType{index} = 4; %

clear index; index=84; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '080721'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=85; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '080721'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 
clear index; index=86; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '080721'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 2; 

clear index; index=87; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=88; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=89; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=90; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=91; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=92; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=93; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=94; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=95; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 
clear index; index=96; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120821'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; 

%% Recorded in RIG 2, Phase 1, Common & Mapping Protocols: Stim type 4: timeStartFromBaseLineList(4) = -0.848;deltaTList(4)=2.048; % in seconds- with fs as 1000 samples/sec you get 2048 samples with this 

clear index; index=97; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=98; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=99; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=100; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=101; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=102; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=103; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=104; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=105; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=106; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=107; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=108; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '060921'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8

clear index; index=109; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=110; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=111; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=112; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=113; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=114; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=115; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=116; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=117; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=118; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=119; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=120; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=121; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms
clear index; index=122; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '070921'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=123; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=124; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=125; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=126; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=127; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=128; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=129; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=130; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=131; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=132; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=133; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=134; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=135; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=136; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '160921'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=137; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=138; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=139; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=140; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=141; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=142; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=143; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=144; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=145; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=146; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=147; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6; 1st correct trial was run on without fixation mode; used fixOn.timeMS as fixate.timeMS; bad trials due to eye outside fixation window, if any, will be removed by findBadTrialsWithEEG later
clear index; index=148; deviceName{index} = 'BP'; subjectNames{index} = 'SA'; expDates{index} = '210921'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8

clear index; index=149; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=150; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=151; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=152; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=153; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=154; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=155; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=156; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=157; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=158; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=159; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=160; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=161; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=162; deviceName{index} = 'BP'; subjectNames{index} = 'SB'; expDates{index} = '280921'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=163; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=164; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=165; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=166; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=167; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=168; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=169; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=170; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=171; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=172; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=173; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=174; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=175; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=176; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '290921'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=177; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=178; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=179; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=180; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=181; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=182; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=183; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=184; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=185; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=186; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=187; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=188; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=189; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=190; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '011021'; protocolNames{index} = 'GRF_015'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=191; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=192; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=193; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=194; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=195; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=196; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=197; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=198; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=199; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=200; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=201; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=202; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=203; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=204; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '051021'; protocolNames{index} = 'GRF_015'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=205; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=206; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=207; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=208; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=209; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=210; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=211; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=212; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=213; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=214; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=215; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=216; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=217; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=218; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '061021'; protocolNames{index} = 'GRF_015'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=219; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=220; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=221; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=222; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=223; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=224; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=225; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=226; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=227; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=228; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=229; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=230; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=231; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=232; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '071021'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

clear index; index=233; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=234; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=235; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=236; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Continuous
clear index; index=237; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 2
clear index; index=238; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 4
clear index; index=239; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 6
clear index; index=240; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_008'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Diagonal 8
clear index; index=241; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_009'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 2
clear index; index=242; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_010'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 4
clear index; index=243; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_011'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 6
clear index; index=244; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_012'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Mapping Horizontal 8
clear index; index=245; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_013'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Horizontal 8; Stim: 3000ms, interStim: 1500ms 
clear index; index=246; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '111021'; protocolNames{index} = 'GRF_014'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Mapping Diagonal   8; Stim: 3000ms, interStim: 1500ms

%% Recorded in RIG 2, Phase 2, Attention Protocol SRC-Long: Stim type 7: timeStartFromBaseLineList(7) = -1.5;deltaTList(7)=3.072; % in seconds- with fs as 1000 samples/sec you get 3072 samples with this 
clear index; index=247; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '101121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available

clear index; index=248; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '161121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=249; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '171121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available
clear index; index=250; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '171121'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=251; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '241121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available

clear index; index=252; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '261121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=253; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '261121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available

clear index; index=254; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '291121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available

clear index; index=255; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '301121'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available

clear index; index=256; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '011221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available
clear index; index=257; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '011221'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available
clear index; index=258; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '011221'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available

clear index; index=259; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '011221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior; Only LL File available

clear index; index=260; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '021221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=261; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '021221'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=262; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '021221'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=263; deviceName{index} = 'BP'; subjectNames{index} = 'SSP'; expDates{index} = '021221'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=264; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '031221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=265; deviceName{index} = 'BP'; subjectNames{index} = 'SK'; expDates{index} = '031221'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=266; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '041221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=267; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '041221'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=268; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '041221'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=269; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '041221'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=270; deviceName{index} = 'BP'; subjectNames{index} = 'AB'; expDates{index} = '041221'; protocolNames{index} = 'SRC_005'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=271; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '081221'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=272; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '081221'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=273; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '081221'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=274; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '070122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=275; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '070122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=276; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '070122'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=277; deviceName{index} = 'BP'; subjectNames{index} = 'DG'; expDates{index} = '070122'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=278; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '100122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=279; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '100122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=280; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '100122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=281; deviceName{index} = 'BP'; subjectNames{index} = 'AD'; expDates{index} = '100122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG

clear index; index=282; deviceName{index} = 'BP'; subjectNames{index} = 'SW'; expDates{index} = '120122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=283; deviceName{index} = 'BP'; subjectNames{index} = 'SW'; expDates{index} = '120122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=284; deviceName{index} = 'BP'; subjectNames{index} = 'SW'; expDates{index} = '120122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=285; deviceName{index} = 'BP'; subjectNames{index} = 'SW'; expDates{index} = '120122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=286; deviceName{index} = 'BP'; subjectNames{index} = 'SW'; expDates{index} = '120122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=287; deviceName{index} = 'BP'; subjectNames{index} = 'SW'; expDates{index} = '120122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=288; deviceName{index} = 'BP'; subjectNames{index} = 'SW'; expDates{index} = '120122'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=289; deviceName{index} = 'BP'; subjectNames{index} = 'AR'; expDates{index} = '190122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=290; deviceName{index} = 'BP'; subjectNames{index} = 'AR'; expDates{index} = '190122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=291; deviceName{index} = 'BP'; subjectNames{index} = 'AR'; expDates{index} = '190122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=292; deviceName{index} = 'BP'; subjectNames{index} = 'AR'; expDates{index} = '190122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=293; deviceName{index} = 'BP'; subjectNames{index} = 'AR'; expDates{index} = '190122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=294; deviceName{index} = 'BP'; subjectNames{index} = 'AR'; expDates{index} = '190122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=295; deviceName{index} = 'BP'; subjectNames{index} = 'SN'; expDates{index} = '210122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=296; deviceName{index} = 'BP'; subjectNames{index} = 'SN'; expDates{index} = '210122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=297; deviceName{index} = 'BP'; subjectNames{index} = 'SN'; expDates{index} = '210122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=298; deviceName{index} = 'BP'; subjectNames{index} = 'SN'; expDates{index} = '210122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=299; deviceName{index} = 'BP'; subjectNames{index} = 'SN'; expDates{index} = '210122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=300; deviceName{index} = 'BP'; subjectNames{index} = 'SN'; expDates{index} = '210122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=301; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '230122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=302; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '230122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Closed
clear index; index=303; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '230122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=304; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '230122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI No InterStim
clear index; index=305; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '230122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Psychophysics
clear index; index=306; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '230122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Attention

clear index; index=307; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=308; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=309; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=310; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=311; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=312; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=313; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=314; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=315; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=316; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=317; deviceName{index} = 'BP'; subjectNames{index} = 'SU'; expDates{index} = '260122'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=318; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Eyes Open
clear index; index=319; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % EyesClosed
clear index; index=320; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Aborted
clear index; index=321; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI Aborted
clear index; index=322; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % SF-ORI
clear index; index=323; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Psychophysics
clear index; index=324; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Attention Aborted
clear index; index=325; deviceName{index} = 'BP'; subjectNames{index} = 'ArD'; expDates{index} = '290122'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Attention

clear index; index=326; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=327; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=328; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=329; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=330; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=331; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=332; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'GRF_007'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=333; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=334; deviceName{index} = 'BP'; subjectNames{index} = 'UG'; expDates{index} = '300122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=335; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=336; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=337; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=338; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=339; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=340; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=341; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=342; deviceName{index} = 'BP'; subjectNames{index} = 'SS'; expDates{index} = '310122'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=343; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=344; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=345; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=346; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=347; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=348; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=349; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=350; deviceName{index} = 'BP'; subjectNames{index} = 'SNS'; expDates{index} = '060222'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=351; deviceName{index} = 'BP'; subjectNames{index} = 'SG'; expDates{index} = '080222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=352; deviceName{index} = 'BP'; subjectNames{index} = 'SG'; expDates{index} = '080222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=353; deviceName{index} = 'BP'; subjectNames{index} = 'SG'; expDates{index} = '080222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=354; deviceName{index} = 'BP'; subjectNames{index} = 'SG'; expDates{index} = '080222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=355; deviceName{index} = 'BP'; subjectNames{index} = 'SG'; expDates{index} = '080222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=356; deviceName{index} = 'BP'; subjectNames{index} = 'SG'; expDates{index} = '080222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=357; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '090222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=358; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '090222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=359; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '090222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=360; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '090222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=361; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '090222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=362; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '090222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=363; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '100222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=364; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '100222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=365; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '100222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=366; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '100222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=367; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '100222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 2; % Behavior + EEG

clear index; index=368; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=369; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=370; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=371; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=372; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120222'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=373; deviceName{index} = 'BP'; subjectNames{index} = 'NN'; expDates{index} = '120222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG

clear index; index=374; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '150222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=375; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '150222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=376; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '150222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=377; deviceName{index} = 'BP'; subjectNames{index} = 'NM'; expDates{index} = '150222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG

clear index; index=378; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '160222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=379; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '160222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=380; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '160222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=381; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '160222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=382; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '160222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=383; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '160222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=384; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=385; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=386; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=387; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=388; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=389; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=390; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'SRC_005'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=391; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'SRC_006'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=392; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '190222'; protocolNames{index} = 'SRC_007'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=393; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '210222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=394; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '210222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=395; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '210222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=396; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '210222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=397; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '210222'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=398; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '210222'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=399; deviceName{index} = 'BP'; subjectNames{index} = 'MK'; expDates{index} = '210222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG

clear index; index=400; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '220222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=401; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '220222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=402; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '220222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=403; deviceName{index} = 'BP'; subjectNames{index} = 'SP'; expDates{index} = '220222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=404; deviceName{index} = 'BP'; subjectNames{index} = 'RG'; expDates{index} = '240222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=405; deviceName{index} = 'BP'; subjectNames{index} = 'RG'; expDates{index} = '240222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=406; deviceName{index} = 'BP'; subjectNames{index} = 'RG'; expDates{index} = '240222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=407; deviceName{index} = 'BP'; subjectNames{index} = 'RG'; expDates{index} = '240222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=408; deviceName{index} = 'BP'; subjectNames{index} = 'RG'; expDates{index} = '240222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=409; deviceName{index} = 'BP'; subjectNames{index} = 'RG'; expDates{index} = '240222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=410; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=411; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=412; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=413; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=414; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=415; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=416; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=417; deviceName{index} = 'BP'; subjectNames{index} = 'IJ'; expDates{index} = '260222'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=418; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=419; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=420; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=421; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=422; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=423; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=424; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=425; deviceName{index} = 'BP'; subjectNames{index} = 'SM'; expDates{index} = '270222'; protocolNames{index} = 'SRC_004'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=426; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '010322'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=427; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '010322'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=428; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '010322'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=429; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '010322'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=430; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '010322'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=431; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '010322'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=432; deviceName{index} = 'BP'; subjectNames{index} = 'LK'; expDates{index} = '010322'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG

clear index; index=433; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '040322'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=434; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '040322'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=435; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '040322'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=436; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '040322'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=437; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '040322'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=438; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '040322'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=439; deviceName{index} = 'BP'; subjectNames{index} = 'SR'; expDates{index} = '040322'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG

clear index; index=440; deviceName{index} = 'BP'; subjectNames{index} = 'SJP'; expDates{index} = '050322'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=441; deviceName{index} = 'BP'; subjectNames{index} = 'SJP'; expDates{index} = '050322'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=442; deviceName{index} = 'BP'; subjectNames{index} = 'SJP'; expDates{index} = '050322'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=443; deviceName{index} = 'BP'; subjectNames{index} = 'SJP'; expDates{index} = '050322'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=444; deviceName{index} = 'BP'; subjectNames{index} = 'SJP'; expDates{index} = '050322'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=445; deviceName{index} = 'BP'; subjectNames{index} = 'SJP'; expDates{index} = '050322'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=446; deviceName{index} = 'BP'; subjectNames{index} = 'SJP'; expDates{index} = '050322'; protocolNames{index} = 'SRC_003'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

clear index; index=447; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=448; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=449; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=450; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=451; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'GRF_005'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=452; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'GRF_006'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=453; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=454; deviceName{index} = 'BP'; subjectNames{index} = 'AM'; expDates{index} = '060322'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG

clear index; index=455; deviceName{index} = 'BP'; subjectNames{index} = 'SC'; expDates{index} = '080322'; protocolNames{index} = 'GRF_001'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=456; deviceName{index} = 'BP'; subjectNames{index} = 'SC'; expDates{index} = '080322'; protocolNames{index} = 'GRF_002'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=457; deviceName{index} = 'BP'; subjectNames{index} = 'SC'; expDates{index} = '080322'; protocolNames{index} = 'GRF_003'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=458; deviceName{index} = 'BP'; subjectNames{index} = 'SC'; expDates{index} = '080322'; protocolNames{index} = 'GRF_004'; capLayout{index} = 'actiCap64'; stimType{index} = 4; % Behavior + EEG
clear index; index=459; deviceName{index} = 'BP'; subjectNames{index} = 'SC'; expDates{index} = '080322'; protocolNames{index} = 'SRC_001'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG
clear index; index=460; deviceName{index} = 'BP'; subjectNames{index} = 'SC'; expDates{index} = '080322'; protocolNames{index} = 'SRC_002'; capLayout{index} = 'actiCap64'; stimType{index} = 7; % Behavior + EEG

end