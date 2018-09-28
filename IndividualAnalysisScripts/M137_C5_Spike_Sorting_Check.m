%% CFA\M137_C5===============================================================================================
%define the directories
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\CFA\M137_C5\'
% get the data for M137_C5
fname='M137_C5_Mech_L6 05mW'
load([homedir fname 'analysis'], 'Triggers', 'Conditions')
load([homedir fname '_all_channels'])
load([homedir 'CrossCoeffData'])
LabelStr='UMS2k'
datastr=fname;
units=sortedData(:,1);
indices=sort([goods 6])
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
sp=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms

%issues: 1. unit 6 is laser contaminated. 2. Units 1 and 2 should be merged
%due to shifted xcorrelation
% GOOD UNITS=?

SpikeSortingChecks

%% CFA\M16_C2===============================================================================================
%define the directories
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\CFA\M16_C2\'

% get the data for M16_C2
fname='M16_C2_Mech+L6_4mW'
load([homedir fname 'analysis'], 'Triggers', 'Conditions')
load([homedir fname '_all_channels'])
load([homedir 'CrossCoeffData'])
units=sortedData(:,1);
indices=sort([goods bads])% ALL UNITS
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
LabelStr='UMS2k'
datastr=fname;
sp=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms

%ISSUES/NOTES
% 'chPack3_cl_1 chPack4_cl_34' are double counted!   possible bursting
% 'chPack4_cl_40 chPack5_cl_26' are double counted!   possible bursting
% no laser artifacts
SpikeSortingChecks

%% CFA\M137_C4===============================================================================================
%define the directories
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\CFA\M137_C4\'

% get the data for M137_C4
fname='M137_C4_Mech+L6 05mW'
load([homedir fname 'analysis'], 'Triggers', 'Conditions')
load([homedir fname '_all_channels'])
load([homedir 'CrossCoeffData'])
units=sortedData(:,1);
indices=sort([goods bads])
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
LabelStr='UMS2k'
datastr=fname;
sp=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms

%ISSUES/NOTES
%this recording has tons of artifacts and cross pack double detection
%but not clear why some recordings are labeled bad....
% bads=
%
%     'chPack1_cl_18' ??????
%     'chPack1_cl_32' light
%     'chPack1_cl_33' ??????
%     'chPack2_cl_13' ??????
%     'chPack2_cl_18' light
%     'chPack3_cl_3' light
%     'chPack4_cl_5' light
%     'chPack5_cl_2' light

%MERGING of real neural signals
% 'chPack1_cl_36 chPack2_cl_10'
% 'chPack1_cl_40 chPack2_cl_56'
names(bads)
%
SpikeSortingChecks

%%

%% Saline\M7_C3===============================================================================================
%define the directories
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\Saline\M7_C3\'

% get the data
fname='M7_C3_Mech_05mW'
load([homedir fname 'analysis'], 'Triggers', 'Conditions')
load([homedir fname '_all_channels'])
load([homedir 'CrossCoeffData'])
units=sortedData(:,1);
indices=sort([goods bads])
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
LabelStr='UMS2k'
datastr=fname;
sp=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms

%ISSUES/NOTES
%  'chPack4_cl_1' why bad?  weird small intervals


%   'chPack1_cl_9 chPack2_cl_122'   %merging issue/bursting
%   'chPack1_cl_112 chPack3_cl_5'     %not clear but needs examination
%   'chPack2_cl_108 chPack3_cl_29' %merging issue/bursting
%   'chPack4_cl_25 chPack5_cl_2' %merging issue/bursting


%possible issues
% 'chPack3_cl_5chPack3_cl_6'
%   'chPack1_cl_9 chPack1_cl_13'
%    'chPack1_cl_9chPack2_cl_109'
%    'chPack1_cl_13chPack2_cl_109'

names(sort(bads))

%%
SpikeSortingChecks

%% CFA\M27_C2_L6===============================================================================================
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\CFA\HL_L6alone\M27_C2_L6_recordings\'

% get the data
fname='M27_C2_Mechmaybe'
load([homedir fname 'analysis'], 'Triggers', 'Conditions')
load([homedir fname '_all_channels'])
load([homedir 'LaserResponse'])
units=sortedData(:,1);

indices=sort([goods bads])
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
LabelStr='UMS2k'
datastr=fname;
sp=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms

% ISSUES/NOTES

%very suspicious
'chPack1_cl_5chPack2_cl_1'
'chPack1_cl_7chPack2_cl_1'
'chPack1_cl_5chPack2_cl_18'

%%

'chPack1_cl_22 chPack2_cl_1'
'chPack1_cl_1 chPack1_cl_5'
'chPack1_cl_1 chPack1_cl_7'
'chPack1_cl_7 chPack2_cl_18'
'chPack2_cl_1chPack2_cl_18'
'chPack2_cl_1chPack2_cl_46'
'chPack1_cl_7chPack2_cl_48'


names(sort(bads))

%%
SpikeSortingChecks

%% CFA\M27_C1_L6===============================================================================================
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\CFA\HL_L6alone\M27_C1_L6_recordings\'

% get the data
fname='M27_C1_L6Mechmaybe'
load([homedir fname 'analysis'], 'Triggers', 'Conditions')
load([homedir fname '_all_channels'])
%load([homedir 'LaserResponse']) we don't have this information yet
units=sortedData(:,1);
bads=[];goods=1:numel(units);
indices=sort([goods bads])
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
LabelStr='UMS2k'
datastr=fname;
sp=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms

% SpikeSortingChecks
