%% let's compare UMS and kilosort...
clear all
close all
cd D:\GD_Rebecca\Work\1_Active_Projects\PAIN\MUA_analysis

homedir='PracticeCheck'
cd(homedir)
load M16_C2_Mech+L6_4mWanalysis.mat Triggers Fs Conditions
stimulus=lazynorm(Triggers.whisker);
ppms=Fs/1000
T=numel(Triggers.whisker)/ppms/1000; %duration of experiment in seconds!

%% find every single trigger for step stimuli
figure
tmech=diff(lazynorm(Triggers.whisker));tmech=find(tmech>.5);
tlight=diff(lazynorm(Triggers.light));tlight=find(tlight>.5);
t=1:numel(stimulus);t=t/ppms;

plot(t,lazynorm(Triggers.whisker),'r',t,lazynorm(Triggers.light),'b')
hold on
plot(tmech/ppms,ones(size(tmech)),'ro',tlight/ppms,ones(size(tlight)),'ob')

%% kilosort
LabelStr='kilosort'
datastr='M16_C2_Mech+L6_4mW'
load M16_C2_Mech+L6_4mW_KS.mat 
units=Spikes(:,1);
sp=cellfun(@double,Spikes(:,2),'UniformOutput', false);
sp=cellfun(@(x) x/ppms,sp,'UniformOutput', false); %ms
SpikeSortingChecks

%% UMS
LabelStr='UMS2k'
load M16_C2_Mech+L6_4mW_all_channels.mat  
datastr='M16_C2_Mech+L6_4mW'
units=sortedData(:,1);
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);
sp=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms
SpikeSortingChecks
