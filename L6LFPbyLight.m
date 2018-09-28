

%% for recording 'M106_C2_Mech+5mW_new' check to see if L6 stimulation paradigm is working
% load data

close all
clear all

cd('F:\Experiments_2018\Sailaja Juxtasomal Data_RM\CFA injected\3_1_2018\M133_C2')
filename='M113_C2_Mech+L6alone edited'
load([filename '.mat'],'chan22')
load([filename 'analysis.mat'],'Triggers','spikeFindingData', 'EEG')


%%  processing LFP, getting light triggers, overplotting the two
LFP=zscore(EEG.data);   %is this always the correct channel?
light=Triggers.light;light=light/max(light);
ppms=spikeFindingData.ppms;

figure
plot(LFP,'k')
hold on
plot(light,'c')
%get light onsets
tL=getTriggersNew(light,100,.9);%this is stupid, this should be saved already somewhere else
plot(tL,light(tL),'o')

%% get triggered signals, overplot and produce figure
figure
timeBefore=-50*ppms;timeAfter=150*ppms;
LFP_l = TriggeredSegments(LFP, tL, timeBefore, timeAfter);
light_l = TriggeredSegments(light, tL, timeBefore, timeAfter);
t=[1:size(LFP_l,1)]/ppms+timeBefore/ppms;
m=mean(LFP_l');
plot(t,LFP_l,'color',[.5 .5 .5])
hold on
plot(t,m,'k','linewidth',2)
hold on
plot(t,mean(light_l'),'c');
xlabel ms
title(filename,'interpreter','none')

%% for recording 'M7_C3_Mech_05mW' check to see if L6 stimulation paradigm is working
% load data

close all
clear all

cd('F:\Experiments_2018\24_8_2018\M16_C2')
filename='M16_C2_Mech+L6_4mW'
load([filename '.mat'],'chan22')
load([filename 'analysis.mat'],'Triggers','spikeFindingData', 'EEG')


%%  processing LFP, getting light triggers, overplotting the two
LFP=zscore(chan22);   %is this always the correct channel?
light=Triggers.light;light=light/max(light);
ppms=spikeFindingData.ppms;

figure
plot(LFP,'k')
hold on
plot(light,'c')
%get light onsets
tL=getTriggersNew(light,100,.9);%this is stupid, this should be saved already somewhere else
plot(tL,light(tL),'o')

%% get triggered signals, overplot and produce figure
figure
timeBefore=-50*ppms;timeAfter=150*ppms;
LFP_l = TriggeredSegments(LFP, tL, timeBefore, timeAfter);
light_l = TriggeredSegments(light, tL, timeBefore, timeAfter);
t=[1:size(LFP_l,1)]/ppms+timeBefore/ppms;
m=mean(LFP_l');
plot(t,LFP_l,'color',[.5 .5 .5])
hold on
plot(t,m,'k','linewidth',2)
hold on
plot(t,mean(light_l'),'c');
xlabel ms
title(filename,'interpreter','none')




