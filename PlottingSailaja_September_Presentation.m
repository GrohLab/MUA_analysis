% plotting for Sailaja talk
% cortex mechanical stimulation
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\CFA\HL_L6alone\M27_C1_L6_recordings\'

% get the data
fname='M27_C1_L6Mechmaybe'
load([homedir fname 'analysis'], 'Triggers', 'Conditions','EEG')
load([homedir 'MechPressure_M27_C1_L6Mechmaybe'])
load([homedir fname '_all_channels'])
%load([homedir 'LaserResponse'])
units=sortedData(:,1);
goods=1:size(units);
indices=sort([goods])
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
LabelStr='UMS2k'
datastr=fname;
Spikes=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms
ppms=fs/1000

%% mechanical triggers
triggers=Conditions{1}.Triggers(2:end)/ppms; %in ms
AllSpikes={cell2mat(Spikes)} %pool all spikes due to shitty spike sorting
N_Sp=numel(AllSpikes);n_trig=numel(triggers); %number of units, number of triggers

label='mech'
timeBefore=1000;timeAfter=7000; %in ms
binsize=100; %ms
mech=lazynorm(Triggers.whisker); %continuous signals of interest
press=lazynorm(double(chan21));
[TriggeredSpikeTimes Xvals Yvals]=TriggeredSpikes(AllSpikes,triggers,timeBefore, timeAfter);  %get all spike times for a single trigger

Xs=ManyTriggeredSegments({mech press},triggers*ppms,timeBefore*ppms, timeAfter*ppms); %This is in sample for continuous
t=[1:size(Xs{1},1)]/ppms-timeBefore;
%f=plotMUA(TriggeredSpikeTimes,timeBefore,timeAfter, Xs,na,label, Names)

close all
figure
plot(mech)
hold on
plot(triggers*ppms,ones(size(triggers)),'ro')
figure
plot(t,Xs{1});ylim([-.1 1.1])

%%
%make histograms
H=[];
bins=[-timeBefore:binsize:timeAfter];
for I=1:N_Sp
    H(I,:)=hist(Xvals{I},bins)/numel(triggers);
end


%%
%means of continuous variables
Xmeans = cellfun(@transpose, Xs,'UniformOutput', false);
Xmeans = cellfun(@mean, Xmeans,'UniformOutput', false);

figure
bar(bins,H/numel(triggers),1)
title(fname,'interpreter','none')
axis tight
ylabel 'mean spike/bin, 100 ms'
xlabel ms
ylim([0 1.2])
hold on
plot(t,Xmeans{1},'k')   
plot(t,Xmeans{2},'r')





%%

gray=[.7 .7 .7]
% get color for each neuron, just for plotting
colors=brighten(distinguishable_colors(N_Sp),-.25);%color for each unit

% plot it all
figure
subplot(6,4, 1:18)
for I=1:N_Sp
    plot(Xvals{I},Yvals{I},'.','color',colors(I,:),'markersize',5)
    hold on
    
    %label for each neuron
    thisname=Names{I}{1}
    thisname_x=timeAfter
    thisname_y=I*n_trig-n_trig/2;
    text(thisname_x,thisname_y,thisname,'interpreter','none')
    
end

title(label,'interpreter','none')
% ylabel 'trials/neuron'
box off
xlim([min(bins) max(bins)])

%% Triggered conditions analysis for the very good recording==========================
%delay===
clear all
close all
homedir='D:\GD_Rebecca\Work\1_Active_Projects\PAIN\Files for RM_22_9_2018\CFA\M16_C2\'
% get the data for M16_C2
fname='M16_C2_Mech+L6_4mW'
load([homedir fname 'analysis'], 'Triggers', 'Conditions','EEG')
load([homedir fname '_all_channels'])
load([homedir 'CrossCoeffData'])
units=sortedData(:,1);
indices=sort([goods])
sp=cellfun(@double,sortedData(:,2),'UniformOutput', false);sp=sp(indices);
names=units(indices);
LabelStr='UMS2k'
datastr=fname;
Spikes=cellfun(@(x) x'*1000,sp,'UniformOutput', false); %ms
ppms=fs/1000


%%  mechanical triggers  %do the stats
DELAY=200;%in ms
%for II=1:numel(Conditions);
triggers=Conditions{1}.Triggers(1:end)/ppms+DELAY; %in ms
AllSpikes=Spikes; %pool all spikes due to shitty spike sorting
N_Sp=numel(AllSpikes);n_trig=numel(triggers); %number of units, number of triggers
%290-300
label='mech'
timeBefore=1000;timeAfter=2000; %in ms
binsize=100; %ms
mech=lazynorm(Triggers.whisker); %continuous signals of interest
[TriggeredSpikeTimes Xvals Yvals]=TriggeredSpikes(AllSpikes,triggers,timeBefore, timeAfter);  %get all spike times for a single trigger

Xs=ManyTriggeredSegments({mech},triggers*ppms,timeBefore*ppms, timeAfter*ppms); %This is in sample for continuous
t=[1:size(Xs{1},1)]/ppms-timeBefore;
%f=plotMUA(TriggeredSpikeTimes,timeBefore,timeAfter, Xs,na,label, Names)

close all
figure
plot(mech)
hold on
plot(triggers*ppms,ones(size(triggers)),'ro')
figure
plot(t,Xs{1});ylim([-.1 1.1])
%
%make histograms
H=[];
bins=[-timeBefore:binsize:timeAfter];
for I=1:N_Sp
    H(I,:)=hist(Xvals{I},bins)/numel(triggers);
end

%%
figure
bar(bins,sum(H),1);hold on,plot(bins,mean(H),'r')
hold on
plot(t,mean(Xs{1}'-2))



%%  mechanical triggers  %do the stats
DELAY=200;%in ms
window=1000
TriggeredSpikeTimes={};
for II=1:numel(Conditions);
    if II==4
        triggers=Conditions{II}.Triggers(1:end)/ppms; %in ms
    else
        triggers=Conditions{II}.Triggers(1:end)/ppms+DELAY; %in s
    end
    AllSpikes=Spikes;
    N_Sp=numel(AllSpikes);n_trig=numel(triggers); %number of units, number of triggers
    %290-300
    label='mech'
    timeBefore=500;timeAfter=window; %in ms
    binsize=100; %ms
    mech=lazynorm(Triggers.whisker); %continuous signals of interest
    laser=lazynorm(Triggers.light); %continuous signals of interest
    eeg=lazynorm(EEG.data); %continuous signals of interest
    
    [TriggeredSpikeTimes{II} Xvals Yvals]=TriggeredSpikes(AllSpikes,triggers,timeBefore, timeAfter);  %get all spike times for a single trigger
    Xs=ManyTriggeredSegments({mech, laser,eeg},triggers*ppms,timeBefore*ppms, timeAfter*ppms); %This is in sample for continuous
    t=[1:size(Xs{1},1)]/ppms-timeBefore;
    %f=plotMUA(TriggeredSpikeTimes,timeBefore,timeAfter, Xs,na,label, Names)
    
    figure
    plot(t,mean(Xs{1}'),'r',t,mean(Xs{2}'),'c',t,mean(Xs{3}')*20-.5,'k')
    set(gcf,'paperpositionmode','auto')
end

%control for mech! 1000 before mech triggers!
triggers=Conditions{1}.Triggers(1:end)/ppms+DELAY; %in ms
[TriggeredSpikeTimes{II+1} Xvals Yvals]=TriggeredSpikes(AllSpikes,triggers,window, 0);  %get all spike times for a single trigger

%control for light! 1000 ms before light riggers
triggers=Conditions{4}.Triggers(1:end)/ppms+DELAY; %in ms
[TriggeredSpikeTimes{II+2} Xvals Yvals]=TriggeredSpikes(AllSpikes,triggers,window, 0);  %get all spike times for a single trigger

%%Condition>Neurons>Trials

for II=1:numel(TriggeredSpikeTimes)
    for JJ=1:numel(TriggeredSpikeTimes{II}) %number of neurons
        Counts{II,JJ}=cellfun(@numel, TriggeredSpikeTimes{II}{JJ})
    end
end

%% Is mechanical statistically sig?
% Is does light stim significantly decrease or increase?
Control={};Exp={};LM1={};LM10={};Light={};ControlL={};
for i=1:N_Sp
    Exp{i}=cell2mat(cellfun(@numel, TriggeredSpikeTimes{1}{i},'UniformOutput', false)) %mech
    LM1{i}=cell2mat(cellfun(@numel, TriggeredSpikeTimes{2}{i},'UniformOutput', false)) %1 Hz
    LM10{i}=cell2mat(cellfun(@numel, TriggeredSpikeTimes{3}{i},'UniformOutput', false)) %10 Hz
    Light{i}=cell2mat(cellfun(@numel, TriggeredSpikeTimes{4}{i},'UniformOutput', false)) %light control
    Control{i}=cell2mat(cellfun(@numel, TriggeredSpikeTimes{5}{i},'UniformOutput', false)) %spontaneous before mech
    ControlL{i}=cell2mat(cellfun(@numel, TriggeredSpikeTimes{6}{i},'UniformOutput', false)) %spontaneous before light
end

%median or mean by cell
M_Control=cell2mat(cellfun(@mean, Control,'UniformOutput', false))
M_Exp=cell2mat(cellfun(@mean,Exp,'UniformOutput', false))
M_LM1=cell2mat(cellfun(@mean,LM1,'UniformOutput', false))
M_LM10=cell2mat(cellfun(@mean,LM10,'UniformOutput', false))
M_Light=cell2mat(cellfun(@mean,Light,'UniformOutput', false))

%does mechanical stim drive response  YES, in 6 out of 12 neurons
P_mech=cell2mat(cellfun(@signrank,Control,Exp,'UniformOutput', false))
%does light drive response alone   NO, in 0 out of 12 neurons
P_light=cell2mat(cellfun(@signrank,ControlL,Light,'UniformOutput', false))

% does light significantly change mech response for 1 hz stim?  NO 0/12
P_LM1=cell2mat(cellfun(@signrank,LM1,Exp,'UniformOutput', false))

% does light significantly change mech response for 1 hz stim?  NO 0/12
P_LM10=cell2mat(cellfun(@signrank,LM10,Exp,'UniformOutput', false))

%does light stim drive response

plot(M_Control,M_Exp,'o')
mechIndices=find(P_mech<.05);
lightIndices=find(P_LM10<.05)
ModSubset=intersect(mechIndices,lightIndices)

figure
subplot(1,2,1)
plot(M_Control,M_Exp,'.k','linewidth',2)
hold on
plot(M_Control(mechIndices),M_Exp(mechIndices),'or','linewidth',2)
plot([0 6],[0 6],':k')
axis square
xlabel 'spontaneous [Hz]'
ylabel 'mech evoked [Hz]'
legend('sorted unit', 'p<0.05')


figure
subplot(1,2,1)
plot(M_Control,M_Exp,'.k','linewidth',2)
hold on
plot(M_Control(mechIndices),M_Exp(mechIndices),'or','linewidth',2)
plot([0 7],[0 7],':k')
axis square
xlabel 'spontaneous [Hz]'
ylabel 'mech evoked [Hz]'
legend('sorted unit', 'p<0.05')

subplot(1,2,2)
plot(M_Exp,M_LM10,'.k','linewidth',2)
hold on
plot(M_Exp(ModSubset),M_LM10(ModSubset),'oc','linewidth',2)
plot([0 7],[0 7],':k')
axis square
xlabel 'mech evoked [Hz]'
ylabel 'mech evoked +L10Hz [Hz]'
legend('sorted unit', 'p<0.05')

%% triggered LFP





%%