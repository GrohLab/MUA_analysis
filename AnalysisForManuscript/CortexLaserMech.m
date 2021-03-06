clear all

file='Z:\Ross\Experiments\3.5mW_Saline_VPL\M160_L6Stim+Mech_Saline'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME 'analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])


%% light driven units in cortex  CFA??? CORTEX
clear all
%at the moment, this has light artifacts Feb
file='Z:\Ross\Experiments\10mW_CFA_Cortex\M167_HL_Cortex'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
load('re-curated_all_channels.mat');
load('re-curatedanalysis.mat','Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])

stimCond=Conditions([1:4]);
ppms=fs/1000;
%% DREADDS ANALYSIS inhibitory  CORTEX CONTROL

clear all
file='Z:\Ross\Experiments\DREADDs\Gi\M34_GiConc'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME '_analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])

%% excitatory   CORTEX CONTROL
clear all
file='Z:\Ross\Experiments\DREADDs\Gq\M30_GqConc'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME '_analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])

%%  CFA??
clear all
file='Z:\Ross\Experiments\10mW_CFA_VPL\M167_L6-Stim10mW_Mech_CFA_ForePaw_Conc'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME 'analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])


%%
clear all
file='Z:\Ross\Experiments\10mW_Saline_VPL\M168_10mW_MechStim_Saline_VPL'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME 'analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])

%% rerun to include all triggers

clear all
file='Z:\Ross\Experiments\2.5mW_CFA_VPL\CFA_VPL_2p5mw'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME 'analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])



%%  not much activity
clear all
file='Z:\Sailaja\Manuscript\Data\16 channel\Poly2A probe\Saline\M7_C3_Kilosorted\Kilosort\Th6 4\M7_C3_Mech_05mW'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME 'analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])

%% looks pretty good!



clear all
file='Z:\Sailaja\Manuscript\Data\16 channel\standard 16 channel probe\M137_C5\Kilosort\Th6 2\M137_C5_Mech_L6 05mW'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
%spikefile=[NAME '_all_channels.mat'];
load('M137_all_channels.mat');
load([NAME 'analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])
%%
clear stimCond
%necessary for Sailaja's data in old format
for i=1:numel(Conditions)
stimCond(i)=Conditions{i};
if size(stimCond(i).Triggers,2)>2,stimCond(i).Triggers=stimCond(i).Triggers';end
end
neworder=[2 3 4 1] %to match ross's order!
stimCond=stimCond(neworder)
%% organized spike trains
% well isolated clusters, only isolated


ppms=fs/1000;
stimCond=Conditions(1:2) % CHANGE THIS THIS!! FOR ROSS DATA!!!
goodsInd = cellfun(@(x) x==1,sortedData(:,3));multiInd=[];
multiInd = cellfun(@(x) x==2,sortedData(:,3));
goods = [find(goodsInd); find(multiInd)];
%Ns = min(structfun(@numel,Triggers));
clusters=sortedData(goods,2);ids=sortedData(goods,1);
Ns = max(max(cell2mat(clusters)))*ppms*1000; %this is a fix for no mech stim, DREADDS exp

% Total duration of the recording
Nt = Ns/fs;  %seconds

if size(clusters{1},2)~=1,clusters=cellfun(@transpose,clusters,'UniformOutput',0);end
maxNonZero = sum(cellfun(@numel, clusters));
biSpikes=spalloc(floor(Nt*1000*ppms),numel(goods),maxNonZero);% preallocate sparse matrix to save space

% ms x num clusters
T=size(biSpikes,1); %recording length in ms
biISIs=biSpikes;

ISIs=cell(size(clusters));
for i=1:numel(clusters)
    sp=[clusters{i};Nt];
    rm=returnmap(sp,0);
    ISIs{i}=rm.d1;
end


%make big matrix of all spike times
for j=1:numel(goods)
    tt=round(clusters{j}*1000*ppms); okays=(tt>0 & tt<=T);%spike times in ms*ppms (samples)
    S=zeros(T,1);S(tt(okays))=1;
    biSpikes(:,j)=sparse(S);
    Sisi=zeros(T,1);Sisi(tt(okays))=ISIs{j}(okays);
    biISIs(:,j)=sparse(Sisi);
    j/numel(goods)
end



%% for all conditions
%maxTrials=max(cellfun(@(x) size(x,1),{stimCond.Triggers}));
%in ms because of division before, though triggered seg takes samples
timeBefore=(-3000*ppms);
timeAfter=round(7000*ppms);
trialLength=-timeBefore+timeAfter+1;

%preallocate omg could be huge!

%firingRates: N x S x T x maxTrialNum
%firingRates=zeros(numel(clusters),numel(stimCond),round(trialLength),maxTrials);
relativeSpikeTimes=struct([]);relativeISIs=struct([]);
catSp=cell(numel(clusters),numel(stimCond));catISIs=catSp;
for n =1:numel(clusters)
    for s=1:numel(stimCond)
        triggers=round(stimCond(s).Triggers(:,1)/fs*1000*ppms); %triggers in ms*ppms
        %triggers=triggers(1:end-1);
        X=TriggeredSegments(biSpikes(:,n),triggers,timeBefore,timeAfter);
        %firingRates(n,s,1:size(X,1),1:size(X,2))=(shiftdim(X,-2));  %some fuckery to deal with 4d array
        sp={};tnum={};
        for ii=1:size(X,2),sp{ii}=find(X(:,ii))+timeBefore; tnum{ii}=ones(size(sp{ii}))*ii;end
        relativeSpikeTimes(n,s).sp=sp;
        relativeSpikeTimes(n,s).tnum=tnum;
        relativeSpikeTimes(n,s).clustername=ids(n);
        relativeSpikeTimes(n,s).condition=stimCond(s).name;
        catSp{n,s}=cell2mat(sp');
        isis={};
        X=TriggeredSegments(biISIs(:,n),triggers,timeBefore,timeAfter);
        for ii=1:size(X,2),isis{ii}=X(sp{ii}-timeBefore,ii);end
        catISIs{n,s}=cell2mat(isis');
        relativeSpikeTimes(n,s).isis=isis;
    end
    n/numel(clusters)
end
%%   histograms

close all
binsize=100*ppms;
bins=[timeBefore:binsize:timeAfter];
H=cellfun(@(c) hist(c/ppms,bins/ppms),catSp,'UniformOutput',0);

f1=figure
caxes=[];
for i=1:numel(stimCond)
    Hcond=cell2mat(H(:,i));
    
    %Hcond(Hcond==0)=nan;
    subplot(2,2,i)
    %plot histogram normalized to counts
    image(bins/ppms,1:numel(clusters),Hcond/numel(relativeSpikeTimes(1,i).isis)/(binsize/fs),'CDataMapping','scaled')
    
    shading flat
    colormap fire
    
    % figure
    % plot(bins/ppms,Hcond)
    %
    % xlabel ms
    ylabel(['counts/bin binsize=' num2str(binsize/ppms) ' ms'])
    title(stimCond(i).name)
    caxes=[caxes;caxis];
    colorbar
end
for i=1:numel(stimCond)
    subplot(2,2,i)
    caxis([0 25])
end
savefig('Histograms')

save TriggeredAnalysis bins ppms H catSp stimCond ids relativeSpikeTimes fs binsize catISIs


%%
close all
binsize=5*ppms;
bins=[timeBefore:binsize:timeAfter];
H=cellfun(@(c) hist(c/ppms,bins/ppms),catSp,'UniformOutput',0);

i=3
   Hcond=cell2mat(H(:,i));
    
   
   
    
figure
plot(bins/ppms,Hcond)
xlim([-200 2000])
xlabel ms
    ylabel(['counts/bin binsize=' num2str(binsize/ppms) ' ms'])
    title(stimCond(i).name)
    caxes=[caxes;caxis];
    colorbar

savefig('Zoom Laser Control')


window=find(bins>(-1*ppms) & bins<(5*ppms))
maxs=[];indices=[];
for i=1:size(Hcond,1)
    m=Hcond(i,window);maxs(i)=max(m);
    indices(i)=find(m==maxs(i),1);
end


figure
plot(bins(indices)/ppms,maxs,'o')
text(bins(indices)/ppms,maxs,ids)
hold on;
%plot(tx/ppms,-X*max(max(Hcond)),'linewidth',3)
axis tight
xlabel ms
ylabel(['counts/bin binsize=' num2str(binsize/ppms) ' ms'])
ti=title(['Possible light artifacts ' spikefile],'Interpreter','none')
savefig('LightArtifacts')





%%
% rasters
figure
lastlevel=0;
for n=1:numel(clusters)
    if mod(n,2)==0, col='b';else col='k';end
    [R lastlevel]=manyRasters(relativeSpikeTimes(n,3).sp,col,.9,ppms,0+lastlevel);
end





%% for laser trials only  NEEDS TO BE FIXED
triggers={lt laser_control}
maxTrials=max(cellfun(@(x) size(x,1),triggers));
%in ms because of division before, though triggered seg takes samples
timeBefore=(-5*ppms);
timeAfter=round(35*ppms);
trialLength=-timeBefore+timeAfter+1;

%preallocate omg could be huge!
% firingRates: N x S x T x maxTrialNum
%firingRates=zeros(numel(clusters),numel(triggers),round(trialLength),maxTrials);
relativeSpikeTimes=struct([]);
catSp=cell(numel(clusters),numel(triggers));
%need to parallelize
for n =1:numel(clusters)
    for s=1:numel(triggers)
        trigs=round(triggers{s}/fs*1000*ppms); %triggers in ms*ppms
        X=TriggeredSegments(biSpikes(:,n),trigs,timeBefore,timeAfter);
        %firingRates(n,s,1:size(X,1),1:size(X,2))=(shiftdim(X,-2));  %some fuckery to deal with 4d array
        sp={};tnum={};
        for ii=1:size(X,2),sp{ii}=find(X(:,ii))+timeBefore; tnum{ii}=ones(size(sp{ii}))*ii;end
        relativeSpikeTimes(n,s).sp=sp;
        relativeSpikeTimes(n,s).tnum=tnum;
        relativeSpikeTimes(n,s).clustername=ids(n);
        relativeSpikeTimes(n,s).condition='laser';
        catSp{n,s}=cell2mat(sp');
    end
end
 histograms
thisCond=1;

close all
binsize=.5*ppms;
bins=[timeBefore:binsize:timeAfter];
H=cellfun(@(c) hist(c/ppms,bins/ppms),catSp,'UniformOutput',0);
Hcond=cell2mat(H(:,thisCond));
figure
%Hcond(Hcond==0)=nan;
subplot(3,1,1:2)
pcolor(bins/ppms,1:numel(clusters),Hcond)
shading flat
colormap fire

X=mean(TriggeredSegments(laser,triggers{thisCond},timeBefore,timeAfter)');
tx=[timeBefore:timeAfter]'
subplot(3,1,3)
plot(tx/ppms,X)
axis tight

laserMed=cell2mat(cellfun(@median,catSp(:,thisCond),'UniformOutput',0))/ppms;
figure
histogram(laserMed,1000)
plot(bins/ppms,Hcond)

xlabel ms
ylabel(['counts/bin binsize=' num2str(binsize/ppms) ' ms'])
ti=title(['Possible light artifacts ' spikefile],'Interpreter','none')
savefig('LightArtifactsHist')


maxs=[];indices=[];
for i=1:size(Hcond,1)
    m=Hcond(i,:);maxs(i)=max(m);
    indices(i)=find(m==maxs(i),1);
end


figure
plot(bins(indices)/ppms,maxs,'o')
text(bins(indices)/ppms,maxs,ids)
hold on;
plot(tx/ppms,-X*max(max(Hcond)),'linewidth',3)
axis tight
xlabel ms
ylabel(['counts/bin binsize=' num2str(binsize/ppms) ' ms'])
ti=title(['Possible light artifacts ' spikefile],'Interpreter','none')
savefig('LightArtifacts')




%stimCond=Conditions([1:4]);

% get single laser pulses

thislaser=Triggers.laser;
thislaser=thislaser-min(thislaser);thislaser=thislaser/max(thislaser);
thislaser(thislaser>.5)=1;
thislaser(thislaser~=1)=0;
laser=thislaser; %every single laser pulse;
lt=find(diff(thislaser)==1); %triggers for all lasers

whisker=Triggers.whisker;
whisker=whisker-min(whisker); whisker=whisker/max(whisker);
whisker(whisker>.9)=1;


for i=1:2
    trig=stimCond(i).Triggers;
    pad=round(ppms*1000);
    for j=1:size(trig,1)
        indices=(trig(j,1)-pad):(trig(j,2)+pad);
        indices=indices(indices>=1 & indices<=(numel(whisker)));
        whisker(indices)=ones(size(indices))';
    end
end

whisker(whisker~=1)=0;
thislaser=thislaser-whisker;
thislaser(thislaser~=1)=0;
laser_control=find(diff(thislaser)==1); %triggers for all lasers

figure
plot(lt/10,laser(lt),'o')
hold on
plot(laser(1:10:end))

save laserinfo lt laser_control
