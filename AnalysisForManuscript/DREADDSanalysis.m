%DREADDS ANALYSIS inhibitory

clear all
file='Z:\Ross\Experiments\DREADDs\Gi\M34_GiConc'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME '_analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])

%% excitatory
clear all
file='Z:\Ross\Experiments\DREADDs\Gq\M30_GqConc'
[FILEPATH,NAME,EXT] = fileparts(file);
cd(FILEPATH)
spikefile=[NAME '_all_channels.mat'];
load(spikefile);
load([NAME '_analysis.mat'],'Conditions','Triggers')
load([NAME '_sampling_frequency.mat'])

%%
%spontaneous responses (before mech trigger)

%mechanical responses (after mech trigger)
stimCond=Conditions([1:2]);
ppms=fs/1000;



%% organized spike trains
% well isolated clusters, only isolated
goodsInd = cellfun(@(x) x==1,sortedData(:,3));multiInd=[];
multiInd = cellfun(@(x) x==2,sortedData(:,3));
goods = [find(goodsInd); find(multiInd)];
clusters=sortedData(goods,2);ids=sortedData(goods,1);


Ns = max(max(cell2mat(clusters)))*ppms*1000; %this is a fix for no mech stim
% Total duration of the recording
Nt = Ns/fs;  %seconds
biSpikes=spalloc(floor(Nt*1000*ppms),numel(goods),numel(cell2mat(clusters)));% preallocate sparse matrix to save space
biISIs=biSpikes;
% ms x num clusters
T=size(biSpikes,1); %recording length in ms


ISIs=cell(size(clusters));
for i=1:numel(clusters)
    sp=[clusters{i};Nt];
    rm=returnmap(sp,0);
    ISIs{i}=rm.d1;
end


%stimCond(end-1).Triggers=stimCond(end-1).Triggers(1:5:end-1,:)
%make big matrix of all spike times
for j=1:numel(goods)
    tt=round(clusters{j}*1000*ppms);tt=tt(tt>0 & tt<=T);%spike times in ms*ppms (samples)
    S=zeros(T,1);S(tt)=1;
    biSpikes(:,j)=sparse(S);
    Sisi=zeros(T,1);Sisi(tt)=ISIs{j};
    biISIs(:,j)=sparse(Sisi);
end


figure
sp=clusters{1};tr=Conditions(1).Triggers(:,1);;tr2=Conditions(2).Triggers(:,1);
plot(sp*1000*ppms,ones(size(sp)),'o')
hold on
plot(tr,ones(size(tr))+.2,'o')
plot(tr2,ones(size(tr2))+.2,'o')



%% for all conditions
maxTrials=max(cellfun(@(x) size(x,1),{stimCond.Triggers}));
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
    subplot(1,2,i)
    image(bins/ppms,1:numel(clusters),Hcond/numel(stimCond(i).Triggers(:,1))/(binsize/fs),'CDataMapping','scaled')
    
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
    subplot(1,2,i)
    caxis([0 50])
end
savefig('BeforeAndAfter_Histogram')


%
% maxs=[];indices=[];
% for i=1:size(Hcond,1)
%     m=Hcond(i,:);maxs(i)=max(m);
%     indices(i)=find(m==maxs(i),1);
% end

save TriggeredAnalysis bins ppms H catSp stimCond ids relativeSpikeTimes fs binsize
%%
% rasters
figure
lastlevel=0;
for n=1:numel(clusters)
    if mod(n,2)==0, col='b';else col='k';end
    [R lastlevel]=manyRasters(relativeSpikeTimes(n,3).sp,col,.9,ppms,0+lastlevel);
end

