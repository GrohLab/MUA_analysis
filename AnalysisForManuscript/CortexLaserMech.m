 
%% light driven units in cortex

%at the moment, this has light artifacts Feb
datadir='Z:\Ross\Experiments\17.12.19\HL Cortex'
cd(datadir)
clear all
spikefile='M167_HL_Cortex_all_channels.mat'
load(spikefile);
load M167_HL_Cortex_sampling_frequency.mat  
load M167_HL_Cortexanalysis Conditions Triggers


stimCond=Conditions([1:4]);
ppms=fs/1000;
%% get single laser pulses
thislaser=Triggers.laser;
thislaser=thislaser-min(thislaser);thislaser=thislaser/max(thislaser);
thislaser(thislaser>.5)=1;
thislaser(thislaser~=1)=0;
laser=thislaser;

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
%thislaser=thislaser-whisker;
thislaser(thislaser~=1)=0;
lt=find(diff(thislaser)==1); %triggers for all lasers

figure
plot(lt/10,Triggers.laser(lt),'o')
hold on
plot(Triggers.laser(1:10:end))


j=numel(stimCond)+1;
stimCond(j).name='laser ind.'
stimCond(j).Triggers=[lt lt+round(ppms*10)];

%% make it simple and find everything 
% well isolated clusters, only isolated
goodsInd = cellfun(@(x) x==1,sortedData(:,3));multiInd=[];
multiInd = cellfun(@(x) x==2,sortedData(:,3));
goods = [find(goodsInd); find(multiInd)];
Ns = min(structfun(@numel,Triggers));
% Total duration of the recording
Nt = Ns/fs;  %seconds
clusters=sortedData(goods,2);ids=sortedData(goods,1);

biSpikes=spalloc(floor(Nt*1000*ppms),numel(goods),numel(cell2mat(clusters)));% preallocate sparse matrix to save space
% ms x num clusters
T=size(biSpikes,1); %recording length in ms




%stimCond(end-1).Triggers=stimCond(end-1).Triggers(1:5:end-1,:)
%make big matrix of all spike times
for j=1:numel(goods)
    tt=round(clusters{j}*1000*ppms);tt=tt(tt>0 & tt<=T);%spike times in ms*ppms (samples)
    S=zeros(T,1);S(tt)=1;
    biSpikes(:,j)=sparse(S);
end


%%  for all condtions
maxTrials=max(cellfun(@(x) size(x,1),{stimCond.Triggers}));
%in ms because of division before, though triggered seg takes samples
timeBefore=(-5*ppms);
timeAfter=round(35*ppms);
trialLength=-timeBefore+timeAfter+1;

%preallocate omg could be huge!
% firingRates: N x S x T x maxTrialNum
firingRates=zeros(numel(clusters),numel(stimCond),round(trialLength),maxTrials);
relativeSpikeTimes=struct([]);
catSp=cell(numel(clusters),numel(stimCond));
for n =1:numel(clusters)
    for s=1:numel(stimCond)
        triggers=round(stimCond(s).Triggers(:,1)/fs*1000*ppms); %triggers in ms*ppms
        X=TriggeredSegments(biSpikes(:,n),triggers,timeBefore,timeAfter);
        firingRates(n,s,1:size(X,1),1:size(X,2))=(shiftdim(X,-2));  %some fuckery to deal with 4d array
        sp={};tnum={};
        for ii=1:size(X,2),sp{ii}=find(X(:,ii))+timeBefore; tnum{ii}=ones(size(sp{ii}))*ii;end
        relativeSpikeTimes(n,s).sp=sp;
        relativeSpikeTimes(n,s).tnum=tnum;
        relativeSpikeTimes(n,s).clustername=ids(n);
        relativeSpikeTimes(n,s).condition=stimCond(s).name;
        catSp{n,s}=cell2mat(sp');
    end
end
%%   histograms
binsize=.05*ppms;
bins=[timeBefore:binsize:timeAfter];
H=cellfun(@(c) hist(c/ppms,bins/ppms),catSp,'UniformOutput',0);
Hcond=cell2mat(H(:,5));
figure
%Hcond(Hcond==0)=nan;
subplot(3,1,1:2)
pcolor(bins/ppms,1:numel(clusters),Hcond)

shading flat
colormap fire

X=mean(TriggeredSegments(laser,stimCond(5).Triggers(:,1),timeBefore,timeAfter)');
tx=[timeBefore:timeAfter]'
subplot(3,1,3)
plot(tx/ppms,X)
axis tight

laserMed=cell2mat(cellfun(@median,catSp(:,5),'UniformOutput',0))/ppms;
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
plot(tx/ppms,-X*100,'linewidth',3)
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

