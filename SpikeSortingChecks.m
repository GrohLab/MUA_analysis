
%% exploratory script for looking at MUA recordings
%% find every single trigger for step stimuli
ppms=fs/1000
stimulus=lazynorm(Triggers.whisker);
T=numel(Triggers.whisker)/ppms/1000; %duration of experiment in seconds!
figure
tmech=diff(lazynorm(Triggers.whisker));tmech=find(tmech>.5);
tlight=diff(lazynorm(Triggers.light));tlight=find(tlight>.5);
tlight_alone=Conditions{4}.Triggers;
t=1:numel(stimulus);t=t/ppms;
plot(t,lazynorm(Triggers.whisker),'r',t,lazynorm(Triggers.light),'b')
hold on
plot(tmech/ppms,ones(size(tmech)),'ro',tlight/ppms,ones(size(tlight)),'ob',tlight_alone/ppms,ones(size(tlight_alone)),'c*')
cd(homedir)
save AllTriggers tmech tlight tlight_alone

%%
gray=[.7 .7 .7]
N=numel(sp)
% get color for each neuron, just for plotting
colors=brighten(distinguishable_colors(N),-.25);%color for each unit

%% rasters by Unit: how do they look over time in comparison to stimulus
t=1:numel(stimulus);t=t/ppms;
figure
subplot(2,2,1:2)
plot(t,lazynorm(stimulus),'color',gray)
hold on
for i=1:numel(sp)
    yval=ones(size(sp{i}))*i;
    plot(sp{i},yval,'.','color',colors(i,:),'markersize',15)
    hold on
    xlabel ms;
end

% firing rates by cell-- is this reasonable?
%calculate firing rates
Rates=cell2mat(cellfun(@numel,sp,'UniformOutput', false))/T;
[rh bins]=hist(Rates,50);
subplot(2,2,3)
bar(bins,rh,1)
ylabel units
xlabel 'Spiking rate Hz'
title([LabelStr '    ' num2str(T) '=time in seconds'])
subplot(2,2,4)
for i=1:N
    stem(i,Rates(i),'linewidth',2,'color',colors(i,:))
    hold on
end
ylabel 'Spiking rate Hz'
%% look at ISIs...are any very precisely locked to stimulation frequency?
Isis=cellfun(@(x) [x(1); diff(x)],sp,'UniformOutput', false); %ISIs
dt=.5%ms
[h,bins]=hist(cell2mat(Isis),0:dt:max(cell2mat(Isis)));%get bins for all
Bins=cell(size(Isis));[Bins{:}]=deal(bins); %bins for all
[H Bins]=cellfun(@hist,Isis,Bins,'UniformOutput', false); %histogram for all

figure
n=ceil(sqrt(numel(Bins)));
for i=1:numel(Bins)
    subplot(n,n,i)
    plot((Bins{i}),H{i},'linewidth',2,'color',colors(i,:))
    title(['unit ' num2str(i) ' ' units{i}], 'interpreter','none')
    xlim([0 10])
    set(gca,'xtick',[0:.25:1 2:10])
    grid on
end

%normalize H
H=cellfun(@(x) x/sum(x),H,'UniformOutput', false)



figure
for i=1:numel(Bins)
    subplot(n,n,i)
    plot(log10(Bins{i}),(H{i}),'linewidth',2,'color',colors(i,:))
    title(['unit ' num2str(i) ' ' units{i}], 'interpreter','none')
    xlim([.01 6])
end

figure
for i=1:numel(Bins)
    subplot(n,n,i)
    plot(log10(Bins{i}),cumsum(H{i}),'linewidth',2,'color',colors(i,:))
    title(['unit ' num2str(i) ' ' units{i}], 'interpreter','none')
    xlim([.01 6])
end


%% cross correlation between spike times
% %this should find 1)double-counted spikes (identified in different
% %groupings, and 2) spikes that are part of a burst but separated as
% %different neurons due to waveform differences...
%
maxlag=50
lag=[-maxlag:maxlag];
if exist([homedir 'CrossCorrelation.mat'])==2
    load CrossCorrelation
    display('loaded cross correlation data')
else   %if it hasn't been run yet, then run i
    Z=zeros(N,round(T*1000));%ms precision
    
    for i=1:N
        index=round(sp{i});index(index==0)=1;
        Z(i,index)=1;
    end
    Z2=Z;
    indices=1:size(Z,2);
    r=xcorr(Z(1,indices),Z(2,indices),maxlag,'coeff');
    R=nan(N,N,maxlag*2+1);
    for i=1:N
        parfor j=1:N
            if j>i
                R(i,j,:)=xcorr(Z(i,indices),Z2(j,indices),maxlag,'coeff');
            end
        end
        display(i/N)
    end
    beep
    save CrossCorrelation R maxlag lag
end
%% find highly correlated units, determine threshold
figure
for i=1:N
    hold on
    plot(lag,squeeze(R(1,:,:)),'linewidth',2)
    hold on
end
thresh=0.02% is the correlation ever greater than this value...
rsum=sum(R>=thresh,3);
[I J]=ind2sub(size(rsum),find(rsum>0))

%% look at highly correlated units-- events and xcorr PLUS WAVEFORMS???
figure
for ii=1:numel(I)
    i=I(ii);j=J(ii);
    
    yvali=ones(size(sp{i}))*1.1;
    yvalj=ones(size(sp{j}));
    subplot(2,2,1:2)
    plot(sp{i},yvali,'o','markersize',5,'linewidth',2,'color',colors(i,:))
    hold on
    plot(sp{j},yvalj,'.','markersize',15,'color',colors(j,:))
    title([num2str(i) ' ' num2str(j)])
    xlabel ms;
    l=legend(units{i}, units{j})
    set(l,'interpreter','none')
    ylim([0 3])
    grid on
    hold off
    subplot(2,2,3)
    plot(lag,squeeze(R(i,j,:)))
    
    subplot(2,2,4)
    title('average waveforms')
    [units{i}, units{j}]
    pause
    hold off
end

%% now look at all triggers times, regardless of conditions FOR LIGHT
timeBefore=5;timeAfter=35;
binsize=.01;
triggers=tlight;

%all light triggers IN TIME for spikes
[TriggeredSpikeTimes Xvals Yvals]=TriggeredSpikes(sp,triggers/ppms,timeBefore, timeAfter);  %get all spike times for a single trigger

%all
LSpikes=cellfun(@(x) cell2mat(x'),TriggeredSpikeTimes,'UniformOutput', false);
bins=[-timeBefore:binsize:timeAfter];
Bins=cell(size(LSpikes));[Bins{:}]=deal(bins); %bins for all
[H Bins]=cellfun(@hist,LSpikes,Bins,'UniformOutput', false); %histogram for all

Xs=ManyTriggeredSegments({lazynorm(Triggers.light),lazynorm(Triggers.whisker)},triggers,timeBefore*ppms, timeAfter*ppms);
trig_t=[1:size(Xs{1},1)]/ppms-timeBefore;
%means of continuos variables
Xmeans = cellfun(@transpose, Xs,'UniformOutput', false);
MeanStimuli = cellfun(@mean, Xmeans,'UniformOutput', false);
Ymax=max(cellfun(@max,H));
%

figure
for i=1:N
    subplot(n,n,i)
    plot(trig_t,MeanStimuli{1}*Ymax,'c','linewidth',2)
    hold on
    plot(trig_t,MeanStimuli{2}*Ymax,'color',gray,'linewidth',2)
    bar(bins,H{i},1,'facecolor',colors(i,:),'edgecolor',colors(i,:))
    title(i)
    axis tight
    %ylim([0 Ymax*1.25])
end
legend('light','mech')

%light only only
xs=cellfun(@sum,Xs,'UniformOutput', false)
lightControl=find(xs{2}<1)

timeBefore=5;timeAfter=25;
binsize=.5;
triggers=tlight(lightControl);

%all light triggers IN TIME for spikes
[TriggeredSpikeTimes Xvals Yvals]=TriggeredSpikes(sp,triggers/ppms,timeBefore, timeAfter);  %get all spike times for a single trigger

%all
LSpikes=cellfun(@(x) cell2mat(x'),TriggeredSpikeTimes,'UniformOutput', false);
bins=[-timeBefore:binsize:timeAfter];
Bins=cell(size(LSpikes));[Bins{:}]=deal(bins); %bins for all
[H Bins]=cellfun(@hist,LSpikes,Bins,'UniformOutput', false); %histogram for all


Xs=ManyTriggeredSegments({lazynorm(Triggers.light),lazynorm(Triggers.whisker)},triggers,timeBefore*ppms, timeAfter*ppms);
trig_t=[1:size(Xs{1},1)]/ppms-timeBefore;
%means of continuos variables
Xmeans = cellfun(@transpose, Xs,'UniformOutput', false);
MeanStimuli = cellfun(@mean, Xmeans,'UniformOutput', false);
Ymax=max(cellfun(@max,H));
%

figure
for i=1:N
   
    subplot(n,n,i)
    plot(trig_t,MeanStimuli{1}*Ymax,'c','linewidth',2)
    hold on
    plot(trig_t,MeanStimuli{2}*Ymax,'color',gray,'linewidth',2)
    bar(bins,H{i},1,'facecolor',colors(i,:),'edgecolor',colors(i,:))
    title(i)
    
    axis tight
    ylim([0 Ymax*1.25])
end
legend('light','mech')



%% now look at all triggers times, regardless of conditions FOR MECH
timeBefore=100;timeAfter=600;
binsize=10;
triggers=tmech;

%all light triggers IN TIME for spikes
[TriggeredSpikeTimes Xvals Yvals]=TriggeredSpikes(sp,triggers/ppms,timeBefore, timeAfter);  %get all spike times for a single trigger

%all
MSpikes=cellfun(@(x) cell2mat(x'),TriggeredSpikeTimes,'UniformOutput', false);
bins=[-timeBefore:binsize:timeAfter];
Bins=cell(size(LSpikes));[Bins{:}]=deal(bins); %bins for all
[H Bins]=cellfun(@hist,LSpikes,Bins,'UniformOutput', false); %histogram for all


Xs=ManyTriggeredSegments({lazynorm(Triggers.light),lazynorm(Triggers.whisker)},triggers,timeBefore*ppms, timeAfter*ppms);
trig_t=[1:size(Xs{1},1)]/ppms-timeBefore;
%means of continuos variables
Xmeans = cellfun(@transpose, Xs,'UniformOutput', false);
MeanStimuli = cellfun(@mean, Xmeans,'UniformOutput', false);
Ymax=max(cellfun(@max,H));
%

figure
for i=1:N
   
    subplot(n,n,i)
    plot(trig_t,MeanStimuli{1}*Ymax,'c','linewidth',2)
    hold on
    plot(trig_t,MeanStimuli{2}*Ymax,'color',gray,'linewidth',2)
    bar(bins,H{i},1,'facecolor',colors(i,:),'edgecolor',colors(i,:))
    title(i)
    
    axis tight
    ylim([0 Ymax*1.25])
end
legend('light','mech')

%mech only only
xs=cellfun(@sum,Xs,'UniformOutput', false)
mechControl=find(xs{1}<10)

binsize=1;
%triggers=tmech(mechControl);

%all light triggers IN TIME for spikes
[TriggeredSpikeTimes Xvals Yvals]=TriggeredSpikes(sp,triggers/ppms,timeBefore, timeAfter);  %get all spike times for a single trigger

%all
MSpikes=cellfun(@(x) cell2mat(x'),TriggeredSpikeTimes,'UniformOutput', false);
bins=[-timeBefore:binsize:timeAfter];
Bins=cell(size(LSpikes));[Bins{:}]=deal(bins); %bins for all
[H Bins]=cellfun(@hist,LSpikes,Bins,'UniformOutput', false); %histogram for all


Xs=ManyTriggeredSegments({lazynorm(Triggers.light),lazynorm(Triggers.whisker)},triggers,timeBefore*ppms, timeAfter*ppms);
trig_t=[1:size(Xs{1},1)]/ppms-timeBefore;
%means of continuos variables
Xmeans = cellfun(@transpose, Xs,'UniformOutput', false);
MeanStimuli = cellfun(@mean, Xmeans,'UniformOutput', false);
%Ymax=max(cellfun(@max,H));
%

figure
for i=1:N
   
    subplot(n,n,i)
    plot(trig_t,MeanStimuli{1}*Ymax,'c','linewidth',2)
    hold on
    plot(trig_t,MeanStimuli{2}*Ymax,'color',gray,'linewidth',2)
    bar(bins,H{i},1,'facecolor',colors(i,:),'edgecolor',colors(i,:))
    title(i)
    
    axis tight
    ylim([0 Ymax*1.25])
end
legend('light','mech')
