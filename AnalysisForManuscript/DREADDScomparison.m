%%  1. load files
clear all
cd 'C:\Users\Rebecca Mease\Documents\GitHub\MUA_analysis\AnalysisForManuscript';
VPL_file_list %define all the files here

% get data; make sure to deal with different sampling frequencies!!!
filestr = '\TriggeredAnalysis.mat'
STRS={'CFA_cortex','control_Gi','control_Gq'}
S={};
for I=1:3
    S{I}{1}=load([cortex_mu{I}  filestr])
end



%% 2. pooled comparison across animals, not paired!

%get spont rate chunks
%convert to rates
%get spont ISIs corresponding to window
%convert to MS

%how many neurons total?

%merge relative spike times structures across recordings
RSP={};PPMS={};
for I=1:numel(S)
    rsp=struct([])
    ppms=[];
    for i=1:numel(S{I})
        rsp=[rsp;S{I}{i}.relativeSpikeTimes];
        ppms=[ppms;S{I}{i}.ppms*ones(size(S{I}{i}.relativeSpikeTimes,1),1)];
    end
    RSP{I}=rsp;
    PPMS{I}=ppms;
end


% for spontaneous activity, we take all counts before the trigger,
% regardless of triggers

%% DREADDS
results=struct([])

% standard order is   [1 Hz 10 Hz laser mech]

windowsize=1
pads=[.2 .2 0 .2]  %don't need to compensate for lag

%Gather everything

I=3   %INHIBITORY
rsp=RSP{I}; %relative spike times
nUnits=size(rsp,1)
ISIs=cell(nUnits,2);C=cell(nUnits,2);%n neuron by 2 conditions (before and after)
Sp=cell(nUnits,2);

Cs={};
ISIS={};


%% just one experimental condition, before and after
for j=1:2
    for n=1:nUnits  %for all units
        ppms=PPMS{I}(n);
        Sp={};Sp2={};Isis={};Isis2={};
        sp=rsp(n,j).sp;
        c=[];c2=[];
        for ii=1:numel(sp) %for each trial
            
            %before
            spikes=rsp(n,j).sp{ii}/(ppms*1000); %one unit, one condition, one trial samples --> spike times
            indices=find(spikes>windowC(1) & spikes<=windowC(2));%which spikes are in control window?
            isis=rsp(n,j).isis{ii}*1000; %ISIs in milliseconds
            Isis{ii}=isis(indices);
            Sp{ii}=spikes(indices);
            c(ii)=numel(indices);
            
            %after
            indices=find(spikes>windowE(1) & spikes<=windowE(2));%which spikes are in exp? window?
            isis=rsp(n,j).isis{ii}*1000; %ISIs in milliseconds
            Isis2{ii}=isis(indices);
            Sp2{ii}=spikes(indices);
            c2(ii)=numel(indices);
            
        end
        C{n,1}=c; C{n,2}=c2;
        ISIs{n,1}=Isis;ISIs{n,2}=Isis2;
    end
    
    Cs{j}=C;
    ISIS{j}=ISIs;
    
end

%%  compare mech responses before and after
figure
%exp_bi
%spont_bi

% first rank-sum (paired by pairs of observations) for each experimental
% condition
close all
fig=figure
f2=figure;

fractions=nan(3,2);%no sig change, smaller, bigger
f3=figure;
Deltas={};
EFs={};CFs={};
categories={}
for I=[2]
    nUnits=size(Cs{I},1)
    p=nan(nUnits,1);
    deltas=nan(nUnits,1);
    
    eF=[];
    cF=[];
    for n=1:nUnits
        p(n) =  signrank(Cs{I}{n,1}',Cs{I}{n,2}'); % sign rank test
        deltas(n)=median(Cs{I}{n,2}')-median(Cs{I}{n,1}');  %experimental - control!
    end
    
    eF=cellfun(@median, Cs{I}(:,2))/windowsize
    cF=cellfun(@median, Cs{I}(:,1))/windowsize
    
    EFs{I}=eF;
    CFs{I}=cF;
    H=p;H(H>.05)=0;H(H~=0)=1;  %each column a different experimental condition
    
    decrease=(deltas<0 & H==1)
    increase=(deltas>0 & H==1)
    nochange=H==0;
    Deltas{I}=(eF-cF);
    %Deltas{I}=Deltas{I}(increase);
    figure(f2)
    
    if I==1, col='r'; else col='b';end
    plot(cF(nochange), eF(nochange),['.'],'color',[.7 .7 .7 ])
    hold on
    plot(cF(increase), eF(increase),[col 'o'])
    plot(cF(decrease), eF(decrease),[col 'o'])

    plot([0:.01:50], [0:.01:50],'-k')
    xlim([0 18])
    ylim([0 18])
    axis square
    ylabel 'Evoked Firing [Hz]'
    xlabel 'Spontaneous Firing [Hz]'
    
    
    median(Cs{I}{n,2}')-median(Cs{I}{n,1}');
    %
    
    fractions(1,I)=numel(find(H==0)); %not sig
    fractions(2,I)=numel(find(deltas<0 & H==1)); % sig decrease
    fractions(3,I)=sum(deltas>0 & H==1); %sig increase
    
    figure(fig)
    subplot(1,2,I)
    mypie=pie(fractions(:,I))
    title(rsp(I).condition)
    
    %legend('not sig. p>.05 W sr','decrease','increase')
    
    figure(f3)
    subplot(1,2,I)
    b=plot([1 2]+I-1,[cF eF],':','color',[.7 .7 .7])
    hold on
    b=plot([1 2]+I-1,[cF(increase) eF(increase)],'-o','color',[.7 .7 .7]*0)
    
    % b=boxplot([cF eF],'Notch','on','plotstyle','traditional','datalim',[0 8],'Widths',0.7,'labels',{'Spont','Mech'})
    
    
    hold on
    hold on
    xlim([.5 2.5]+I-1)
    box off
    ylim([0 50])
    % b=boxplot([cF eF],'Notch','on','plotstyle','traditional','datalim',[0 20],'Widths',0.7,'labels',{'Spont','Mech'})
    categories{I}=[nochange decrease increase];
end

%%
figure
M=squishcell(Deltas);
%violinplot(M,{'Saline','CFA'})
ylim([0 15])
b=boxplot(M,'Notch','on','plotstyle','traditional','datalim',[0 8],'Widths',0.7,'labels',{'control','CNO'})

%paired because we do this cell by cell
p=signrank(Deltas{1},Deltas{2})
p1=signrank(EFs{1},EFs{2})
p2=signrank(CFs{1},CFs{2})



%% NOW SPONTANEOUS RESPONSES ONLY



%%  compare mech responses before and after
figure
%exp_bi
%spont_bi
I=1;
% first rank-sum (paired by pairs of observations) for each experimental
% condition
close all
fig=figure
f2=figure;

fractions=nan(3,2);%no sig change, smaller, bigger
f3=figure;
Deltas={};
EFs={};CFs={};
categories={}

nUnits=size(Cs{I},1)
p=nan(nUnits,1);
deltas=nan(nUnits,1);

eF=[];
cF=[];
for n=1:nUnits
    exp=Cs{2}{n,1};control=Cs{1}{n,1};
    exp=exp(1:30);
    control=control(1:30);
    p(n) =  signrank(exp, control); % sign rank test
    deltas(n)=median(exp)-median(control);  %experimental - control!
    eF(n)=mean(exp);
    cF(n)=mean(control);
end


EFs{I}=eF;
CFs{I}=cF;
H=p;H(H>.05)=0;H(H~=0)=1;  %each column a different experimental condition

decrease=(deltas<0 & H==1)
increase=(deltas>0 & H==1)
nochange=H==0;
Deltas{I}=(eF-cF);
%Deltas{I}=Deltas{I}(increase);
figure(f2)
subplot(1,2,1)
if I==2, col='r'; else col='b';end
plot(cF(nochange), eF(nochange),['.'],'color',[.7 .7 .7 ])
hold on
plot(cF(increase), eF(increase),[col 'o'])
plot(cF(decrease), eF(decrease),[col 'o'])

plot([0:.01:50], [0:.01:50],'-k')
xlim([0 20])
ylim([0 20])
axis square
ylabel 'post-CNO Firing [Hz]'
xlabel 'pre-CNO  Firing [Hz]'


median(Cs{I}{n,2}')-median(Cs{I}{n,1}');
%
subplot(1,2,2)
fractions(1,I)=numel(find(H==0)); %not sig
fractions(2,I)=numel(find(deltas<0 & H==1)); % sig decrease
fractions(3,I)=sum(deltas>0 & H==1); %sig increase
mypie=pie(fractions(:,I))
title(rsp(I).condition)

%legend('not sig. p>.05 W sr','decrease','increase')

figure(f3)
subplot(1,2,I)
b=plot([1 2]+I-1,[cF eF],':','color',[.7 .7 .7])
hold on
b=plot([1 2]+I-1,[cF(increase) eF(increase)],'-o','color',[.7 .7 .7]*0)

b=boxplot([cF' eF'],'Notch','on','plotstyle','traditional','datalim',[0 20],'Widths',0.7,'labels',{'Spont','Mech'})


hold on
hold on
xlim([.5 2.5]+I-1)
box off
ylim([0 50])
% b=boxplot([cF eF],'Notch','on','plotstyle','traditional','datalim',[0 20],'Widths',0.7,'labels',{'Spont','Mech'})
categories{I}=[nochange decrease increase];

p=signrank(Deltas{1},Deltas{2})
p1=signrank(EFs{1},CFs{1})  %this is meaningless though


