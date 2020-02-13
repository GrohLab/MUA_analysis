
%%  1. load files
clear all
cd 'C:\Users\Rebecca Mease\Documents\GitHub\MUA_analysis\AnalysisForManuscript';
VPL_file_list %define all the files here

% get data; make sure to deal with different sampling frequencies!!!
filestr = '\TriggeredAnalysis.mat'
STRS={'CFA_cortex','control_Gi','control_Gq'}
S={};
S{1}{1}=load([cortex_mu{2}  filestr])
S{2}{1}=load([cortex_mu{1}  filestr])




%% 2. pooled comparison across animals, not paired!
S{2}{1}.relativeSpikeTimes(:,1)=S{2}{1}.relativeSpikeTimes(:,4);
S{2}{1}.catSp(:,1)=S{2}{1}.catSp(:,4);
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






%% spontaneous window for CFA and Saline
window=[-1 0]; %window of interest, in seconds;

Cs={}; Trials={};SP={};
for I=1:2
    rsp=RSP{I};
    nUnits=size(rsp,1)
    nCond=size(rsp,2);
    ISIs=cell(size(rsp));C=cell(size(rsp));%counts
    Sp=cell(size(rsp));%counts
    for j=1:1
        for n=1:nUnits
            ppms=PPMS{I}(n);
            Sp={};Isis={};
            sp=rsp(n,j).sp;
            c=[];
            for ii=1:numel(sp)
                spikes=rsp(n,j).sp{ii}/(ppms*1000); %one unit, one condition, one trial samples --> spike times
                indices=find(spikes>window(1) & spikes<=window(2));%which spikes are in desired window?
                isis=rsp(n,j).isis{ii}*1000; %ISIs in milliseconds
                Isis{ii}=isis(indices);
                Sp{ii}=spikes(indices);
                c(ii)=numel(indices);
            end
            C{n,j}=c;ISIs{n,j}=Isis;
        end
    end
    %Trials{I}=cell2mat(cellfun(@numel,SP,'UniformOutput',0));
    ISIS{I}=ISIs;
    Cs{I}=C;
end


Saline={}
I=1;
nUnits=size(Cs{I},1)
nCond=size(Cs{I},2);
for n=1:nUnits
    Saline{n}=cell2mat(Cs{I}(n,1:nCond));
end
MMsaline=[cellfun(@median,Saline)' cellfun(@mean,Saline)']

CFA={}
I=2;
nUnits=size(Cs{I},1)
nCond=size(Cs{I},2);
for n=1:nUnits
    CFA{n}=cell2mat(Cs{I}(n,1:nCond));
end
MMcfa=[cellfun(@median,CFA)' cellfun(@mean,CFA)']

% MAKES THE FIGURE


close all
figure
subplot(2,1,2)
data=nan(max([size(MMcfa,1) size(MMsaline,1)]),2);
data(1:size(MMsaline,1),1)=MMsaline(:,2);
data(1:size(MMcfa,1),2)=MMcfa(:,2);

p=ranksum(MMsaline(:,1),MMcfa(:,1))
p=ranksum(MMsaline(:,2),MMcfa(:,2))


b=boxplot(data,'Notch','on','plotstyle','traditional','datalim',[0 15],'orientation','horizontal','Widths',0.7,'labels',{'Saline','CFA'})

sact=(sum(MMsaline(:,2)>1)/size(MMsaline(:,2),1))
fract=[1-sact sact]
sact=sum(MMcfa(:,2)>1)/size(MMcfa(:,2),1)
fract=[fract;[1-sact sact]]

xlabel 'Spontaneous Firing Rate [Hz]'

%subplot(1,3,1);pie(fract(1,:)); title(['Saline n= ' num2str(size(MMsaline(:,1),1)) ' units' ])
%subplot(1,3,2);pie(fract(2,:)); title(['CFA n= ' num2str(size(MMcfa(:,1),1)) ' units' ])


bins=[0:.1:40];
h1=hist(MMsaline(:,2),bins)
h2=hist(MMcfa(:,2),bins)

h1=(h1)/sum(h1);
h2=(h2)/sum(h2);

h1=cumsum(h1)/sum(h1);
h2=cumsum(h2)/sum(h2);


subplot(2,1,1)
plot(bins,h1,bins,h2)
xlim([0 35])

%% interspike intervals  CDF compare saline and cfa
X={}

for I=1:2;
    isis=ISIS{I};
    for n=1:size(isis,1)
        x=isis(n,:);   %can pool across conditions
        xx=[]
        for i=1:numel(x)
            
            new=squishcell(x{i});new=new(~isnan(new));
            xx=[xx; new(:)];
        end
        X{n,I}=xx;
    end
end

bins=[1:100000];
H=cellfun(@(c) cumsum(hist(c,bins))/numel(c),X,'UniformOutput',0);

[h,p] =kstest2(cell2mat(X(:,2)),cell2mat(X(:,1)))  %should be sampled evenly from all neurons....
Hsaline=cell2mat(H(:,1));Hcfa=cell2mat(H(:,2));

figure
plot(log10(bins),nanmean(Hsaline),':k','linewidth', 2)
hold on
plot(log10(bins),nanmean(Hcfa),':r','linewidth', 2)
ylabel('CDF(ISI)')
xlabel('Log_{10} Interspike interval [ms]')
%%    now compare mechanical response
% standard order is   [1 Hz 10 Hz laser mech]
% we do mar
windowsize=1
windowC=[-windowsize 0]; %window of interest, in seconds;
pad=.2
windowE=[pad windowsize+pad]; %window of interest, in seconds;
Cs={}; Trials={};
ISIS={};
for I=1:2  %for saline and CFA
    rsp=RSP{I}; %relative spike times
    nUnits=size(rsp,1)
    ISIs=cell(nUnits,2);C=cell(nUnits,2);%n neuron by 2 conditions (before and after)
    Sp=cell(nUnits,2);
    j=1  %just one stimulus
    for n=1:nUnits  %for all units
        ppms=PPMS{I}(n);
        Sp={};Sp2={};Isis={};Isis2={};
        sp=rsp(n,j).sp;
        c=[];c2=[];
        for ii=1:numel(sp) %for each trial
            spikes=rsp(n,j).sp{ii}/(ppms*1000); %one unit, one condition, one trial samples --> spike times
            indices=find(spikes>windowC(1) & spikes<=windowC(2));%which spikes are in control window?
            isis=rsp(n,j).isis{ii}*1000; %ISIs in milliseconds
            Isis{ii}=isis(indices);
            Sp{ii}=spikes(indices);
            c(ii)=numel(indices);
            
            indices=find(spikes>windowE(1) & spikes<=windowE(2));%which spikes are in exp? window?
            isis=rsp(n,j).isis{ii}*1000; %ISIs in milliseconds
            Isis2{ii}=isis(indices);
            Sp2{ii}=spikes(indices);
            c2(ii)=numel(indices);
            
        end
        C{n,1}=c; C{n,2}=c2;
        ISIs{n,1}=Isis;ISIs{n,2}=Isis2;
    end
    
    Trials{I}=cell2mat(cellfun(@numel,SP,'UniformOutput',0));
    %median spike count per neuron per condition
    Cs{I}=C;
    ISIS{I}=ISIs;
end


%%  ISIs for mechanical response
X={}

for I=1:2
    isis=ISIS{I};  %saline/CFA
    for j=1:size(isis,2)
        for n=1:size(isis,1) %for spont and mech (before and after stim)
            x=isis(n,j);   %one neuron, one condition, all trials
            xx=[]
            for i=1:numel(x)
                new=squishcell(x{i});new=new(~isnan(new));
                xx=[xx; new(:)];
            end
            X{I}{n,j}=xx;  %gather everything
        end
    end
end



%%
bins=[1:100000];
Hall_saline=cellfun(@(c) cumsum(hist(c,bins))/numel(c),X{1},'UniformOutput',0);
Hall_cfa=cellfun(@(c) cumsum(hist(c,bins))/numel(c),X{2},'UniformOutput',0);


%test difference of ISIs for mechanical response
[h,p] =kstest2(cell2mat(X{1}(:,2)),cell2mat(X{2}(:,2)))  %should be sampled evenly from all neurons....

Hsaline=cell2mat(Hall_saline(:,2));Hcfa=cell2mat(Hall_cfa(:,2)); %get cdfs

figure
hold on
plot(log10(bins),nanmean(Hsaline),'k','linewidth', 2)
hold on
plot(log10(bins),nanmean(Hcfa),'r','linewidth', 2)
ylabel('CDF(ISI)')
xlabel('Log_{10} Interspike interval [ms]')

title 'mech evoked'

%%  now do statistical tests

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
EFs={};
for I=[2 1]
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
    
    EFS{I}=eF;
    H=p;H(H>.05)=0;H(H~=0)=1;  %each column a different experimental condition
    
    decrease=(deltas<0 & H==1)
    increase=(deltas>0 & H==1)
    nochange=H==0;
    Deltas{I}=(eF-cF);
    Deltas{I}=Deltas{I}(increase);
    figure(f2)
    
    if I==2, col='r'; else col='b';end
    plot(cF(nochange), eF(nochange),['.'],'color',[.7 .7 .7 ])
    hold on
    plot(cF(increase), eF(increase),[col 'o'])
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
    title(STRS{I})
    
    %legend('not sig. p>.05 W sr','decrease','increase')
    
    figure(f3)
    subplot(1,2,I)
    b=plot([1 2]+I-1,[cF eF],':','color',[.7 .7 .7])
    hold on
    b=plot([1 2]+I-1,[cF(increase) eF(increase)],'-o','color',[.7 .7 .7]*0)
    
    b=boxplot([cF eF],'Notch','on','plotstyle','traditional','datalim',[0 8],'Widths',0.7,'labels',{'Spont','Mech'})
    
    
    hold on
    hold on
    xlim([0 3]+I-1)
    box off
    ylim([0 50])
    % b=boxplot([cF eF],'Notch','on','plotstyle','traditional','datalim',[0 20],'Widths',0.7,'labels',{'Spont','Mech'})
    
end

figure
M=squishcell(Deltas);
%violinplot(M,{'Saline','CFA'})
ylim([0 15])
b=boxplot(M,'Notch','on','plotstyle','traditional','datalim',[0 8],'Widths',0.7,'labels',{'Saline','CFA'})
p=ranksum(Deltas{1},Deltas{2})
p1=ranksum(EFS{1},EFS{2})




%% now we look at the effect of light alone on spontaneous activity



%%    now compare mechanical response
% standard order is   [1 Hz 10 Hz laser mech]
% we do mar
windowsize=1
windowC=[-windowsize 0]; %window of interest, in seconds;
pad=.2
windowE=[pad windowsize+pad]; %window of interest, in seconds;

Cs={}; Trials={};
for I=1:2  %for saline and CFA
    rsp=RSP{I}; %relative spike times
    nUnits=size(rsp,1)
    ISIs=cell(size(rsp));C=cell(nUnits,2);%n neuron by 2 conditions (before and after)
    Sp=cell(nUnits,2);
    j=4  %just one stimulus
    for n=1:nUnits  %for all units
        ppms=PPMS{I}(n);
        Sp={};Sp2={};Isis={};Isis2={};
        sp=rsp(n,j).sp;
        c=[];c2=[];
        for ii=1:numel(sp) %for each trial
            spikes=rsp(n,j).sp{ii}/(ppms*1000); %one unit, one condition, one trial samples --> spike times
            indices=find(spikes>windowC(1) & spikes<=windowC(2));%which spikes are in control window?
            isis=rsp(n,j).isis{ii}/ppms; %ISIs;
            Isis{ii}=isis(indices);
            Sp{ii}=spikes(indices);
            c(ii)=numel(indices);
            
            indices=find(spikes>windowE(1) & spikes<=windowE(2));%which spikes are in exp? window?
            isis=rsp(n,j).isis{ii}/ppms; %ISIs;
            Isis2{ii}=isis(indices);
            Sp2{ii}=spikes(indices);
            c2(ii)=numel(indices);
            
        end
        C{n,1}=c; C{n,2}=c2;
        ISIs{n,1}=Isis;ISIs{n,2}=Isis2;
    end
    
    Trials{I}=cell2mat(cellfun(@numel,SP,'UniformOutput',0));
    %median spike count per neuron per condition
    Cs{I}=C;
end




%%  now do statistical tests

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
EFs={};
for I=[2 1]
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
    
    EFS{I}=eF;
    H=p;H(H>.05)=0;H(H~=0)=1;  %each column a different experimental condition
    
    decrease=(deltas<0 & H==1)
    increase=(deltas>0 & H==1)
    nochange=H==0;
    Deltas{I}=(eF-cF);
    Deltas{I}=Deltas{I}(increase);
    figure(f2)
    
    if I==2, col='r'; else col='b';end
    plot(cF(nochange), eF(nochange),['.'],'color',[.7 .7 .7 ])
    hold on
    plot(cF(increase), eF(increase),[col 'o'])
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
    title(STRS{I})
    
    %legend('not sig. p>.05 W sr','decrease','increase')
    
    figure(f3)
    subplot(1,2,I)
    b=plot([1 2]+I-1,[cF eF],':','color',[.7 .7 .7])
    hold on
    b=plot([1 2]+I-1,[cF(increase) eF(increase)],'-o','color',[.7 .7 .7]*0)
    
    b=boxplot([cF eF],'Notch','on','plotstyle','traditional','datalim',[0 8],'Widths',0.7,'labels',{'Spont','Mech'})
    
    
    hold on
    hold on
    xlim([0 3]+I-1)
    box off
    ylim([0 50])
    % b=boxplot([cF eF],'Notch','on','plotstyle','traditional','datalim',[0 20],'Widths',0.7,'labels',{'Spont','Mech'})
    
end

figure
M=squishcell(Deltas);
%violinplot(M,{'Saline','CFA'})
ylim([0 15])
b=boxplot(M,'Notch','on','plotstyle','traditional','datalim',[0 8],'Widths',0.7,'labels',{'Saline','CFA'})
p=ranksum(Deltas{1},Deltas{2})
p1=ranksum(EFS{1},EFS{2})



%%    now compare LIGHT response, condition 3
% standard order is   [1 Hz 10 Hz laser mech]

windowsize=1.5
windowC=[-windowsize 0]; %window of interest, in seconds;
pad=0   %don't need to compensate for lag
windowE=[pad windowsize+pad]; %window of interest, in seconds;
Cs={}; Trials={};
ISIS={};
for I=1:2  %for saline and CFA
    rsp=RSP{I}; %relative spike times
    nUnits=size(rsp,1)
    ISIs=cell(nUnits,2);C=cell(nUnits,2);%n neuron by 2 conditions (before and after)
    Sp=cell(nUnits,2);
    j=3  %just one stimulus
    for n=1:nUnits  %for all units
        ppms=PPMS{I}(n);
        Sp={};Sp2={};Isis={};Isis2={};
        sp=rsp(n,j).sp;
        c=[];c2=[];
        for ii=1:numel(sp) %for each trial
            spikes=rsp(n,j).sp{ii}/(ppms*1000); %one unit, one condition, one trial samples --> spike times
            indices=find(spikes>windowC(1) & spikes<=windowC(2));%which spikes are in control window?
            isis=rsp(n,j).isis{ii}*1000; %ISIs in milliseconds
            Isis{ii}=isis(indices);
            Sp{ii}=spikes(indices);
            c(ii)=numel(indices);
            
            indices=find(spikes>windowE(1) & spikes<=windowE(2));%which spikes are in exp? window?
            isis=rsp(n,j).isis{ii}*1000; %ISIs in milliseconds
            Isis2{ii}=isis(indices);
            Sp2{ii}=spikes(indices);
            c2(ii)=numel(indices);
            
        end
        C{n,1}=c; C{n,2}=c2;
        ISIs{n,1}=Isis;ISIs{n,2}=Isis2;
    end
    
    Trials{I}=cell2mat(cellfun(@numel,SP,'UniformOutput',0));
    %median spike count per neuron per condition
    Cs{I}=C;
    ISIS{I}=ISIs;
end


%%  ISIs for light response
X={}

for I=1:2
    isis=ISIS{I};  %saline/CFA
    for j=1:size(isis,2)
        for n=1:size(isis,1) %for spont and mech (before and after stim)
            x=isis(n,j);   %one neuron, one condition, all trials
            xx=[]
            for i=1:numel(x)
                new=squishcell(x{i});new=new(~isnan(new));
                xx=[xx; new(:)];
            end
            X{I}{n,j}=xx;  %gather everything
        end
    end
end



%%
bins=[1:100000];
Hall_saline=cellfun(@(c) cumsum(hist(c,bins))/numel(c),X{1},'UniformOutput',0);
Hall_cfa=cellfun(@(c) cumsum(hist(c,bins))/numel(c),X{2},'UniformOutput',0);


%test difference of ISIs for light response
[h,p] =kstest2(cell2mat(X{1}(:,2)),cell2mat(X{2}(:,2)))  %should be sampled evenly from all neurons....

Hsaline=cell2mat(Hall_saline(:,2));Hcfa=cell2mat(Hall_cfa(:,2)); %get cdfs

figure
hold on
plot(log10(bins),nanmean(Hsaline),'k','linewidth', 2)
hold on
plot(log10(bins),nanmean(Hcfa),'r','linewidth', 2)
ylabel('CDF(ISI)')
xlabel('Log_{10} Interspike interval [ms]')

title 'light evoked'

%%  now do statistical tests

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
EFs={};
for I=[2 1]
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
    
    EFS{I}=eF;
    H=p;H(H>.05)=0;H(H~=0)=1;  %each column a different experimental condition
    
    decrease=(deltas<0 & H==1)
    increase=(deltas>0 & H==1)
    nochange=H==0;
    Deltas{I}=(eF-cF);
    Deltas{I}=Deltas{I}(increase);
    figure(f2)
    
    if I==2, col='r'; else col='b';end
    plot(cF(nochange), eF(nochange),['.'],'color',[.7 .7 .7 ])
    hold on
    plot(cF(increase), eF(increase),[col 'o'])
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
    title(STRS{I})
    
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
    
end



%%
figure
M=squishcell(Deltas);
%violinplot(M,{'Saline','CFA'})
ylim([0 15])
b=boxplot(M,'Notch','on','plotstyle','traditional','datalim',[0 8],'Widths',0.7,'labels',{'Saline','CFA'})
p=ranksum(Deltas{1},Deltas{2})
p1=ranksum(EFS{1},EFS{2})

